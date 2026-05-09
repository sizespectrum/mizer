#' Animate size spectra or rates through time
#'
#' Creates an interactive plotly animation in which a play button steps through
#' time, drawing one line per species at each frame.
#'
#' The function dispatches on the class of `x`:
#'
#' * **`MizerSim`** — animates the community abundance spectra (number density
#'   or biomass density vs body size). Resource, background species, and a
#'   community total can be added via the `resource`, `background`, and `total`
#'   arguments. The `power` argument controls whether the y-axis shows number
#'   density (`power = 0`), biomass density (`power = 1`, default), or biomass
#'   density in logarithmic size bins (`power = 2`). Both axes are log10 by
#'   default and can each be switched to linear with `log_x = FALSE` or
#'   `log_y = FALSE`.
#'
#' * **`ArrayTimeBySpeciesBySize`** — animates any per-species, size-resolved
#'   rate returned by a `MizerSim` accessor, such as [getFMort()],
#'   [getFeedingLevel()], or [getPredMort()]. Both axes are log10 by default
#'   and can each be switched to linear with `log_x = FALSE` or `log_y = FALSE`.
#'   Species colours follow `params@linecolour`. Background species and a
#'   species total can be added via the `background` and `total` arguments.
#'
#' `animateSpectra()` is retained as a backward-compatible alias.
#'
#' @param x A `MizerSim` or `ArrayTimeBySpeciesBySize` object.
#' @param species Name or vector of names of the species to be plotted. By
#'   default all species are plotted.
#' @param time_range The time range to animate over. Either a vector of time
#'   values to include, or a length-two vector giving the min and max of the
#'   range. Default is the entire time range of `x`.
#' @param log_x If `TRUE` (default), use a log10 x-axis for body size.
#' @param log_y If `TRUE` (default), use a log10 y-axis.
#' @param total A boolean value that determines whether the total over all
#'   selected species is plotted as an additional trace called `"Total"`.
#'   Default is `FALSE`.
#' @param background A boolean value that determines whether background species
#'   are included. Ignored if the model does not contain background species.
#'   Default is `TRUE`.
#' @param wlim A numeric vector of length two providing lower and upper limits
#'   for the body-size (x) axis. Use `NA` to refer to the existing minimum or
#'   maximum.
#' @param ylim A numeric vector of length two providing lower and upper limits
#'   for the value (y) axis. Use `NA` to refer to the existing minimum or
#'   maximum. Values below `1e-20` are always cut off on log scales.
#' @param ... Additional arguments passed to the method.
#'
#' @return A plotly object with one animated line trace per plotted group. Use
#'   the play button or the slider to step through time.
#' @export
#' @family plotting functions
#' @examples
#' \donttest{
#' # Animate biomass density spectra, showing only sizes above 0.1 g
#' animate(NS_sim, power = 2, wlim = c(0.1, NA), time_range = 1997:2007)
#'
#' # Animate fishing mortality through time
#' animate(getFMort(NS_sim))
#'
#' # Animate feeding level for two species only
#' animate(getFeedingLevel(NS_sim), species = c("Cod", "Herring"))
#' }
animate <- function(x, ...) UseMethod("animate")

#' @rdname animate
#' @param power The abundance is plotted as the number density times the weight
#'   raised to \code{power}. The default \code{power = 1} gives the biomass
#'   density, whereas \code{power = 2} gives the biomass density with respect
#'   to logarithmic size bins. Only applies to `MizerSim`.
#' @param resource A boolean value that determines whether resource is included.
#'   If `TRUE`, the resource spectrum is plotted as an additional trace called
#'   `"Resource"`. Default is `TRUE`. Only applies to `MizerSim`.
#' @export
animate.MizerSim <- function(x, species = NULL, time_range = NULL,
                              log_x = TRUE, log_y = TRUE,
                              wlim = c(NA, NA), ylim = c(NA, NA),
                              power = 1, total = FALSE, resource = TRUE,
                              background = TRUE, ...) {
    sim <- x
    assert_that(is.flag(total), is.flag(resource), is.flag(background),
                is.number(power),
                length(wlim) == 2,
                length(ylim) == 2)

    species <- valid_species_arg(sim, species)
    if (is.null(time_range)) {
        time_range  <- as.numeric(dimnames(sim@n)$time)
    }
    time_elements <- get_time_elements(sim, time_range)

    nf <- melt(sim@n[time_elements,
                     as.character(dimnames(sim@n)$sp) %in% species,
                               , drop = FALSE])
    names(nf)[names(nf) == "sp"] <- "Species"
    nf$legend_name <- as.character(nf$Species)

    # Add resource ----
    if (resource) {
        nf_pp <- melt(sim@n_pp[time_elements, , drop = FALSE])
        nf_pp$Species <- "Resource"
        nf_pp$legend_name <- "Resource"
        nf <- rbind(nf, nf_pp)
    }
    # Add background ----
    # Keep each background species as its own trace (avoids oscillation from
    # interleaved data points) but label them all as "Background" in the legend.
    if (background && any(sim@params@species_params$is_background)) {
        back_n <- sim@n[time_elements, sim@params@species_params$is_background, , drop = FALSE]
        nf_back <- melt(back_n)
        names(nf_back)[names(nf_back) == "sp"] <- "Species"
        nf_back$legend_name <- "Background"
        nf <- rbind(nf, nf_back)
    }
    # Add total ----
    if (total) {
        # Calculate total community abundance
        fish_idx <- (length(sim@params@w_full) -
                         length(sim@params@w) + 1):length(sim@params@w_full)
        total_n <- sim@n_pp
        total_n[, fish_idx] <- total_n[, fish_idx] +
            rowSums(aperm(sim@n, c(1, 3, 2)), dims = 2)
        nf_total <- melt(total_n[time_elements, , drop = FALSE])
        nf_total$Species <- "Total"
        nf_total$legend_name <- "Total"
        nf <- rbind(nf, nf_total)
    }

    # Deal with power argument ----
    if (power %in% c(0, 1, 2)) {
        y_label <- c("Number density [1/g]", "Biomass density",
                     "Biomass density [g]")[power + 1]
    } else {
        y_label <- paste0("Number density * w^", power)
    }
    nf <- mutate(nf, value = value * w^power)

    # Impose limits ----
    if (is.na(wlim[1])) wlim[1] <- min(sim@params@w) / 100
    if (is.na(wlim[2])) wlim[2] <- max(sim@params@w_full)
    if (is.na(ylim[1])) ylim[1] <- 10^-20
    if (is.na(ylim[2])) ylim[2] <- 10^20
    nf <- nf %>%
        filter(value >= ylim[1],
               value <= ylim[2],
               w >= wlim[1],
               w <= wlim[2])

    animate_plotly(nf, sim@params, log_x, log_y, y_label)
}

# Build a plotly animation from a prepared long-format data frame.
# df must have columns: Species, legend_name, w, time, value.
# Traces are ordered by legend_name first (following params@linecolour), then
# by individual Species within each legend group — so background species always
# appear together and share a single legend entry.
animate_plotly <- function(df, params, log_x, log_y, y_label) {
    legend_name_order <- intersect(names(params@linecolour),
                                   unique(df$legend_name))
    sp_order <- unlist(lapply(legend_name_order, function(ln) {
        intersect(names(params@linecolour),
                  unique(df$Species[df$legend_name == ln]))
    }))
    p <- plotly::plot_ly()
    shown_legend_names <- character(0)
    for (sp in sp_order) {
        df_sp <- df[df$Species == sp, , drop = FALSE]
        ln <- unique(df_sp$legend_name)
        col <- params@linecolour[[ln]]
        showlegend <- !(ln %in% shown_legend_names)
        shown_legend_names <- c(shown_legend_names, ln)
        p <- plotly::add_lines(
            p,
            data = df_sp,
            x = ~w, y = ~value,
            frame = ~time,
            name = ln,
            legendgroup = ln,
            line = list(color = col, simplify = FALSE),
            showlegend = showlegend
        )
    }
    plotly::layout(p,
                   xaxis = list(type = if (log_x) "log" else "-",
                                exponentformat = "power",
                                title = "Size [g]"),
                   yaxis = list(type = if (log_y) "log" else "-",
                                exponentformat = "power",
                                title = y_label),
                   legend = list(traceorder = "normal"))
}

#' @rdname animate
#' @export
animateSpectra <- function(sim, ...) animate(sim, ...)
