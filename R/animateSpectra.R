#' Animate size-dependent quantities through time
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
#'   quantity returned by a `MizerSim` accessor, such as [getFMort()],
#'   [getFeedingLevel()], or [getPredMort()]. Both axes are log10 by default
#'   and can each be switched to linear with `log_x = FALSE` or `log_y = FALSE`.
#'   Background species and a species total can be added via the `background`
#'   and `total` arguments.
#'
#' Species linecolours and linetypes follow `params@linecolour` and
#' `params@linetype`.
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
#'   maximum. Limits are applied as Plotly axis ranges, so points outside the
#'   limits are clipped by the viewport rather than removed from the animation
#'   frames.
#' @param frame_duration Duration in milliseconds for which each saved frame is
#'   displayed. Default is 500.
#' @param transition_duration Duration in milliseconds of the interpolation
#'   between frames. Use `transition_duration = 0` to step directly from one
#'   saved frame to the next. Default is `frame_duration`.
#' @param easing The Plotly easing function to use when interpolating between
#'   frames. Default is `"linear"`. Available options are `"linear"`, `"quad"`,
#'   `"cubic"`, `"sin"`, `"exp"`, `"circle"`, `"elastic"`, `"back"`,
#'   `"bounce"`, and each of those with suffix `"-in"`, `"-out"`, or
#'   `"-in-out"` appended, for example `"cubic-in-out"`.
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
                              background = TRUE,
                              frame_duration = 500,
                              transition_duration = frame_duration,
                              easing = "linear", ...) {
    sim <- x
    assert_that(is.flag(total), is.flag(resource), is.flag(background),
                is.number(power),
                is.number(frame_duration), frame_duration >= 0,
                is.number(transition_duration), transition_duration >= 0,
                is.string(easing),
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

    animate_plotly(nf, sim@params, log_x, log_y, y_label, wlim, ylim,
                   frame_duration, transition_duration, easing)
}

# Build a plotly animation from a prepared long-format data frame.
# df must have columns: Species, legend_name, w, time, value.
# Traces are ordered by legend_name first (following params@linecolour), then
# by individual Species within each legend group — so background species always
# appear together and share a single legend entry.
animate_plotly <- function(df, params, log_x, log_y, y_label,
                           wlim = c(NA, NA), ylim = c(NA, NA),
                           frame_duration = 500, transition_duration = 500,
                           easing = "linear") {
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
    p <- plotly::layout(p,
                        xaxis = plotly_axis(df$w, wlim, log_x, "Size [g]"),
                        yaxis = plotly_axis(df$value, ylim, log_y, y_label),
                        legend = list(traceorder = "normal"))
    plotly::animation_opts(p, frame = frame_duration,
                           transition = transition_duration,
                           easing = easing)
}

plotly_axis <- function(values, limits, log_axis, title) {
    axis <- list(type = if (log_axis) "log" else "-",
                 exponentformat = "power",
                 title = title)
    range <- plotly_axis_range(values, limits, log_axis)
    if (!is.null(range)) axis$range <- range
    axis
}

plotly_axis_range <- function(values, limits, log_axis) {
    if (all(is.na(limits))) return(NULL)

    finite_values <- values[is.finite(values)]
    if (log_axis) finite_values <- finite_values[finite_values > 0]

    if (is.na(limits[1]) && length(finite_values) > 0) {
        limits[1] <- min(finite_values)
    }
    if (is.na(limits[2]) && length(finite_values) > 0) {
        limits[2] <- max(finite_values)
    }
    if (any(is.na(limits))) return(NULL)

    if (log_axis) {
        if (limits[1] <= 0 && length(finite_values) > 0) {
            limits[1] <- min(finite_values)
        }
        if (limits[2] <= 0 && length(finite_values) > 0) {
            limits[2] <- max(finite_values)
        }
        if (any(limits <= 0)) return(NULL)
        limits <- log10(limits)
    }

    if (!all(is.finite(limits))) return(NULL)
    if (limits[1] == limits[2]) {
        delta <- if (log_axis) 0.5 else max(abs(limits[1]) * 0.05, 1)
        limits <- limits + c(-1, 1) * delta
    }
    limits
}

#' @rdname animate
#' @export
animateSpectra <- function(sim, ...) animate(sim, ...)
