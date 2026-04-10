#' Animation of the abundance spectra
#'
#' `r lifecycle::badge("experimental")`
#'
#' @param sim A MizerSim object
#' @param species Name or vector of names of the species to be plotted. By
#'   default all species are plotted.
#' @param time_range The time range to animate over. Either a vector of values
#'   or a vector of min and max time. Default is the entire time range of the
#'   simulation.
#' @param wlim A numeric vector of length two providing lower and upper limits
#'   for the w axis. Use NA to refer to the existing minimum or maximum.
#' @param ylim A numeric vector of length two providing lower and upper limits
#'   for the y axis. Use NA to refer to the existing minimum or maximum. Any
#'   values below 1e-20 are always cut off.
#' @param power The abundance is plotted as the number density times the weight
#' raised to \code{power}. The default \code{power = 1} gives the biomass
#' density, whereas \code{power = 2} gives the biomass density with respect
#' to logarithmic size bins.
#' @param total A boolean value that determines whether the total over all
#'   species in the system is plotted as an additional trace called `"Total"`.
#'   Default is FALSE.
#' @param resource A boolean value that determines whether resource is included.
#'   If `TRUE`, the resource spectrum is plotted as an additional trace called
#'   `"Resource"`. Default is TRUE.
#' @param background A boolean value that determines whether background species
#'   are included. Ignored if the model does not contain background species.
#'   Default is TRUE.
#' @param ... Additional arguments passed to the method.
#'
#' @return A plotly object with one animated line trace per plotted group. The
#'   y-axis title depends on `power`.
#' @export
#' @family plotting functions
#' @examples
#' \donttest{
#' animateSpectra(NS_sim, power = 2, wlim = c(0.1, NA), time_range = 1997:2007)
#' }
animateSpectra <- function(sim, species, time_range,
                           wlim,
                           ylim,
                           power,
                           total,
                           resource,
                           background, ...)
    UseMethod("animateSpectra")

#' @export
animateSpectra.MizerSim <- function(sim, species = NULL, time_range = NULL,
                                    wlim = c(NA, NA), ylim = c(NA, NA),
                                    power = 1, total = FALSE, resource = TRUE,
                                    background = TRUE, ...) {
    assert_that(is.flag(total), is.flag(resource), is.flag(background),
                is.number(power),
                length(wlim) == 2,
                length(ylim) == 2)

    species <- valid_species_arg(sim, species)
    if (missing(time_range)) {
        time_range  <- as.numeric(dimnames(sim@n)$time)
    }
    time_elements <- get_time_elements(sim, time_range)

    nf <- melt(sim@n[time_elements,
                     as.character(dimnames(sim@n)$sp) %in% species,
                               , drop = FALSE])
    # legend_name drives the legend entry and colour lookup; sp drives trace splitting.
    # For regular species these are identical.
    nf$legend_name <- as.character(nf$sp)

    # Add resource ----
    if (resource) {
        nf_pp <- melt(sim@n_pp[time_elements, , drop = FALSE])
        nf_pp$sp <- "Resource"
        nf_pp$legend_name <- "Resource"
        nf <- rbind(nf, nf_pp)
    }
    # Add background ----
    # Keep each background species as its own trace (avoids oscillation from
    # interleaved data points) but label them all as "Background" in the legend.
    if (background && any(sim@params@species_params$is_background)) {
        back_n <- sim@n[time_elements, sim@params@species_params$is_background, , drop = FALSE]
        nf_back <- melt(back_n)
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
        nf_total$sp <- "Total"
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

    # Determine trace order: follow linecolour order for legend_names, then
    # within each legend_name group follow linecolour order for individual sp.
    legend_name_order <- intersect(names(sim@params@linecolour),
                                   unique(nf$legend_name))
    sp_order <- unlist(lapply(legend_name_order, function(ln) {
        intersect(names(sim@params@linecolour),
                  unique(nf$sp[nf$legend_name == ln]))
    }))

    # Build one trace per sp; background species share a legend group so only
    # the first one gets a visible legend entry.
    p <- plotly::plot_ly()
    shown_legend_names <- character(0)
    for (sp_lev in sp_order) {
        df_sp <- nf[nf$sp == sp_lev, , drop = FALSE]
        ln <- unique(df_sp$legend_name)
        col <- sim@params@linecolour[[ln]]
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
                   xaxis = list(type = "log", exponentformat = "power",
                                title = "Size [g]"),
                   yaxis = list(type = "log", exponentformat = "power",
                                title = y_label),
                   legend = list(traceorder = "normal"))
}
