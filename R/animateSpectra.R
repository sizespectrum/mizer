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
#'   species in the system is plotted as well. Default is FALSE.
#' @param resource A boolean value that determines whether resource is included.
#'   Default is TRUE.
#' 
#' @return A plotly object
#' @export
#' @family plotting functions
#' @examples
#' \donttest{
#' animateSpectra(NS_sim, power = 2, wlim = c(0.1, NA))
#' }
animateSpectra <- function(sim,
                           species = NULL,
                           time_range,
                           wlim = c(NA, NA),
                           ylim = c(NA, NA),
                           power = 1,
                           total = FALSE,
                           resource = TRUE) {

    species <- valid_species_arg(sim, species)
    if (missing(time_range)) {
        time_range  <- as.numeric(dimnames(sim@n)$time)
    }
    time_elements <- get_time_elements(sim, time_range)
    nf <- melt(sim@n[time_elements,
                     as.character(dimnames(sim@n)$sp) %in% species,
                               , drop = FALSE])

    # Add resource ----
    if (resource) {
        nf_pp <- melt(sim@n_pp[time_elements, , drop = FALSE])
        nf_pp$sp <- "Resource"
        nf <- rbind(nf, nf_pp)
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

    nf %>%
        plotly::plot_ly() %>%
        plotly::add_lines(x = ~w, y = ~value,
                          color = ~sp, colors = sim@params@linecolour,
                          frame = ~time,
                          line = list(simplify = FALSE)) %>%
        plotly::layout(xaxis = list(type = "log", exponentformat = "power",
                                    title = "Size [g]"),
                       yaxis = list(type = "log", exponentformat = "power",
                                    title = y_label))
}
