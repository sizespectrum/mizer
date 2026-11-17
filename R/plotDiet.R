# Diet plots and their interactive compatibility wrapper.

#' Plot diet, resolved by prey species, as function of predator at size.
#'
#' `r lifecycle::badge("experimental")` Plots the proportions with which each
#' prey species contributes to the total biomass consumed by the specified
#' predator species, as a function of the predator's size. These proportions are
#' obtained with `getDiet()`.
#'
#' Prey species that contribute less than 1 permille to the diet are suppressed
#' in the plot. The plot only extends to predator sizes where the predator has a
#' meaningful abundance (defined as having a biomass density greater than 0.1%
#' of its maximum biomass density).
#'
#' If more than one predator species is selected, then the plot contains one
#' facet for each species.
#'
#' @param object An object of class \linkS4class{MizerSim} or
#'   \linkS4class{MizerParams}.
#' @param species The species to be selected. Optional. By default all target
#'   species are selected. A vector of species names, or a numeric vector with
#'   the species indices, or a logical vector indicating for each species
#'   whether it is to be selected (TRUE) or not.
#' @param wlim A numeric vector of length two providing lower and upper limits
#'   for the weight (x) axis. Use `NA` to auto-scale to the data range.
#' @param llim A numeric vector of length two providing lower and upper limits
#'   for the length (x) axis when `size_axis = "l"`. Use `NA` to auto-scale to
#'   the data range.
#' @param size_axis Whether to plot size as weight (`"w"`, default) or length
#'   (`"l"`), using the allometric weight-length relationship.
#' @param return_data A boolean value that determines whether the formatted data
#'   used for the plot is returned instead of the plot itself. Default is FALSE.
#' @param log_x If `TRUE` (default), use a log10 x-axis.
#' @param log_y If `TRUE`, use a log10 y-axis. Default is `FALSE`.
#' @param log Character string specifying which axes should use log10 scales,
#'   in the same form as the base [plot()] argument. For example, `"x"`,
#'   `"y"`, `"xy"` or `""`. If supplied, this overrides `log_x` and `log_y`.
#' @param ... Further arguments used by only some of the methods:
#'
#'   **For `MizerSim` methods:**
#'   * `time_range`: The time range (either a vector of values, a vector
#'     of min and max time, or a single value) over which to average the diet.
#'     The consumption rates are averaged over this range and then normalised
#'     to proportions. Default is the final time step.
#'
#' @return A ggplot2 object, unless `return_data = TRUE`, in which case a data
#'   frame with the four variables 'Predator', 'w' (or 'l' if
#'   `size_axis = "l"`), 'Proportion', 'Prey' is returned.
#' @export
#' @seealso [getDiet()]
#' @family plotting functions
#' @examples
#' \donttest{
#' plotDiet(NS_params, species = "Cod")
#' plotDiet(NS_params, species = 5:9)
#'
#' # Returning the data frame
#' fr <- plotDiet(NS_params, species = "Cod", return_data = TRUE)
#' str(fr)
#' }
#' @rdname plotDiet
#' @export
plotDiet <- function(object, species = NULL, wlim = c(NA, NA),
                     llim = c(NA, NA), size_axis = c("w", "l"),
                     return_data = FALSE, log_x = TRUE, log_y = FALSE,
                     log = NULL, ...) {
    UseMethod("plotDiet")
}

#' @rdname plotDiet
#' @usage NULL
#' @export
plotDiet.MizerSim <- function(object, species = NULL,
                              wlim = c(NA, NA), llim = c(NA, NA),
                              size_axis = c("w", "l"),
                              return_data = FALSE,
                              log_x = TRUE, log_y = FALSE, log = NULL,
                              time_range, ...) {
    log_axes <- parsePlotLog(log, log_x = log_x, log_y = log_y)
    size_axis <- plot_size_axis(size_axis)
    assert_that(is.flag(return_data),
                length(wlim) == 2,
                length(llim) == 2)
    if (missing(time_range)) {
        time_range <- max(as.numeric(dimnames(object@n)$time))
    }
    params <- validParams(object@params)
    # Work with the consumption rates (proportion = FALSE) so that the diet
    # averaged over a time range is the mean rate at which biomass of each prey
    # is consumed. Only afterwards do we normalise to proportions. Averaging the
    # per-time-step proportions directly would be wrong, since each step's
    # proportions are normalised independently.
    diet <- getDiet(object, time_range = time_range, drop = FALSE,
                    proportion = FALSE)
    if (length(dim(diet)) == 4) {
        diet <- apply(diet, c(2, 3, 4), mean)
    }
    # Normalise the time-averaged rates to proportions, as getDiet() does.
    total <- rowSums(diet, dims = 2)
    diet <- sweep(diet, c(1, 2), total, "/")
    diet[is.nan(diet)] <- 0
    # Use the abundance averaged over the same time range to restrict the plot
    # to each predator's meaningful size range.
    time_elements <- get_time_elements(object, time_range)
    n <- apply(object@n[time_elements, , , drop = FALSE], c(2, 3), mean)
    plot_diet(params, n = n, diet = diet, species = species,
              log_x = log_axes$log_x, log_y = log_axes$log_y,
              wlim = wlim, llim = llim, size_axis = size_axis,
              return_data = return_data)
}

#' @rdname plotDiet
#' @usage NULL
#' @export
plotDiet.MizerParams <- function(object, species = NULL,
                                 wlim = c(NA, NA), llim = c(NA, NA),
                                 size_axis = c("w", "l"),
                                 return_data = FALSE,
                                 log_x = TRUE, log_y = FALSE, log = NULL, ...) {
    log_axes <- parsePlotLog(log, log_x = log_x, log_y = log_y)
    size_axis <- plot_size_axis(size_axis)
    assert_that(is.flag(return_data),
                length(wlim) == 2,
                length(llim) == 2)
    params <- validParams(object)
    diet <- getDiet(params)
    plot_diet(params, n = params@initial_n, diet = diet, species = species,
              log_x = log_axes$log_x, log_y = log_axes$log_y,
              wlim = wlim, llim = llim, size_axis = size_axis,
              return_data = return_data)
}

#' Build the diet-composition plot
#'
#' Internal worker shared by the [plotDiet()] methods. It melts the diet array
#' into a data frame, restricts to meaningful size ranges, converts to a length
#' axis if requested and draws a stacked-area plot of prey proportions.
#'
#' @param params A MizerParams object.
#' @param n Array of species abundances (species by size).
#' @param diet Array of diet proportions (predator by size by prey).
#' @param species The predator species to be plotted.
#' @param log_x,log_y Logical flags for log10 axes.
#' @param wlim,llim Numeric vectors of length two giving the weight and length
#'   limits.
#' @param size_axis Either `"w"` (weight) or `"l"` (length).
#' @param return_data If `TRUE`, return the plotting data frame instead of the
#'   plot.
#' @return A `mizer_plot` (ggplot2) object, or the plotting data frame if
#'   `return_data = TRUE`.
#' @keywords internal
plot_diet <- function(params, n, diet, species, log_x, log_y, wlim, llim,
                      size_axis, return_data) {
    size_axis <- plot_size_axis(size_axis)
    species <- valid_species_arg(params, species, return.logical = TRUE)
    diet <- diet[species, , , drop = FALSE]
    names(dimnames(diet)) <- c("Predator", "w", "Prey")
    plot_dat <- melt(diet, value.name = "Proportion")
    prey <- dimnames(diet)$Prey
    # the plot looks better upsided down
    plot_dat$Prey <- factor(plot_dat$Prey, levels = rev(prey))

    plot_dat <- plot_dat[plot_dat$Proportion > 0.001, ]

    # Restrict plot to relevant size ranges where abundance is meaningful
    # For each predator species, find the maximum size where density is meaningful
    predator_names <- dimnames(diet)$Predator
    for (pred in predator_names) {
        pred_idx <- which(dimnames(n)[[1]] == pred)
        if (length(pred_idx) > 0) {
            density <- n[pred_idx, ] * params@w^2
            max_density <- max(density)
            if (max_density > 0) {
                # Find sizes where density > 1e-3 * max_density
                meaningful_idx <- which(density > 1e-3 * max_density)
                if (length(meaningful_idx) > 0) {
                    w_max_meaningful <- params@w[max(meaningful_idx)]
                    # Filter plot data for this predator
                    plot_dat <- plot_dat[!(plot_dat$Predator == pred & plot_dat$w > w_max_meaningful), ]
                }
            }
        }
    }

    plot_dat <- convert_plot_size_axis(plot_dat, params, size_axis,
                                       species_col = "Predator")
    if (return_data) return(plot_dat)
    x_var <- plot_size_x_var(size_axis)

    legend_levels <-
        intersect(names(params@linecolour), plot_dat$Prey)
    p <- ggplot(plot_dat) +
        geom_area(aes(x = .data[[x_var]], y = Proportion, fill = Prey)) +
        scale_x_continuous(trans = if (log_x) "log10" else "identity") +
        coord_cartesian(xlim = plot_size_xlim(wlim, size_axis, llim)) +
        scale_y_continuous(trans = if (log_y) "log10" else "identity") +
        labs(x = plot_size_xlab(size_axis), y = "Proportion") +
        scale_fill_manual(values = params@linecolour[legend_levels],
                          limits = legend_levels)
    if (sum(species) > 1) {
        p <- p + facet_wrap(vars(Predator))
    }
    make_mizer_plot(p, plot_size_tooltip(size_axis,
                                         after = c("Proportion", "Prey")))
}

#' @rdname plotDiet
#' @usage NULL
#' @return `plotlyDiet()` returns a plotly object.
#' @export
plotlyDiet <- function(object, species = NULL,
                       wlim = c(NA, NA), llim = c(NA, NA),
                       size_axis = c("w", "l"),
                       log_x = TRUE, log_y = FALSE, log = NULL, ...) {
    plotly_forward("plotDiet")
}
