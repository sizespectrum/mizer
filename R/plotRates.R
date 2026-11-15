# Feeding-level and mortality plots built from rate arrays.

# Prepare the rate's arithmetic mean over the requested saved times. Keep the
# dedicated plots' labels and average representation, including when a custom
# getter returns a plain matrix instead of a wrapped array.
prepare_plot_rate <- function(object, getter, time_range, value_name, units,
                              type = NULL, params = object@params) {
    if (missing(time_range)) {
        time_range <- max(as.numeric(dimnames(object@n)$time))
    }
    rate <- getter(object, time_range = time_range, drop = FALSE)
    if (length(dim(rate)) == 3) {
        rate <- apply(rate, c(2, 3), mean)
    }
    ArraySpeciesBySize(rate, value_name = value_name, units = units,
                       type = type, params = params,
                       representation = "average")
}

#' Plot the feeding level of species by size
#'
#' After running a projection, plot the feeding level of each species by size.
#' The feeding level is averaged over the specified time range (a single value
#' for the time range can be used).
#'
#' When called with a \linkS4class{MizerSim} object, the feeding level is averaged
#' over the specified time range (a single value for the time range can be used
#' to plot a single time step). When called with a \linkS4class{MizerParams}
#' object the initial feeding level is plotted.
#'
#' If `include_critical = TRUE` then the critical feeding level (the feeding
#' level at which the intake just covers the metabolic cost) is also plotted,
#' with a thinner line. This line should always stay below the line of the
#' actual feeding level, because the species would stop growing at any point
#' where the feeding level drops to the critical feeding level.
#'
#' @param object An object of class \linkS4class{MizerSim} or
#'   \linkS4class{MizerParams}.
#' @param species The species to be selected. Optional. By default all target
#'   species are selected. A vector of species names, or a numeric vector with
#'   the species indices, or a logical vector indicating for each species
#'   whether it is to be selected (TRUE) or not.
#' @param all.sizes If TRUE, then feeding level is plotted also for sizes
#'   outside a species' size range. Default FALSE.
#' @param highlight Name or vector of names of the species to be highlighted.
#' @param include_critical If TRUE, then the critical feeding level is also
#'   plotted. Default FALSE.
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
#'   \describe{
#'     \item{`time_range`}{The time range (either a vector of values, a vector
#'       of min and max time, or a single value) to average the feeding level
#'       over. Default is the final time step.}
#'   }
#'
#' @return A ggplot2 object, unless `return_data = TRUE`, in which case a data
#'   frame with the variables 'w' (or 'l' if `size_axis = "l"`), 'value' and
#'   'Species' is returned. If also `include_critical = TRUE` then the data
#'   frame contains a fourth variable 'Type' that distinguishes between
#'   'actual' and 'critical' feeding level.
#' @export
#' @family plotting functions
#' @seealso [plotting_functions], [getFeedingLevel()]
#' @examples
#' \donttest{
#' params <-  NS_params
#' sim <- project(params, effort=1, t_max=20, t_save = 2, progress_bar = FALSE)
#' plotFeedingLevel(sim)
#' plotFeedingLevel(sim, time_range = 10:20, species = c("Cod", "Herring"),
#'                  include_critical = TRUE)
#'
#' # Returning the data frame
#' fr <- plotFeedingLevel(sim, return_data = TRUE)
#' str(fr)
#' }
#' @rdname plotFeedingLevel
#' @export
plotFeedingLevel <- function(object, species = NULL, all.sizes = FALSE,
                             highlight = NULL, include_critical = FALSE,
                             wlim = c(NA, NA), llim = c(NA, NA),
                             size_axis = c("w", "l"), return_data = FALSE,
                             log_x = TRUE, log_y = FALSE, log = NULL, ...) {
    UseMethod("plotFeedingLevel")
}

#' @rdname plotFeedingLevel
#' @usage NULL
#' @export
plotFeedingLevel.MizerSim <- function(object, species = NULL,
            all.sizes = FALSE,
            highlight = NULL, include_critical = FALSE,
            wlim = c(NA, NA), llim = c(NA, NA),
            size_axis = c("w", "l"),
            return_data = FALSE,
            log_x = TRUE, log_y = FALSE, log = NULL,
            time_range, ...) {
    params <- validParams(object@params)
    feed <- prepare_plot_rate(object, getFeedingLevel, time_range,
                              value_name = "Feeding level", units = "",
                              type = "proportion", params = params)
    plot_feeding_level(params, feed, species = species,
                       highlight = highlight, all.sizes = all.sizes,
                       include_critical = include_critical,
                       wlim = wlim, llim = llim,
                       size_axis = size_axis,
                       return_data = return_data,
                       log_x = log_x, log_y = log_y, log = log, ...)
}

#' @rdname plotFeedingLevel
#' @usage NULL
#' @export
plotFeedingLevel.MizerParams <- function(object, species = NULL,
            all.sizes = FALSE,
            highlight = NULL, include_critical = FALSE,
            wlim = c(NA, NA), llim = c(NA, NA),
            size_axis = c("w", "l"),
            return_data = FALSE,
            log_x = TRUE, log_y = FALSE, log = NULL, ...) {
    params <- validParams(object)
    feed <- getFeedingLevel(params)
    plot_feeding_level(params, feed, species = species,
                       highlight = highlight, all.sizes = all.sizes,
                       include_critical = include_critical,
                       wlim = wlim, llim = llim,
                       size_axis = size_axis,
                       return_data = return_data,
                       log_x = log_x, log_y = log_y, log = log, ...)
}

#' Build the feeding-level plot
#'
#' Internal worker shared by the [plotFeedingLevel()] methods. It assembles the
#' plotting data, optionally adds the critical feeding level, restricts to each
#' species' size range, converts to a length axis if requested and draws the
#' plot.
#'
#' @param params A MizerParams object.
#' @param feed Array of feeding levels (species by size).
#' @param species The species to be plotted.
#' @param highlight Name or vector of names of species to be highlighted.
#' @param all.sizes If `FALSE`, feeding levels outside each species' size range
#'   are removed.
#' @param include_critical Whether to also plot the critical feeding level.
#' @param log_x,log_y Logical flags for log10 axes.
#' @param log Optional base-R log argument string or boolean.
#' @param wlim,llim Numeric vectors of length two giving the weight and length
#'   limits.
#' @param size_axis Either `"w"` (weight) or `"l"` (length).
#' @param return_data If `TRUE`, return the plotting data frame instead of the
#'   plot.
#' @param ... Additional arguments passed to [plot.ArraySpeciesBySize()].
#' @return A `mizer_plot` (ggplot2) object, or the plotting data frame if
#'   `return_data = TRUE`.
#' @keywords internal
plot_feeding_level <- function(params, feed, species, highlight,
                               all.sizes, include_critical,
                               wlim, llim, size_axis, return_data,
                               log_x = TRUE, log_y = FALSE, log = NULL, ...) {
    size_axis <- plot_size_axis(size_axis)
    log_axes <- parsePlotLog(log, log_x = log_x, log_y = log_y)
    log_x <- log_axes$log_x
    log_y <- log_axes$log_y

    if (!include_critical) {
        return(plot(feed, species = species, all.sizes = all.sizes,
                    highlight = highlight, wlim = wlim, llim = llim,
                    size_axis = size_axis, return_data = return_data,
                    log_x = log_x, log_y = log_y, ...))
    }
    crit <- getCriticalFeedingLevel(params)
    df_actual <- plot(feed, species = species, all.sizes = all.sizes,
                      wlim = wlim, llim = llim, size_axis = size_axis,
                      log_x = log_x, log_y = log_y,
                      return_data = TRUE, ...)
    df_crit <- plot(crit, species = species, all.sizes = all.sizes,
                    wlim = wlim, llim = llim, size_axis = size_axis,
                    log_x = log_x, log_y = log_y,
                    return_data = TRUE, ...)
    names(df_crit)[2] <- names(df_actual)[2]
    df_actual$Type <- "actual"
    df_crit$Type <- "critical"
    if (return_data) {
        return(rbind(df_actual, df_crit))
    }
    p <- plot(feed, species = species, all.sizes = all.sizes,
              highlight = highlight, wlim = wlim, llim = llim,
              size_axis = size_axis, log_x = log_x, log_y = log_y, ...)
    p <- addPlot(p, crit, species = species, all.sizes = all.sizes,
                 wlim = wlim, llim = llim, size_axis = size_axis,
                 alpha = 0.5)
    if (!log_y) {
        ylim <- list(...)$ylim %||% c(NA, NA)
        feeding_level_ylim <- array_ylim(feed, ylim, log_y,
                                         c(df_actual[[2]], df_crit[[2]]))
        p$coordinates <- coord_cartesian(xlim = p$coordinates$limits$x,
                                         ylim = feeding_level_ylim)
    }
    p
}

#' @rdname plotFeedingLevel
#' @usage NULL
#' @export
plotlyFeedingLevel <- function(object, ...) {
    plotHover(plotFeedingLevel(object, ...))
}

#' Plot predation mortality rate of each species against size
#'
#' After running a projection, plot the predation mortality rate of each species
#' by size. The mortality rate is averaged over the specified time range (a
#' single value for the time range can be used to plot a single time step).
#'
#' @param object An object of class \linkS4class{MizerSim} or
#'   \linkS4class{MizerParams}.
#' @param species The species to be selected. Optional. By default all target
#'   species are selected. A vector of species names, or a numeric vector with
#'   the species indices, or a logical vector indicating for each species
#'   whether it is to be selected (TRUE) or not.
#' @param all.sizes If TRUE, then predation mortality is plotted also for sizes
#'   outside a species' size range. Default FALSE.
#' @param highlight Name or vector of names of the species to be highlighted.
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
#'   \describe{
#'     \item{`time_range`}{The time range (either a vector of values, a vector
#'       of min and max time, or a single value) to average the predation
#'       mortality over. Default is the final time step.}
#'   }
#'
#' @return  A ggplot2 object, unless `return_data = TRUE`, in which case a data
#'   frame with the three variables 'w' (or 'l' if `size_axis = "l"`), 'value',
#'   'Species' is returned.
#' @export
#' @family plotting functions
#' @seealso [plotting_functions],  [getPredMort()]
#' @examples
#' \donttest{
#' params <-  NS_params
#' sim <- project(params, effort=1, t_max=20, t_save = 2, progress_bar = FALSE)
#' plotPredMort(sim)
#' plotPredMort(sim, time_range = 10:20)
#'
#' # Returning the data frame
#' fr <- plotPredMort(sim, return_data = TRUE)
#' str(fr)
#' }
#' @rdname plotPredMort
#' @export
plotPredMort <- function(object, species = NULL, all.sizes = FALSE,
                         highlight = NULL, wlim = c(NA, NA), llim = c(NA, NA),
                         size_axis = c("w", "l"), return_data = FALSE,
                         log_x = TRUE, log_y = FALSE, log = NULL, ...) {
    UseMethod("plotPredMort")
}

#' @rdname plotPredMort
#' @usage NULL
#' @export
plotPredMort.MizerSim <- function(object, species = NULL,
                         all.sizes = FALSE,
                         highlight = NULL,
                         wlim = c(NA, NA), llim = c(NA, NA),
                         size_axis = c("w", "l"),
                         return_data = FALSE,
                         log_x = TRUE, log_y = FALSE, log = NULL,
                         time_range, ...) {
    size_axis <- plot_size_axis(size_axis)
    pred_mort <- prepare_plot_rate(object, getPredMort, time_range,
                                    value_name = "Predation mortality",
                                    units = "1/year")
    plot(pred_mort, species = species, all.sizes = all.sizes,
         highlight = highlight, wlim = wlim, llim = llim,
         log_x = log_x, log_y = log_y, log = log,
         size_axis = size_axis,
         return_data = return_data,
         ylim = c(0, NA))
}

#' @rdname plotPredMort
#' @usage NULL
#' @export
plotPredMort.MizerParams <- function(object, species = NULL,
                         all.sizes = FALSE,
                         highlight = NULL,
                         wlim = c(NA, NA), llim = c(NA, NA),
                         size_axis = c("w", "l"),
                         return_data = FALSE,
                         log_x = TRUE, log_y = FALSE, log = NULL,
                         ...) {
    size_axis <- plot_size_axis(size_axis)
    plot(getPredMort(validParams(object)), species = species,
         all.sizes = all.sizes, highlight = highlight,
         log_x = log_x, log_y = log_y, log = log,
         wlim = wlim, llim = llim, size_axis = size_axis,
         return_data = return_data, ylim = c(0, NA))
}

#' Alias for `plotPredMort()`
#'
#' `r lifecycle::badge("superseded")`
#' An alias provided for backward compatibility with mizer version <= 1.0
#' @inherit plotPredMort
#' @export
#' @concept deprecated
plotM2 <- plotPredMort

#' @rdname plotPredMort
#' @usage NULL
#' @export
plotlyPredMort <- function(object, species = NULL,
                           time_range,
                           highlight = NULL,
                           wlim = c(NA, NA), llim = c(NA, NA),
                           size_axis = c("w", "l"),
                           log_x = TRUE, log_y = FALSE, log = NULL, ...) {
    size_axis <- plot_size_axis(size_axis)
    argg <- as.list(environment())
    plotHover(do.call("plotPredMort", argg))
}

#' Plot total fishing mortality of each species by size
#'
#' After running a projection, plot the total fishing mortality of each species
#' by size. The total fishing mortality is averaged over the specified time
#' range (a single value for the time range can be used to plot a single time
#' step).
#'
#' @param object An object of class \linkS4class{MizerSim} or
#'   \linkS4class{MizerParams}.
#' @param species The species to be selected. Optional. By default all target
#'   species are selected. A vector of species names, or a numeric vector with
#'   the species indices, or a logical vector indicating for each species
#'   whether it is to be selected (TRUE) or not.
#' @param all.sizes If TRUE, then fishing mortality is plotted also for sizes
#'   outside a species' size range. Default FALSE.
#' @param highlight Name or vector of names of the species to be highlighted.
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
#'   \describe{
#'     \item{`time_range`}{The time range (either a vector of values, a vector
#'       of min and max time, or a single value) to average the fishing
#'       mortality over. Default is the final time step.}
#'   }
#' @return A ggplot2 object, unless `return_data = TRUE`, in which case a data
#'   frame with the three variables 'w' (or 'l' if `size_axis = "l"`), 'value',
#'   'Species' is returned.
#' @export
#' @family plotting functions
#' @seealso [plotting_functions], [getFMort()]
#' @examples
#' \donttest{
#' params <-  NS_params
#' sim <- project(params, effort=1, t_max=20, t_save = 2, progress_bar = FALSE)
#' plotFMort(sim)
#' plotFMort(sim, highlight = c("Cod", "Haddock"))
#'
#' # Returning the data frame
#' fr <- plotFMort(sim, return_data = TRUE)
#' str(fr)
#' }
#' @rdname plotFMort
#' @export
plotFMort <- function(object, species = NULL, all.sizes = FALSE,
                      highlight = NULL, wlim = c(NA, NA), llim = c(NA, NA),
                      size_axis = c("w", "l"), return_data = FALSE,
                      log_x = TRUE, log_y = FALSE, log = NULL, ...) {
    UseMethod("plotFMort")
}

#' @rdname plotFMort
#' @usage NULL
#' @export
plotFMort.MizerSim <- function(object, species = NULL,
                      all.sizes = FALSE,
                      highlight = NULL,
                      wlim = c(NA, NA), llim = c(NA, NA),
                      size_axis = c("w", "l"),
                      return_data = FALSE,
                      log_x = TRUE, log_y = FALSE, log = NULL,
                      time_range, ...) {
    size_axis <- plot_size_axis(size_axis)
    f <- prepare_plot_rate(object, getFMort, time_range,
                           value_name = "Fishing mortality", units = "1/year")
    plot(f, species = species, all.sizes = all.sizes,
         highlight = highlight, wlim = wlim, llim = llim,
         log_x = log_x, log_y = log_y, log = log,
         size_axis = size_axis,
         return_data = return_data)
}

#' @rdname plotFMort
#' @usage NULL
#' @export
plotFMort.MizerParams <- function(object, species = NULL,
                      all.sizes = FALSE,
                      highlight = NULL,
                      wlim = c(NA, NA), llim = c(NA, NA),
                      size_axis = c("w", "l"),
                      return_data = FALSE,
                      log_x = TRUE, log_y = FALSE, log = NULL,
                      ...) {
    size_axis <- plot_size_axis(size_axis)
    plot(getFMort(validParams(object)), species = species,
         all.sizes = all.sizes, highlight = highlight,
         log_x = log_x, log_y = log_y, log = log,
         wlim = wlim, llim = llim, size_axis = size_axis,
         return_data = return_data)
}

#' @rdname plotFMort
#' @usage NULL
#' @export
plotlyFMort <- function(object, species = NULL,
                        time_range,
                        highlight = NULL,
                        wlim = c(NA, NA), llim = c(NA, NA),
                        size_axis = c("w", "l"),
                        log_x = TRUE, log_y = FALSE, log = NULL, ...) {
    size_axis <- plot_size_axis(size_axis)
    argg <- as.list(environment())
    plotHover(do.call("plotFMort", argg))
}
