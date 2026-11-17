# Yield time series, gear-resolved yield and interactive wrappers.

#' Plot the total yield of species through time
#'
#' After running a projection, the total yield of each species across all
#' fishing gears can be plotted against time. The yield is obtained with
#' [getYield()].
#'
#' @param object An object of class \linkS4class{MizerSim}
#' @param sim2 An optional second object of class \linkS4class{MizerSim}. If
#'   this is provided its yields will be shown on the same plot in bolder lines.
#' @param species The species to be selected. Optional. By default all target
#'   species are selected. A vector of species names, or a numeric vector with
#'   the species indices, or a logical vector indicating for each species
#'   whether it is to be selected (TRUE) or not.
#' @param total A boolean value that determines whether the total yield from all
#'   species is plotted as well. Default is FALSE.
#' @param log_x If `TRUE`, use a log10 x-axis. Default is `FALSE`.
#' @param log_y If `TRUE`, use a log10 y-axis. Default is `TRUE`.
#' @param log Character string specifying which axes should use log10 scales,
#'   in the same form as the base [plot()] argument. For example, `"x"`,
#'   `"y"`, `"xy"` or `""`. If supplied, this overrides `log_x` and `log_y`.
#'   For backward compatibility, `TRUE` and `FALSE` are interpreted as setting
#'   only `log_y`.
#' @param ylim A numeric vector of length two providing lower and upper limits
#'   for the y axis. Use `NA` to refer to the existing minimum or maximum.
#' @param tlim A numeric vector of length two providing lower and upper limits
#'   for the time axis, e.g. `c(1980, 2000)`. Use `NA` to apply no limit at
#'   that end. Default is `c(NA, NA)`.
#' @param highlight Name or vector of names of the species to be highlighted.
#' @param return_data A boolean value that determines whether the formatted data
#'   used for the plot is returned instead of the plot itself. Default is FALSE.
#' @param ... Arguments passed to [getYield()].
#'
#' @return A ggplot2 object, unless `return_data = TRUE`, in which case a data
#'   frame with the three variables 'Year', 'Yield', 'Species' is returned. If
#'   `sim2` is supplied, the data frame also has a 'Simulation' variable.
#' @export
#' @family plotting functions
#' @seealso [plotting_functions],  [getYield()]
#' @examples
#' \donttest{
#' params <- NS_params
#' sim <- project(params, effort = 1, t_max = 20, t_save = 0.2, progress_bar = FALSE)
#' plotYield(sim)
#' plotYield(sim, species = c("Cod", "Herring"), total = TRUE)
#'
#' # Comparing with yield from twice the effort
#' sim2 <- project(params, effort=2, t_max=20, t_save = 0.2, progress_bar = FALSE)
#' plotYield(sim, sim2, species = c("Cod", "Herring"), log = FALSE)
#'
#' # Returning the data frame
#' fr <- plotYield(sim, return_data = TRUE)
#' str(fr)
#' }
#' @rdname plotYield
#' @export
plotYield <- function(object, sim2 = NULL, species = NULL, total = FALSE,
                      log_x = FALSE, log_y = TRUE, log = NULL,
                      ylim = c(NA, NA), tlim = c(NA, NA), highlight = NULL,
                      return_data = FALSE, ...) {
    UseMethod("plotYield")
}

#' @rdname plotYield
#' @usage NULL
#' @export
plotYield.MizerSim <- function(object, sim2 = NULL,
                      species = NULL,
                      total = FALSE,
                      log_x = FALSE, log_y = TRUE, log = NULL,
                      ylim = c(NA, NA), tlim = c(NA, NA),
                      highlight = NULL, return_data = FALSE,
                      ...) {
    log_axes <- parsePlotLog(log, log_x = log_x, log_y = log_y)
    assert_that(is(object, "MizerSim"),
                is.flag(total),
                is.flag(log_axes$log_x),
                is.flag(log_axes$log_y),
                is.flag(return_data))
    if (is.null(sim2)) {
        y <- getYield(object, ...)
        attr(y, "value_name") <- "Yield"
        ans <- plot(y, species = species,
                    tlim = tlim, ylim = ylim,
                    total = total, highlight = highlight,
                    log_x = log_axes$log_x, log_y = log_axes$log_y,
                    return_data = return_data)
        if (return_data) ans$Legend <- NULL
        return(ans)
    }
    lifecycle::deprecate_warn("2.6.0", "plotYield(sim2)", "plot2()",
                              details = "Use plot2(getYield(sim1), getYield(sim2)) instead.")
    params <- object@params
    # We need to combine two plots
    if (!all(dimnames(object@n)$time == dimnames(sim2@n)$time)) {
        stop("The two simulations do not have the same times")
    }
    ym <- plotYield(object, species = species,
                    tlim = tlim, total = total,
                    log_x = log_axes$log_x, log_y = log_axes$log_y,
                    ylim = ylim,
                    highlight = highlight, return_data = TRUE, ...)
    ym2 <- plotYield(sim2, species = species,
                     tlim = tlim, total = total,
                     log_x = log_axes$log_x, log_y = log_axes$log_y,
                     ylim = ylim,
                     highlight = highlight, return_data = TRUE, ...)
    ym$Simulation <- rep(1, nrow(ym))
    ym2$Simulation <- rep(2, nrow(ym2))
    ym <- rbind(ym, ym2)

    if (return_data) return(ym)

    plotDataFrame(ym, params,
                  ylab = "Yield [g/year]",
                  xtrans = ifelse(log_axes$log_x, "log10", "identity"),
                  ytrans = ifelse(log_axes$log_y, "log10", "identity"),
                  xlim = tlim, ylim = ylim,
                  highlight = highlight, wrap_var = "Simulation")
}

#' @rdname plotYield
#' @usage NULL
#' @export
plotlyYield <- function(object, sim2 = NULL,
                        species = NULL,
                        total = FALSE,
                        log_x = FALSE, log_y = TRUE, log = NULL,
                        ylim = c(NA, NA), tlim = c(NA, NA),
                        highlight = NULL,
                        ...) {
    plotly_forward("plotYield", tooltip = c("Species", "Year", "Yield"))
}

#' Plot the total yield of each species by gear through time
#'
#' After running a projection, the total yield of each species by fishing gear
#' can be plotted against time.
#'
#' This plot is pretty easy to do by hand. It just
#' gets the biomass using the [getYieldGear()] method and plots using
#' the ggplot2 package. You can then fiddle about with colours and linetypes
#' etc. Just look at the source code for details.
#'
#' @param object An object of class \linkS4class{MizerSim}
#' @param species The species to be selected. Optional. By default all target
#'   species are selected. A vector of species names, or a numeric vector with
#'   the species indices, or a logical vector indicating for each species
#'   whether it is to be selected (TRUE) or not.
#' @param gears A vector of gear names to be included in the plot. Default is
#'  all gears.
#' @param total A boolean value that determines whether the total yield from all
#'   species is plotted as well. Default is FALSE.
#' @param log_x If `TRUE`, use a log10 x-axis. Default is `FALSE`.
#' @param log_y If `TRUE`, use a log10 y-axis. Default is `TRUE`.
#' @param log Character string specifying which axes should use log10 scales,
#'   in the same form as the base [plot()] argument. For example, `"x"`,
#'   `"y"`, `"xy"` or `""`. If supplied, this overrides `log_x` and `log_y`.
#'   For backward compatibility, `TRUE` and `FALSE` are interpreted as setting
#'   only `log_y`.
#' @param ylim A numeric vector of length two providing lower and upper limits
#'   for the y axis. Use `NA` to refer to the existing minimum or maximum.
#' @param tlim A numeric vector of length two providing lower and upper limits
#'   for the time axis, e.g. `c(1980, 2000)`. Use `NA` to apply no limit at
#'   that end. Default is `c(NA, NA)`.
#' @param highlight Name or vector of names of the species to be highlighted.
#' @param return_data A boolean value that determines whether the formatted data
#'   used for the plot is returned instead of the plot itself. Default is FALSE.
#' @param ... Arguments passed to [getYieldGear()].
#'
#' @return A ggplot2 object, unless `return_data = TRUE`, in which case a data
#'   frame with the four variables 'Year', 'Yield', 'Species' and 'Gear' is
#'   returned.
#' @export
#' @family plotting functions
#' @seealso [plotting_functions],  [getYieldGear()]
#' @examples
#' \donttest{
#' params <-  NS_params
#' sim <- project(params, effort=1, t_max=20, t_save = 0.2, progress_bar = FALSE)
#' plotYieldGear(sim)
#' plotYieldGear(sim, species = c("Cod", "Herring"), total = TRUE)
#'
#' # Returning the data frame
#' fr <- plotYieldGear(sim, return_data = TRUE)
#' str(fr)
#' }
#' @rdname plotYieldGear
#' @export
plotYieldGear <- function(object, species = NULL, gears = NULL,
                          total = FALSE,
                          log_x = FALSE, log_y = TRUE, log = NULL,
                          ylim = c(NA, NA), tlim = c(NA, NA),
                          highlight = NULL, return_data = FALSE, ...) {
    UseMethod("plotYieldGear")
}

#' @rdname plotYieldGear
#' @usage NULL
#' @export
plotYieldGear.MizerSim <- function(object,
                          species = NULL,
                          gears = NULL,
                          total = FALSE,
                          log_x = FALSE, log_y = TRUE, log = NULL,
                          ylim = c(NA, NA), tlim = c(NA, NA),
                          highlight = NULL, return_data = FALSE,
                          ...) {
    log_axes <- parsePlotLog(log, log_x = log_x, log_y = log_y)
    assert_that(is(object, "MizerSim"),
                is.flag(total),
                is.flag(log_axes$log_x),
                is.flag(log_axes$log_y),
                is.flag(return_data))
    params <- object@params
    species <- valid_species_arg(object, species, error_on_empty = TRUE)
    gears <- valid_gears_arg(object, gears, error_on_empty = TRUE)

    y <- getYieldGear(object, ...)
    y_total <- rowSums(y, dims = 2)
    y <- y[, dimnames(y)$gear %in% gears, dimnames(y)$sp %in% species, drop = FALSE]
    names(dimnames(y))[names(dimnames(y)) == "sp"] <- "Species"
    ym <- reshape2::melt(y)
    if (total) {
        yt <- reshape2::melt(y_total)
        yt$Species <- "Total"
        ym <- rbind(ym, yt)
    }
    ym <- subset(ym, ym$value > 0)

    ym <- ym[, c(1, 4, 3, 2)]
    names(ym) <- c("Year", "Yield", "Species", "Gear")

    if (return_data) return(ym)

    # This does not use `plotDataFrame()` because it uses Gear to set
    # the linetype
    # Need to keep species in order for legend
    species_levels <- intersect(names(params@linecolour), ym$Species)
    ym$Species <- factor(ym$Species, levels = species_levels)
    ym$Legend <- ym$Species
    linesize <- make_linesize(species_levels, highlight)
    make_mizer_plot(
        ggplot(ym) +
            geom_line(aes(x = Year, y = Yield, colour = Species,
                          linetype = Gear, linewidth = Legend)) +
            scale_y_continuous(trans = if (log_axes$log_y) "log10" else "identity",
                               name = "Yield [g]") +
            scale_x_continuous(trans = if (log_axes$log_x) "log10" else "identity") +
            coord_cartesian(xlim = tlim, ylim = ylim) +
            scale_colour_manual(values = params@linecolour[species_levels]) +
            scale_discrete_manual("linewidth", values = linesize),
        c("Species", "Gear", "Year", "Yield")
    )
}

#' @rdname plotYieldGear
#' @usage NULL
#' @export
plotlyYieldGear <- function(object, species = NULL,
                            gears = NULL,
                            total = FALSE,
                            log_x = FALSE, log_y = TRUE, log = NULL,
                            ylim = c(NA, NA), tlim = c(NA, NA),
                            highlight = NULL, ...) {
    plotly_forward("plotYieldGear",
                   tooltip = c("Species", "Gear", "Year", "Yield"))
}
