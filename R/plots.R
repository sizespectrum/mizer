# Plotting functions ----

# Copyright 2012 Finlay Scott and Julia Blanchard.
# Copyright 2019 Gustav Delius and Richard Southwell.
# Development has received funding from the European Commission's Horizon 2020
# Research and Innovation Programme under Grant Agreement No. 634495
# for the project MINOUW (http://minouw-project.eu/).
# Distributed under the GPL 3 or later
# Maintainer: Gustav Delius, University of York, <gustav.delius@york.ac.uk>

#' Description of the plotting functions
#'
#' Mizer provides a range of plotting functions for visualising the results
#' of running a simulation, stored in a MizerSim object, or the initial state
#' stored in a MizerParams object.
#'
#' The quickest way to make a standard plot is often to call `plot()` directly.
#' mizer provides `plot()` methods for `MizerSim` and `MizerParams` objects, and
#' also for the array classes returned by many summary and rate functions:
#'
#' * `plot(<MizerSim>)` produces a five-panel summary plot with
#'   [plotFeedingLevel()], [plotBiomass()], [plotPredMort()], [plotFMort()] and
#'   [plotSpectra()].
#'
#' * `plot(<MizerParams>)` produces a three-panel summary plot with
#'   [plotFeedingLevel()], [plotPredMort()] and [plotSpectra()] from the initial
#'   state.
#'
#' * `plot(<ArrayTimeBySpecies>)` plots any time-by-species array, such as
#'   those returned by [getBiomass()], [getSSB()], [getYield()] and [getN()] on
#'   a `MizerSim`, as lines of value against time.
#'
#' * `plot(<ArraySpeciesBySize>)` plots any species-by-size array, such as
#'   those returned by [getEncounter()], [getFeedingLevel()], [getPredMort()],
#'   [getFMort()] and [getSearchVolume()], as lines of value against body size.
#'
#' * `plot(<ArrayTimeBySpeciesBySize>)` plots a time slice from a
#'   time-by-species-by-size array, such as those returned by [getFMort()] or
#'   [getPredMort()] on a `MizerSim`.
#'
#' The same array objects can be passed to [ggplotly()] to produce interactive
#' versions, for example `ggplotly(getBiomass(sim))` or
#' `ggplotly(getEncounter(params))`. To add another compatible array to an
#' existing ggplot, use [addPlot()]. To compare two compatible mizer arrays
#' directly, use [plot2()]. To plot cumulative distributions over body size,
#' use [plotCDF()]. To visualise how spectra or rates change through time, use
#' [animate()] on a `MizerSim` or an
#' `ArrayTimeBySpeciesBySize` object.
#'
#' The named plotting functions give more specialised control. This table shows
#' the available named plotting functions.
#' \tabular{ll}{
#'   Plot \tab Description \cr
#'   [plotBiomass()] \tab Plots the total biomass of each species through time. A time range to be plotted can be specified. The size range of the community can be specified in the same way as for [getBiomass()]. \cr
#'   [plotYield()] \tab Plots the total yield of each species across all fishing gears against time. \cr
#'   [plotYieldGear()] \tab Plots the total yield of each species by gear against time. \cr
#'   [plotSpectra()] \tab Plots the abundance (biomass or numbers) spectra of each species and the background community. It is possible to specify a minimum size which is useful for truncating the plot. \cr
#'   [plotCDF()] \tab Plots cumulative distributions of abundance or biomass over size. \cr
#'   [plotCDF2()] \tab Compares cumulative distributions from two simulations or parameter objects in one plot. \cr
#'   [plotSpectra2()] \tab Compares the spectra from two simulations or parameter objects in one plot. \cr
#'   [plotFeedingLevel()] \tab Plots the feeding level of each species against size. \cr
#'   [plotPredMort()] \tab Plots the predation mortality of each species against size. \cr
#'   [plotFMort()] \tab Plots the total fishing mortality of each species against size. \cr
#'   [plotGrowthCurves()] \tab Plots the size as a function of age. \cr
#'   [plotDiet()] \tab Plots the diet composition at size for a given predator species. \cr
#'   [plotBiomassObservedVsModel()] \tab Compares observed biomass with model biomass. \cr
#'   [plotYieldObservedVsModel()] \tab Compares observed yield with model yield. \cr
#'   [animate()] \tab Animates spectra or rate arrays through time. The older [animateSpectra()] name is retained as an alias. \cr
#' }
#'
#' The static plotting functions use ggplot2 and return a ggplot object. This
#' means that you can manipulate the plot further after its creation using the
#' ggplot grammar of graphics. The named high-level plot functions have plotly
#' counterparts, for example [plotlyBiomass()] or [plotlySpectra()], for
#' interactive exploration. Generic and compositional plotting APIs, such as
#' [plot()], [plot2()], [plotRelative()] and [addPlot()], do not have separate
#' plotly wrappers. Use [ggplotly()] on the ggplot object they return.
#'
#' While most plot functions take their data from a MizerSim object, some of
#' those that make plots representing data at a single time can also take their
#' data from the initial values in a MizerParams object.
#'
#' Where plots show results for species, the line colour and line type for each
#' species are specified by the `linecolour` and `linetype` slots in
#' the MizerParams object. These were either taken from a default palette
#' hard-coded into [emptyParams()] or they were specified by the user
#' in the species parameters dataframe used to set up the MizerParams object.
#' The `linecolour` and `linetype` slots hold named vectors, named by
#' the species. They can be overwritten by the user at any time.
#'
#' Most plots allow the user to select to show only a subset of species,
#' specified as a vector in the `species` argument to the plot function.
#'
#' The ordering of the species in the legend is the same as the ordering in
#' the species parameter data frame.
#'
#' @seealso [summary_functions], [indicator_functions]
#' @family plotting functions
#' @name plotting_functions
#' @examples
#' \donttest{
#' sim <- NS_sim
#'
#' # Generic plot methods
#' plot(sim)
#' plot(getBiomass(sim), species = c("Cod", "Herring"))
#' ggplotly(getBiomass(sim))
#'
#' # Named plot functions
#' plotFeedingLevel(sim)
#'
#' # Plotting only a subset of species
#' plotFeedingLevel(sim, species = c("Cod", "Herring"))
#'
#' # Adding another compatible array to an existing plot
#' p <- plot(getBiomass(sim), species = "Cod")
#' addPlot(p, getBiomass(sim), species = "Herring", linetype = "dashed")
#'
#' # Specifying new colours and linetypes for some species
#' sim@params@linetype["Cod"] <- "dashed"
#' sim@params@linecolour["Cod"] <- "red"
#' plotFeedingLevel(sim, species = c("Cod", "Herring"))
#'
#' # Manipulating the plot
#' library(ggplot2)
#' p <- plotFeedingLevel(sim)
#' p <- p + geom_hline(aes(yintercept = 0.7))
#' p <- p + theme_bw()
#' p
#' }
NULL

# Hackiness to get past the 'no visible binding ... ' warning when running check
utils::globalVariables(c("time", "value", "Species", "w", "gear", "Age",
                         "x", "y", "Year", "Yield", "Biomass", "Size",
                         "Proportion", "Prey", "Legend", "Type", "Gear",
                         "Predator", "weight", "a", "b", "age", "w_max",
                         "Model", "rel_diff"))

#' Make a plot from a data frame
#'
#' This is used internally by most plotting functions.
#'
#' @param frame A data frame with at least three variables.
#'   The first three variables are used, in that order, as:
#'   1. Variable to be plotted on x-axis
#'   2. Variable to be plotted on y-axis
#'   3. Grouping variable
#' @param params A MizerParams object, which is used for the line colours and
#'   line types.
#' @param style The style of the plot. Available options are `"line"` for
#'   `geom_line()` and `"area"` for `geom_area()`. Default is `"line"`.
#' @param legend_var The name of the variable that should be used in the legend
#'   and to determine the line style. If NULL then the grouping variable is
#'   used for this purpose.
#' @param wrap_var Optional. The name of the variable that should be used for
#'  creating wrapped facets.
#' @param wrap_scale Optional. Used to pass the scales argument to facet_wrap().
#' @param xlab Label for the x-axis
#' @param ylab Label for the y-axis
#' @param xtrans Transformation for the x-axis. Often "log10" may be useful
#'   instead of the default of "identity".
#' @param ytrans Transformation for the y-axis.
#' @param xlim A numeric vector of length two providing lower and upper limits
#'   for the x axis. Use NA to refer to the existing minimum or maximum.
#' @param ylim A numeric vector of length two providing lower and upper limits
#'   for the y axis. Use NA to refer to the existing minimum or maximum.
#' @param y_ticks The approximate number of ticks desired on the y axis
#' @param highlight Name or vector of names of the species to be highlighted.
#' @return A ggplot2 object
#' @keywords internal
#' @export
plotDataFrame <- function(frame, params, style = "line", xlab = waiver(),
                          ylab = waiver(), xtrans = "identity", ytrans = "identity",
                          xlim = c(NA, NA), ylim = c(NA, NA),
                          y_ticks = 6, highlight = NULL, legend_var = NULL,
                          wrap_var = NULL, wrap_scale = NULL) {
    assert_that(is.data.frame(frame),
                is(params, "MizerParams"))
    if (ncol(frame) < 3) {
        stop("The data frame needs to have at least 3 variables.")
    }

    var_names <- names(frame)
    x_var <- var_names[[1]]
    y_var <- var_names[[2]]
    group_var <- var_names[[3]]
    if (is.null(legend_var)) {
        frame$Legend <- frame[[group_var]]
        legend_var <- "Legend"
    } else {
        if (!(legend_var %in% var_names)) {
            stop("The `legend_var` argument must be the name of a variable ",
                 "in the data frame.")
        }
    }

    # Need to keep species in order for legend
    legend_levels <-
        intersect(names(params@linecolour), frame[[legend_var]])
    frame[[legend_var]] <- factor(frame[[legend_var]], levels = legend_levels)

    if (sum(is.na(frame$Legend))) {
        warning("missing legend in params@linecolour, some groups won't be displayed")
    }

    linecolour <- params@linecolour[legend_levels]
    linetype <- params@linetype[legend_levels]
    linesize <- make_linesize(legend_levels, highlight)

    xbreaks <- waiver()
    if (xtrans == "log10") xbreaks <- log_breaks()
    ybreaks <- waiver()
    if (ytrans == "log10") ybreaks <- log_breaks(n = y_ticks)

    # Set up axis limits. NA values mean auto-scale to data range.
    # The reason why below `group = species` is included in `ggplot()`
    # rather than in `geom_line` is because that puts it first in the
    # plotly tooltips, due to a bug in plotly.
    p <- ggplot(frame, aes(group = .data[[group_var]])) +
        scale_y_continuous(trans = ytrans, breaks = ybreaks,
                           labels = prettyNum, name = ylab,
                           limits = ylim) +
        scale_x_continuous(trans = xtrans, breaks = xbreaks, name = xlab,
                           limits = xlim)

    switch(style,
           "line" = {
               p <- p +
                   geom_line(aes(x = .data[[x_var]], y = .data[[y_var]],
                                 colour = .data[[legend_var]],
                                 linetype = .data[[legend_var]],
                                 linewidth = .data[[legend_var]])) +
                   scale_colour_manual(values = linecolour) +
                   scale_linetype_manual(values = linetype) +
                   scale_discrete_manual("linewidth", values = linesize)
           },
           "area" = {
               p <- p +
                   geom_area(aes(x = .data[[x_var]], y = .data[[y_var]],
                                 fill = .data[[legend_var]])) +
                   scale_fill_manual(values = linecolour)
           },
           stop("unknown style selected")
    )

    if (!is.null(wrap_var)) {
        if (!(wrap_var %in% var_names)) {
            stop("The `wrap_var` argument must be the name of a variable ",
                 "in the data frame.")
        }
        p <- p + facet_wrap(wrap_var, scales = wrap_scale)
    }

    make_mizer_plot(p, mizer_tooltip_vars(frame, group_var, x_var, y_var,
                                          legend_var))
}

make_linesize <- function(levels, highlight) {
    linesize <- rep(0.8, length(levels))
    names(linesize) <- levels
    linesize[highlight] <- 1.6
    linesize
}

make_mizer_plot <- function(plot, tooltip) {
    attr(plot, "mizer_tooltip") <- tooltip
    class(plot) <- unique(c("mizer_plot", class(plot)))
    plot
}

mizer_tooltip_vars <- function(frame, group_var, x_var, y_var,
                               legend_var = NULL, extra = NULL) {
    tooltip <- c(group_var, x_var, y_var)
    if (!is.null(legend_var) && legend_var %in% names(frame) &&
            !identical(legend_var, group_var) &&
            any(as.character(frame[[legend_var]]) !=
                    as.character(frame[[group_var]]), na.rm = TRUE)) {
        tooltip <- c(tooltip, legend_var)
    }
    unique(c(tooltip, extra))
}

#' @exportS3Method plotly::ggplotly
ggplotly.mizer_plot <- function(p = ggplot2::last_plot(), ...,
                                tooltip = attr(p, "mizer_tooltip") %||%
                                    "all") {
    class(p) <- setdiff(class(p), "mizer_plot")
    ggplotly(p, ..., tooltip = tooltip)
}

plotComparisonDataFrame <- function(frame1, frame2, params,
                                    name1 = "First", name2 = "Second",
                                    xlab = waiver(), ylab = waiver(),
                                    xtrans = "identity", ytrans = "identity",
                                    xlim = c(NA, NA), ylim = c(NA, NA),
                                    y_ticks = 6, legend_var = "Legend") {
    assert_that(is.data.frame(frame1),
                is.data.frame(frame2),
                is(params, "MizerParams"))

    names(frame2)[seq_len(min(3, ncol(frame2)))] <-
        names(frame1)[seq_len(min(3, ncol(frame1)))]
    frame1$Model <- name1
    frame2$Model <- name2
    frame <- rbind(frame1, frame2)
    frame$Model <- factor(frame$Model, levels = c(name1, name2))

    var_names <- names(frame)
    x_var <- var_names[[1]]
    y_var <- var_names[[2]]
    group_var <- var_names[[3]]
    if (!(legend_var %in% var_names)) {
        stop("The `legend_var` argument must be the name of a variable ",
             "in the data frame.")
    }

    legend_levels <- intersect(names(params@linecolour), frame[[legend_var]])
    frame[[legend_var]] <- factor(frame[[legend_var]], levels = legend_levels)
    if (sum(is.na(frame[[legend_var]]))) {
        warning("missing legend in params@linecolour, some groups won't be displayed")
    }
    linecolour <- params@linecolour[legend_levels]

    xbreaks <- waiver()
    if (xtrans == "log10") xbreaks <- log_breaks()
    ybreaks <- waiver()
    if (ytrans == "log10") ybreaks <- log_breaks(n = y_ticks)

    p <- ggplot(frame,
                aes(group = interaction(.data[[group_var]], .data[["Model"]]))) +
        scale_y_continuous(trans = ytrans, breaks = ybreaks,
                           labels = prettyNum, name = ylab,
                           limits = ylim) +
        scale_x_continuous(trans = xtrans, breaks = xbreaks, name = xlab,
                           limits = xlim) +
        geom_line(aes(x = .data[[x_var]], y = .data[[y_var]],
                      colour = .data[[legend_var]],
                      linetype = .data[["Model"]])) +
        scale_colour_manual(values = linecolour) +
        scale_linetype_discrete(drop = FALSE)
    make_mizer_plot(p, mizer_tooltip_vars(frame, group_var, x_var, y_var,
                                          legend_var, extra = "Model"))
}

plotRelativeDataFrame <- function(frame1, frame2, params,
                                  xlab = waiver(),
                                  xtrans = "identity",
                                  xlim = c(NA, NA),
                                  ylim = c(NA, NA),
                                  legend_var = "Legend") {
    assert_that(is.data.frame(frame1),
                is.data.frame(frame2),
                is(params, "MizerParams"))

    names(frame2)[seq_len(min(3, ncol(frame2)))] <-
        names(frame1)[seq_len(min(3, ncol(frame1)))]
    var_names <- names(frame1)
    x_var <- var_names[[1]]
    y_var <- var_names[[2]]
    group_var <- var_names[[3]]
    if (!(legend_var %in% var_names)) {
        stop("The `legend_var` argument must be the name of a variable ",
             "in the data frame.")
    }

    by_vars <- c(x_var, group_var, legend_var)
    frame <- dplyr::inner_join(frame1, frame2, by = by_vars,
                               suffix = c(".x", ".y"))
    frame$rel_diff <- relative_difference(frame[[paste0(y_var, ".x")]],
                                          frame[[paste0(y_var, ".y")]])
    frame <- frame[is.finite(frame$rel_diff), ]

    legend_levels <- intersect(names(params@linecolour), frame[[legend_var]])
    frame[[legend_var]] <- factor(frame[[legend_var]], levels = legend_levels)
    if (sum(is.na(frame[[legend_var]]))) {
        warning("missing legend in params@linecolour, some groups won't be displayed")
    }
    linecolour <- params@linecolour[legend_levels]

    xbreaks <- waiver()
    if (xtrans == "log10") xbreaks <- log_breaks()

    p <- ggplot(frame, aes(group = .data[[group_var]])) +
        scale_y_continuous(name = "Relative difference", limits = ylim) +
        scale_x_continuous(trans = xtrans, breaks = xbreaks, name = xlab,
                           limits = xlim) +
        geom_hline(yintercept = 0, linetype = 1,
                   colour = "dark grey", linewidth = 0.75) +
        geom_line(aes(x = .data[[x_var]], y = .data[["rel_diff"]],
                      colour = .data[[legend_var]])) +
        scale_colour_manual(values = linecolour)
    make_mizer_plot(p, mizer_tooltip_vars(frame, group_var, x_var, "rel_diff",
                                          legend_var))
}

relative_difference <- function(first, second) {
    2 * (second - first) / (first + second)
}

#' Helper function to produce nice breaks on logarithmic axes
#'
#' This is needed when the logarithmic y-axis spans less than one order of
#' magnitude, in which case the ggplot2 default produces no ticks.
#'
#' Thanks to Heather Turner at
#' \url{https://stackoverflow.com/questions/14255533/pretty-ticks-for-log-normal-scale-using-ggplot2-dynamic-not-manual}
#'
#' @param n Approximate number of ticks
#'
#' @return A function that can be used as the break argument in calls to
#'   scale_y_continuous() or scale_x_continuous()
#' @export
#' @keywords internal
#' @family helper
log_breaks <- function(n = 6) {
    n <- max(1, n)  # Because n=0 could lead to R crash
    function(x) {
        grDevices::axisTicks(log10(range(x, na.rm = TRUE)),
                             log = TRUE, nint = n)
    }
}


#' Plot the biomass of species through time
#'
#' After running a projection, the biomass of each species can be plotted
#' against time. The biomass is calculated within user defined size limits
#' (see [getBiomass()]).
#'
#' @param object An object of class \linkS4class{MizerSim}
#' @inheritParams valid_species_arg
#' @param start_time The first time to be plotted. Default (`NULL`) is the
#'   beginning of the time series.
#' @param end_time The last time to be plotted. Default (`NULL`) is the end of
#'   the time series.
#' @param y_ticks The approximate number of ticks desired on the y axis.
#' @param ylim A numeric vector of length two providing lower and upper limits
#'   for the y axis. Use `NA` to refer to the existing minimum or maximum. Any
#'   values below 1e-20 are always cut off.
#' @param total A boolean value that determines whether the total biomass from
#'   all species is plotted as well. Default is FALSE.
#' @param background A boolean value that determines whether background species
#'   are included. Ignored if the model does not contain background species.
#'   Default is TRUE.
#' @param highlight Name or vector of names of the species to be highlighted.
#' @param log If `TRUE` (default), use a log10 y-axis.
#' @param return_data A boolean value that determines whether the formatted data
#'   used for the plot is returned instead of the plot itself. Default is FALSE.
#' @param use_cutoff If TRUE, the `biomass_cutoff` column in the species
#'   parameters is used as the minimum weight for each species.
#' @inheritParams get_size_range_array
#'
#' @return A ggplot2 object, unless `return_data = TRUE`, in which case a data
#'   frame with the four variables 'Year', 'Biomass', 'Species', 'Legend' is
#'   returned.
#' @export
#' @family plotting functions
#' @seealso [plotting_functions], [getBiomass()]
#' @examples
#' \donttest{
#' plotBiomass(NS_sim)
#' plotBiomass(NS_sim, species = c("Sandeel", "Herring"), total = TRUE)
#' plotBiomass(NS_sim, start_time = 1980, end_time = 1990)
#'
#' # Returning the data frame
#' fr <- plotBiomass(NS_sim, return_data = TRUE)
#' str(fr)
#' }
#' @rdname plotBiomass
#' @export
plotBiomass <- function(object, ...) {
    UseMethod("plotBiomass")
}

#' @rdname plotBiomass
#' @export
plotBiomass.MizerSim <- function(object, species = NULL,
                        start_time = NULL, end_time = NULL,
                        y_ticks = 6, ylim = c(NA, NA),
                        total = FALSE, background = TRUE,
                        highlight = NULL, log = TRUE,
                        return_data = FALSE,
                        use_cutoff = FALSE,
                        min_w = min(object@params@w),
                        max_w = max(object@params@w),
                        min_l = NULL, max_l = NULL, ...) {
    bm <- getBiomass(object, use_cutoff = use_cutoff,
                     min_w = min_w, max_w = max_w,
                     min_l = min_l, max_l = max_l)
    plot(bm, species = species,
         start_time = start_time, end_time = end_time,
         y_ticks = y_ticks, ylim = ylim,
         total = total, background = background,
         highlight = highlight, log_y = log,
         return_data = return_data)
}

#' @rdname plotBiomass
#' @export
plotlyBiomass <- function(object,
             species = NULL,
             start_time = NULL,
             end_time = NULL,
             y_ticks = 6,
             ylim = c(NA, NA),
             total = FALSE,
             background = TRUE,
             highlight = NULL,
             log = TRUE,
             use_cutoff = FALSE,
             min_w = min(object@params@w),
             max_w = max(object@params@w),
             min_l = NULL,
             max_l = NULL) {
    argg <- as.list(environment())
    ggplotly(do.call("plotBiomass", argg),
             tooltip = c("Species", "Year", "Biomass"))
}


#' Plot the total yield of species through time
#'
#' After running a projection, the total yield of each species across all
#' fishing gears can be plotted against time. The yield is obtained with
#' [getYield()].
#'
#' @param object An object of class \linkS4class{MizerSim}
#' @param sim2 An optional second object of class \linkS4class{MizerSim}. If
#'   this is provided its yields will be shown on the same plot in bolder lines.
#' @inheritParams plotSpectra
#' @param log Boolean whether yield should be plotted on a logarithmic axis.
#'   Defaults to true.
#'
#' @return A ggplot2 object, unless `return_data = TRUE`, in which case a data
#'   frame with the three variables 'Year', 'Yield', 'Species' is returned.
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
plotYield <- function(object, ...) {
    UseMethod("plotYield")
}

#' @rdname plotYield
#' @export
plotYield.MizerSim <- function(object, sim2,
                      species = NULL,
                      total = FALSE, log = TRUE,
                      highlight = NULL, return_data = FALSE,
                      ...) {
    assert_that(is(object, "MizerSim"),
                is.flag(total),
                is.flag(log),
                is.flag(return_data))
    params <- object@params
    species <- valid_species_arg(object, species, error_on_empty = TRUE)
    if (missing(sim2)) {
        y <- getYield(object, ...)
        y_total <- rowSums(y)
        y <- y[, (as.character(dimnames(y)[[2]]) %in% species),
               drop = FALSE]
        if (total) {
            # Include total
            y <- cbind(y, "Total" = y_total)
        }
        plot_dat <- reshape2::melt(y, varnames = c("Year", "Species"),
                                   value.name = "Yield")
        plot_dat <- subset(plot_dat, plot_dat$Yield > 0)
        # plotDataFrame() needs the columns in a particular order
        plot_dat <- plot_dat[, c(1, 3, 2)]

        if (nrow(plot_dat) == 0) {
            warning("There is no yield to include.")
        }
        if (return_data) return(plot_dat)

        plotDataFrame(plot_dat, params,
                      ylab = "Yield [g/year]",
                      ytrans = ifelse(log, "log10", "identity"),
                      highlight = highlight)
    } else {
        # We need to combine two plots
        if (!all(dimnames(object@n)$time == dimnames(sim2@n)$time)) {
            stop("The two simulations do not have the same times")
        }
        ym <- plotYield(object, species = species,
                            total = total, log = log,
                            highlight = highlight, return_data = TRUE, ...)
        ym2 <- plotYield(sim2, species = species,
                            total = total, log = log,
                            highlight = highlight, return_data = TRUE, ...)
        ym$Simulation <- rep(1, nrow(ym)) # We don't use recycling because that
                                          # fails when there are zero rows.
        ym2$Simulation <- rep(2, nrow(ym2))
        ym <- rbind(ym, ym2)

        if (return_data) return(ym)

        plotDataFrame(ym, params,
                      ylab = "Yield [g/year]",
                      ytrans = ifelse(log, "log10", "identity"),
                      highlight = highlight, wrap_var = "Simulation")
    }
}

#' @rdname plotYield
#' @export
plotlyYield <- function(object, sim2,
                        species = NULL,
                        total = FALSE, log = TRUE,
                        highlight = NULL, ...) {
    argg <- as.list(environment())
    ggplotly(do.call("plotYield", argg),
             tooltip = c("Species", "Year", "Yield"))
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
#' @inheritParams plotSpectra
#' @param gears A vector of gear names to be included in the plot. Default is
#'  all gears.
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
plotYieldGear <- function(object, ...) {
    UseMethod("plotYieldGear")
}

#' @rdname plotYieldGear
#' @export
plotYieldGear.MizerSim <- function(object,
                          species = NULL,
                          gears = NULL,
                          total = FALSE,
                          highlight = NULL, return_data = FALSE,
                          ...) {
    assert_that(is(object, "MizerSim"),
                is.flag(total),
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
    linesize <- make_linesize(species_levels, highlight)
    ggplot(ym) +
        geom_line(aes(x = Year, y = Yield, colour = Species,
                      linetype = Gear, linewidth = Species)) +
        scale_y_continuous(trans = "log10", name = "Yield [g]") +
        scale_colour_manual(values = params@linecolour[species_levels]) +
        scale_discrete_manual("linewidth", values = linesize)
}

#' @rdname plotYieldGear
#' @export
plotlyYieldGear <- function(object, species = NULL,
                            total = FALSE, highlight = NULL, ...) {
    argg <- as.list(environment())
    ggplotly(do.call("plotYieldGear", argg),
             tooltip = c("Species", "Year", "Yield"))
}

#' Plot the abundance spectra
#'
#' Plots the number density multiplied by a power of the weight, with the power
#' specified by the `power` argument.
#'
#' When called with a \linkS4class{MizerSim} object, the abundance is averaged
#' over the specified time range (a single value for the time range can be used
#' to plot a single time step). When called with a \linkS4class{MizerParams}
#' object the initial abundance is plotted.
#'
#' @param object An object of class \linkS4class{MizerSim} or
#'   \linkS4class{MizerParams}.
#' @inheritParams valid_species_arg
#' @param time_range The time range (either a vector of values, a vector of min
#'   and max time, or a single value) to average the abundances over. Default is
#'   the final time step. Ignored when called with a \linkS4class{MizerParams}
#'   object.
#' @param geometric_mean `r lifecycle::badge("experimental")`
#'   If TRUE then the average of the abundances over the
#'   time range is a geometric mean instead of the default arithmetic mean.
#' @param wlim A numeric vector of length two providing lower and upper limits
#'   for the w axis. Use NA for the default: the lower default is
#'   `min(params@w) / 100` when `resource = TRUE` (to show some resource below
#'   the fish grid) or `min(params@w)` when `resource = FALSE`; the upper
#'   default is `max(params@w_full)`. Data is filtered to this range and the
#'   axis limits are set accordingly.
#' @param ylim A numeric vector of length two providing lower and upper limits
#'   for the y axis. Use NA to auto-scale to the data range. Values below 1e-20
#'   are always filtered out from the data regardless of `ylim[1]`. Data above
#'   `ylim[2]` is filtered and the upper axis limit is set accordingly.
#' @param power The abundance is plotted as the number density times the weight
#' raised to `power`. The default \code{power = 1} gives the biomass
#' density, whereas \code{power = 2} gives the biomass density with respect
#' to logarithmic size bins.
#' @param biomass `r lifecycle::badge("deprecated")`
#'  Only used if `power` argument is missing. Then
#'   \code{biomass = TRUE} is equivalent to \code{power=1} and
#'   \code{biomass = FALSE} is equivalent to \code{power=0}
#' @param total A boolean value that determines whether the total over all
#'   species in the system is plotted as well. Note that even if the plot
#'   only shows a selection of species, the total is including all species.
#'   Default is FALSE.
#' @param resource A boolean value that determines whether resource is included.
#'   Default is TRUE.
#' @param background A boolean value that determines whether background species
#'   are included. Ignored if the model does not contain background species.
#'   Default is TRUE.
#' @param highlight Name or vector of names of the species to be highlighted.
#' @param log_x If `TRUE` (default), use a log10 x-axis.
#' @param log_y If `TRUE` (default), use a log10 y-axis.
#' @param log Character string specifying which axes should use log10 scales,
#'   in the same form as the base [plot()] argument. For example, `"x"`,
#'   `"y"`, `"xy"` or `""`. If supplied, this overrides `log_x` and `log_y`.
#' @param return_data A boolean value that determines whether the formatted data
#' used for the plot is returned instead of the plot itself. Default value is FALSE
#' @param ... Other arguments (currently unused)
#'
#' @return A ggplot2 object, unless `return_data = TRUE`, in which case a data
#'   frame with the four variables 'w', 'value', 'Species', 'Legend' is
#'   returned.
#' @export
#' @family plotting functions
#' @seealso [plotting_functions]
#' @examples
#' \donttest{
#' params <-  NS_params
#' sim <- project(params, effort=1, t_max=20, t_save = 2, progress_bar = FALSE)
#' plotSpectra(sim)
#' plotSpectra(sim, wlim = c(1e-6, NA))
#' plotSpectra(sim, time_range = 10:20)
#' plotSpectra(sim, time_range = 10:20, power = 0)
#' plotSpectra(sim, species = c("Cod", "Herring"), power = 1)
#'
#' # Returning the data frame
#' fr <- plotSpectra(sim, return_data = TRUE)
#' str(fr)
#' }
#' @rdname plotSpectra
#' @export
plotSpectra <- function(object, ...) {
    UseMethod("plotSpectra")
}

#' @rdname plotSpectra
#' @export
plotSpectra.MizerSim <- function(object, species = NULL,
                        time_range,
                        geometric_mean = FALSE,
                        wlim = c(NA, NA), ylim = c(NA, NA),
                        power = 1, biomass = TRUE,
                        total = FALSE, resource = TRUE,
                        background = TRUE,
                        highlight = NULL, log_x = TRUE, log_y = TRUE,
                        log = NULL, return_data = FALSE, ...) {
    # to deal with old-type biomass argument
    if (missing(power)) {
        power <- as.numeric(biomass)
    }
    log_axes <- parsePlotLog(log, log_x = log_x, log_y = log_y)
    log_x <- log_axes$log_x
    log_y <- log_axes$log_y

    assert_that(is.flag(total), is.flag(resource),
                is.flag(background),
                is.number(power),
                length(wlim) == 2,
                length(ylim) == 2)
    species <- valid_species_arg(object, species)
    if (length(species) == 0 && !total && !resource) {
        stop("There is nothing to plot as no valid species have been selected.")
    }
    if (missing(time_range)) {
        time_range  <- max(as.numeric(dimnames(object@n)$time))
    }
    time_elements <- get_time_elements(object, time_range)
    mean_fn <- mean
    if (geometric_mean) {
        mean_fn <- function(x) {
            exp(mean(log(x)))
        }
    }
    n <- apply(object@n[time_elements, , , drop = FALSE], c(2, 3), mean_fn)
    n_pp <- apply(object@n_pp[time_elements, , drop = FALSE], 2, mean_fn)
    plot_spectra(object@params, n = n, n_pp = n_pp,
                 species = species, wlim = wlim, ylim = ylim,
                 power = power, total = total, resource = resource,
                 background = background, highlight = highlight,
                 log_x = log_x, log_y = log_y,
                 return_data = return_data)
}

#' @rdname plotSpectra
#' @export
plotSpectra.MizerParams <- function(object, species = NULL,
                        wlim = c(NA, NA), ylim = c(NA, NA),
                        power = 1, biomass = TRUE,
                        total = FALSE, resource = TRUE,
                        background = TRUE,
                        highlight = NULL, log_x = TRUE, log_y = TRUE,
                        log = NULL, return_data = FALSE, ...) {
    # to deal with old-type biomass argument
    if (missing(power)) {
        power <- as.numeric(biomass)
    }
    log_axes <- parsePlotLog(log, log_x = log_x, log_y = log_y)
    log_x <- log_axes$log_x
    log_y <- log_axes$log_y

    assert_that(is.flag(total), is.flag(resource),
                is.flag(background),
                is.number(power),
                length(wlim) == 2,
                length(ylim) == 2)
    species <- valid_species_arg(object, species)
    if (length(species) == 0 && !total && !resource) {
        stop("There is nothing to plot as no valid species have been selected.")
    }
    plot_spectra(object, n = object@initial_n,
                 n_pp = object@initial_n_pp,
                 species = species, wlim = wlim, ylim = ylim,
                 power = power, total = total, resource = resource,
                 background = background, highlight = highlight,
                 log_x = log_x, log_y = log_y,
                 return_data = return_data)
}


plot_spectra <- function(params, n, n_pp,
                         species, wlim, ylim, power,
                         total, resource, background,
                         highlight, log_x, log_y, return_data) {
    params <- validParams(params)
    if (is.na(wlim[1])) {
        wlim[1] <- if (resource) min(params@w) / 100 else min(params@w)
    }
    if (is.na(wlim[2])) {
        wlim[2] <- max(params@w_full)
    }

    if (total) {
        # Calculate total community abundance
        fish_idx <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)
        total_n <- n_pp
        total_n[fish_idx] <- total_n[fish_idx] + colSums(n)
        total_n <- total_n * params@w_full^power
    }
    species <- valid_species_arg(params, species)
    # Deal with power argument
    y_label <- spectra_y_label(power)
    n <- sweep(n, 2, params@w^power, "*")
    # Select only the desired species
    spec_n <- n[as.character(dimnames(n)[[1]]) %in% species, , drop = FALSE]
    # Make data.frame for plot
    plot_dat <- data.frame(w = rep(params@w,
                                   each = dim(spec_n)[[1]]),
                           value = c(spec_n),
                           Species = dimnames(spec_n)[[1]],
                           Legend = dimnames(spec_n)[[1]])
    if (resource) {
        resource_sel <- (params@w_full >= wlim[1]) &
                        (params@w_full <= wlim[2])
        # Do we have any resource to plot?
        if (sum(resource_sel) > 0) {
            w_resource <- params@w_full[resource_sel]
            plank_n <- n_pp[resource_sel] * w_resource^power
            plot_dat <- rbind(plot_dat,
                              data.frame(w = w_resource,
                                         value = c(plank_n),
                                         Species = "Resource",
                                         Legend = "Resource")
            )
        }
    }
    if (total) {
        plot_dat <- rbind(plot_dat,
                          data.frame(w = params@w_full,
                                     value = c(total_n),
                                     Species = "Total",
                                     Legend = "Total")
                          )
    }
    if (background && any(params@species_params$is_background)) {
        back_n <- n[params@species_params$is_background, , drop = FALSE]
        plot_dat <-
            rbind(plot_dat,
                  data.frame(w = rep(params@w,
                                     each = dim(back_n)[[1]]),
                             value = c(back_n),
                             Species = as.factor(dimnames(back_n)[[1]]),
                             Legend = "Background")
            )
    }
    # lop off 0s and apply wlim
    plot_dat <- plot_dat[(plot_dat$value > 0) &
                             (plot_dat$w >= wlim[1]) &
                             (plot_dat$w <= wlim[2]), ]
    # Impose ylim
    if (!is.na(ylim[2])) {
        plot_dat <- plot_dat[plot_dat$value <= ylim[2], ]
    }
    filter_min <- if (is.na(ylim[1])) 1e-20 else ylim[1]
    plot_dat <- plot_dat[plot_dat$value > filter_min, ]

    if (return_data) return(plot_dat)

    plotDataFrame(plot_dat, params, xlab = "Size [g]", ylab = y_label,
                  xtrans = if (log_x) "log10" else "identity",
                  ytrans = if (log_y) "log10" else "identity",
                  xlim = wlim, ylim = ylim,
                  highlight = highlight, legend_var = "Legend")
}

#' Plot cumulative abundance or biomass distributions
#'
#' `plotCDF()` plots the cumulative distribution over body size from small to
#' large sizes. It uses the same spectra data preparation as [plotSpectra()].
#' The density is first multiplied by `w^power`, then integrated over size.
#' With `normalise = TRUE`, each curve is divided by its final value so that it
#' ends at 1.
#'
#' @inheritParams plotSpectra
#' @param resource A boolean value that determines whether resource is included.
#'   Default is FALSE.
#' @param normalise If `TRUE` (default), plot the cumulative proportion. If
#'   `FALSE`, plot the cumulative abundance, biomass, or other unnormalised
#'   integral.
#' @param log_x If `TRUE` (default), use a log10 x-axis.
#' @param log_y If `TRUE`, use a log10 y-axis. Default is `FALSE`.
#' @param log Character string specifying which axes should use a log10 scale,
#'   in the same form as the base [plot()] argument. If supplied, this overrides
#'   `log_x` and `log_y`.
#'
#' @return A ggplot2 object, unless `return_data = TRUE`, in which case a data
#'   frame with the four variables 'w', 'value', 'Species', 'Legend' is
#'   returned.
#' @export
#' @family plotting functions
#' @seealso [plotSpectra()]
#' @examples
#' \donttest{
#' plotCDF(NS_params, species = c("Cod", "Herring"))
#' plotCDF(NS_sim, power = 0, normalise = FALSE)
#' }
plotCDF <- function(object, ...) {
    UseMethod("plotCDF")
}

#' @rdname plotCDF
#' @export
plotCDF.MizerSim <- function(object, species = NULL,
                             time_range,
                             geometric_mean = FALSE,
                             wlim = c(NA, NA), ylim = c(NA, NA),
                             power = 1, biomass = TRUE,
                             total = FALSE, resource = FALSE,
                             background = TRUE,
                             highlight = NULL, normalise = TRUE,
                             log_x = TRUE, log_y = FALSE, log = NULL,
                             return_data = FALSE, ...) {
    if (missing(power)) {
        power <- as.numeric(biomass)
    }
    log_axes <- parsePlotLog(log, log_x = log_x, log_y = log_y)
    assert_that(is.flag(total), is.flag(resource),
                is.flag(background), is.flag(normalise),
                is.number(power),
                length(wlim) == 2,
                length(ylim) == 2)

    args <- list(object = object, species = species,
                 geometric_mean = geometric_mean,
                 wlim = wlim, ylim = c(NA, NA),
                 power = power, total = total,
                 resource = resource, background = background,
                 return_data = TRUE)
    if (!missing(time_range)) {
        args$time_range <- time_range
    }
    plot_dat <- do.call(plotSpectra, args)
    plot_cdf(plot_dat, object@params, power = power, normalise = normalise,
             log_x = log_axes$log_x, log_y = log_axes$log_y,
             wlim = wlim, ylim = ylim,
             highlight = highlight, return_data = return_data)
}

#' @rdname plotCDF
#' @export
plotCDF.MizerParams <- function(object, species = NULL,
                                wlim = c(NA, NA), ylim = c(NA, NA),
                                power = 1, biomass = TRUE,
                                total = FALSE, resource = FALSE,
                                background = TRUE,
                                highlight = NULL, normalise = TRUE,
                                log_x = TRUE, log_y = FALSE, log = NULL,
                                return_data = FALSE, ...) {
    if (missing(power)) {
        power <- as.numeric(biomass)
    }
    log_axes <- parsePlotLog(log, log_x = log_x, log_y = log_y)
    assert_that(is.flag(total), is.flag(resource),
                is.flag(background), is.flag(normalise),
                is.number(power),
                length(wlim) == 2,
                length(ylim) == 2)

    plot_dat <- plotSpectra(object, species = species,
                            wlim = wlim, ylim = c(NA, NA),
                            power = power, total = total,
                            resource = resource, background = background,
                            return_data = TRUE)
    plot_cdf(plot_dat, object, power = power, normalise = normalise,
             log_x = log_axes$log_x, log_y = log_axes$log_y,
             wlim = wlim, ylim = ylim,
             highlight = highlight, return_data = return_data)
}

plot_cdf <- function(plot_dat, params, power, normalise, log_x, log_y, wlim,
                     ylim, highlight, return_data) {
    cdf_dat <- prepare_spectra_cdf_data(plot_dat, params,
                                        normalise = normalise)
    if (return_data) return(cdf_dat)

    plotDataFrame(cdf_dat, validParams(params),
                  xlab = "Size [g]", ylab = cdf_y_label(power, normalise),
                  xtrans = if (log_x) "log10" else "identity",
                  ytrans = if (log_y) "log10" else "identity",
                  xlim = wlim, ylim = ylim,
                  highlight = highlight, legend_var = "Legend")
}

prepare_spectra_cdf_data <- function(plot_dat, params, normalise = TRUE) {
    params <- validParams(params)
    plot_dat <- plot_dat[order(plot_dat$Species, plot_dat$w), ]
    plot_dat$value <- plot_dat$value * spectra_bin_width(plot_dat$w, params)
    plot_dat$value <- stats::ave(plot_dat$value, plot_dat$Species,
                                 FUN = cumsum)
    if (normalise) {
        totals <- stats::ave(plot_dat$value, plot_dat$Species, FUN = max)
        plot_dat$value <- plot_dat$value / totals
    }
    plot_dat
}

spectra_bin_width <- function(w, params) {
    idx <- match(w, params@w_full)
    if (anyNA(idx)) {
        missing <- which(is.na(idx))
        for (i in missing) {
            idx[i] <- which.min(abs(params@w_full - w[i]))
        }
        if (!isTRUE(all.equal(w, params@w_full[idx], scale = 1))) {
            stop("Could not determine size-bin widths for the spectra data.")
        }
    }
    params@dw_full[idx]
}

cdf_y_label <- function(power, normalise) {
    if (normalise) {
        if (power == 0) {
            return("Cumulative proportion of abundance")
        }
        if (power == 1) {
            return("Cumulative proportion of biomass")
        }
        return("Cumulative proportion")
    }
    if (power == 0) {
        return("Cumulative abundance")
    }
    if (power == 1) {
        return("Cumulative biomass [g]")
    }
    paste0("Cumulative number density * w^", power)
}


#' Compare two cumulative abundance or biomass distributions
#'
#' `plotCDF2()` compares cumulative distributions from two `MizerParams` or
#' `MizerSim` objects in a single plot. Colours identify species or groups and
#' linetype identifies the object.
#'
#' @param object1 First `MizerParams` or `MizerSim` object.
#' @param object2 Second `MizerParams` or `MizerSim` object.
#' @param name1,name2 Labels for the two objects, used in the linetype legend.
#' @inheritParams plotCDF
#' @param ... Arguments passed to [plotCDF()] for preparing the cumulative
#'   distribution data, for example `species`, `time_range`, `wlim`,
#'   `resource`, `background` or `total`.
#'
#' @return A ggplot2 object.
#' @export
#' @family plotting functions
#'
#' @examples
#' \donttest{
#' sim1 <- project(NS_params, t_max = 10, progress_bar = FALSE)
#' sim2 <- project(NS_params, effort = 0.5, t_max = 10, progress_bar = FALSE)
#' plotCDF2(sim1, sim2, "Original", "Effort = 0.5")
#' }
plotCDF2 <- function(object1, object2, name1 = "First", name2 = "Second",
                     power = 1, normalise = TRUE, log_x = TRUE, log_y = FALSE,
                     log = NULL, resource = FALSE, ...) {
    log_axes <- parsePlotLog(log, log_x = log_x, log_y = log_y)
    assert_that(is.number(power), is.flag(normalise))

    args <- list(...)
    wlim <- args$wlim %||% c(NA, NA)
    ylim <- args$ylim %||% c(NA, NA)

    cf1 <- plotCDF(object1, power = power, normalise = normalise,
                   return_data = TRUE, ...)
    cf2 <- plotCDF(object2, power = power, normalise = normalise,
                   return_data = TRUE, ...)
    params <- if (is(object1, "MizerSim")) object1@params else object1

    plotComparisonDataFrame(cf1, cf2, validParams(params),
                            name1 = name1, name2 = name2,
                            xlab = "Size [g]",
                            ylab = cdf_y_label(power, normalise),
                            xtrans = if (log_axes$log_x) "log10" else "identity",
                            ytrans = if (log_axes$log_y) "log10" else "identity",
                            xlim = wlim, ylim = ylim,
                            legend_var = "Legend")
}

#' Compare two size spectra in the same plot
#'
#' `plotSpectra2()` compares the abundance spectra from two `MizerParams` or
#' `MizerSim` objects in a single plot. Colours identify species or groups and
#' linetype identifies the object.
#'
#' @param object1 First `MizerParams` or `MizerSim` object.
#' @param object2 Second `MizerParams` or `MizerSim` object.
#' @param name1,name2 Labels for the two objects, used in the linetype legend.
#' @inheritParams plotSpectra
#' @param log_x If `TRUE` (default), use a log10 x-axis.
#' @param log_y If `TRUE` (default), use a log10 y-axis.
#' @param log Character string specifying which axes should use log10 scales,
#'   in the same form as the base [plot()] argument. For example, `"x"`,
#'   `"y"`, `"xy"` or `""`. If supplied, this overrides `log_x` and `log_y`.
#' @param ... Arguments passed to [plotSpectra()] for preparing the spectra
#'   data, for example `species`, `time_range`, `wlim`, `ylim`, `resource`,
#'   `background` or `total`.
#'
#' @return A ggplot2 object.
#' @export
#' @family plotting functions
#'
#' @examples
#' \donttest{
#' sim1 <- project(NS_params, t_max = 10, progress_bar = FALSE)
#' sim2 <- project(NS_params, effort = 0.5, t_max = 10, progress_bar = FALSE)
#' plotSpectra2(sim1, sim2, "Original", "Effort = 0.5")
#' }
plotSpectra2 <- function(object1, object2, name1 = "First", name2 = "Second",
                         power = 1, log_x = TRUE, log_y = TRUE,
                         log = NULL, ...) {
    log_axes <- parsePlotLog(log, log_x = log_x, log_y = log_y)
    log_x <- log_axes$log_x
    log_y <- log_axes$log_y

    args <- list(...)
    wlim <- args$wlim %||% c(NA, NA)
    ylim <- args$ylim %||% c(NA, NA)

    sf1 <- plotSpectra(object1, power = power, return_data = TRUE, ...)
    sf2 <- plotSpectra(object2, power = power, return_data = TRUE, ...)
    params <- if (is(object1, "MizerSim")) object1@params else object1

    plotComparisonDataFrame(sf1, sf2, validParams(params),
                            name1 = name1, name2 = name2,
                            xlab = "Size [g]",
                            ylab = spectra_y_label(power),
                            xtrans = if (log_x) "log10" else "identity",
                            ytrans = if (log_y) "log10" else "identity",
                            xlim = wlim, ylim = ylim,
                            legend_var = "Legend")
}

spectra_y_label <- function(power) {
    if (power %in% c(0, 1, 2)) {
        return(c("Number density [1/g]", "Biomass density",
                 "Biomass density [g]")[power + 1])
    }
    paste0("Number density * w^", power)
}

#' @rdname plotSpectra2
#' @return `plotlySpectra2()` returns a plotly object.
#' @export
plotlySpectra2 <- function(object1, object2, name1 = "First",
                           name2 = "Second", power = 1,
                           log_x = TRUE, log_y = TRUE, log = NULL, ...) {
    ggplotly(plotSpectra2(object1, object2, name1 = name1, name2 = name2,
                          power = power, log_x = log_x, log_y = log_y,
                          log = log, ...),
             tooltip = c("Species", "w", "value", "Model"))
}

#' Plot the relative difference between two spectra
#'
#' `plotSpectraRelative()` plots the difference between the spectra relative to
#' their average. If we denote the number density from the first object as
#' \eqn{N_1(w)} and that from the second object as \eqn{N_2(w)}, then this plot
#' shows
#' \deqn{2 (N_2(w) - N_1(w)) / (N_2(w) + N_1(w)).}
#'
#' The individual spectra are calculated by [plotSpectra()], to which all
#' additional arguments are passed. For example, you can determine a time range
#' over which to average simulation results via `time_range`. See
#' [plotSpectra()] for more options.
#'
#' Note that it does not matter whether the relative difference is calculated
#' for number density, biomass density, or biomass density in log weight,
#' because the factors of \eqn{w} by which the densities differ cancel out in
#' the relative difference.
#'
#' @param object1 First `MizerParams` or `MizerSim` object.
#' @param object2 Second `MizerParams` or `MizerSim` object.
#' @param log_x If `TRUE` (default), use a log10 x-axis.
#' @param ylim A numeric vector of length two providing lower and upper limits
#'   for the relative difference (y) axis. Use `NA` to refer to the existing
#'   minimum or maximum.
#' @param ... Arguments passed to [plotSpectra()] for preparing the spectra
#'   data, for example `species`, `time_range`, `wlim`, `resource`,
#'   `background` or `total`.
#'
#' @return A ggplot2 object.
#' @export
#' @family plotting functions
#'
#' @examples
#' \donttest{
#' sim1 <- project(NS_params, t_max = 10, progress_bar = FALSE)
#' sim2 <- project(NS_params, effort = 0.5, t_max = 10, progress_bar = FALSE)
#' plotSpectraRelative(sim1, sim2)
#' }
plotSpectraRelative <- function(object1, object2, log_x = TRUE,
                                ylim = c(NA, NA), ...) {
    args <- list(...)
    wlim <- args$wlim %||% c(NA, NA)

    sf1 <- plotSpectra(object1, return_data = TRUE, ...)
    sf2 <- plotSpectra(object2, return_data = TRUE, ...)
    params <- if (is(object1, "MizerSim")) object1@params else object1

    plotRelativeDataFrame(sf1, sf2, validParams(params),
                          xlab = "Size [g]",
                          xtrans = if (log_x) "log10" else "identity",
                          xlim = wlim, ylim = ylim,
                          legend_var = "Legend")
}

#' @rdname plotSpectraRelative
#' @export
plotlySpectraRelative <- function(object1, object2, log_x = TRUE,
                                  ylim = c(NA, NA), ...) {
    ggplotly(plotSpectraRelative(object1, object2, log_x = log_x,
                                  ylim = ylim, ...),
             tooltip = c("Legend", "w", "rel_diff"))
}

#' @rdname plotCDF
#' @return `plotlyCDF()` returns a plotly object.
#' @export
plotlyCDF <- function(object, species = NULL,
                      time_range, geometric_mean = FALSE,
                      wlim = c(NA, NA), ylim = c(NA, NA),
                      power = 1, biomass = TRUE,
                      total = FALSE, resource = FALSE,
                      background = TRUE,
                      highlight = NULL, normalise = TRUE,
                      log_x = TRUE, log = NULL, ...) {
    args <- list(object = object, species = species,
                 geometric_mean = geometric_mean,
                 wlim = wlim, ylim = ylim,
                 biomass = biomass, total = total,
                 resource = resource, background = background,
                 highlight = highlight, normalise = normalise,
                 log_x = log_x, log = log, ...)
    if (!missing(time_range)) {
        args$time_range <- time_range
    }
    if (!missing(power)) {
        args$power <- power
    }
    ggplotly(do.call("plotCDF", args),
             tooltip = c("Species", "w", "value"))
}

#' @rdname plotCDF2
#' @return `plotlyCDF2()` returns a plotly object.
#' @export
plotlyCDF2 <- function(object1, object2, name1 = "First", name2 = "Second",
                       power = 1, normalise = TRUE,
                       log_x = TRUE, log = NULL, resource = FALSE, ...) {
    args <- list(object1 = object1, object2 = object2,
                 name1 = name1, name2 = name2,
                 normalise = normalise, log_x = log_x, log = log, ...)
    if (!missing(power)) {
        args$power <- power
    }
    ggplotly(do.call("plotCDF2", args),
             tooltip = c("Species", "w", "value", "Model"))
}

#' @rdname plotSpectra
#' @export
plotlySpectra <- function(object, species = NULL,
                        time_range, geometric_mean = FALSE,
                        wlim = c(NA, NA), ylim = c(NA, NA),
                        power = 1, biomass = TRUE,
                        total = FALSE, resource = TRUE,
                        background = TRUE,
                        highlight = NULL, log_x = TRUE, log_y = TRUE,
                        log = NULL, ...) {
    argg <- as.list(environment())
    ggplotly(do.call("plotSpectra", argg),
             tooltip = c("Species", "w", "value"))
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
#' @inheritParams plotSpectra
#' @param all.sizes If TRUE, then feeding level is plotted also for sizes
#'   outside a species' size range. Default FALSE.
#' @param include_critical If TRUE, then the critical feeding level is also
#'   plotted. Default FALSE.
#'
#' @return A ggplot2 object, unless `return_data = TRUE`, in which case a data
#'   frame with the variables 'w', 'value' and 'Species' is returned. If also
#'   `include_critical = TRUE` then the data frame contains a fourth variable
#'   'Type' that distinguishes between 'actual' and 'critical' feeding level.
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
plotFeedingLevel <- function(object, ...) {
    UseMethod("plotFeedingLevel")
}

#' @rdname plotFeedingLevel
#' @export
plotFeedingLevel.MizerSim <- function(object, species = NULL,
            time_range, highlight = NULL,
            all.sizes = FALSE, include_critical = FALSE,
            return_data = FALSE, ...) {
    assert_that(is.flag(all.sizes),
                is.flag(include_critical),
                is.flag(return_data))
    if (missing(time_range)) {
        time_range  <- max(as.numeric(dimnames(object@n)$time))
    }
    params <- validParams(object@params)
    feed <- getFeedingLevel(object, time_range = time_range, drop = FALSE)
    # If a time range was returned, average over it
    if (length(dim(feed)) == 3) {
        feed <- apply(feed, c(2, 3), mean)
    }
    plot_feeding_level(params, feed, species = species,
                       highlight = highlight, all.sizes = all.sizes,
                       include_critical = include_critical,
                       return_data = return_data)
}

#' @rdname plotFeedingLevel
#' @export
plotFeedingLevel.MizerParams <- function(object, species = NULL,
            highlight = NULL,
            all.sizes = FALSE, include_critical = FALSE,
            return_data = FALSE, ...) {
    assert_that(is.flag(all.sizes),
                is.flag(include_critical),
                is.flag(return_data))
    params <- validParams(object)
    feed <- getFeedingLevel(params, drop = FALSE)
    plot_feeding_level(params, feed, species = species,
                       highlight = highlight, all.sizes = all.sizes,
                       include_critical = include_critical,
                       return_data = return_data)
}

plot_feeding_level <- function(params, feed, species, highlight,
                               all.sizes, include_critical,
                               return_data) {

    # selector for desired species
    sel_sp <- valid_species_arg(params, species, return.logical = TRUE,
                                error_on_empty = TRUE)
    species <- dimnames(params@initial_n)$sp[sel_sp]
    feed <- feed[sel_sp, , drop = FALSE]

    plot_dat <- data.frame(w = rep(params@w, each = length(species)),
                           value = c(feed),
                           Species = species)

    if (include_critical) {
        feed_crit <- getCriticalFeedingLevel(params)[sel_sp, , drop = FALSE]
        plot_dat_crit <- data.frame(
            w = rep(params@w, each = length(species)),
            value = c(feed_crit),
            Species = species)
        plot_dat$Type <- "actual"
        plot_dat_crit$Type <- "critical"
        plot_dat <- rbind(plot_dat, plot_dat_crit)
    }

    if (!all.sizes) {
        # Remove feeding level for sizes outside a species' size range
        for (sp in species) {
            plot_dat$value[plot_dat$Species == sp &
                               (plot_dat$w < params@species_params[sp, "w_min"] |
                                    plot_dat$w > params@species_params[sp, "w_max"])] <- NA
        }
        plot_dat <- plot_dat[complete.cases(plot_dat), ]
    }

    if (return_data) return(plot_dat)

    # Need to keep species in order for legend
    legend_levels <-
        intersect(names(params@linecolour), plot_dat$Species)
    plot_dat$Legend <- factor(plot_dat$Species, levels = legend_levels)
    linesize <- make_linesize(legend_levels, highlight)

    # We do not use `plotDataFrame()` to create the plot because it would not
    # handle the alpha transparency for the critical feeding level.

    # The reason why below `group = species` is included in `ggplot()`
    # rather than in `geom_line` is because that puts it first in the
    # plotly tooltips, due to a bug in plotly.
    if (include_critical) {
        plot_dat$Species <- interaction(plot_dat$Species, plot_dat$Type)
        p <- ggplot(plot_dat, aes(group = Species,
                                  alpha = Type)) +
            scale_discrete_manual("alpha", name = "Feeding Level",
                                  values = c(actual = 1, critical = 0.5))
    } else {
        p <- ggplot(plot_dat, aes(group = Species))
    }
    p + geom_line(aes(x = w, y = value,
                      colour = Legend, linetype = Legend, linewidth = Legend)) +
        scale_x_continuous(name = "Size [g]", trans = "log10") +
        scale_y_continuous(name = "Feeding Level") +
        coord_cartesian(ylim = c(0, 1)) +
        scale_colour_manual(values = params@linecolour[legend_levels]) +
        scale_linetype_manual(values = params@linetype[legend_levels]) +
        scale_discrete_manual("linewidth", values = linesize)
}

#' @rdname plotFeedingLevel
#' @export
plotlyFeedingLevel <- function(object,
                             species = NULL,
                             time_range,
                             highlight = NULL,
                             include_critical = FALSE, ...) {
    argg <- as.list(environment())
    p <- ggplotly(do.call("plotFeedingLevel", argg),
                  tooltip = c("Species", "w", "value"))

    # When critical feeding level is included, ggplotly creates traces split by the
    # interaction of Species and Type, which produces a very long combined legend.
    # The code below reshapes the legend to mirror the ggplot output:
    # - Species appear once under a "Legend" group
    # - A separate "Feeding Level" group shows "actual" and "critical"
    if (isTRUE(include_critical)) {
        species_seen <- character(0)
        for (i in seq_along(p$x$data)) {
            tr <- p$x$data[[i]]
            # Only adjust line traces with names like "(actual, Cod)"
            if (!identical(tr$type, "scatter") ||
                is.null(tr$mode) || !grepl("lines", tr$mode)) next
            nm <- tr$name
            if (!is.null(nm) && grepl("^\\(", nm)) {
                nm_clean <- gsub("^\\(|\\)$", "", nm)
                parts <- strsplit(nm_clean, ",\\s*")[[1]]
                if (length(parts) >= 2) {
                    typ <- parts[[1]]
                    sp <- paste(parts[-1], collapse = ",")
                    # Rename the trace to species name and group under "Legend"
                    p$x$data[[i]]$name <- sp
                    p$x$data[[i]]$legendgroup <- "Legend"
                    # Add group title once
                    if (length(species_seen) == 0) {
                        p$x$data[[i]]$legendgrouptitle <- list(text = "Legend")
                    }
                    # Hide duplicate legend entries for "critical" traces
                    if (identical(typ, "critical")) {
                        p$x$data[[i]]$showlegend <- FALSE
                    }
                    species_seen <- unique(c(species_seen, sp))
                }
            }
        }

        # Add two legend-only traces to show the "Feeding Level" group
        p <- plotly::add_trace(
            p,
            x = c(0, 1), y = c(0, 1),
            type = "scatter", mode = "lines",
            name = "actual",
            legendgroup = "Feeding Level",
            legendgrouptitle = list(text = "Feeding Level"),
            line = list(color = "blue"),
            opacity = 1,
            hoverinfo = "skip",
            visible = "legendonly",
            showlegend = TRUE,
            inherit = FALSE
        )
        p <- plotly::add_trace(
            p,
            x = c(0, 1), y = c(0, 1),
            type = "scatter", mode = "lines",
            name = "critical",
            legendgroup = "Feeding Level",
            line = list(color = "blue"),
            opacity = 0.5,
            hoverinfo = "skip",
            visible = "legendonly",
            showlegend = TRUE,
            inherit = FALSE
        )
    }
    p
}

#' Plot predation mortality rate of each species against size
#'
#' After running a projection, plot the predation mortality rate of each species
#' by size. The mortality rate is averaged over the specified time range (a
#' single value for the time range can be used to plot a single time step).
#'
#' @inheritParams plotSpectra
#' @param all.sizes If TRUE, then predation mortality is plotted also for sizes
#'   outside a species' size range. Default FALSE.
#'
#' @return  A ggplot2 object, unless `return_data = TRUE`, in which case a data
#'   frame with the three variables 'w', 'value', 'Species' is returned.
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
plotPredMort <- function(object, ...) {
    UseMethod("plotPredMort")
}

#' @rdname plotPredMort
#' @export
plotPredMort.MizerSim <- function(object, species = NULL,
                         time_range, all.sizes = FALSE,
                         highlight = NULL, return_data = FALSE,
                         ...) {
    if (missing(time_range)) {
        time_range <- max(as.numeric(dimnames(object@n)$time))
    }
    pred_mort <- getPredMort(object, time_range = time_range, drop = FALSE)
    if (length(dim(pred_mort)) == 3) {
        pred_mort <- apply(pred_mort, c(2, 3), mean)
    }
    pred_mort <- ArraySpeciesBySize(pred_mort,
                                    value_name = "Predation mortality",
                                    units = "1/year",
                                    params = object@params)
    plot(pred_mort, species = species, all.sizes = all.sizes,
         highlight = highlight, return_data = return_data,
         ylim = c(0, NA))
}

#' @rdname plotPredMort
#' @export
plotPredMort.MizerParams <- function(object, species = NULL,
                         all.sizes = FALSE,
                         highlight = NULL, return_data = FALSE,
                         ...) {
    plot(getPredMort(validParams(object)), species = species,
         all.sizes = all.sizes, highlight = highlight,
         return_data = return_data, ylim = c(0, NA))
}

#' Alias for `plotPredMort()`
#'
#' `r lifecycle::badge("deprecated")`
#' An alias provided for backward compatibility with mizer version <= 1.0
#' @inherit plotPredMort
#' @export
#' @concept deprecated
plotM2 <- plotPredMort

#' @rdname plotPredMort
#' @export
plotlyPredMort <- function(object, species = NULL,
                           time_range,
                           highlight = NULL, ...) {
    argg <- as.list(environment())
    ggplotly(do.call("plotPredMort", argg),
             tooltip = c("Species", "w", "value"))
}

#' Plot total fishing mortality of each species by size
#'
#' After running a projection, plot the total fishing mortality of each species
#' by size. The total fishing mortality is averaged over the specified time
#' range (a single value for the time range can be used to plot a single time
#' step).
#'
#' @inheritParams plotSpectra
#' @param all.sizes If TRUE, then fishing mortality is plotted also for sizes
#'   outside a species' size range. Default FALSE.
#' @return A ggplot2 object, unless `return_data = TRUE`, in which case a data
#'   frame with the three variables 'w', 'value', 'Species' is returned.
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
plotFMort <- function(object, ...) {
    UseMethod("plotFMort")
}

#' @rdname plotFMort
#' @export
plotFMort.MizerSim <- function(object, species = NULL,
                      time_range, all.sizes = FALSE,
                      highlight = NULL, return_data = FALSE,
                      ...) {
    if (missing(time_range)) {
        time_range <- max(as.numeric(dimnames(object@n)$time))
    }
    f <- getFMort(object, time_range = time_range, drop = FALSE)
    if (length(dim(f)) == 3) {
        f <- apply(f, c(2, 3), mean)
    }
    f <- ArraySpeciesBySize(f, value_name = "Fishing mortality",
                            units = "1/year", params = object@params)
    plot(f, species = species, all.sizes = all.sizes,
         highlight = highlight, return_data = return_data)
}

#' @rdname plotFMort
#' @export
plotFMort.MizerParams <- function(object, species = NULL,
                      all.sizes = FALSE,
                      highlight = NULL, return_data = FALSE,
                      ...) {
    plot(getFMort(validParams(object)), species = species,
         all.sizes = all.sizes, highlight = highlight,
         return_data = return_data)
}

#' @rdname plotFMort
#' @export
plotlyFMort <- function(object, species = NULL,
                        time_range,
                        highlight = NULL, ...) {
    argg <- as.list(environment())
    ggplotly(do.call("plotFMort", argg),
             tooltip = c("Species", "w", "value"))
}


#' Plot growth curves
#'
#' The growth curves represent the average age of all the living fish of a
#' species as a function of their size. So it would be natural to plot size
#' on the x-axis. But to follow the usual convention from age-based models, we
#' plot size on the y-axis and age on the x-axis.
#'
#' In each panel for a single species, a horizontal line is included that
#' indicate the maturity size of the species and a vertical line indicating its
#' maturity age.
#'
#' If size at age data is passed via the `size_at_age` argument, this is plotted
#' on top of the growth curve. When comparing this to the growth curves, you
#' need to remember that the growth curves should only represent the average
#' age at each size. So a scatter in the x-direction around the curve is to be
#' expected.
#'
#' For legacy reasons, if the species parameters contain the variables `a` and
#' `b` for length to weight conversion and the von Bertalanffy parameter `k_vb`,
#' `w_inf` (and optionally `t0`), then the von Bertalanffy growth curve is
#' superimposed in black. This was implemented before we understood that the von
#' Bertalanffy curves (which approximates the average length at each age) should
#' not be compared to the mizer growth curves (which approximate the average age
#' at each length).
#'
#' @inheritParams getGrowthCurves
#' @inheritParams plotSpectra
#' @param species_panel If TRUE (default), and `percentage = FALSE`, display all
#'   species as facets. Otherwise puts all species into a single panel.
#' @param size_at_age A data frame with observed size at age data to be plotted
#'   on top of growth curve graphs. Should contain columns `species` (species
#'   name as used in the model), `age` (in years) and either `weight` (in grams)
#'   or `length` (in cm). If both `weight` and `length` are provided, only
#'   `weight` is used.
#' @return A ggplot2 object
#' @export
#' @family plotting functions
#' @seealso [plotting_functions]
#' @examples
#' \donttest{
#' params <-  NS_params
#' sim <- project(params, effort=1, t_max=20, t_save = 2, progress_bar = FALSE)
#' plotGrowthCurves(sim, percentage = TRUE)
#' plotGrowthCurves(sim, species = "Cod", max_age = 24)
#' plotGrowthCurves(sim, species_panel = TRUE)
#'
#' # Returning the data frame
#' fr <- plotGrowthCurves(sim, return_data = TRUE)
#' str(fr)
#' }
#' @rdname plotGrowthCurves
#' @export
plotGrowthCurves <- function(object, ...) {
    UseMethod("plotGrowthCurves")
}

#' @rdname plotGrowthCurves
#' @export
plotGrowthCurves.MizerSim <- function(object, species = NULL,
                             max_age = 20, percentage = FALSE,
                             species_panel = FALSE, highlight = NULL,
                             size_at_age = NULL,
                             return_data = FALSE, ...) {
    assert_that(is.flag(percentage),
                is.flag(species_panel),
                is.flag(return_data),
                is.number(max_age))
    params <- object@params
    params <- setInitialValues(params, object)
    plot_growth_curves(params, species = species,
                       max_age = max_age, percentage = percentage,
                       species_panel = species_panel, highlight = highlight,
                       size_at_age = size_at_age,
                       return_data = return_data)
}

#' @rdname plotGrowthCurves
#' @export
plotGrowthCurves.MizerParams <- function(object, species = NULL,
                             max_age = 20, percentage = FALSE,
                             species_panel = FALSE, highlight = NULL,
                             size_at_age = NULL,
                             return_data = FALSE, ...) {
    assert_that(is.flag(percentage),
                is.flag(species_panel),
                is.flag(return_data),
                is.number(max_age))
    params <- validParams(object)
    plot_growth_curves(params, species = species,
                       max_age = max_age, percentage = percentage,
                       species_panel = species_panel, highlight = highlight,
                       size_at_age = size_at_age,
                       return_data = return_data)
}

plot_growth_curves <- function(params, species,
                               max_age, percentage,
                               species_panel, highlight,
                               size_at_age,
                               return_data) {
    sp <- params@species_params
    sp <- set_species_param_default(sp, "age_mat", age_mat_vB(params))

    # size at age
    if (!is.null(size_at_age)) {
        if (!"species" %in% names(size_at_age)) {
            stop("The size at age data frame needs to have a 'species' column.")
        }
        if (!"age" %in% names(size_at_age)) {
            stop("The size at age data frame needs to have an 'age' column.")
        }
        if (!any(c("weight", "length") %in% names(size_at_age))) {
            stop("The size at age data frame needs to have either a 'length' or a 'weight' column.")
        }
        if (!"weight" %in% names(size_at_age)) {
            sp <- set_species_param_default(sp, "a", 0.004, message = "Using a = 0.004 for missing weight-length conversion parameters.")
            sp <- set_species_param_default(sp, "b", 3, message = "Using b = 3 for missing weight-length conversion parameters.")
            size_at_age <- left_join(size_at_age, select(sp, species, a, b),
                                     by = "species")
            size_at_age$weight <- size_at_age$a * size_at_age$length ^ size_at_age$b
        }
    }

    species <- valid_species_arg(params, species)
    # needed later to not confuse variable and column name
    selected_species <- species

    sp_sel <- sp$species %in% species
    ws <- getGrowthCurves(params, species, max_age, percentage)
    plot_dat <- reshape2::melt(ws)
    plot_dat$Species <- factor(plot_dat$Species, sp$species)
    plot_dat$Legend <- "model"

    # creating some VB
    if (all(c("a", "b", "k_vb", "w_inf") %in% names(sp))) {
        if ("t0" %in% names(sp)) {
            t0 <- sp$t0
        } else {
            t0 <- 0
        }
        VBdf <- data.frame("species" = sp$species,
                           "w_inf" = sp$w_inf,
                           "a" = sp$a,
                           "b" = sp$b,
                           "k_vb" = sp$k_vb,
                           "t0" = t0)
        VBdf$L_inf <- (VBdf$w_inf / VBdf$a) ^ (1 / VBdf$b)
        plot_dat2 <- plot_dat
        plot_dat2$value <-
            apply(plot_dat, 1,
                  function(x) {
                      sel <- VBdf$species == x[1]
                      length <- VBdf$L_inf[sel] *
                          (1 - exp(-VBdf$k_vb[sel] *
                                       (as.numeric(x[2]) - VBdf$t0[sel])))
                      VBdf$a[sel] * length ^ VBdf$b[sel]
                  })
        plot_dat2$Legend <- "von Bertalanffy"
        plot_dat <- rbind(plot_dat, plot_dat2)
    }
    if (return_data) return(plot_dat)

    p <- ggplot(filter(plot_dat, Legend == "model")) +
        geom_line(aes(x = Age, y = value,
                      colour = Species, linetype = Species,
                      linewidth = Species))
    y_label <- if (percentage)
        "Percent of maximum size"
    else "Size [g]"
    # Need to keep species in order for legend
    legend_levels <-
        intersect(c(dimnames(params@initial_n)$sp,
                    "Background", "Resource", "Total"),
                  plot_dat$Species)
    plot_dat$Species <- factor(plot_dat$Species, levels = legend_levels)
    linesize <- make_linesize(legend_levels, highlight)
    p <- p + scale_x_continuous(name = "Age [Years]") +
        scale_y_continuous(name = y_label) +
        scale_colour_manual(values = params@linecolour[legend_levels]) +
        scale_linetype_manual(values = params@linetype[legend_levels]) +
        scale_discrete_manual("linewidth", values = linesize)

    # starting cases now
    if (!percentage)  {
        if (length(species) == 1) {
            idx <- which(sp$species == species)
            w_mat <- sp$w_mat[idx]
            age_mat <- sp$age_mat[idx]
            p <- p + geom_hline(yintercept = w_mat, linetype = "dashed",
                                colour = "grey") +
                geom_vline(xintercept = age_mat, linetype = "dashed",
                           colour = "grey") +
                annotate("text", 0, w_mat, vjust = -1, label = "Maturity")
            if ("von Bertalanffy" %in% plot_dat$Legend) {
                p <- p + geom_line(data = filter(plot_dat,
                                                 Legend == "von Bertalanffy"),
                                   aes(x = Age, y = value))
            }
            if (!is.null(size_at_age)) {
                size_at_age <- filter(size_at_age, species == selected_species)
                p <- p + geom_point(aes(x = age, y = weight),
                                    data = size_at_age,
                                    alpha = 0.2)
            }

        } else if (species_panel) { # need to add either no panel if no param
            # for VB or create a panel without VB
            p <- ggplot(plot_dat) +
                geom_line(aes(x = Age, y = value, colour = Legend)) +
                scale_x_continuous(name = "Age [years]") +
                scale_y_continuous(name = "Size [g]") +
                geom_hline(aes(yintercept = w_mat),
                           data = tibble(Species = factor(legend_levels),
                                         w_mat = sp$w_mat[sp_sel]),
                           linetype = "dashed",
                           colour = "grey") +
                geom_vline(aes(xintercept = age_mat),
                           data = tibble(Species = factor(legend_levels),
                                         age_mat = sp$age_mat[sp_sel]),
                           linetype = "dashed",
                           colour = "grey") +
                geom_hline(aes(yintercept = w_max),
                           data = tibble(Species = factor(legend_levels),
                                         w_max = sp$w_max[sp_sel]),
                           linetype = "solid",
                           colour = "grey") +
                facet_wrap(~Species, scales = "free_y")

        }

    }
    return(p)

}

#' @rdname plotGrowthCurves
#' @export
plotlyGrowthCurves <- function(object, species = NULL,
                               max_age = 20,
                               percentage = FALSE,
                               species_panel = FALSE,
                               highlight = NULL) {
    argg <- as.list(environment())
    ggplotly(do.call("plotGrowthCurves", argg),
             tooltip = c("Species", "Age", "value"))
}


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
#' @inheritParams plotSpectra
#'
#' @return A ggplot2 object, unless `return_data = TRUE`, in which case a data
#'   frame with the four variables 'Predator', 'w', 'Proportion', 'Prey' is
#'   returned.
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
plotDiet <- function(object, ...) {
    UseMethod("plotDiet")
}

#' @rdname plotDiet
#' @export
plotDiet.MizerSim <- function(object, species = NULL,
                              return_data = FALSE, ...) {
    plotDiet(object@params, species = species, return_data = return_data, ...)
}

#' @rdname plotDiet
#' @export
plotDiet.MizerParams <- function(object, species = NULL, return_data = FALSE, ...) {
    assert_that(is.flag(return_data))
    params <- validParams(object)
    diet <- getDiet(params)
    plot_diet(params, n = params@initial_n, diet = diet, species = species,
              return_data = return_data)
}

plot_diet <- function(params, n, diet, species, return_data) {
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

    if (return_data) return(plot_dat)

    legend_levels <-
        intersect(names(params@linecolour), plot_dat$Prey)
    p <- ggplot(plot_dat) +
        geom_area(aes(x = w, y = Proportion, fill = Prey)) +
        scale_x_log10() +
        labs(x = "Size [g]", y = "Proportion") +
        scale_fill_manual(values = params@linecolour[legend_levels],
                          limits = legend_levels)
    if (sum(species) > 1) {
        p <- p + facet_wrap(vars(Predator))
    }
    p
}

#' @rdname plotDiet
#' @return `plotlyDiet()` returns a plotly object.
#' @export
plotlyDiet <- function(object, species = NULL, ...) {
    ggplotly(plotDiet(object, species = species, ...),
             tooltip = c("Predator", "w", "Proportion", "Prey"))
}


#### plot ####
#' Summary plot for `MizerSim` objects
#'
#' After running a projection, produces 5 plots in the same window: feeding
#' level, abundance spectra, predation mortality and fishing mortality of each
#' species by size; and biomass of each species through time. This method just
#' uses the other plotting functions and puts them all in one window.
#'
#' @param x An object of class \linkS4class{MizerSim}
#' @param ...  For additional arguments see the documentation for
#'   [plotBiomass()],
#'   [plotFeedingLevel()],[plotSpectra()],[plotPredMort()]
#'   and [plotFMort()].
#' @return A viewport object
#' @export
#' @family plotting functions
#' @seealso [plotting_functions]
#' @rdname plotMizerSim
#' @name plotMizerSim
#' @examples
#' \donttest{
#' params <-  NS_params
#' sim <- project(params, effort=1, t_max=20, t_save = 2, progress_bar = FALSE)
#' plot(sim)
#' }
plot.MizerSim <- function(x, ...) {
    p1 <- plotFeedingLevel(x, ...)
    p2 <- plotSpectra(x, ...)
    p3 <- plotBiomass(x, y_ticks = 3, ...)
    p4 <- plotPredMort(x, ...)
    p5 <- plotFMort(x, ...)
    grid::grid.newpage()
    glayout <- grid::grid.layout(3, 2) # widths and heights arguments
    vp <- grid::viewport(layout = glayout)
    grid::pushViewport(vp)
    vplayout <- function(x, y) {
        grid::viewport(layout.pos.row = x, layout.pos.col = y)
    }
    print(p1 + theme(legend.position = "none"), vp = vplayout(1, 1))
    print(p3 + theme(legend.position = "none"), vp = vplayout(1, 2))
    print(p4 + theme(legend.position = "none"), vp = vplayout(2, 1))
    print(p5 + theme(legend.position = "none"), vp = vplayout(2, 2))
    print(p2 + theme(legend.position = "right",
                     legend.key.size = unit(0.1, "cm")),
          vp = vplayout(3, 1:2))
}

#' Summary plot for `MizerParams` objects
#'
#' Produces 3 plots in the same window: abundance spectra, feeding
#' level and predation mortality of each species through time. This method just
#' uses the other plotting functions and puts them all in one window.
#'
#' @param x An object of class \linkS4class{MizerParams}
#' @param ...  For additional arguments see the documentation for
#'   [plotFeedingLevel()],[plotSpectra()],[plotPredMort()]
#' @return A viewport object
#' @export
#' @family plotting functions
#' @seealso [plotting_functions]
#' @name plotMizerParams
#' @examples
#' \donttest{
#' params <-  NS_params
#' plot(params)
#' }
plot.MizerParams <- function(x, ...) {
    params <- validParams(x)
    p11 <- plotFeedingLevel(params, ...)
    p2 <- plotSpectra(params, ...)
    p12 <- plotPredMort(params, ...)
    grid::grid.newpage()
    glayout <- grid::grid.layout(2, 2) # widths and heights arguments
    vp <- grid::viewport(layout = glayout)
    grid::pushViewport(vp)
    vplayout <- function(x, y) {
        grid::viewport(layout.pos.row = x, layout.pos.col = y)
    }
    print(p11 + theme(legend.position = "none"), vp = vplayout(1, 1))
    print(p12 + theme(legend.position = "none"), vp = vplayout(1, 2))
    print(p2 + theme(legend.position = "right",
                     legend.key.size = unit(0.1, "cm")),
          vp = vplayout(2, 1:2))
}
