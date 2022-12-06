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
#' Every plotting function exists in two versions, `plotSomething` and 
#' `plotlySomething`. The plotly version is more interactive but not
#' suitable for inclusion in documents.
#'
#' This table shows the available plotting functions.
#' \tabular{ll}{
#'   Plot \tab Description \cr
#'   [plotBiomass()] \tab Plots the total biomass of each species through time. A time range to be plotted can be specified. The size range of the community can be specified in the same way as for [getBiomass()]. \cr
#'   [plotSpectra()] \tab Plots the abundance (biomass or numbers) spectra of each species and the background community. It is possible to specify a minimum size which is useful for truncating the plot. \cr
#'   [plotFeedingLevel()] \tab Plots the feeding level of each species against size. \cr
#'   [plotPredMort()] \tab Plots the predation mortality of each species against size. \cr
#'   [plotFMort()] \tab Plots the total fishing mortality of each species against size. \cr
#'   [plotYield()] \tab Plots the total yield of each species across all fishing gears against time. \cr
#'   [plotYieldGear()] \tab Plots the total yield of each species by gear against time. \cr
#'   [plotDiet()] \tab Plots the diet composition at size for a given predator species. \cr
#'   [plotGrowthCurves()] \tab Plots the size as a function of age. \cr
#'   [plot()] \tab Produces 5 plots ([plotFeedingLevel()], [plotBiomass()], [plotPredMort()], [plotFMort()] and [plotSpectra()]) in the same window. \cr
#' }
#' 
#' These functions use the ggplot2 package and return the plot as a ggplot
#' object. This means that you can manipulate the plot further after its 
#' creation using the ggplot grammar of graphics. The corresponding function
#' names with `plot` replaced by `plotly` produce interactive plots
#' with the help of the plotly package.
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
#' # Some example plots
#' plotFeedingLevel(sim)
#' 
#' # Plotting only a subset of species
#' plotFeedingLevel(sim, species = c("Cod", "Herring"))
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
                         "Proportion", "Prey", "Legend", "Type", "Gear"))

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
#' @param style The style of the plot. Availalble options are "line' for geom_line
#' and "area" for geom_area. Default is "line".
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
#' @param y_ticks The approximate number of ticks desired on the y axis
#' @param highlight Name or vector of names of the species to be highlighted.
#' @keywords internal
#' @export
plotDataFrame <- function(frame, params, style = "line", xlab = waiver(),
                          ylab = waiver(), xtrans = "identity", ytrans = "identity",
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
    linesize <- rep_len(0.8, length(legend_levels))
    names(linesize) <- legend_levels
    linesize[highlight] <- 1.6
    
    xbreaks <- waiver()
    if (xtrans == "log10") xbreaks <- log_breaks()
    ybreaks <- waiver()
    if (ytrans == "log10") ybreaks <- log_breaks(n = y_ticks)
    
    # The reason why below `group = species` is included in `ggplot()`
    # rather than in `geom_line` is because that puts it first in the
    # plotly tooltips, due to a bug in plotly.
    p <- ggplot(frame, aes(group = .data[[group_var]])) +
        scale_y_continuous(trans = ytrans, breaks = ybreaks,
                           labels = prettyNum, name = ylab) +
        scale_x_continuous(trans = xtrans, name = xlab)
    
    switch(style,
           "line" = {p <- p +
               geom_line(aes(x = .data[[x_var]], y = .data[[y_var]],
                             colour = .data[[legend_var]],
                             linetype = .data[[legend_var]],
                             linewidth = .data[[legend_var]])) +
               scale_colour_manual(values = linecolour) +
               scale_linetype_manual(values = linetype) +
               scale_discrete_manual("linewidth", values = linesize)
           },
           "area" = {p <- p +
               geom_area(aes(x = .data[[x_var]], y = .data[[y_var]],
                             fill = .data[[legend_var]])) +
               scale_fill_manual(values = linecolour)
           },
           {"unknown style selected"}
    )
    
    if (!is.null(wrap_var)) {
        if (!(wrap_var %in% var_names)) {
            stop("The `wrap_var` argument must be the name of a variable ",
                 "in the data frame.")
        }
        p <- p + facet_wrap(wrap_var, scales = wrap_scale)
    }
    
    p
}

#' Helper function to produce nice breaks on logarithmic axes
#'
#' This is needed when the logarithmic y-axis spans less than one order of
#' magnitude, in which case the ggplot2 default produces no ticks.
#' 
#' Thanks to Heather Turner at
#' https://stackoverflow.com/questions/14255533/pretty-ticks-for-log-normal-scale-using-ggplot2-dynamic-not-manual
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
#' (min_w, max_w, min_l, max_l, see [getBiomass()]). 
#' 
#' @param sim An object of class \linkS4class{MizerSim}
#' @inheritParams valid_species_arg
#' @inheritParams plotDataFrame
#' @param start_time The first time to be plotted. Default is the beginning
#'   of the time series.
#' @param end_time The last time to be plotted. Default is the end of the
#'   time series.
#' @param ylim A numeric vector of length two providing lower and upper limits
#'   for the y axis. Use NA to refer to the existing minimum or maximum. Any
#'   values below 1e-20 are always cut off.
#' @param total A boolean value that determines whether the total biomass from
#'   all species is plotted as well. Default is FALSE.
#' @inheritParams plotSpectra
#' @inheritDotParams get_size_range_array -params
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
plotBiomass <- function(sim, species = NULL, 
                        start_time, end_time, 
                        y_ticks = 6, ylim = c(NA, NA), 
                        total = FALSE, background = TRUE, 
                        highlight = NULL, return_data = FALSE,
                        ...) {
    assert_that(is(sim, "MizerSim"),
                is.flag(total),
                is.flag(background),
                is.flag(return_data),
                length(ylim) == 2)
    params <- sim@params
    species <- valid_species_arg(sim, species, error_on_empty = TRUE)
    if (missing(start_time)) start_time <- 
            as.numeric(dimnames(sim@n)[[1]][1])
    if (missing(end_time)) end_time <- 
            as.numeric(dimnames(sim@n)[[1]][dim(sim@n)[1]])
    if (start_time >= end_time) {
        stop("start_time must be less than end_time")
    }
    # First we get the data frame for all species, including the background,
    # for all times but only the desired size range, by passing any size range
    # arguments on to getBiomass()
    bm <- getBiomass(sim, ...)
    # Select time range
    bm <- bm[(as.numeric(dimnames(bm)[[1]]) >= start_time) &
               (as.numeric(dimnames(bm)[[1]]) <= end_time), , drop = FALSE]

    # Include total
    if (total) {
        bm <- cbind(bm, Total = rowSums(bm))
    }
    
    bm <- reshape2::melt(bm)
    
    # Implement ylim and a minimal cutoff and bring columns in desired order
    min_value <- 1e-20
    bm <- bm[bm$value >= min_value &
                 (is.na(ylim[1]) | bm$value >= ylim[1]) &
                 (is.na(ylim[2]) | bm$value <= ylim[2]), c(1, 3, 2)]
    names(bm) <- c("Year", "Biomass", "Species")
    
    # Select species
    plot_dat <- bm[bm$Species %in% c("Total", species), ]
    plot_dat$Legend <- plot_dat$Species
    
    if (background) {
        # Add background species in light grey
        bkgrd_sp <- dimnames(sim@n)$sp[is.na(sim@params@A)]
        if (length(bkgrd_sp) > 0) {
            bm_bkgrd <- bm[bm$Species %in% bkgrd_sp, ]
            bm_bkgrd$Legend <- "Background"
            plot_dat <- rbind(plot_dat, bm_bkgrd)
        }
    }
    
    if (return_data) return(plot_dat) 
    
    plotDataFrame(plot_dat, params, xlab = "Year", ylab = "Biomass [g]",
                  ytrans = "log10", 
                  y_ticks = y_ticks, highlight = highlight,
                  legend_var = "Legend")
}

#' @rdname plotBiomass
#' @export
plotlyBiomass <- function(sim,
             species = NULL,
             start_time,
             end_time,
             y_ticks = 6,
             ylim = c(NA, NA),
             total = FALSE,
             background = TRUE,
             highlight = NULL,
             ...) {
    argg <- c(as.list(environment()), list(...))
    ggplotly(do.call("plotBiomass", argg),
             tooltip = c("Species", "Year", "Biomass"))
}


#' Plot the total yield of species through time
#'
#' After running a projection, the total yield of each species across all 
#' fishing gears can be plotted against time. The yield is obtained with
#' [getYield()].
#' 
#' @param sim An object of class \linkS4class{MizerSim}
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
plotYield <- function(sim, sim2,
                      species = NULL,
                      total = FALSE, log = TRUE,
                      highlight = NULL, return_data = FALSE,
                      ...) {
    assert_that(is(sim, "MizerSim"),
                is.flag(total),
                is.flag(log),
                is.flag(return_data))
    params <- sim@params
    species <- valid_species_arg(sim, species, error_on_empty = TRUE)
    if (missing(sim2)) {
        y <- getYield(sim, ...)
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
        if (!all(dimnames(sim@n)$time == dimnames(sim2@n)$time)) {
            stop("The two simulations do not have the same times")
        }
        ym <- plotYield(sim, species = species,
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
plotlyYield <- function(sim, sim2,
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
#' @param sim An object of class \linkS4class{MizerSim}
#' @inheritParams plotSpectra
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
plotYieldGear <- function(sim,
                          species = NULL,
                          total = FALSE,
                          highlight = NULL, return_data = FALSE,
                          ...) {
    assert_that(is(sim, "MizerSim"),
                is.flag(total),
                is.flag(return_data))
    params <- sim@params
    species <- valid_species_arg(sim, species, error_on_empty = TRUE)
    
    y <- getYieldGear(sim, ...)
    y_total <- rowSums(y, dims = 2)
    y <- y[, , dimnames(y)$sp %in% species, drop = FALSE]
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
    linesize <- rep(0.8, length(species_levels))
    names(linesize) <- names(params@linetype[species_levels])
    linesize[highlight] <- 1.6
    ggplot(ym) +
        geom_line(aes(x = Year, y = Yield, colour = Species,
                      linetype = Gear, linewidth = Species)) +
        scale_y_continuous(trans = "log10", name = "Yield [g]") +
        scale_colour_manual(values = params@linecolour[species_levels]) +
        scale_discrete_manual("linewidth", values = linesize)
}

#' @rdname plotYieldGear
#' @export
plotlyYieldGear <- function(sim, species = NULL,
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
#' @param wlim A numeric vector of length two providing lower and upper limits
#'   for the w axis. Use NA to refer to the existing minimum or maximum.
#' @param ylim A numeric vector of length two providing lower and upper limits
#'   for the y axis. Use NA to refer to the existing minimum or maximum. Any
#'   values below 1e-20 are always cut off.
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
plotSpectra <- function(object, species = NULL,
                        time_range,
                        wlim = c(NA, NA), ylim = c(NA, NA),
                        power = 1, biomass = TRUE,
                        total = FALSE, resource = TRUE, 
                        background = TRUE,
                        highlight = NULL, return_data = FALSE, ...) {
    # to deal with old-type biomass argument
    if (missing(power)) {
        power <- as.numeric(biomass)
    }
    assert_that(is.flag(total), is.flag(resource),
                is.flag(background),
                is.number(power), 
                length(wlim) == 2,
                length(ylim) == 2)
    species <- valid_species_arg(object, species)
    if (length(species) == 0 && !total && !resource) {
        stop("There is nothing to plot as no valid species have been selected.")
    }
    if (is(object, "MizerSim")) {
        if (missing(time_range)) {
            time_range  <- max(as.numeric(dimnames(object@n)$time))
        }
        time_elements <- get_time_elements(object, time_range)
        n <- apply(object@n[time_elements, , , drop = FALSE], c(2, 3), mean)
        n_pp <- apply(object@n_pp[time_elements, , drop = FALSE], 2, mean)
        ps <- plot_spectra(object@params, n = n, n_pp = n_pp,
                           species = species, wlim = wlim, ylim = ylim,
                           power = power, total = total, resource = resource,
                           background = background, highlight = highlight, 
                           return_data = return_data)
        return(ps)
    } else if (is(object, "MizerParams")) {
        ps <- plot_spectra(object, n = object@initial_n,
                           n_pp = object@initial_n_pp,
                           species = species, wlim = wlim, ylim = ylim,
                           power = power, total = total, resource = resource,
                           background = background, highlight = highlight, 
                           return_data = return_data)
        return(ps)
    } else {
        stop("First argument of `plotSpectra()` needs to be a MizerSim or ",
             "a MizerParams object.")
    }
}


plot_spectra <- function(params, n, n_pp,
                         species, wlim, ylim, power,
                         total, resource, background, 
                         highlight, return_data) {
    params <- validParams(params)
    if (is.na(wlim[1])) {
        wlim[1] <- min(params@w) / 100
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
    if (power %in% c(0, 1, 2)) {
        y_label <- c("Number density [1/g]", "Biomass density",
                    "Biomass density [g]")[power + 1]
    } else {
        y_label <- paste0("Number density * w^", power)
    }
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
    if (background && anyNA(params@A)) {
        back_n <- n[is.na(params@A), , drop = FALSE]
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
    if (is.na(ylim[1])) {
        ylim[1] <- 1e-20
    }
    plot_dat <- plot_dat[plot_dat$value > ylim[1], ]
    
    if (return_data) return(plot_dat) 
    
    plotDataFrame(plot_dat, params, xlab = "Size [g]", ylab = y_label,
                  xtrans = "log10", ytrans = "log10", 
                  highlight = highlight, legend_var = "Legend")
}

#' @rdname plotSpectra
#' @export
plotlySpectra <- function(object, species = NULL,
                        time_range,
                        wlim = c(NA, NA), ylim = c(NA, NA),
                        power = 1, biomass = TRUE,
                        total = FALSE, resource = TRUE, 
                        background = TRUE,
                        highlight = NULL, ...) {
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
plotFeedingLevel <- function(object, species = NULL,
            time_range, highlight = NULL,
            all.sizes = FALSE, include_critical = FALSE,
            return_data = FALSE, ...) {
    assert_that(is.flag(all.sizes),
                is.flag(include_critical),
                is.flag(return_data))
    if (is(object, "MizerSim")) {
        if (missing(time_range)) {
            time_range  <- max(as.numeric(dimnames(object@n)$time))
        }
        params <- validParams(object@params)
        feed <- getFeedingLevel(object, time_range = time_range, drop = FALSE)
    } else if (is(object, "MizerParams")) {
        params <- validParams(object)
        feed <- getFeedingLevel(params, drop = FALSE)
    }
    # If a time range was returned, average over it
    if (length(dim(feed)) == 3) {
        feed <- apply(feed, c(2, 3), mean)
    }
        
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
                                    plot_dat$w > params@species_params[sp, "w_inf"])] <- NA
        }
        plot_dat <- plot_dat[complete.cases(plot_dat), ]
    }
    
    if (return_data) return(plot_dat)
    
    # Need to keep species in order for legend
    legend_levels <- 
        intersect(names(params@linecolour), plot_dat$Species)
    plot_dat$Legend <- factor(plot_dat$Species, levels = legend_levels)
    linesize <- rep(0.8, length(legend_levels))
    names(linesize) <- names(params@linetype[legend_levels])
    linesize[highlight] <- 1.6
    
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
        scale_y_continuous(name = "Feeding Level", limits = c(0, 1)) +
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
                             include_critical, ...) {
    argg <- as.list(environment())
    ggplotly(do.call("plotFeedingLevel", argg),
             tooltip = c("Species", "w", "value"))
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
plotPredMort <- function(object, species = NULL,
                         time_range, all.sizes = FALSE,
                         highlight = NULL, return_data = FALSE,
                         ...) {
    assert_that(is.flag(all.sizes),
                is.flag(return_data))
    if (is(object, "MizerSim")) {
        if (missing(time_range)) {
            time_range  <- max(as.numeric(dimnames(object@n)$time))
        }
        params <- object@params
    } else {
        params <- validParams(object)
    }
    pred_mort <- getPredMort(object, time_range = time_range, drop = FALSE)
    # If a time range was returned, average over it
    if (length(dim(pred_mort)) == 3) {
        pred_mort <- apply(pred_mort, c(2, 3), mean)
    }
    
    species <- valid_species_arg(params, species, error_on_empty = TRUE)
    pred_mort <- pred_mort[as.character(dimnames(pred_mort)[[1]]) %in% species, , drop = FALSE]
    plot_dat <- data.frame(
        w = rep(params@w, each = length(species)),
        value = c(pred_mort),
        Species = species)
    
    if (!all.sizes) {
        # Remove feeding level for sizes outside a species' size range
        for (sp in species) {
            plot_dat$value[plot_dat$Species == sp &
                               (plot_dat$w < params@species_params[sp, "w_min"] |
                                    plot_dat$w > params@species_params[sp, "w_inf"])] <- NA
        }
        plot_dat <- plot_dat[complete.cases(plot_dat), ]
    }

    if (return_data) return(plot_dat)
    
    p <- plotDataFrame(plot_dat, params, xlab = "Size [g]", xtrans = "log10",
                       highlight = highlight)
    suppressMessages(
        p <- p + scale_y_continuous(labels = prettyNum, 
                                    name = "Predation mortality [1/year]",
                                    limits = c(0, max(plot_dat$value))))
    p
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
plotFMort <- function(object, species = NULL,
                      time_range, all.sizes = FALSE,
                      highlight = NULL, return_data = FALSE,
                      ...) {
    assert_that(is.flag(all.sizes),
                is.flag(return_data))
    if (is(object, "MizerSim")) {
        if (missing(time_range)) {
            time_range  <- max(as.numeric(dimnames(object@n)$time))
        }
        params <- object@params
    } else {
        params <- validParams(object)
    }
    f <- getFMort(object, time_range = time_range, drop = FALSE)
    # If a time range was returned, average over it
    if (length(dim(f)) == 3) {
        f <- apply(f, c(2, 3), mean)
    }
    species <- valid_species_arg(params, species, error_on_empty = TRUE)
    f <- f[as.character(dimnames(f)[[1]]) %in% species, , drop = FALSE]
    plot_dat <- data.frame(w = rep(params@w, each = length(species)),
                           value = c(f),
                           Species = species)
    
    if (!all.sizes) {
        # Remove feeding level for sizes outside a species' size range
        for (sp in species) {
            plot_dat$value[plot_dat$Species == sp &
                               (plot_dat$w < params@species_params[sp, "w_min"] |
                                    plot_dat$w > params@species_params[sp, "w_inf"])] <- NA
        }
        plot_dat <- plot_dat[complete.cases(plot_dat), ]
    }
    
    if (return_data) return(plot_dat)
    
    plotDataFrame(plot_dat, params, xlab = "Size [g]", xtrans = "log10",
                  ylab = "Fishing mortality [1/Year]", highlight = highlight)
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


#' Plot growth curves giving weight as a function of age
#' 
#' When the growth curve for only a single species is plotted, horizontal
#' lines are included that indicate the maturity size and the maximum size for 
#' that species. If furthermore the species parameters contain the variables
#' a and b for length to weight conversion and the von Bertalanffy parameter
#' k_vb (and optionally t0), then the von Bertalanffy growth curve is
#' superimposed in black.
#' 
#' @inheritParams getGrowthCurves
#' @inheritParams plotSpectra
#' @param species_panel `r lifecycle::badge("experimental")`
#'   If TRUE, display all species with their Von Bertalanffy curves as facets 
#'   (need species and percentage to be set to default). Default FALSE.
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
plotGrowthCurves <- function(object, species = NULL, 
                             max_age = 20, percentage = FALSE, 
                             species_panel = FALSE, highlight = NULL,
                             return_data = FALSE, ...) {
    assert_that(is.flag(percentage),
                is.flag(species_panel),
                is.flag(return_data),
                is.number(max_age))
    if (is(object, "MizerSim")) {
        params <- object@params
        params <- setInitialValues(params, object)
    } else if (is(object, "MizerParams")) {
        params <- validParams(object)
    }
    species <- valid_species_arg(params, species, error_on_empty = TRUE)
    sp_sel <- params@species_params$species %in% species
    ws <- getGrowthCurves(params, species, max_age, percentage)
    plot_dat <- reshape2::melt(ws)
    plot_dat$Species <- factor(plot_dat$Species, params@species_params$species)
    plot_dat$Legend <- "model"
    
    # creating some VB
    if (all(c("a", "b", "k_vb") %in% names(params@species_params))) {
        if ("t0" %in% names(params@species_params)) {
            t0 <- params@species_params$t0
        } else {
            t0 <- 0
        }
        VBdf <- data.frame("species" = params@species_params$species, 
                           "w_inf" = params@species_params$w_inf, 
                           "a" = params@species_params$a, 
                           "b" = params@species_params$b, 
                           "k_vb" = params@species_params$k_vb, 
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
                      colour = Species, linetype = Species, linewidth = Species))
    y_label <- if (percentage) 
        "Percent of maximum size"
    else "Size [g]"
    # Need to keep species in order for legend
    legend_levels <- 
        intersect(c(dimnames(params@initial_n)$sp,
                    "Background", "Resource", "Total"),
                  plot_dat$Species)
    plot_dat$Species <- factor(plot_dat$Species, levels = legend_levels)
    linesize <- rep(0.8, length(legend_levels))
    names(linesize) <- names(params@linetype[legend_levels])
    linesize[highlight] <- 1.6
    p <- p + scale_x_continuous(name = "Age [Years]") + 
        scale_y_continuous(name = y_label) + 
        scale_colour_manual(values = params@linecolour[legend_levels]) + 
        scale_linetype_manual(values = params@linetype[legend_levels]) + 
        scale_discrete_manual("linewidth", values = linesize)
    
    # starting cases now
    if (!percentage)  {
        if (length(species) == 1) {
            idx <- which(params@species_params$species == species)
            w_inf <- params@species_params$w_inf[idx]
            p <- p + geom_hline(yintercept = w_inf, colour = "grey") + 
                annotate("text", 0, w_inf, vjust = -1, label = "Maximum")
            w_mat <- params@species_params$w_mat[idx]
            p <- p + geom_hline(yintercept = w_mat, linetype = "dashed", 
                                colour = "grey") + 
                annotate("text", 0, w_mat, vjust = -1, label = "Maturity")
            if ("von Bertalanffy" %in% plot_dat$Legend) 
                p <- p + geom_line(data = filter(plot_dat, Legend == "von Bertalanffy"), 
                                   aes(x = Age, y = value))
            
        } else if (species_panel) { # need to add either no panel if no param 
                                    # for VB or create a panel without VB
            p <- ggplot(plot_dat) +
                geom_line(aes(x = Age, y = value , colour = Legend)) +
                scale_x_continuous(name = "Age [years]") +
                scale_y_continuous(name = "Size [g]") +
                geom_hline(aes(yintercept = w_mat),
                           data = tibble(Species = factor(legend_levels),
                                         w_mat = params@species_params$w_mat[sp_sel]),
                           linetype = "dashed",
                           colour = "grey") +
                geom_hline(aes(yintercept = w_inf),
                           data = tibble(Species = factor(legend_levels),
                                         w_inf = params@species_params$w_inf[sp_sel]),
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
#' in the plot.
#'
#' If more than one predator species is selected, then the plot contains one
#' facet for each species.
#'
#' @inheritParams plotSpectra
#'
#' @return A ggplot2 object, unless `return_data = TRUE`, in which case a data
#'   frame with the four variables 'Predator', 'w', 'Proportion', 'Prey' is returned.
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
plotDiet <- function(object, species = NULL, return_data = FALSE) {
    assert_that(is.flag(return_data))
    params <- validParams(object)
    species <- valid_species_arg(object, species, return.logical = TRUE)
    diet <- getDiet(params)[species, , , drop = FALSE]
    names(dimnames(diet)) <- c("Predator", "w", "Prey")
    plot_dat <- melt(diet, value.name = "Proportion")
    prey <- dimnames(diet)$Prey
    # the plot looks better upsided down
    plot_dat$Prey <- factor(plot_dat$Prey, levels = rev(prey))
    
    plot_dat <- plot_dat[plot_dat$Proportion > 0.001, ]
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


#### plot ####
#' Summary plot for `MizerSim` objects
#' 
#' After running a projection, produces 5 plots in the same window: feeding
#' level, abundance spectra, predation mortality and fishing mortality of each
#' species by size; and biomass of each species through time. This method just
#' uses the other plotting functions and puts them all in one window.
#' 
#' @param x An object of class \linkS4class{MizerSim}
#' @param y Not used
#' @param ...  For additional arguments see the documentation for
#'   [plotBiomass()],
#'   [plotFeedingLevel()],[plotSpectra()],[plotPredMort()]
#'   and [plotFMort()].
#' @return A viewport object
#' @export
#' @family plotting functions
#' @seealso [plotting_functions]
#' @rdname plotMizerSim
#' @examples
#' \donttest{
#' params <-  NS_params
#' sim <- project(params, effort=1, t_max=20, t_save = 2, progress_bar = FALSE)
#' plot(sim)
#' plot(sim, time_range = 10:20) # change time period for size-based plots
#' plot(sim, min_w = 10, max_w = 1000) # change size range for biomass plot
#' }
setMethod("plot", signature(x = "MizerSim", y = "missing"),
          function(x, ...) {
              p1 <- plotFeedingLevel(x, ...)
              p2 <- plotSpectra(x, ...)
              p3 <- plotBiomass(x, y_ticks = 3, ...)
              p4 <- plotPredMort(x, ...)
              p5 <- plotFMort(x, ...)
              grid::grid.newpage()
              glayout <- grid::grid.layout(3, 2) # widths and heights arguments
              vp <- grid::viewport(layout = glayout)
              grid::pushViewport(vp)
              vplayout <- function(x, y)
                  grid::viewport(layout.pos.row = x, layout.pos.col = y)
              print(p1 + theme(legend.position = "none"), vp = vplayout(1, 1))
              print(p3 + theme(legend.position = "none"), vp = vplayout(1, 2))
              print(p4 + theme(legend.position = "none"), vp = vplayout(2, 1))
              print(p5 + theme(legend.position = "none"), vp = vplayout(2, 2))
              print(p2 + theme(legend.position = "right",
                               legend.key.size = unit(0.1, "cm")),
                    vp = vplayout(3, 1:2))
          }
)

#' Summary plot for `MizerParams` objects
#' 
#' Produces 3 plots in the same window: abundance spectra, feeding
#' level and predation mortality of each species through time. This method just
#' uses the other plotting functions and puts them all in one window.
#' 
#' @return A viewport object
#' @export
#' @family plotting functions
#' @seealso [plotting_functions]
#' @rdname plotMizerSim
#' @examples
#' \donttest{
#' params <-  NS_params
#' plot(params)
#' plot(params, min_w = 10, max_w = 1000) # change size range for biomass plot
#' }
setMethod("plot", signature(x = "MizerParams", y = "missing"),
          function(x, ...) {
              params <- validParams(x)
              p11 <- plotFeedingLevel(params, ...)
              p2 <- plotSpectra(params, ...)
              p12 <- plotPredMort(params, ...)
              grid::grid.newpage()
              glayout <- grid::grid.layout(2, 2) # widths and heights arguments
              vp <- grid::viewport(layout = glayout)
              grid::pushViewport(vp)
              vplayout <- function(x, y)
                  grid::viewport(layout.pos.row = x, layout.pos.col = y)
              print(p11 + theme(legend.position = "none"), vp = vplayout(1, 1))
              print(p12 + theme(legend.position = "none"), vp = vplayout(1, 2))
              print(p2 + theme(legend.position = "right",
                               legend.key.size = unit(0.1, "cm")),
                    vp = vplayout(2, 1:2))
          }
)
