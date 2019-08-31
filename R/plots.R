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
#' Every plotting function exists in two versions, \code{plotSomething} and 
#' \code{plotlySomething}. The plotly version is more interactive but not
#' suitable for inclusion in documents.
#'
#' This table shows the available plotting functions.
#' \tabular{ll}{
#'   Plot \tab Description \cr
#'   \code{\link{plotBiomass}} \tab Plots the total biomass of each species through time. A time range to be plotted can be specified. The size range of the community can be specified in the same way as for \code{\link{getBiomass}}. \cr
#'   \code{\link{plotSpectra}} \tab Plots the abundance (biomass or numbers) spectra of each species and the background community. It is possible to specify a minimum size which is useful for truncating the plot. \cr
#'   \code{\link{plotFeedingLevel}} \tab Plots the feeding level of each species against size. \cr
#'   \code{\link{plotPredMort}} \tab Plots the predation mortality of each species against size. \cr
#'   \code{\link{plotFMort}} \tab Plots the total fishing mortality of each species against size. \cr
#'   \code{\link{plotYield}} \tab Plots the total yield of each species across all fishing gears against time. \cr
#'   \code{\link{plotYieldGear}} \tab Plots the total yield of each species by gear against time. \cr
#'   \code{\link{plot}} \tab Produces 5 plots (\code{\link{plotFeedingLevel}}, \code{\link{plotBiomass}}, \code{\link{plotPredMort}}, \code{\link{plotFMort}} and \code{\link{plotSpectra}}) in the same window as a summary. \cr
#' }
#' 
#' These functions use the ggplot2 package and return the plot as a ggplot
#' object. This means that you can manipulate the plot further after its 
#' creation using the ggplot grammar of graphics. The corresponding function
#' names with \code{plot} replaced by \code{plotly} produce interactive plots
#' with the help of the plotly package.
#' 
#' While most plot functions take their data from a MizerSim object, some of
#' those that make plots representing data at a single time can also take their
#' data from the initial values in a MizerParams object.
#' 
#' Where plots show results for species, the line colour and line type for each 
#' species are specified by the \code{linecolour} and \code{linetype} slots in
#' the MizerParams object. These were either taken from a default palette
#' hard-coded into \code{\link{emptyParams}} or they were specified by the user
#' in the species parameters dataframe used to set up the MizerParams object.
#' The \code{linecolour} and \code{linetype} slots hold named vectors, named by
#' the species. They can be overwritten by the user at any time.
#' 
#' Most plots allow the user to select to show only a subset of species,
#' specified as a vector in the \code{species} argument to the plot function.
#' 
#' The ordering of the species in the legend is the same as the ordering in
#' the species parameter data frame.
#' 
#' Sometimes one wants to show two plots side-by-side with the same axes and
#' the same legend. This is made possible for some of the plots via the
#' \code{\link{displayFrames}} function.
#' 
#' @seealso \code{\link{summary_functions}}
#' @family plotting functions
#' @name plotting_functions
#' @examples
#' # Set up example MizerParams and MizerSim objects
#' data(NS_species_params_gears)
#' data(inter)
#' params <- suppressMessages(set_multispecies_model(NS_species_params_gears, inter))
#' sim <- project(params, effort=1, t_max=20, t_save = 2, progress_bar = FALSE)
#' 
#' # Some example plots
#' plotFeedingLevel(sim)
#' 
#' # Plotting only a subset of species
#' plotFeedingLevel(sim, species = c("Cod", "Herring"))
#' 
#' # Specifying new colours and linetypes for some species
#' sim@params@linetype["Cod"] <- "solid"
#' sim@params@linecolour["Cod"] <- "red"
#' plotFeedingLevel(sim, species = c("Cod", "Herring"))
#' 
#' # Manipulating the plot
#' library(ggplot2)
#' p <- plotFeedingLevel(sim)
#' p <- p + geom_hline(aes(yintercept = 0.7))
#' p <- p + theme_bw()
#' p
NULL

# Hackiness to get past the 'no visible binding ... ' warning when running check
utils::globalVariables(c("time", "value", "Species", "w", "gear", "Age",
                         "x", "y", "Year", "Yield", "Biomass", "Size",
                         "Proportion", "Prey"))

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
log_breaks <- function(n = 6){
    n <- max(1, n)  # Because n=0 could lead to R crash
    function(x) {
        grDevices::axisTicks(log10(range(x, na.rm = TRUE)),
                             log = TRUE, nint = n)
    }
}


#' Display frames
#' 
#' Takes two data frames with plotting data and displays them side-by-side,
#' using the same axes and legend.
#' 
#' The two data frames each need to have the same three variables. The first
#' variable will go on the x-axis, the third on the y-axis with a logarithmic
#' scale. The second variable should be the species and will be used to group
#' the data and display with the linetype and linecolour specified by the
#' \code{linetype} and \code{linecolour} slots of the \code{params} object.
#' 
#' The recommended way is to obtain the data frames using one of the supplied
#' functions, e.g., \code{\link{getBiomassFrame}}, \code{\link{getSSBFrame}}.
#' 
#' @param f1 Data frame for left plot
#' @param f2 Data frame for right plot
#' @param params A MizerParams object
#' @param xlab Label for x-axis. Defaults to first variable name.
#' @param ylab Label for y-axis. Defaults to third variable name.
#' @param y_ticks The approximate number of ticks desired on the y axis
#' 
#' @return ggplot2 object
#' @export
#' @family frame functions
#' @family plotting functions
#' @seealso \code{\link{plotting_functions}},
#'   \code{\link{getBiomassFrame}}, \code{\link{getSSBFrame}}
#' @examples 
#' # Set up example MizerParams and MizerSim objects
#' data(NS_species_params_gears)
#' data(inter)
#' params <- suppressMessages(set_multispecies_model(NS_species_params_gears, inter))
#' sim0 <- project(params, effort=0, t_max=20, progress_bar = FALSE)
#' sim1 <- project(params, effort=1, t_max=20, progress_bar = FALSE)
#' 
#' # Display biomass from each simulation next to each other
#' displayFrames(getBiomassFrame(sim0), getBiomassFrame(sim1), params)
displayFrames <- function(f1, f2, params, 
                           xlab = NA, ylab = NA,
                           y_ticks = 6) {
    var_names <- names(f1)
    if (!(length(var_names) == 3)) {
        stop("A frame needs to have three variables.")
    }
    if (!all(names(f2) == var_names)) {
        stop("Both frames need to have the same variable names.")
    }
    f <- rbind(cbind(f1, Simulation = 1), cbind(f2, Simulation = 2))
    
    if (is.na(xlab)) {
        xlab <- var_names[1]
    }
    if (is.na(ylab)) {
        ylab <- var_names[3]
    }
    ytrans <- "log10"
    breaks <- log_breaks(n = y_ticks)
    
    p <- ggplot(f, aes_string(x = names(f)[1], y = names(f)[3],
                              colour = names(f)[2], linetype = names(f)[2])) +
        scale_y_continuous(trans = ytrans, breaks = breaks,
                           labels = prettyNum, name = ylab) +
        scale_x_continuous(name = xlab) +
        geom_line() +
        facet_wrap(~ Simulation) +
        scale_colour_manual(values = params@linecolour) +
        scale_linetype_manual(values = params@linetype)
    return(p)
}


#' Get data frame of spawning stock biomass of species through time, 
#' ready for ggplot2
#'
#' After running a projection, the spawning stock biomass of each species can be
#' plotted against time.
#' 
#' @param sim An object of class \linkS4class{MizerSim}
#' @param species Name or vector of names of the species to be plotted. By
#'   default all foreground species are plotted.
#' @param start_time The first time to be plotted. Default is the beginning
#'   of the time series.
#' @param end_time The last time to be plotted. Default is the end of the
#'   time series.
#' @param ylim A numeric vector of length two providing limits of for the
#'   y axis. Use NA to refer to the existing minimum or maximum. Any values
#'   below 1e-20 are always cut off.
#' @param total A boolean value that determines whether the total SSB from
#'   all species is plotted as well. Default is FALSE
#'   
#' @return A data frame that can be used in \code{\link{displayFrames}}
#' @export
#' @family frame functions
#' @seealso \code{\link{getSSB}}
getSSBFrame <- function(sim,
            species = dimnames(sim@n)$sp[!is.na(sim@params@A)],
            start_time = as.numeric(dimnames(sim@n)[[1]][1]),
            end_time = as.numeric(dimnames(sim@n)[[1]][dim(sim@n)[1]]),
            ylim = c(NA, NA), total = FALSE){
    b <- getSSB(sim)
    if (start_time >= end_time) {
        stop("start_time must be less than end_time")
    }
    # Select time range
    b <- b[(as.numeric(dimnames(b)[[1]]) >= start_time) &
               (as.numeric(dimnames(b)[[1]]) <= end_time), , drop = FALSE]
    b_total <- rowSums(b)
    # Include total
    if (total) {
        b <- cbind(b, Total = b_total)
    }
    bm <- reshape2::melt(b)
    # Implement ylim and a minimal cutoff
    min_value <- 1e-20
    bm <- bm[bm$value >= min_value &
                 (is.na(ylim[1]) | bm$value >= ylim[1]) &
                 (is.na(ylim[2]) | bm$value <= ylim[2]), ]
    names(bm) <- c("Year", "Species", "SSB")
    # Force Species column to be a factor (otherwise if numeric labels are
    # used they may be interpreted as integer and hence continuous)
    bm$Species <- as.factor(bm$Species)
    # Select species
    bm <- bm[bm$Species %in% species, ]
    return(bm)
}


#' Get data frame of biomass of species through time, ready for ggplot2
#'
#' After running a projection, the biomass of each species can be plotted
#' against time. The biomass is calculated within user defined size limits 
#' (min_w, max_w, min_l, max_l, see \code{\link{getBiomass}}). This function
#' returns a dataframe that can be displayed with 
#' \code{\link{displayFrames}}.
#' 
#' @param sim An object of class \linkS4class{MizerSim}
#' @param species Name or vector of names of the species to be plotted. By
#'   default all species are plotted.
#' @param start_time The first time to be plotted. Default is the beginning
#'   of the time series.
#' @param end_time The last time to be plotted. Default is the end of the
#'   time series.
#' @param ylim A numeric vector of length two providing lower and upper limits
#'   for the y axis. Use NA to refer to the existing minimum or maximum. Any
#'   values below 1e-20 are always cut off.
#' @param total A boolean value that determines whether the total biomass from
#'   all species is plotted as well. Default is FALSE.
#' @inheritDotParams get_size_range_array -params
#'   
#' @return A data frame that can be used in \code{\link{displayFrames}}
#' @export
#' @family frame functions
#' @seealso \code{\link{getBiomass}}, \code{\link{displayFrames}}
#' @examples 
#' # Set up example MizerParams and MizerSim objects
#' data(NS_species_params_gears)
#' data(inter)
#' params <- suppressMessages(set_multispecies_model(NS_species_params_gears, inter))
#' sim0 <- project(params, effort=0, t_max=20, progress_bar = FALSE)
#' sim1 <- project(params, effort=1, t_max=20, progress_bar = FALSE)
#' 
#' # Display biomass from each simulation next to each other
#' displayFrames(getBiomassFrame(sim0), getBiomassFrame(sim1), params)
getBiomassFrame <- function(sim,
            species = dimnames(sim@n)$sp[!is.na(sim@params@A)],
            start_time = as.numeric(dimnames(sim@n)[[1]][1]),
            end_time = as.numeric(dimnames(sim@n)[[1]][dim(sim@n)[1]]),
            ylim = c(NA, NA), total = FALSE, ...) {
    b <- getBiomass(sim, ...)
    if (start_time >= end_time) {
        stop("start_time must be less than end_time")
    }
    # Select time range
    b <- b[(as.numeric(dimnames(b)[[1]]) >= start_time) &
               (as.numeric(dimnames(b)[[1]]) <= end_time), , drop = FALSE]
    b_total <- rowSums(b)
    # Include total
    if (total) {
        b <- cbind(b, Total = b_total)
        species <- c("Total", species)
    }
    bm <- reshape2::melt(b)
    
    # Implement ylim and a minimal cutoff
    min_value <- 1e-20
    bm <- bm[bm$value >= min_value &
                 (is.na(ylim[1]) | bm$value >= ylim[1]) &
                 (is.na(ylim[2]) | bm$value <= ylim[2]), ]
    names(bm) <- c("Year", "Species", "Biomass")
    
    # Force Species column to be a factor (otherwise if numeric labels are
    # used they may be interpreted as integer and hence continuous).
    # Need to keep species in order for legend.
    species_levels <- c(dimnames(sim@n)$sp, "Background", "Plankton", "Total")
    bm$Species <- factor(bm$Species, levels = species_levels)
    
    # Select species
    bm <- bm[bm$Species %in% species, ]

    return(bm)
}


#' Plot the biomass of species through time
#'
#' After running a projection, the biomass of each species can be plotted
#' against time. The biomass is calculated within user defined size limits 
#' (min_w, max_w, min_l, max_l, see \code{\link{getBiomass}}). 
#' 
#' @inheritParams getBiomassFrame
#' @inheritParams plotSpectra
#' @param y_ticks The approximate number of ticks desired on the y axis
#' @inheritDotParams get_size_range_array -params
#'   
#' @return A ggplot2 object
#' @export
#' @family plotting functions
#' @seealso \code{\link{plotting_functions}}, \code{\link{getBiomass}}
#' @examples
#' # Set up example MizerParams and MizerSim objects
#' data(NS_species_params_gears)
#' data(inter)
#' params <- suppressMessages(set_multispecies_model(NS_species_params_gears, inter))
#' sim <- project(params, effort = 1, t_max = 20, t_save = 0.2, progress_bar = FALSE)
#' 
#' plotBiomass(sim)
#' plotBiomass(sim, species = c("Cod", "Haddock"), total = TRUE)
#' plotBiomass(sim, min_w = 10, max_w = 1000)
#' plotBiomass(sim, start_time = 10, end_time = 15)
#' plotBiomass(sim, y_ticks = 3)
#' 
plotBiomass <- function(sim,
            species = dimnames(sim@n)$sp[!is.na(sim@params@A)],
            start_time = as.numeric(dimnames(sim@n)[[1]][1]),
            end_time = as.numeric(dimnames(sim@n)[[1]][dim(sim@n)[1]]),
            y_ticks = 6,
            ylim = c(NA, NA),
            total = FALSE, background = TRUE, highlight = NULL, ...) {
    # First we get the data frame for all species, including the background
    bm <- getBiomassFrame(sim, species = dimnames(sim@n)$sp,
                          start_time = start_time,
                          end_time = end_time,
                          ylim = ylim, total = total, ...)
    # Select species
    spec_bm <- bm[bm$Species %in% c("Total", species), ]
    x_label <- "Year"
    y_label <- "Biomass [g]"
    p <- ggplot(spec_bm, aes(x = Year, y = Biomass)) +
        scale_y_continuous(trans = "log10", breaks = log_breaks(n = y_ticks),
                           labels = prettyNum, name = y_label) +
        scale_x_continuous(name = x_label) +
        scale_colour_manual(values = sim@params@linecolour) +
        scale_linetype_manual(values = sim@params@linetype)

    if (background) {
        # Add background species in light grey
        back_sp <- dimnames(sim@n)$sp[is.na(sim@params@A)]
        back_bm <- bm[bm$Species %in% back_sp, ]
        if (nrow(back_bm) > 0) {
            p <- p + geom_line(aes(group = Species), data = back_bm,
                               colour = sim@params@linecolour["Background"],
                               linetype = sim@params@linetype["Background"])
        }
    }
    
    linesize <- rep(0.8, length(sim@params@linetype))
    names(linesize) <- names(sim@params@linetype)
    linesize[highlight] <- 1.6
    p <- p + scale_size_manual(values = linesize) +
        geom_line(aes(colour = Species, linetype = Species, size = Species))

    return(p)
}

#' Plot the biomass of species against time with plotly
#' 
#' @inherit plotBiomass params return description details seealso
#' @inheritDotParams get_size_range_array -params
#' @export
#' @family plotting functions
plotlyBiomass <- function(sim,
             species = dimnames(sim@n)$sp[!is.na(sim@params@A)],
             start_time = as.numeric(dimnames(sim@n)[[1]][1]),
             end_time = as.numeric(dimnames(sim@n)[[1]][dim(sim@n)[1]]),
             y_ticks = 6,
             ylim = c(NA, NA),
             total = FALSE,
             background = TRUE,
             highlight = NULL,
             ...) {
    argg <- c(as.list(environment()), list(...))
    ggplotly(do.call("plotBiomass", argg))
    # ggplotly(do.callplotBiomass(sim, species, start_time, end_time, y_ticks,
    #                      ylim, total, background, ...))
}


#' Plot the total yield of species through time
#'
#' After running a projection, the total yield of each species across all 
#' fishing gears can be plotted against time. The yield is obtained with
#' \code{\link{getYield}}.
#' 
#' @param sim An object of class \linkS4class{MizerSim}
#' @param sim2 An optional second object of class \linkS4class{MizerSim}. If
#'   this is provided its yields will be shown on the same plot in bolder lines.
#' @inheritParams plotSpectra
#' @param log Boolean whether yield should be plotted on a logarithmic axis. 
#'   Defaults to true.
#'
#' @return A ggplot2 object
#' @export
#' @family plotting functions
#' @seealso \code{\link{plotting_functions}},  \code{\link{getYield}}
#' @examples
#' 
#' data(NS_species_params_gears)
#' data(inter)
#' params <- suppressMessages(set_multispecies_model(NS_species_params_gears, inter))
#' sim <- project(params, effort=1, t_max=20, t_save = 0.2, progress_bar = FALSE)
#' plotYield(sim)
#' plotYield(sim, species = c("Cod", "Herring"), total = TRUE)
#' 
#' # Comparing with yield from twice the effort
#' sim2 <- project(params, effort=2, t_max=20, t_save = 0.2, progress_bar = FALSE)
#' plotYield(sim, sim2, species = c("Cod", "Herring"), log = FALSE)
#' 
plotYield <- function(sim, sim2,
                      species = dimnames(sim@n)$sp,
                      total = FALSE, log = TRUE,
                      highlight = NULL, ...){
    # Need to keep species in order for legend
    species_levels <- c(dimnames(sim@n)$sp, "Background", "Plankton", "Total")
    if (missing(sim2)) {
        y <- getYield(sim, ...)
        y_total <- rowSums(y)
        y <- y[, (as.character(dimnames(y)[[2]]) %in% species),
               drop = FALSE]
        if (total) {
            # Include total
            y <- cbind(y, "Total" = y_total)
        }
        ym <- reshape2::melt(y, varnames = c("Year", "Species"),
                             value.name = "Yield")
        ym$Species <- factor(ym$Species, levels = species_levels)
        ym <- subset(ym, ym$Yield > 0)
        p <- ggplot(ym) +
                geom_line(aes(x = Year, y = Yield,
                              colour = Species, linetype = Species,
                              size = Species))

        if (log) {
            p <- p + scale_y_continuous(trans = "log10", name = "Yield [g/year]",
                                        breaks = log_breaks(),
                                        labels = prettyNum)
        } else {
            p <- p + scale_y_continuous(name = "Yield [g/year]")
        }
        linesize <- rep(0.8, length(sim@params@linetype))
        names(linesize) <- names(sim@params@linetype)
        linesize[highlight] <- 1.6
        p <- p +
            scale_colour_manual(values = sim@params@linecolour) +
            scale_linetype_manual(values = sim@params@linetype) +
            scale_size_manual(values = linesize)
        return(p)
    } else {
        if (!all(dimnames(sim@n)$time == dimnames(sim2@n)$time)) {
            stop("The two simulations do not have the same times")
        }
        y <- getYield(sim, ...)
        y2 <- getYield(sim2, ...)
        y_total <- rowSums(y)
        y <- y[, (as.character(dimnames(y)[[2]]) %in% species) & colSums(y) > 0,
               drop = FALSE]
        y2_total <- rowSums(y2)
        y2 <- y2[, (as.character(dimnames(y2)[[2]]) %in% species),
                 drop = FALSE]
        if (total) {
            # Include total
            y <- cbind(y, Total = y_total)
            y2 <- cbind(y2, Total = y2_total)
        }
        ym <- reshape2::melt(y, varnames = c("Year", "Species"),
                             value.name = "Yield")
        ym2 <- reshape2::melt(y2, varnames = c("Year", "Species"),
                              value.name = "Yield")
        ym$Simulation <- 1
        ym2$Simulation <- 2
        ym <- rbind(ym, ym2)
        ym$Species <- factor(ym$Species, levels = species_levels)
        ym$Simulation <- as.factor(ym$Simulation)
        ym <- subset(ym, ym$Yield > 0)
        p <- ggplot(ym) +
                geom_line(aes(x = Year, y = Yield, colour = Species,
                              linetype = Species))

        if (log) {
            p <- p + scale_y_continuous(trans = "log10", name = "Yield [g/year]")
        } else {
            p <- p + scale_y_continuous(name = "Yield [g/year]")
        }
        p <- p + facet_wrap(~ Simulation)
        return(p)
    }
}

#' Plot the total yield of species through time with plotly
#' @inherit plotYield params return description details seealso
#' @export
#' @family plotting functions
plotlyYield <- function(sim, sim2,
                        species = dimnames(sim@n)$sp,
                        total = FALSE, log = TRUE,
                        highlight = NULL, ...) {
    argg <- as.list(environment())
    ggplotly(do.call("plotYield", argg))
}


#' Plot the total yield of each species by gear through time
#'
#' After running a projection, the total yield of each species by fishing gear
#' can be plotted against time. 
#' 
#' This plot is pretty easy to do by hand. It just
#' gets the biomass using the \code{\link{getYieldGear}} method and plots using
#' the ggplot2 package. You can then fiddle about with colours and linetypes
#' etc. Just look at the source code for details.
#' 
#' @param sim An object of class \linkS4class{MizerSim}
#' @inheritParams plotSpectra
#'
#' @return A ggplot2 object
#' @export
#' @family plotting functions
#' @seealso \code{\link{plotting_functions}},  \code{\link{getYieldGear}}
#' @examples
#' 
#' data(NS_species_params_gears)
#' data(inter)
#' params <- suppressMessages(set_multispecies_model(NS_species_params_gears, inter))
#' sim <- project(params, effort=1, t_max=20, t_save = 0.2, progress_bar = FALSE)
#' plotYieldGear(sim)
#' plotYieldGear(sim, species = c("Cod", "Herring"), total = TRUE)
#' 
plotYieldGear <- function(sim,
                          species = dimnames(sim@n)$sp,
                          total = FALSE,
                          highlight = NULL, ...){
    # Need to keep species in order for legend
    species_levels <- c(dimnames(sim@n)$sp, "Background", "Plankton", "Total")
    
    y <- getYieldGear(sim, ...)
    y_total <- rowSums(y, dims = 2)
    y <- y[, , dimnames(y)$sp %in% species, drop = FALSE]
    names(dimnames(y))[names(dimnames(y)) == "sp"] <- "Species"
    ym <- reshape2::melt(y)
    ym$Species <- factor(ym$Species, levels = species_levels)
    if (total) {
        yt <- reshape2::melt(y_total)
        yt$Species <- "Total"
        ym <- rbind(ym, yt)
    }
    ym <- subset(ym, ym$value > 0)
    p <- ggplot(ym) +
            geom_line(aes(x = time, y = value, colour = Species, 
                          linetype = gear, size = Species))

    linesize <- rep(0.8, length(sim@params@linetype))
    names(linesize) <- names(sim@params@linetype)
    linesize[highlight] <- 1.6
    p <- p + scale_y_continuous(trans = "log10", name = "Yield [g]") +
        scale_x_continuous(name = "Year") +
        scale_colour_manual(values = sim@params@linecolour) +
        scale_size_manual(values = linesize)
    return(p)
}

#' Plot the total yield of each species by gear through time with plotly
#' @inherit plotYieldGear params return description details seealso
#' @export
#' @family plotting functions
plotlyYieldGear <- function(sim,
                            species = dimnames(sim@n)$sp,
                            total = FALSE, highlight = NULL, ...) {
    argg <- as.list(environment())
    ggplotly(do.call("plotYieldGear", argg))
}

#' Plot the abundance spectra
#' 
#' Plots the number density multiplied by a power of the weight, with the power
#' specified by the \code{power} argument.
#'
#' When called with a \linkS4class{MizerSim} object, the abundance is averaged
#' over the specified time range (a single value for the time range can be used
#' to plot a single time step). When called with a \linkS4class{MizerParams}
#' object the initial abundance is plotted.
#' 
#' @param object An object of class \linkS4class{MizerSim} or \linkS4class{MizerParams}.
#' @param species Name or vector of names of the species to be plotted. By
#'   default all species are plotted.
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
#' raised to \code{power}. The default \code{power = 1} gives the biomass 
#' density, whereas \code{power = 2} gives the biomass density with respect
#' to logarithmic size bins.
#' @param biomass Obsolete. Only used if \code{power} argument is missing. Then
#'   \code{biomass = TRUE} is equivalent to \code{power=1} and 
#'   \code{biomass = FALSE} is equivalent to \code{power=0}
#' @param total A boolean value that determines whether the total over all
#'   species in the system is plotted as well. Default is FALSE
#' @param plankton A boolean value that determines whether plankton is included.
#'   Default is TRUE.
#' @param background A boolean value that determines whether background species
#'   are included. Ignored if the model does not contain background species.
#'   Default is TRUE.
#' @param highlight Name or vector of names of the species to be highlighted.
#' @param ... Other arguments (currently unused)
#'   
#' @return A ggplot2 object
#' @export
#' @family plotting functions
#' @seealso \code{\link{plotting_functions}}
#' @examples
#' 
#' data(NS_species_params_gears)
#' data(inter)
#' params <- suppressMessages(set_multispecies_model(NS_species_params_gears, inter))
#' sim <- project(params, effort=1, t_max=20, t_save = 2, progress_bar = FALSE)
#' plotSpectra(sim)
#' plotSpectra(sim, wlim = c(1e-6, NA))
#' plotSpectra(sim, time_range = 10:20)
#' plotSpectra(sim, time_range = 10:20, power = 0)
#' plotSpectra(sim, species = c("Cod", "Herring"), power = 1)
#' 
plotSpectra <- function(object, species = NULL,
                        time_range,
                        wlim = c(NA, NA), ylim = c(NA, NA),
                        power = 1, biomass = TRUE,
                        total = FALSE, plankton = TRUE, 
                        background = TRUE,
                        highlight = NULL, ...) {
    # to deal with old-type biomass argument
    if (missing(power)) {
        power <- as.numeric(biomass)
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
                           power = power,
                           total = total, plankton = plankton,
                           background = background, highlight = highlight)
        return(ps)
    } else {
        ps <- plot_spectra(object, n = object@initial_n,
                           n_pp = object@initial_n_pp,
                           species = species, wlim = wlim, ylim = ylim,
                           power = power,
                           total = total, plankton = plankton,
                           background = background, highlight = highlight)
        return(ps)
    }
}


plot_spectra <- function(params, n, n_pp,
                         species, wlim, ylim, power,
                         total, plankton, background, highlight) {
    if (is.na(wlim[1])) {
        wlim[1] <- min(params@w) / 100
    }
    if (is.na(wlim[2])) {
        wlim[2] <- max(params@w_full)
    }
    # Need to keep species in order for legend
    species_levels <- c(dimnames(params@initial_n)$sp,
                        "Background", "Plankton", "Total")
    if (total) {
        # Calculate total community abundance
        fish_idx <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)
        total_n <- n_pp
        total_n[fish_idx] <- total_n[fish_idx] + colSums(n)
        total_n <- total_n * params@w_full^power
    }
    # Set species if missing to list of all non-background species
    if (is.null(species)) {
        species <- dimnames(params@initial_n)$sp[!is.na(params@A)]
    }
    species <- as.character(species)
    invalid_species <- 
        !(species %in% as.character(dimnames(params@initial_n)[[1]]))
    if (any(invalid_species)) {
        warning(paste("The following species do not exist in the model and are ignored:",
                      species[invalid_species]))
    }
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
    plot_dat <- data.frame(value = c(spec_n),
                           # ordering of factor is important for legend
                           Species = factor(dimnames(spec_n)[[1]],
                                            levels = species_levels),
                           w = rep(params@w,
                                   each = dim(spec_n)[[1]]))
    if (plankton) {
        plankton_sel <- (params@w_full >= wlim[1]) & 
                        (params@w_full <= wlim[2])
        # Do we have any plankton to plot?
        if (sum(plankton_sel) > 0) {
            w_plankton <- params@w_full[plankton_sel]
            plank_n <- n_pp[plankton_sel] * w_plankton^power
            plot_dat <- rbind(plot_dat,
                              data.frame(value = c(plank_n),
                                         Species = "Plankton",
                                         w = w_plankton))
        }
    }
    if (total) {
        plot_dat <- rbind(plot_dat,
                          data.frame(value = c(total_n),
                                     Species = "Total",
                                     w = params@w_full))
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
    # Create plot
    p <- ggplot(plot_dat, aes(x = w, y = value)) +
        scale_x_continuous(name = "Size [g]", trans = "log10",
                           breaks = log_breaks()) +
        scale_y_continuous(name = y_label, trans = "log10",
                           breaks = log_breaks()) +
        scale_colour_manual(values = params@linecolour) +
        scale_linetype_manual(values = params@linetype)
    if (background) {
        back_n <- n[is.na(params@A), , drop = FALSE]
        plot_back <- data.frame(value = c(back_n),
                                Species = as.factor(dimnames(back_n)[[1]]),
                                w = rep(params@w,
                                        each = dim(back_n)[[1]]))
        # lop off 0s and apply wlim
        plot_back <- plot_back[(plot_back$value > 0) & 
                                   (plot_back$w >= wlim[1]) &
                                   (plot_back$w <= wlim[2]), ]
        # Impose ylim
        if (!is.na(ylim[2])) {
            plot_back <- plot_back[plot_back$value <= ylim[2], ]
        }
        plot_back <- plot_back[plot_back$value > ylim[1], ]
        if (nrow(plot_back) > 0) {
            # Add background species
            p <- p +
                geom_line(aes(group = Species),
                          colour = params@linecolour["Background"],
                          linetype = params@linetype["Background"],
                          data = plot_back)
        }
    }
    linesize <- rep(0.8, length(params@linetype))
    names(linesize) <- names(params@linetype)
    linesize[highlight] <- 1.6
    p <- p + scale_size_manual(values = linesize) + 
        geom_line(aes(colour = Species, linetype = Species, size = Species))
    return(p)
}

#' Plotly plot of the abundance spectra
#' @inherit plotSpectra params return description details seealso
#' @export
#' @family plotting functions
plotlySpectra <- function(object, species = NULL,
                        time_range,
                        wlim = c(NA, NA), ylim = c(NA, NA),
                        power = 1, biomass = TRUE,
                        total = FALSE, plankton = TRUE, 
                        background = TRUE,
                        highlight = NULL, ...) {
    argg <- as.list(environment())
    ggplotly(do.call("plotSpectra", argg))
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
#' @inheritParams plotSpectra
#' @param all.sizes If TRUE, then feeding level is plotted also for sizes
#'   outside a species' size range. Default FALSE.
#'
#' @return A ggplot2 object
#' @export
#' @family plotting functions
#' @seealso \code{\link{plotting_functions}}, \code{\link{getFeedingLevel}}
#' @examples
#' 
#' data(NS_species_params_gears)
#' data(inter)
#' params <- suppressMessages(set_multispecies_model(NS_species_params_gears, inter))
#' sim <- project(params, effort=1, t_max=20, t_save = 2, progress_bar = FALSE)
#' plotFeedingLevel(sim)
#' plotFeedingLevel(sim, time_range = 10:20, species = c("Cod", "Herring"))
#' 
plotFeedingLevel <- function(object,
            species = NULL,
            time_range,
            highlight = NULL,
            all.sizes = FALSE, ...) {
    if (is(object, "MizerSim")) {
        if (missing(time_range)) {
            time_range  <- max(as.numeric(dimnames(object@n)$time))
        }
        params <- object@params
    } else {
        assert_that(is(object, "MizerParams"))
        params <- object
    }
    feed <- getFeedingLevel(params, time_range = time_range, drop = FALSE)
    # If a time range was returned, average over it
    if (length(dim(feed)) == 3) {
        feed <- apply(feed, c(2, 3), mean)
    }
        
    # selector for desired species
    # Set species if missing to list of all non-background species
    if (is.null(species)) {
        species <- dimnames(params@initial_n)$sp[!is.na(params@A)]
    }
    sel_sp <- as.character(dimnames(feed)$sp) %in% species
    feed <- feed[sel_sp, , drop = FALSE]
    plot_dat <- data.frame(value = c(feed),
                           # ggplot orders the legend according to the ordering
                           # of the factors, hence we need the levels argument
                           Species = factor(dimnames(feed)$sp, 
                                            levels = dimnames(feed)$sp),
                           w = rep(params@w, each = length(species)))
    
    if (!all.sizes) {
        # Remove feeding level for sizes outside a species' size range
        for (sp in species) {
            plot_dat$value[plot_dat$Species == sp &
                           (plot_dat$w < params@species_params[sp, "w_min"] |
                            plot_dat$w > params@species_params[sp, "w_inf"])] <- NA
        }
        plot_dat <- plot_dat[complete.cases(plot_dat), ]
    }
    
    p <- ggplot(plot_dat) +
            geom_line(aes(x = w, y = value, colour = Species, 
                          linetype = Species, size = Species))

    linesize <- rep(0.8, length(params@linetype))
    names(linesize) <- names(params@linetype)
    linesize[highlight] <- 1.6
    p <- p +
        scale_x_continuous(name = "Size [g]", trans = "log10") +
        scale_y_continuous(name = "Feeding Level", limits = c(0, 1)) +
        scale_colour_manual(values = params@linecolour) +
        scale_linetype_manual(values = params@linetype) +
        scale_size_manual(values = linesize)
    return(p)
}

#' Plot the feeding level of species by size with plotly
#' 
#' @inherit plotFeedingLevel params return description details seealso
#' @export
#' @family plotting functions
plotlyFeedingLevel <- function(object,
                             species = NULL,
                             time_range,
                             highlight = NULL, ...) {
    argg <- as.list(environment())
    ggplotly(do.call("plotFeedingLevel", argg))
}
    

#' Plot predation mortality rate of each species against size
#' 
#' After running a projection, plot the predation mortality rate of each species
#' by size. The mortality rate is averaged over the specified time range (a
#' single value for the time range can be used to plot a single time step).
#' 
#' @inheritParams plotSpectra
#'
#' @return A ggplot2 object
#' @export
#' @family plotting functions
#' @seealso \code{\link{plotting_functions}},  \code{\link{getPredMort}}
#' @examples
#' 
#' data(NS_species_params_gears)
#' data(inter)
#' params <- suppressMessages(set_multispecies_model(NS_species_params_gears, inter))
#' sim <- project(params, effort=1, t_max=20, t_save = 2, progress_bar = FALSE)
#' plotPredMort(sim)
#' plotPredMort(sim, time_range = 10:20)
#' 
plotPredMort <- function(object, species = NULL,
                         time_range,
                         highlight = NULL, ...) {
    if (is(object, "MizerSim")) {
        if (missing(time_range)) {
            time_range  <- max(as.numeric(dimnames(object@n)$time))
        }
        params <- object@params
    } else {
        assert_that(is(object, "MizerParams"))
        params <- object
    }
    m2 <- getPredMort(object, time_range = time_range, drop = FALSE)
    # If a time range was returned, average over it
    if (length(dim(m2)) == 3) {
        m2 <- apply(m2, c(2, 3), mean)
    }
    
    # selector for desired species
    # Set species if missing to list of all non-background species
    if (is.null(species)) {
        species <- dimnames(params@initial_n)$sp[!is.na(params@A)]
    }
    # Need to keep species in order for legend
    species_levels <- c(as.character(params@species_params$species), 
                        "Background", "Plankton", "Total")
    m2 <- m2[as.character(dimnames(m2)[[1]]) %in% species, , drop = FALSE]
    plot_dat <- data.frame(value = c(m2),
                           Species = factor(dimnames(m2)[[1]],
                                            levels = species_levels),
                           w = rep(params@w, each = length(species)))
    p <- ggplot(plot_dat) +
            geom_line(aes(x = w, y = value, colour = Species, 
                          linetype = Species, size = Species))

    linesize <- rep(0.8, length(params@linetype))
    names(linesize) <- names(params@linetype)
    linesize[highlight] <- 1.6
    p <- p +
        scale_x_continuous(name = "Size [g]", trans = "log10") +
        scale_y_continuous(name = "Predation mortality [1/year]",
                           limits = c(0, max(plot_dat$value))) +
        scale_colour_manual(values = params@linecolour) +
        scale_linetype_manual(values = params@linetype) +
        scale_size_manual(values = linesize)
    return(p)
}

#' Alias for plotPredMort
#' 
#' An alias provided for backward compatibility with mizer version <= 1.0
#' @inherit plotPredMort
#' @export
plotM2 <- plotPredMort

#' Plot predation mortality rate of each species against size with plotly
#' @inherit plotPredMort params return description details seealso
#' @export
#' @family plotting functions
plotlyPredMort <- function(object, species = NULL,
                           time_range,
                           highlight = NULL, ...) {
    argg <- as.list(environment())
    ggplotly(do.call("plotPredMort", argg))
}

#' Plot total fishing mortality of each species by size
#' 
#' After running a projection, plot the total fishing mortality of each species
#' by size. The total fishing mortality is averaged over the specified time
#' range (a single value for the time range can be used to plot a single time
#' step).
#' 
#' @inheritParams plotSpectra
#'
#' @return A ggplot2 object
#' @export
#' @family plotting functions
#' @seealso \code{\link{plotting_functions}}, \code{\link{getFMort}}
#' @examples
#' 
#' data(NS_species_params_gears)
#' data(inter)
#' params <- suppressMessages(set_multispecies_model(NS_species_params_gears, inter))
#' sim <- project(params, effort=1, t_max=20, t_save = 2, progress_bar = FALSE)
#' plotFMort(sim)
#' plotFMort(sim, highlight = c("Cod", "Haddock"))
#' 
plotFMort <- function(object, species = NULL,
                      time_range,
                      highlight = NULL, ...) {
    if (is(object, "MizerSim")) {
        if (missing(time_range)) {
            time_range  <- max(as.numeric(dimnames(object@n)$time))
        }
        params <- object@params
    } else {
        assert_that(is(object, "MizerParams"))
        params <- object
    }
    f <- getFMort(object, time_range = time_range, drop = FALSE)
    # If a time range was returned, average over it
    if (length(dim(f)) == 3) {
        f <- apply(f, c(2, 3), mean)
    }
    # selector for desired species
    # Set species if missing to list of all non-background species
    if (is.null(species)) {
        species <- dimnames(params@initial_n)$sp[!is.na(params@A)]
    }
    # Need to keep species in order for legend
    species_levels <- c(as.character(params@species_params$species), 
                        "Background", "Plankton", "Total")
    f <- f[as.character(dimnames(f)[[1]]) %in% species, , drop = FALSE]
    plot_dat <- data.frame(value = c(f),
                           Species = factor(dimnames(f)[[1]],
                                            levels = species_levels),
                           w = rep(params@w, each = length(species)))
    
    linesize <- rep(0.8, length(params@linetype))
    names(linesize) <- names(params@linetype)
    linesize[highlight] <- 1.6
    p <- ggplot(plot_dat) +
            geom_line(aes(x = w, y = value, colour = Species, 
                          linetype = Species, size = Species))

    p <- p +
        scale_x_continuous(name = "Size [g]", trans = "log10") +
        scale_y_continuous(name = "Fishing mortality [1/Year]",
                           limits = c(0, max(plot_dat$value))) +
        scale_colour_manual(values = params@linecolour) +
        scale_linetype_manual(values = params@linetype) + 
        scale_size_manual(values = linesize)
    return(p)
}

#' Plot total fishing mortality of each species by size with plotly
#' @inherit plotPredMort params return description details seealso
#' @export
#' @family plotting functions
plotlyFMort <- function(object, species = NULL,
                        time_range,
                        highlight = NULL, ...) {
    argg <- as.list(environment())
    ggplotly(do.call("plotFMort", argg))
}


#' Plot growth curves giving weight as a function of age
#' 
#' If given a \linkS4class{MizerSim} object, uses the growth rates at the final
#' time of a simulation to calculate the size at age. If given a
#' \linkS4class{MizerParams} object, uses the initial growth rates instead.
#' 
#' When the growth curve for only a single species is plotted, horizontal
#' lines are included that indicate the maturity size and the maximum size for 
#' that species. If furthermore the species parameters contain the variables
#' a and b for length to weight conversion and the von Bertalanffy parameter
#' k_vb, then the von Bertalanffy growth curve is superimposed in black.
#' 
#' @inheritParams getGrowthCurves
#' @inheritParams plotSpectra
#' 
#' @return A ggplot2 object
#' @export
#' @family plotting functions
#' @seealso \code{\link{plotting_functions}}
#' @examples
#' 
#' data(NS_species_params_gears)
#' data(inter)
#' params <- suppressMessages(set_multispecies_model(NS_species_params_gears, inter))
#' sim <- project(params, effort=1, t_max=20, t_save = 2, progress_bar = FALSE)
#' plotGrowthCurves(sim, percentage = TRUE)
#' plotGrowthCurves(sim, species = "Cod", max_age = 24)
#' 
plotGrowthCurves <- function(object, 
                             species,
                             max_age = 20,
                             percentage = FALSE,
                             highlight = NULL) {
    if (is(object, "MizerSim")) {
        params <- object@params
        t <- dim(object@n)[1]
        params@initial_n <- object@n[t, , ]
        params@initial_n_pp <- object@n_pp[t, ]
    } else if (is(object, "MizerParams")) {
        params <- object
    }
    if (missing(species)) {
        species <- params@species_params$species
    }
    ws <- getGrowthCurves(params, species, max_age, percentage)
    plot_dat <- reshape2::melt(ws)
    # Need to keep species in order for legend
    plot_dat$Species <- factor(plot_dat$Species, params@species_params$species)
    p <- ggplot(plot_dat) +
            geom_line(aes(x = Age, y = value,
                          colour = Species, linetype = Species,
                          size = Species))

    y_label <- if (percentage) "Percent of maximum size" else "Size [g]"
    linesize <- rep(0.8, length(params@linetype))
    names(linesize) <- names(params@linetype)
    linesize[highlight] <- 1.6
    p <- p +
        scale_x_continuous(name = "Age [Years]") +
        scale_y_continuous(name = y_label) +
        scale_colour_manual(values = params@linecolour) +
        scale_linetype_manual(values = params@linetype) +
        scale_size_manual(values = linesize)
    
    # Extra stuff for single-species case
    if (length(species) == 1 && !percentage) {
        idx <- which(params@species_params$species == species)
        w_inf <- params@species_params$w_inf[idx]
        # set w_inf to w at next grid point, because that is when growth rate
        # becomes zero
        w_inf <- params@w[min(sum(w_inf > params@w) + 1, length(params@w))]
        p <- p + geom_hline(yintercept = w_inf) +
            annotate("text", 0, w_inf, vjust = -1, label = "Maximum")
        w_mat <- params@species_params$w_mat[idx]
        p <- p + geom_hline(yintercept = w_mat) +
            annotate("text", 0, w_mat, vjust = -1, label = "Maturity")
        if (all(c("a", "b", "k_vb") %in% names(params@species_params))) {
            age <- as.numeric(dimnames(ws)$Age)
            a <- params@species_params$a[idx]
            b <- params@species_params$b[idx]
            k_vb <- params@species_params$k_vb[idx]
            t0 <- params@species_params$t0[idx]
            if (is.null(t0)) {t0 <- 0}
            L_inf <- (w_inf/a)^(1/b)
            vb <- a * (L_inf * (1 - exp(-k_vb * (age - t0))))^b
            dat <- data.frame(x = age, y = vb)
            p <- p + geom_line(data = dat, aes(x = x, y = y))
        }
    }
    return(p)
}

#' Plot growth curves giving weight as a function of age with plotly
#' @inherit plotGrowthCurves params return description details seealso
#' @export
#' @family plotting functions
plotlyGrowthCurves <- function(object, species,
                             max_age = 20,
                             percentage = FALSE,
                             highlight = NULL) {
    argg <- as.list(environment())
    ggplotly(do.call("plotGrowthCurves", argg))
}


#' Plot diet
#' 
#' @inheritParams plotSpectra
#'
#' @return A ggplot2 object
#' @export
#' @family plotting functions
plotDiet <- function(object, species) {
    params <- object
    if (is.integer(species)) {
        species <- params@species_params$species[species]
    }
    diet <- getDiet(params)[params@species_params$species == species, , ]
    prey <- dimnames(diet)$prey
    prey <- factor(prey, levels = rev(prey))
    plot_dat <- data.frame(
        Proportion = c(diet),
        w = params@w,
        Prey = rep(prey, each = length(params@w)))
    plot_dat <- plot_dat[plot_dat$Proportion > 0, ]
    ggplot(plot_dat) +
        geom_area(aes(x = w, y = Proportion, fill = Prey)) +
        scale_x_log10() +
        labs(x = "Size [g]") +
        scale_fill_manual(values = params@linecolour)
}


#### plot ####
#' Summary plot for \code{MizerSim} objects
#' 
#' After running a projection, produces 5 plots in the same window: feeding
#' level, abundance spectra, predation mortality and fishing mortality of each
#' species by size; and biomass of each species through time. This method just
#' uses the other plotting functions and puts them all in one window.
#' 
#' @param x An object of class \linkS4class{MizerSim}
#' @param y Not used
#' @param ...  For additional arguments see the documentation for
#'   \code{\link{plotBiomass}},
#'   \code{\link{plotFeedingLevel}},\code{\link{plotSpectra}},\code{\link{plotPredMort}}
#'   and \code{\link{plotFMort}}.
#' @return A viewport object
#' @export
#' @family plotting functions
#' @seealso \code{\link{plotting_functions}}
#' @rdname plotMizerSim
#' @examples
#' 
#' data(NS_species_params_gears)
#' data(inter)
#' params <- suppressMessages(set_multispecies_model(NS_species_params_gears, inter))
#' sim <- project(params, effort=1, t_max=20, t_save = 2, progress_bar = FALSE)
#' plot(sim)
#' plot(sim, time_range = 10:20) # change time period for size-based plots
#' plot(sim, min_w = 10, max_w = 1000) # change size range for biomass plot
#' 
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

#' Summary plot for \code{MizerParams} objects
#' 
#' Produces 3 plots in the same window: abundance spectra, feeding
#' level and predation mortality of each species through time. This method just
#' uses the other plotting functions and puts them all in one window.
#' 
#' @return A viewport object
#' @export
#' @family plotting functions
#' @seealso \code{\link{plotting_functions}}
#' @rdname plotMizerSim
#' @examples
#' 
#' data(NS_species_params_gears)
#' data(inter)
#' params <- suppressMessages(set_multispecies_model(NS_species_params_gears, inter))
#' plot(params)
#' plot(params, min_w = 10, max_w = 1000) # change size range for biomass plot
#' 
setMethod("plot", signature(x = "MizerParams", y = "missing"),
          function(x, ...) {
              p11 <- plotFeedingLevel(x, ...)
              p2 <- plotSpectra(x, ...)
              p12 <- plotPredMort(x, ...)
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
