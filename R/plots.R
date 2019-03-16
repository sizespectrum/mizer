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
#'   \code{\link{plotBiomass}} \tab Plots the total biomass of each species through time. A time range to be plotted can be specified. The size range of the community can be specified in the same way as the method getBiomass(). \cr
#'   \code{\link{plotSpectra}} \tab Plots the abundance (biomass or numbers) spectra of each species and the background community. It is possible to specify a minimum size which is useful for truncating the plot. \cr
#'   \code{\link{plotFeedingLevel}} \tab Plots the feeding level of each species against size. \cr
#'   \code{\link{plotM2}} \tab Plots the predation mortality of each species against size. \cr
#'   \code{\link{plotFMort}} \tab Plots the total fishing mortality of each species against size. \cr
#'   \code{\link{plotYield}} \tab Plots the total yield of each species across all fishing gears against time. \cr
#'   \code{\link{plotYieldGear}} \tab Plots the total yield of each species by gear against time. \cr
#'   \code{\link{plot}} \tab Produces 5 plots (plotFeedingLevel(), plotBiomass(), plotM2(), plotFMort() and plotSpectra()) in the same window as a summary. \cr
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
#' species are specified by the `linecolour` and `linetype` slots in the 
#' MizerParams object. These were either taken from a default palette
#' hard-coded into \code{\link{emptyParams}} or they were specified by the user
#' in the species parameters dataframe used to set up the MizerParams object.
#' The `linecolour` and `linetype` slots hold named vectors, named by the
#' species. They can be overwritten by the user at any time.
#' 
#' Most plots allow the user to select to show only a subset of species,
#' specified as a vector in the `species` argument to the plot function.
#' 
#' The ordering of the species in the legend is the same as the ordering in
#' the species parameter data frame.
#' 
#' Sometimes one wants to show two plots side-by-side with the same axes and
#' the same legend. This is made possible for some of the plots via the
#' \code{\link{display_frames}} function.
#' 
#' @seealso \code{\link{summary_functions}}
#' @name plotting_functions
#' @examples
#' # Set up example MizerParams and MizerSim objects
#' data(NS_species_params_gears)
#' data(inter)
#' params <- suppressMessages(MizerParams(NS_species_params_gears, inter))
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
                         "x", "y", "Year", "Yield"))

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
#' `linetype` and `linecolour` slots of the `params` object.
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
#' @seealso \code{\link{getBiomassFrame}}, \code{\link{getSSBFrame}}
#' @examples 
#' # Set up example MizerParams and MizerSim objects
#' data(NS_species_params_gears)
#' data(inter)
#' params <- suppressMessages(MizerParams(NS_species_params_gears, inter))
#' sim0 <- project(params, effort=0, t_max=20, progress_bar = FALSE)
#' sim1 <- project(params, effort=1, t_max=20, progress_bar = FALSE)
#' 
#' # Display biomass from each simulation next to each other
#' display_frames(getBiomassFrame(sim0), getBiomassFrame(sim1), params)
display_frames <- function(f1, f2, params, 
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
#' @return A data frame that can be used in \code{\link{display_frames}}
#' @export
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
                 (is.na(ylim[2]) | bm$value <= ylim[1]), ]
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
#' \code{\link{display_frames}}.
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
#'   all species is plotted as well. Default is FALSE
#' @param ... Other arguments to pass to \code{getBiomass} method, for example
#'   \code{min_w} and \code{max_w}
#'   
#' @return A data frame that can be used in \code{\link{display_frames}}
#' @export
#' @seealso \code{\link{getBiomass}}, \code{\link{display_frames}},
#'   \code{\link{getSSBFrame}}
#' @examples 
#' # Set up example MizerParams and MizerSim objects
#' data(NS_species_params_gears)
#' data(inter)
#' params <- suppressMessages(MizerParams(NS_species_params_gears, inter))
#' sim0 <- project(params, effort=0, t_max=20, progress_bar = FALSE)
#' sim1 <- project(params, effort=1, t_max=20, progress_bar = FALSE)
#' 
#' # Display biomass from each simulation next to each other
#' display_frames(getBiomassFrame(sim0), getBiomassFrame(sim1), params)
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
                 (is.na(ylim[2]) | bm$value <= ylim[1]), ]
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
#' @param y_ticks The approximate number of ticks desired on the y axis
#' @param print_it Display the plot, or just return the ggplot2 object. Default
#'   value is FALSE
#' @param background A boolean value that determines whether background species
#'   are included. Ignored if the model does not contain background species.
#'   Default is TRUE.
#' @param ... Other arguments to pass to \code{getBiomass} method, for example
#'   \code{min_w} and \code{max_w}
#'   
#' @return A ggplot2 object
#' @export
#' @seealso \code{\link{getBiomass}}
#' @examples
#' # Set up example MizerParams and MizerSim objects
#' data(NS_species_params_gears)
#' data(inter)
#' params <- suppressMessages(MizerParams(NS_species_params_gears, inter))
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
            y_ticks = 6, print_it = FALSE,
            ylim = c(NA, NA),
            total = FALSE, background = TRUE, ...) {
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
                               colour = "lightgrey")
        }
    }

    if ( (length(species) + total) > 12) {
        p <- p + geom_line(aes(group = Species))
    } else {
        p <- p +
            geom_line(aes(colour = Species, linetype = Species))
    }
    if (print_it) {
        print(p)
    }
    return(p)
}

#' Plot the biomass of species against time with plotly
#' 
#' @inherit plotBiomass params return description details seealso
#' @export
plotlyBiomass <- function(sim,
             species = dimnames(sim@n)$sp[!is.na(sim@params@A)],
             start_time = as.numeric(dimnames(sim@n)[[1]][1]),
             end_time = as.numeric(dimnames(sim@n)[[1]][dim(sim@n)[1]]),
             y_ticks = 6,
             print_it = FALSE,
             ylim = c(NA, NA),
             total = FALSE,
             background = TRUE,
             ...) {
    argg <- c(as.list(environment()), list(...))
    ggplotly(do.call("plotBiomass", argg))
    # ggplotly(do.callplotBiomass(sim, species, start_time, end_time, y_ticks,
    #                      print_it, ylim, total, background, ...))
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
#' @param species Name or vector of names of the species to be plotted. By
#'   default all species contained in \code{sim} are plotted.
#' @param print_it Display the plot, or just return the ggplot2 object.
#'   Defaults to FALSE.
#' @param total A boolean value that determines whether the total yield from
#'   all species in the system is plotted as well. Default is FALSE.
#' @param log Boolean whether yield should be plotted on a logarithmic axis. 
#'   Defaults to true.
#' @param ... Other arguments to pass to \code{\link{getYield}} method.
#'
#' @return A ggplot2 object
#' @export
#' @seealso \code{\link{getYield}}
#' @examples
#' 
#' data(NS_species_params_gears)
#' data(inter)
#' params <- suppressMessages(MizerParams(NS_species_params_gears, inter))
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
                      print_it = FALSE, total = FALSE, log = TRUE, ...){
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
        if (length(species) > 12) {
            p <- ggplot(ym) +
                geom_line(aes(x = Year, y = Yield, group = Species))
        } else {
            p <- ggplot(ym) +
                geom_line(aes(x = Year, y = Yield,
                              colour = Species, linetype = Species))
        }
        if (log) {
            p <- p + scale_y_continuous(trans = "log10", name = "Yield [g/year]",
                                        breaks = log_breaks(),
                                        labels = prettyNum)
        } else {
            p <- p + scale_y_continuous(name = "Yield [g/year]")
        }
        p <- p +
            scale_colour_manual(values = sim@params@linecolour) +
            scale_linetype_manual(values = sim@params@linetype)
        if (print_it) {
            print(p)
        }
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
        if (nlevels(ym$Species) > 12) {
            p <- ggplot(ym) +
                geom_line(aes(x = Year, y = Yield, group = Species))
        } else {
            p <- ggplot(ym) +
                geom_line(aes(x = Year, y = Yield, colour = Species,
                              linetype = Species))
        }
        if (log) {
            p <- p + scale_y_continuous(trans = "log10", name = "Yield [g/year]")
        } else {
            p <- p + scale_y_continuous(name = "Yield [g/year]")
        }
        p <- p + facet_wrap(~ Simulation)
        if (print_it) {
            print(p)
        }
        return(p)
    }
}

#' Plot the total yield of species through time with plotly
#' @inherit plotYield params return description details seealso
#' @export
plotlyYield <- function(sim, sim2,
                      species = dimnames(sim@n)$sp,
                      print_it = FALSE, total = FALSE, log = TRUE, ...) {
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
#' @param species Name or vector of names of the species to be plotted. By
#'   default all species are plotted.
#' @param print_it Display the plot, or just return the ggplot2 object. 
#'   Defaults to FALSE
#' @param total A boolean value that determines whether the total yield
#'   per gear over all species in the system is plotted as well. Default is FALSE
#' @param ... Other arguments to pass to \code{\link{getYieldGear}} method
#'
#' @return A ggplot2 object
#' @export
#' @seealso \code{\link{getYieldGear}}
#' @examples
#' 
#' data(NS_species_params_gears)
#' data(inter)
#' params <- suppressMessages(MizerParams(NS_species_params_gears, inter))
#' sim <- project(params, effort=1, t_max=20, t_save = 0.2, progress_bar = FALSE)
#' plotYieldGear(sim)
#' plotYieldGear(sim, species = c("Cod", "Herring"), total = TRUE)
#' 
plotYieldGear <- function(sim,
                          species = dimnames(sim@n)$sp,
                          print_it = FALSE, total = FALSE, ...){
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
    if (length(species) > 12) {
        p <- ggplot(ym) + geom_line(aes(x = time, y = value, group = Species))
    } else {
        p <- ggplot(ym) +
            geom_line(aes(x = time, y = value, colour = Species, linetype = gear))
    }
    p <- p + scale_y_continuous(trans = "log10", name = "Yield [g]") +
        scale_x_continuous(name = "Year") +
        scale_colour_manual(values = sim@params@linecolour)
    if (print_it) {
        print(p)
    }
    return(p)
}

#' Plot the total yield of each species by gear through time with plotly
#' @inherit plotYieldGear params return description details seealso
#' @export
plotlyYieldGear <- function(sim,
                          species = dimnames(sim@n)$sp,
                          print_it = FALSE, total = FALSE, ...) {
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
#' @param min_w Minimum weight to be plotted (useful for truncating the
#'   plankton spectrum). Default value is a hundredth of the minimum size
#'   value of the community.
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
#' @param print_it Display the plot, or just return the ggplot2 object.
#'   Defaults to FALSE
#' @param total A boolean value that determines whether the total over all
#'   species in the system is plotted as well. Default is FALSE
#' @param plankton A boolean value that determines whether plankton is included.
#'   Default is TRUE.
#' @param background A boolean value that determines whether background species
#'   are included. Ignored if the model does not contain background species.
#'   Default is TRUE.
#' @param ... Other arguments (currently unused)
#'   
#' @return A ggplot2 object
#' @export
#' @examples
#' 
#' data(NS_species_params_gears)
#' data(inter)
#' params <- suppressMessages(MizerParams(NS_species_params_gears, inter))
#' sim <- project(params, effort=1, t_max=20, t_save = 2, progress_bar = FALSE)
#' plotSpectra(sim)
#' plotSpectra(sim, min_w = 1e-6)
#' plotSpectra(sim, time_range = 10:20)
#' plotSpectra(sim, time_range = 10:20, power = 0)
#' plotSpectra(sim, species = c("Cod", "Herring"), power = 1)
#' 
plotSpectra <- function(object, species = NULL,
                        time_range,
                        min_w, ylim = c(NA, NA),
                        power = 1, biomass = TRUE, print_it = FALSE,
                        total = FALSE, plankton = TRUE, 
                        background = TRUE, ...) {
    if (is(object, "MizerSim")) {
        if (missing(time_range)){
            time_range  <- max(as.numeric(dimnames(object@n)$time))
        }
        if (missing(min_w)){
            min_w <- min(object@params@w) / 100
        }
        # to deal with old-type biomass argument
        if (missing(power)) {
            power <- as.numeric(biomass)
        }
        time_elements <- get_time_elements(object,time_range)
        n <- apply(object@n[time_elements, , ,drop = FALSE], c(2, 3), mean)
        n_pp <- apply(object@n_pp[time_elements,,drop = FALSE], 2, mean)
        ps <- plot_spectra(object@params, n = n, n_pp = n_pp,
                           species = species, min_w = min_w, ylim = ylim,
                           power = power, print_it = print_it,
                           total = total, plankton = plankton,
                           background = background)
        return(ps)
    } else {
        if (missing(power)) {
            power <- as.numeric(biomass)
        }
        if (missing(min_w)) {
            min_w <- min(object@w) / 100
        }
        ps <- plot_spectra(object, n = object@initial_n,
                           n_pp = object@initial_n_pp,
                           species = species, min_w = min_w, ylim = ylim,
                           power = power, print_it = print_it,
                           total = total, plankton = plankton,
                           background = background)
        return(ps)
    }
}


plot_spectra <- function(params, n, n_pp,
                         species, min_w, ylim, power, print_it,
                         total, plankton, background) {
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
    # Deal with power argument
    if (power %in% c(0, 1, 2)) {
        y_label <- c("Number density [1/g]", "Biomass density",
                    "Biomass density [g]")[power + 1]
    } else {
        y_label <- paste0("Number density * w^", power)
    }
    n <- sweep(n, 2, params@w^power, "*")
    # Select only the desired species and background species
    spec_n <- n[as.character(dimnames(n)[[1]]) %in% species, , drop = FALSE]
    # Make data.frame for plot
    plot_dat <- data.frame(value = c(spec_n),
                           # ordering of factor is important for legend
                           Species = factor(dimnames(spec_n)[[1]],
                                            levels = species_levels),
                           w = rep(params@w,
                                   each = dim(spec_n)[[1]]))
    if (plankton) {
        # Decide where to cut off plankton
        max_w <- min(params@species_params$w_mat)
        if (is.na(max_w)) {
            max_w <- Inf
        }
        plankton_sel <- params@w_full >= min_w &
            params@w_full < max_w
        w_plankton <- params@w_full[plankton_sel]
        plank_n <- n_pp[plankton_sel] * w_plankton^power
        plot_dat <- rbind(plot_dat,
                          data.frame(value = c(plank_n),
                                     Species = "Plankton",
                                     w = w_plankton))
    }
    if (total) {
        plot_dat <- rbind(plot_dat,
                          data.frame(value = c(total_n),
                                     Species = "Total",
                                     w = params@w_full))
    }
    # lop off 0s and apply min_w
    plot_dat <- plot_dat[(plot_dat$value > 0) & (plot_dat$w >= min_w), ]
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
        # lop off 0s and apply min_w
        plot_back <- plot_back[(plot_back$value > 0) & (plot_back$w >= min_w), ]
        # Impose ylim
        if (!is.na(ylim[1])) {
            plot_back <- plot_back[plot_back$value < ylim[1], ]
        }
        plot_back <- plot_back[plot_back$value > ylim[2], ]
        if (nrow(plot_back) > 0) {
            # Add background species in grey
            p <- p +
                geom_line(aes(group = Species), colour = "grey",
                          data = plot_back)
        }
    }
    if ( (length(species) + plankton + total) > 13) {
        p <- p + geom_line(aes(group = Species))
    } else {
        p <- p + geom_line(aes(colour = Species, linetype = Species))
    }
    if (print_it)
        print(p)
    return(p)
}

#' Plotly plot of the abundance spectra
#' @inherit plotSpectra params return description details seealso
#' @export
plotlySpectra <- function(object, species = NULL,
                        time_range,
                        min_w, ylim = c(NA, NA),
                        power = 1, biomass = TRUE, print_it = FALSE,
                        total = FALSE, plankton = TRUE, 
                        background = TRUE, ...) {
    argg <- as.list(environment())
    ggplotly(do.call("plotSpectra", argg))
}


#' Plot the feeding level of species by size
#' 
#' After running a projection, plot the feeding level of each species by size. 
#' The feeding level is averaged over the specified time range (a single value
#' for the time range can be used).
#' 
#' @param sim An object of class \linkS4class{MizerSim}.
#' @param species Name or vector of names of the species to be plotted. By
#'   default all species are plotted.
#' @param time_range The time range (either a vector of values, a vector of min
#'   and max time, or a single value) to average the abundances over. Default is
#'   the final time step.
#' @param print_it Display the plot, or just return the ggplot2 object.
#'   Defaults to FALSE
#' @param ... Other arguments to pass to \code{getFeedingLevel} method
#'
#' @return A ggplot2 object
#' @export
#' @seealso \code{\link{getFeedingLevel}}
#' @examples
#' 
#' data(NS_species_params_gears)
#' data(inter)
#' params <- suppressMessages(MizerParams(NS_species_params_gears, inter))
#' sim <- project(params, effort=1, t_max=20, t_save = 2, progress_bar = FALSE)
#' plotFeedingLevel(sim)
#' plotFeedingLevel(sim, time_range = 10:20, species = c("Cod", "Herring"))
#' 
plotFeedingLevel <- function(sim,
            species = dimnames(sim@n)$sp,
            time_range = max(as.numeric(dimnames(sim@n)$time)),
            print_it = FALSE, ...) {
    feed_time <- getFeedingLevel(sim, time_range = time_range,
                                 drop = FALSE, ...)
    feed <- apply(feed_time, c(2, 3), mean)
    # selector for desired species
    sel_sp <- as.character(dimnames(feed)$sp) %in% species
    feed <- feed[sel_sp, , drop = FALSE]
    plot_dat <- data.frame(value = c(feed),
                           # ggplot orders the legend according to the ordering
                           # of the factors, hence we need the levels argument
                           Species = factor(dimnames(feed)$sp, 
                                            levels = dimnames(feed)$sp),
                           w = rep(sim@params@w, each = length(species)))
    if (length(species) > 12) {
        p <- ggplot(plot_dat) +
            geom_line(aes(x = w, y = value, group = Species))
    } else {
        p <- ggplot(plot_dat) +
            geom_line(aes(x = w, y = value, colour = Species, linetype = Species))
    }
    p <- p +
        scale_x_continuous(name = "Size [g]", trans = "log10") +
        scale_y_continuous(name = "Feeding Level", limits = c(0, 1)) +
        scale_colour_manual(values = sim@params@linecolour) +
        scale_linetype_manual(values = sim@params@linetype)
    if (print_it) {
        print(p)
    }
    return(p)
}

#' Plot the feeding level of species by size with plotly
#' 
#' @inherit plotFeedingLevel params return description details seealso
#' @export
plotlyFeedingLevel <- function(sim,
                             species = dimnames(sim@n)$sp,
                             time_range = max(as.numeric(dimnames(sim@n)$time)),
                             print_it = FALSE, ...) {
    argg <- as.list(environment())
    ggplotly(do.call("plotFeedingLevel", argg))
}
    

#' Plot predation mortality rate of each species against size
#' 
#' After running a projection, plot the predation mortality rate of each species
#' by size. The mortality rate is averaged over the specified time range (a
#' single value for the time range can be used to plot a single time step).
#' 
#' @param sim An object of class \linkS4class{MizerSim}
#' @param species Name or vector of names of the species to be plotted. By
#'   default all species are plotted.
#' @param time_range The time range (either a vector of values, a vector of min
#'   and max time, or a single value) to average the abundances over. Default is
#'   the final time step.
#' @param print_it Display the plot, or just return the ggplot2 object.
#'   Defaults to FALSE
#' @param ... Other arguments to pass to \code{getM2} method.
#'
#' @return A ggplot2 object
#' @export
#' @seealso \code{\link{getM2}}
#' @examples
#' 
#' data(NS_species_params_gears)
#' data(inter)
#' params <- suppressMessages(MizerParams(NS_species_params_gears, inter))
#' sim <- project(params, effort=1, t_max=20, t_save = 2, progress_bar = FALSE)
#' plotM2(sim)
#' plotM2(sim, time_range = 10:20)
#' 
plotM2 <- function(sim, species = dimnames(sim@n)$sp,
                   time_range = max(as.numeric(dimnames(sim@n)$time)),
                   print_it = FALSE, ...) {
    # Need to keep species in order for legend
    species_levels <- c(dimnames(sim@n)$sp, "Background", "Plankton", "Total")
    
    m2_time <- getM2(sim, time_range = time_range, drop = FALSE, ...)
    m2 <- apply(m2_time, c(2, 3), mean)
    m2 <- m2[as.character(dimnames(m2)[[1]]) %in% species, , 
             drop = FALSE]
    plot_dat <- data.frame(value = c(m2),
                           Species = factor(dimnames(m2)[[1]],
                                            levels = species_levels),
                           w = rep(sim@params@w, each = length(species)))
    if (length(species) > 12) {
        p <- ggplot(plot_dat) +
            geom_line(aes(x = w, y = value, group = Species))
    } else {
        p <- ggplot(plot_dat) +
            geom_line(aes(x = w, y = value, colour = Species, linetype = Species))
    }
    p <- p +
        scale_x_continuous(name = "Size [g]", trans = "log10") +
        scale_y_continuous(name = "Predation mortality [1/year]",
                           limits = c(0, max(plot_dat$value))) +
        scale_colour_manual(values = sim@params@linecolour) +
        scale_linetype_manual(values = sim@params@linetype)
    if (print_it) {
        print(p)
    }
    return(p)
}

#' Plot predation mortality rate of each species against size with plotly
#' @inherit plotM2 params return description details seealso
#' @export
plotlyM2 <- function(sim, species = dimnames(sim@n)$sp,
                   time_range = max(as.numeric(dimnames(sim@n)$time)),
                   print_it = FALSE, ...) {
    argg <- as.list(environment())
    ggplotly(do.call("plotM2", argg))
}

#' Plot total fishing mortality of each species by size
#' 
#' After running a projection, plot the total fishing mortality of each species
#' by size. The total fishing mortality is averaged over the specified time
#' range (a single value for the time range can be used to plot a single time
#' step).
#' 
#' @param sim An object of class \linkS4class{MizerSim}.
#' @param species Name or vector of names of the species to be plotted. By
#'   default all species are plotted.
#' @param time_range The time range (either a vector of values, a vector of min
#'   and max time, or a single value) to average the abundances over. Default is
#'   the final time step.
#' @param print_it Display the plot, or just return the ggplot2 object.
#'   Defaults to FALSE
#' @param ... Other arguments to pass to \code{getFMort} method
#'
#' @return A ggplot2 object
#' @export
#' @seealso \code{\link{getFMort}}
#' @examples
#' 
#' data(NS_species_params_gears)
#' data(inter)
#' params <- suppressMessages(MizerParams(NS_species_params_gears, inter))
#' sim <- project(params, effort=1, t_max=20, t_save = 2, progress_bar = FALSE)
#' plotFMort(sim)
#' plotFMort(sim, time_range = 10:20)
#' 
plotFMort <- function(sim, species = dimnames(sim@n)$sp,
                      time_range = max(as.numeric(dimnames(sim@n)$time)),
                      print_it = FALSE, ...) {
    # Need to keep species in order for legend
    species_levels <- c(dimnames(sim@n)$sp, "Background", "Plankton", "Total")
    
    f_time <- getFMort(sim, time_range = time_range, drop = FALSE, ...)
    f <- apply(f_time, c(2, 3), mean)
    f <- f[as.character(dimnames(f)[[1]]) %in% species, , drop = FALSE]
    plot_dat <- data.frame(value = c(f),
                           Species = factor(dimnames(f)[[1]],
                                            levels = species_levels),
                           w = rep(sim@params@w, each = length(species)))
    if (length(species) > 12) {
        p <- ggplot(plot_dat) + geom_line(aes(x = w, y = value, group = Species))
    } else {
        p <- ggplot(plot_dat) +
            geom_line(aes(x = w, y = value, colour = Species, linetype = Species))
    }
    p <- p +
        scale_x_continuous(name = "Size [g]", trans = "log10") +
        scale_y_continuous(name = "Fishing mortality [1/Year]",
                           limits = c(0, max(plot_dat$value))) +
        scale_colour_manual(values = sim@params@linecolour) +
        scale_linetype_manual(values = sim@params@linetype)
    if (print_it) {
        print(p)
    }
    return(p)
}

#' Plot total fishing mortality of each species by size with plotly
#' @inherit plotM2 params return description details seealso
#' @export
plotlyFMort <- function(sim, species = dimnames(sim@n)$sp,
                      time_range = max(as.numeric(dimnames(sim@n)$time)),
                      print_it = FALSE, ...) {
    argg <- as.list(environment())
    ggplotly(do.call("plotFMort", argg))
}

#' Get growth curves giving weight as a function of age
#' 
#' If given a \linkS4class{MizerSim} object, uses the growth rates at the final
#' time of a simulation to calculate the size at age. If given a
#' \linkS4class{MizerParams} object, uses the initial growth rates instead.
#' 
#' @param object MizerSim or MizerParams object
#' @param max_age The age up to which the weight is to be plotted. Default is 20
#'
#' @return An array (species x age) containing the weight in grams.
#' @export
#' @examples
#' 
#' data(NS_species_params_gears)
#' data(inter)
#' params <- suppressMessages(MizerParams(NS_species_params_gears, inter))
#' getGrowthCurves(params)
#' sim <- project(params, effort=1, t_max=20, t_save = 2, progress_bar = FALSE)
#' getGrowthCurves(sim, max_age = 24)
#' 
getGrowthCurves <- function(object, max_age = 24) {
    if (is(object, "MizerSim")) {
        params <- object@params
        t <- dim(object@n)[1]
        n <- object@n[t, , ]
        n_pp <- object@n_pp[t, ]
    } else if (is(object, "MizerParams")) {
        params <- object
        n <- object@initial_n
        n_pp <- object@initial_n_pp
    }
    species <- dimnames(n)$sp
    age <- seq(0, max_age, length.out = 50)
    ws <- array(dim = c(length(species), length(age)),
                dimnames = list(Species = species, Age = age))
    g <- getEGrowth(params, n, n_pp)
    for (i in 1:length(species)) {
        g_fn <- stats::approxfun(params@w, g[i, ])
        myodefun <- function(t, state, parameters){
            return(list(g_fn(state)))
        }
        ws[i, ] <- deSolve::ode(y = params@species_params$w_min[i], 
                                times = age, func = myodefun)[, 2]
    }
    return(ws)
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
#' @param object MizerSim or MizerParams object
#' @param species Name or vector of names of the species to be plotted. By
#'   default all species are plotted.
#' @param max_age The age up to which the weight is to be plotted. Default is 20
#' @param percentage Boolean value. If TRUE, the size is shown as a percentage
#'   of the maximal size.
#' @param print_it Display the plot, or just return the ggplot2 object.
#'   Defaults to FALSE
#' @param ... Other arguments (unused)
#' 
#' @return A ggplot2 object
#' @export
#' @examples
#' 
#' data(NS_species_params_gears)
#' data(inter)
#' params <- suppressMessages(MizerParams(NS_species_params_gears, inter))
#' sim <- project(params, effort=1, t_max=20, t_save = 2, progress_bar = FALSE)
#' plotGrowthCurves(sim, percentage = TRUE)
#' plotGrowthCurves(sim, species = "Cod", max_age = 24)
#' 
plotGrowthCurves <- function(object, species,
                             max_age = 20,
                             percentage = FALSE,
                             print_it = FALSE) {
    if (is(object, "MizerSim")) {
        params <- object@params
        t <- dim(object@n)[1]
        n <- object@n[t, , ]
        n_pp <- object@n_pp[t, ]
    } else if (is(object, "MizerParams")) {
        params <- object
        n <- object@initial_n
        n_pp <- object@initial_n_pp
    }
    if (missing(species)) {
        species <- dimnames(n)$sp
    }
    # reorder list of species to coincide with order in params
    idx <- which(dimnames(n)$sp %in% species)
    species <- dimnames(n)$sp[idx]
    age <- seq(0, max_age, length.out = 50)
    ws <- array(dim = c(length(species), length(age)),
                dimnames = list(Species = species, Age = age))
    g <- getEGrowth(params, n, n_pp)
    for (j in 1:length(species)) {
        i <- idx[j]
        g_fn <- stats::approxfun(params@w, g[i, ])
        myodefun <- function(t, state, parameters){
            return(list(g_fn(state)))
        }
        ws[j, ] <- deSolve::ode(y = params@species_params$w_min[i], 
                                times = age, func = myodefun)[, 2]
        if (percentage) {
            ws[j, ] <- ws[j, ] / params@species_params$w_inf[i] * 100
        }
    }
    plot_dat <- reshape2::melt(ws)
    # Need to keep species in order for legend
    plot_dat$Species <- factor(plot_dat$Species, dimnames(n)$sp)
    if (length(species) > 12) {
        p <- ggplot(plot_dat) +
            geom_line(aes(x = Age, y = value, group = Species))
    } else {
        p <- ggplot(plot_dat) +
            geom_line(aes(x = Age, y = value,
                          colour = Species, linetype = Species))
    }
    y_label <- if (percentage) "Percent of maximum size" else "Size [g]"
    p <- p +
        scale_x_continuous(name = "Age [Years]") +
        scale_y_continuous(name = y_label) +
        scale_colour_manual(values = params@linecolour) +
        scale_linetype_manual(values = params@linetype)
    
    # Extra stuff for single-species case
    if (length(species) == 1 && !percentage) {
        w_inf <- params@species_params$w_inf[idx[1]]
        p <- p + geom_hline(yintercept = w_inf) +
            annotate("text", 0, w_inf, vjust = -1, label = "Maximum")
        w_mat <- params@species_params$w_mat[idx[1]]
        p <- p + geom_hline(yintercept = w_mat) +
            annotate("text", 0, w_mat, vjust = -1, label = "Maturity")
        if (all(c("a", "b", "k_vb") %in% names(params@species_params))) {
            a <- params@species_params$a[idx[1]]
            b <- params@species_params$b[idx[1]]
            k_vb <- params@species_params$k_vb[idx[1]]
            t0 <- params@species_params$t0[idx[1]]
            if (is.null(t0)) t0 <- 0
            L_inf <- (w_inf/a)^(1/b)
            vb <- a * (L_inf * (1 - exp(-k_vb * (age - t0))))^b
            dat <- data.frame(x = age, y = vb)
            p <- p + geom_line(data = dat, aes(x = x, y = y))
        }
    }
    
    if (print_it) {
        print(p)
    }
    return(p)
}

#' Plot growth curves giving weight as a function of age with plotly
#' @inherit plotM2 params return description details seealso
#' @export
plotlyGrowthCurves <- function(object, species,
                             max_age = 20,
                             percentage = FALSE,
                             print_it = FALSE) {
    argg <- as.list(environment())
    ggplotly(do.call("plotGrowthCurves", argg))
}


#### plot ####
#' Summary plot for \code{MizerSim} objects
#' 
#' After running a projection, produces 5 plots in the same window: feeding
#' level, abundance spectra, predation mortality and fishing mortality of each
#' species by size; and biomass of each species through time. This method just
#' uses the other plotting methods and puts them all in one window.
#' 
#' @param x An object of class \linkS4class{MizerSim}
#' @param y Not used
#' @param ...  For additional arguments see the documentation for
#'   \code{\link{plotBiomass}},
#'   \code{\link{plotFeedingLevel}},\code{\link{plotSpectra}},\code{\link{plotM2}}
#'   and \code{\link{plotFMort}}.
#' @return A viewport object
#' @export
#' @rdname plotMizerSim
#' @examples
#' 
#' data(NS_species_params_gears)
#' data(inter)
#' params <- suppressMessages(MizerParams(NS_species_params_gears, inter))
#' sim <- project(params, effort=1, t_max=20, t_save = 2, progress_bar = FALSE)
#' plot(sim)
#' plot(sim, time_range = 10:20) # change time period for size-based plots
#' plot(sim, min_w = 10, max_w = 1000) # change size range for biomass plot
#' 
setMethod("plot", signature(x = "MizerSim", y = "missing"),
          function(x, ...) {
              p1 <- plotFeedingLevel(x, print_it = FALSE, ...)
              p2 <- plotSpectra(x, print_it = FALSE, ...)
              p3 <- plotBiomass(x, y_ticks = 3, print_it = FALSE, ...)
              p4 <- plotM2(x, print_it = FALSE, ...)
              p5 <- plotFMort(x, print_it = FALSE, ...)
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
