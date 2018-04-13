# Plotting methods for the MizerSim class

# Copyright 2012 Finlay Scott and Julia Blanchard.
# Copyright 2018 Gustav Delius and Richard Southwell.
# Development has received funding from the European Commission’s Horizon 2020 
# Research and Innovation Programme under Grant Agreement No. 634495 
# for the project MINOUW (http://minouw-project.eu/).
# Distributed under the GPL 3 or later 
# Maintainer: Gustav Delius, University of York, <gustav.delius@york.ac.uk>

# Hackiness to get past the 'no visible binding ... ' warning when running check
 utils::globalVariables(c("time", "value", "Species", "w", "gear"))
 
 #' Helper function to produce nice breaks on logarithmic axes
 #'
 #' This is needed when the logarithmic y-axis spans less than one order of
 #' magnitude, in which case the ggplot2 default produces no ticks.
 #' Thanks to Heather Turner at
 #' https://stackoverflow.com/questions/14255533/pretty-ticks-for-log-normal-scale-using-ggplot2-dynamic-not-manual
 #'
 #' @param n Approximate number of ticks
 #'
 #' @return A function that can be used as the break argument in calls to
 #'     scale_y_continuous() or scale_x_continuous()
log_breaks <- function(n = 6){
   n = max(1, n)  # Because n=0 could lead to R crash
   function(x) {
     grDevices::axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, nint = n)
   }
}

#' Plot the biomass of species through time
#'
#' After running a projection, the biomass of each species can be plotted
#' against time. The biomass is calculated within user defined size limits (see
#' \code{\link{getBiomass}}). 
#' 
#' This plot is pretty easy to do by hand. It just
#' gets the biomass using the \code{\link{getBiomass}} method and plots using
#' the ggplot2 package. You can then fiddle about with colours and linetypes
#' etc. Just look at the source code for details.
#' 
#' @param sim An object of class \code{MizerSim}.
#' @param species Name or vector of names of the species to be plotted. By
#'   default all species are plotted.
#' @param start_time The first time to be plotted. Default is the beginning
#'   of the time series.
#' @param end_time The last time to be plotted. Default is the end of the
#'   time series.
#' @param y_ticks The approximate number of ticks desired on the y axis
#' @param print_it Display the plot, or just return the ggplot2 object. Default
#'   value is TRUE
#' @param total A boolean value that determines whether the total biomass from
#'   all species is plotted as well. Default is FALSE
#' @param ... Other arguments to pass to \code{getBiomass} method, for example
#'   \code{min_w} and \code{max_w}
#'   
#' @return A ggplot2 object
#' @export
#' @seealso \code{\link{getBiomass}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=20, t_save = 0.2)
#' plotBiomass(sim)
#' plotBiomass(sim, species = c("Cod", "Herring"), total = TRUE)
#' plotBiomass(sim, min_w = 10, max_w = 1000)
#' plotBiomass(sim, start_time = 10, end_time = 15)
#' plotBiomass(sim, y_ticks = 3)
#' }
setGeneric('plotBiomass', function(sim, ...)
    standardGeneric('plotBiomass'))

#' Plot the biomass using a \code{MizerSim} object.
#' @rdname plotBiomass
setMethod('plotBiomass', signature(sim='MizerSim'),
    function(sim, species = as.character(sim@params@species_params$species),
             start_time = as.numeric(dimnames(sim@n)[[1]][1]), 
             end_time = as.numeric(dimnames(sim@n)[[1]][dim(sim@n)[1]]),
             y_ticks = 6, print_it = TRUE, total = FALSE, ...){
        b <- getBiomass(sim, ...)
        if(start_time >= end_time){
            stop("start_time must be less than end_time")
        }
        # Select time range
        b <- b[(as.numeric(dimnames(b)[[1]]) >= start_time) & 
                   (as.numeric(dimnames(b)[[1]]) <= end_time), , drop = FALSE]
        b_total = rowSums(b)
        # Select species
        b <- b[, as.character(dimnames(b)[[2]]) %in% species, drop = FALSE]
        # Include total
        if (total) {
            b <- cbind(b, Total = b_total)
        }
        names(dimnames(b)) <- c("time", "Species")
        bm <- reshape2::melt(b)
        # Force Species column to be a character (if numbers used - may be
        # interpreted as integer and hence continuous)
        bm$Species <- as.character(bm$Species)
        # Due to log10, need to set a minimum value, seems like a feature in ggplot
        min_value <- 1e-30
        bm <- bm[bm$value >= min_value,]
        x_label <- "Year"
        y_label <- "Biomass [g]"
        if (length(species) > 12) {
            p <- ggplot(bm) + geom_line(aes(x=time, y=value, group=Species)) 
        } else {
            p <- ggplot(bm) + 
                geom_line(aes(x=time, y=value, colour=Species, linetype=Species))
        }
        p <- p  + 
            scale_y_continuous(trans="log10", breaks=log_breaks(n=y_ticks), 
                               labels = prettyNum, name=y_label) + 
            scale_x_continuous(name=x_label)
        if (print_it) {
            print(p)
        }
        return(p)
    }
)

#' Plot the total yield of species through time
#'
#' After running a projection, the total yield of each species across all 
#' fishing gears can be plotted against time. 
#' 
#' This plot is pretty easy to do by
#' hand. It just gets the biomass using the \code{\link{getYield}} method and
#' plots using the ggplot2 package. You can then fiddle about with colours and
#' linetypes etc. Just look at the source code for details.
#' 
#' @param sim An object of class \code{MizerSim}
#' @param species Name or vector of names of the species to be plotted. By
#'   default all species are plotted.
#' @param print_it Display the plot, or just return the ggplot2 object.
#'   Defaults to TRUE
#' @param total A boolean value that determines whether the total Yield from
#'   all species in the system is plotted as well. Default is FALSE
#' @param ... Other arguments to pass to \code{getYield} method
#'
#' @return A ggplot2 object
#' @export
#' @seealso \code{\link{getYield}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=20, t_save = 0.2)
#' plotYield(sim)
#' plotYield(sim, species = c("Cod", "Herring"), total = TRUE)
#' }
setGeneric('plotYield', function(sim, ...)
    standardGeneric('plotYield'))

#' Plot the yield using a \code{MizerSim} object.
#' @rdname plotYield
setMethod('plotYield', signature(sim='MizerSim'),
    function(sim, species = as.character(sim@params@species_params$species),
             print_it = TRUE, total = FALSE, ...){
        y <- getYield(sim, ...)
        y_total <- rowSums(y)
        y <- y[, (as.character(dimnames(y)[[2]]) %in% species) & colSums(y)>0, 
               drop=FALSE]
        if (total) {
            # Include total
            y <- cbind(y, Total = y_total)
        }
        names(dimnames(y)) <- c("time", "Species")
        ym <- reshape2::melt(y)
        ym$Species <- as.character(ym$Species)
        if (dim(y)[2] > 12) {
            p <- ggplot(ym) + 
                geom_line(aes(x=time,y=value, group=Species))
        } else {
            p <- ggplot(ym) + 
                geom_line(aes(x=time,y=value, colour=Species, linetype=Species))
        }
        p <- p + scale_y_continuous(trans="log10", name="Yield [g]") + 
            scale_x_continuous(name="Year")
    if (print_it)
        print(p)
	return(p)
    }
)

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
#' @param sim An object of class \code{MizerSim}
#' @param species Name or vector of names of the species to be plotted. By
#'   default all species are plotted.
#' @param print_it Display the plot, or just return the ggplot2 object. 
#'   Defaults to TRUE
#' @param total A boolean value that determines whether the total yield
#'   per gear over all species in the system is plotted as well. Default is FALSE
#' @param ... Other arguments to pass to \code{getYieldGear} method
#'
#' @return A ggplot2 object
#' @export
#' @seealso \code{\link{getYieldGear}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=20, t_save = 0.2)
#' plotYieldGear(sim)
#' plotYieldGear(sim, species = c("Cod", "Herring"), total = TRUE)
#' }
setGeneric('plotYieldGear', function(sim, ...)
    standardGeneric('plotYieldGear'))

#' Plot the yield of each gear using a \code{MizerSim} object.
#' @rdname plotYieldGear
setMethod('plotYieldGear', signature(sim='MizerSim'),
    function(sim, species = as.character(sim@params@species_params$species),
             print_it=TRUE, total = FALSE, ...){
	y <- getYieldGear(sim, ...)
	y_total <- rowSums(y, dims = 2)
	y <- y[, , as.character(dimnames(y)[[3]]) %in% species,
	       drop=FALSE]
	names(dimnames(y))[names(dimnames(y))=="sp"] <- "Species"
	ym <- reshape2::melt(y)
	if (total) {
	    yt <- reshape2::melt(y_total)
	    yt$Species <- "Total"
	    ym <- rbind(ym, yt)
	}
	ym <- subset(ym, ym$value > 0)
    if (length(species) > 12) {
        p <- ggplot(ym) + geom_line(aes(x=time,y=value, group=Species))
    } else {
        p <- ggplot(ym) + 
            geom_line(aes(x=time,y=value, colour=Species, linetype=gear))
    }
	p <- p + scale_y_continuous(trans="log10", name="Yield [g]") + 
	    scale_x_continuous(name="Year")
    if (print_it) {
        print(p)
    }
	return(p)
    }
)

#' Plot the abundance spectra of species and the background population
#' 
#' After running a projection, the spectra of the abundance of each species and
#' the background population can be plotted. The abundance is averaged over the
#' specified time range (a single value for the time range can be used to plot a
#' single time step). The abundance can be in terms of numbers or biomass,
#' depending on the \code{biomass} argument.
#' 
#' @param sim An object of class \code{MizerSim}.
#' @param species Name or vector of names of the species to be plotted. By
#'   default all species are plotted.
#' @param time_range The time range (either a vector of values, a vector of min
#'   and max time, or a single value) to average the abundances over. Default is
#'   the final time step.
#' @param min_w Minimum weight to be plotted (useful for truncating the
#'   background spectrum). Default value is a hundredth of the minimum size
#'   value of the community.
#' @param biomass A boolean value. Should the biomass spectrum (TRUE) be plotted
#'   or the abundance in numbers (FALSE). Default is TRUE.
#' @param print_it Display the plot, or just return the ggplot2 object.
#'   Defaults to TRUE
#' @param total A boolean value that determines whether the total over all
#'   species in the system is plotted as well. Default is FALSE
#' @param ... Other arguments (currently unused)
#'   
#' @return A ggplot2 object
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=20, t_save = 2)
#' plotSpectra(sim)
#' plotSpectra(sim, min_w = 1e-6)
#' plotSpectra(sim, time_range = 10:20)
#' plotSpectra(sim, time_range = 10:20, biomass = FALSE)
#' plotSpectra(sim, species = c("Cod", "Herring"), total = TRUE)
#' }
setGeneric('plotSpectra', function(sim, ...)
    standardGeneric('plotSpectra'))

#' Plot the abundance spectra using a \code{MizerSim} object.
#' @rdname plotSpectra
setMethod('plotSpectra', signature(sim='MizerSim'),
    function(sim, species = as.character(sim@params@species_params$species),
             time_range = max(as.numeric(dimnames(sim@n)$time)), 
             min_w = min(sim@params@w)/100, biomass = TRUE, print_it = TRUE, 
             total = FALSE, ...){
        time_elements <- get_time_elements(sim,time_range)
        spec_n <- apply(sim@n[time_elements, , ,drop=FALSE], c(2,3), mean)
        background_n <- apply(sim@n_pp[time_elements,,drop=FALSE],2,mean)
        if (total) {
            # Calculate total community abundance
            fish_idx <- (length(sim@params@w_full)-length(sim@params@w)+1):
                length(sim@params@w_full)
            total_n <- background_n
            total_n[fish_idx] <- total_n[fish_idx] + colSums(spec_n)
            if (biomass) {
                total_n <- total_n * sim@params@w_full
            }
        }
        # Select only the desired species
        spec_n <- spec_n[as.character(dimnames(spec_n)[[1]]) %in% species, ,
                         drop = FALSE]
        y_label = "Abundance density [1/g]"
        if (biomass){
            spec_n <- sweep(spec_n,2,sim@params@w,"*")
            background_n <- background_n * sim@params@w_full
            y_label = "Biomass density"
        }
        # Make data.frame for plot
        plot_dat <- data.frame(value = c(spec_n), 
                               Species = dimnames(spec_n)[[1]], 
                               w = rep(sim@params@w, 
                                       each = length(species)))
        plot_dat <- rbind(plot_dat, 
                          data.frame(value = c(background_n), 
                                     Species = "Plankton", 
                                     w = sim@params@w_full))
        if (total) {
            plot_dat <- rbind(plot_dat, 
                              data.frame(value = c(total_n), 
                                         Species = "Total", 
                                         w = sim@params@w_full))
        }
        # lop off 0s in background and apply min_w
        plot_dat <- plot_dat[(plot_dat$value > 0) & (plot_dat$w >= min_w),]
        if (length(species) > 12) {
            p <- ggplot(plot_dat) + 
                geom_line(aes(x=w, y = value, group = Species))
        } else {
            p <- ggplot(plot_dat) + 
                geom_line(aes(x=w, y = value, colour = Species, linetype=Species)) 
        }
        p <- p + 
            scale_x_continuous(name = "Size [g]", trans="log10") + 
            scale_y_continuous(name = y_label, trans="log10")
        if (print_it)
            print(p)
        return(p)
    }
)


#' Plot the feeding level of species by size
#' 
#' After running a projection, plot the feeding level of each species by size. 
#' The feeding level is averaged over the specified time range (a single value
#' for the time range can be used).
#' 
#' @param sim An object of class \code{MizerSim}.
#' @param species Name or vector of names of the species to be plotted. By
#'   default all species are plotted.
#' @param time_range The time range (either a vector of values, a vector of min
#'   and max time, or a single value) to average the abundances over. Default is
#'   the final time step.
#' @param print_it Display the plot, or just return the ggplot2 object.
#'   Defaults to TRUE
#' @param ... Other arguments to pass to \code{getFeedingLevel} method
#'
#' @return A ggplot2 object
#' @export
#' @seealso \code{\link{getFeedingLevel}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=20, t_save = 2)
#' plotFeedingLevel(sim)
#' plotFeedingLevel(sim, time_range = 10:20)
#' }
setGeneric('plotFeedingLevel', function(sim, ...)
    standardGeneric('plotFeedingLevel'))

#' Plot the feeding level using a \code{MizerSim} object.
#' @rdname plotFeedingLevel
setMethod('plotFeedingLevel', signature(sim='MizerSim'),
    function(sim, species = as.character(sim@params@species_params$species),
             time_range = max(as.numeric(dimnames(sim@n)$time)), 
             print_it = TRUE, ...) {
        feed_time <- getFeedingLevel(sim, time_range=time_range, 
                                     drop=FALSE, ...)
        feed <- apply(feed_time, c(2,3), mean)
        feed <- feed[as.character(dimnames(feed)[[1]]) %in% species, , 
                     drop = FALSE]
        plot_dat <- data.frame(value = c(feed), 
                               Species = dimnames(feed)[[1]], 
                               w = rep(sim@params@w, each=length(species)))
        if (length(species) > 12) {
            p <- ggplot(plot_dat) + 
                geom_line(aes(x=w, y = value, group = Species))
        } else {
            p <- ggplot(plot_dat) + 
                geom_line(aes(x=w, y = value, colour = Species, linetype=Species))
        }
        p <- p + 
            scale_x_continuous(name = "Size [g]", trans="log10") + 
            scale_y_continuous(name = "Feeding Level", limits=c(0,1))
        if (print_it) {
            print(p)
        }
        return(p)
    }
)

#' Plot predation mortality rate of each species by size
#' 
#' After running a projection, plot the predation mortality rate of each species
#' by size. The mortality rate is averaged over the specified time range (a
#' single value for the time range can be used to plot a single time step).
#' 
#' @param sim An object of class \code{MizerSim}
#' @param species Name or vector of names of the species to be plotted. By
#'   default all species are plotted.
#' @param time_range The time range (either a vector of values, a vector of min
#'   and max time, or a single value) to average the abundances over. Default is
#'   the final time step.
#' @param print_it Display the plot, or just return the ggplot2 object.
#'   Defaults to TRUE
#' @param ... Other arguments to pass to \code{getM2} method.
#'
#' @return A ggplot2 object
#' @export
#' @seealso \code{\link{getM2}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=20, t_save = 2)
#' plotM2(sim)
#' plotM2(sim, time_range = 10:20)
#' }
setGeneric('plotM2', function(sim, ...)
    standardGeneric('plotM2'))

#' Plot M2 using a \code{MizerSim} object.
#' @rdname plotM2
setMethod('plotM2', signature(sim='MizerSim'),
    function(sim, species = as.character(sim@params@species_params$species),
             time_range = max(as.numeric(dimnames(sim@n)$time)), 
             print_it = TRUE, ...) {
	m2_time <- getM2(sim, time_range=time_range, drop=FALSE, ...)
	m2 <- apply(m2_time, c(2,3), mean)
	m2 <- m2[as.character(dimnames(m2)[[1]]) %in% species, , 
	             drop = FALSE]
	plot_dat <- data.frame(value = c(m2), 
	                       Species = dimnames(m2)[[1]], 
	                       w = rep(sim@params@w, each=length(species)))
    if (length(species) > 12) {
        p <- ggplot(plot_dat) + 
            geom_line(aes(x=w, y = value, group = Species))
    } else {
        p <- ggplot(plot_dat) + 
            geom_line(aes(x=w, y = value, colour = Species, linetype=Species))
    }
	p <- p + 
	    scale_x_continuous(name = "Size [g]", trans="log10") + 
	    scale_y_continuous(name = "Mortality rate [1/year]", 
	                       limits=c(0,max(plot_dat$value)))
    if (print_it) {
        print(p)
    }
    return(p)
    }
)

#' Plot total fishing mortality of each species by size
#' 
#' After running a projection, plot the total fishing mortality of each species
#' by size. The total fishing mortality is averaged over the specified time
#' range (a single value for the time range can be used to plot a single time
#' step).
#' 
#' @param sim An object of class \code{MizerSim}.
#' @param species Name or vector of names of the species to be plotted. By
#'   default all species are plotted.
#' @param time_range The time range (either a vector of values, a vector of min
#'   and max time, or a single value) to average the abundances over. Default is
#'   the final time step.
#' @param print_it Display the plot, or just return the ggplot2 object.
#'   Defaults to TRUE
#' @param ... Other arguments to pass to \code{getFMort} method
#'
#' @return A ggplot2 object
#' @export
#' @seealso \code{\link{getFMort}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=20, t_save = 2)
#' plotFMort(sim)
#' plotFMort(sim, time_range = 10:20)
#' }
setGeneric('plotFMort', function(sim, ...)
    standardGeneric('plotFMort'))

#' Plot total fishing mortality using a \code{MizerSim} object.
#' @rdname plotFMort
setMethod('plotFMort', signature(sim='MizerSim'),
    function(sim, species = as.character(sim@params@species_params$species),
             time_range = max(as.numeric(dimnames(sim@n)$time)), 
             print_it = TRUE, ...){
	f_time <- getFMort(sim, time_range=time_range, drop=FALSE, ...)
	f <- apply(f_time, c(2,3), mean)
	f <- f[as.character(dimnames(f)[[1]]) %in% species, , drop = FALSE]
	plot_dat <- data.frame(value = c(f), 
	                       Species = dimnames(f)[[1]], 
	                       w = rep(sim@params@w, each=length(species)))
    if (length(species) > 12) {
        p <- ggplot(plot_dat) + geom_line(aes(x=w, y = value, group = Species))
    } else {
        p <- ggplot(plot_dat) + 
            geom_line(aes(x=w, y = value, colour = Species, linetype=Species)) + 
            scale_x_continuous(name = "Size", trans="log10") + 
            scale_y_continuous(name = "Total fishing mortality", 
                               limits=c(0,max(plot_dat$value)))
    }
	p <- p + 
	    scale_x_continuous(name = "Size [g]", trans="log10") + 
	    scale_y_continuous(name = "Total fishing mortality [1/Year]", 
	                       limits=c(0,max(plot_dat$value)))
    if (print_it) {
        print(p)
    }
	return(p)
    }
)


#' Summary plot for \code{MizerSim} objects
#' 
#' After running a projection, produces 5 plots in the same window: feeding
#' level, abundance spectra, predation mortality and fishing mortality of each
#' species by size; and biomass of each species through time. This method just
#' uses the other plotting methods and puts them all in one window.
#' 
#' @param x An object of class \code{MizerSim}
#' @param y Not used
#' @param ...  For additional arguments see the documentation for
#'   \code{\link{plotBiomass}},
#'   \code{\link{plotFeedingLevel}},\code{\link{plotSpectra}},\code{\link{plotM2}}
#'   and \code{\link{plotFMort}}.
#' @return A viewport object
#' @export
#' @rdname plotMizerSim
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=20, t_save = 2)
#' plot(sim)
#' plot(sim, time_range = 10:20) # change time period for size-based plots
#' plot(sim, min_w = 10, max_w = 1000) # change size range for biomass plot
#' }
setMethod("plot", signature(x="MizerSim", y="missing"),
    function(x, ...){
	p1 <- plotFeedingLevel(x,print_it = FALSE,...)
	p2 <- plotSpectra(x,print_it = FALSE,...)
	p3 <- plotBiomass(x, y_ticks = 3, print_it = FALSE,...)
	p4 <- plotM2(x,print_it = FALSE,...)
	p5 <- plotFMort(x,print_it = FALSE,...)
	grid::grid.newpage()
	glayout <- grid::grid.layout(3,2) # widths and heights arguments
	vp <- grid::viewport(layout = glayout)
	grid::pushViewport(vp)
	vplayout <- function(x,y)
	  grid::viewport(layout.pos.row=x, layout.pos.col = y)
	print(p1+ theme(legend.position="none"), vp = vplayout(1,1))
	print(p3+ theme(legend.position="none"), vp = vplayout(1,2))
	print(p4+ theme(legend.position="none"), vp = vplayout(2,1))
	print(p5+ theme(legend.position="none"), vp = vplayout(2,2))
	print(p2+ theme(legend.position="right", legend.key.size=unit(0.1,"cm")), 
	      vp = vplayout(3,1:2))
    }
)


#' Plot growth curves giving weight as a function of age
#' 
#' Uses the growth rates at the final time of a simulation to calculate
#' the size at age.
#' 
#' When the growth curve for only a single species is plotted, horizontal
#' lines are included that indicate the maturity size and the maximum size for 
#' that species. If furthermore the species parameters contain the variables
#' a and b for length to weight conversion and the von Bertalanffy parameter
#' k_vb, then the von Bertalanffy growth curve is superimposed in black.
#' 
#' @param sim MizerSim object
#' @param species Name or vector of names of the species to be plotted. By
#'   default all species are plotted.
#' @param max_age The age up to which the weight is to be plotted. Default is 20
#' @param percentage Boolean value. If TRUE, the size is shown as a percentage
#'   of the maximal size.
#' @param print_it Display the plot, or just return the ggplot2 object.
#'   Defaults to TRUE
#' @param ... Other arguments (unused)
#' 
#' @return A ggplot2 object
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=20, t_save = 2)
#' plotGrowthCurves(sim, percentage = TRUE)
#' plotGrowthCurves(sim, species = "Cod", max_age = 24)
#' }
setGeneric('plotGrowthCurves', function(sim, ...)
    standardGeneric('plotGrowthCurves'))

#' Plot growth curves using a \code{MizerSim} object.
#' @rdname plotGrowthCurves
setMethod('plotGrowthCurves', signature(sim='MizerSim'),
    function(sim, species = as.character(sim@params@species_params$species),
             max_age = 20, percentage = FALSE, print_it = TRUE) {
        # reorder list of species to coincide with order in sim
        idx <- which(sim@params@species_params$species %in% species)
        species <- sim@params@species_params$species[idx]
        age <- seq(0, max_age, length.out = 50)
        ws <- array(dim = c(length(species), length(age)), 
                    dimnames = list(Species = species, Age = age))
        g <- getEGrowth(sim@params, sim@n[dim(sim@n)[1], , ], sim@n_pp[dim(sim@n)[1], ])
        for (j in 1:length(species)) {
            i <- idx[j]
            g_fn <- stats::approxfun(sim@params@w, g[i, ])
            myodefun <- function(t, state, parameters){
                return(list(g_fn(state)))
            }
            ws[j, ] <- deSolve::ode(y = sim@params@species_params$w_min[i], 
                           times = age, func = myodefun)[,2]
            if (percentage) {
                ws[j, ] <- ws[j, ] / sim@params@species_params$w_inf[i] * 100
            }
        }	
        plot_dat <- reshape2::melt(ws)
        plot_dat$Species <- as.character(plot_dat$Species)
        if (length(species) > 12) {
            p <- ggplot(plot_dat) + 
                geom_line(aes(x=Age, y = value, group = Species))
        } else {
            p <- ggplot(plot_dat) + 
                geom_line(aes(x = Age, y = value, 
                              colour = Species, linetype=Species))
        }
        y_label <- if (percentage) "Percent of maximum size" else "Size [g]"
        p <- p + 
            scale_x_continuous(name = "Age [Years]") + 
            scale_y_continuous(name = y_label)
        
        # Extra stuff for single-species case
        if (length(species) == 1 && !percentage) {
            w_inf <- sim@params@species_params$w_inf[idx[1]]
            p <- p + geom_hline(yintercept = w_inf) +
                annotate("text", 0, w_inf, vjust = -1, label = "Maximum")
            w_mat <- sim@params@species_params$w_mat[idx[1]]
            p <- p + geom_hline(yintercept = w_mat) +
                annotate("text", 0, w_mat, vjust = -1, label = "Maturity")
            if (all(c("a", "b", "k_vb") %in% names(sim@params@species_params))) {
                a <- sim@params@species_params$a[idx[1]]
                b <- sim@params@species_params$b[idx[1]]
                k_vb <- sim@params@species_params$k_vb[idx[1]]
                L_inf <- (w_inf/a)^(1/b)
                vb <- a * (L_inf * (1 - exp(-k_vb * age)))^b
                dat <- data.frame(x = age, y = vb)
                p <- p + geom_line(data = dat, aes(x = x, y = y))
            }
        }
        
        if (print_it) {
            print(p)
        }
        return(p)
    }
)