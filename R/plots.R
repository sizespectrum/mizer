# Plotting methods for the MizerSim class

# Copyright 2012 Finlay Scott and Julia Blanchard. 
# Distributed under the GPL 3 or later 
# Maintainer: Gustav Delius, University of York, <gustav.delius@york.ac.uk>

# Hackiness to get past the 'no visible binding ... ' warning when running check
 utils::globalVariables(c("time", "value", "Species", "w", "gear"))

# Biomass through time
# Pretty easy - user could do it themselves by hand for fine tuning if necessary
 
#' Plot the biomass of each species through time
#'
#' After running a projection, the biomass of each species can be plotted
#' against time. The biomass is calculated within user defined size limits (see
#' \code{\link{getBiomass}}). This plot is pretty easy to do by hand. It just
#' gets the biomass using the \code{\link{getBiomass}} method and plots using
#' the ggplot2 package. You can then fiddle about with colours and linetypes
#' etc. Just look at the source code for details.
#' 
#' @param object An object of class \code{MizerSim}.
#' @param start_time The first time step to be plotted. Default is the beginning
#'   of the time series.
#' @param end_time The first time step to be plotted. Default is the end of the
#'   time series.
#' @param print_it Display the plot, or just return the ggplot2 object. Default
#'   value is TRUE
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
#' sim <- project(params, effort=1, t_max=20, t_save = 2)
#' plotBiomass(sim)
#' plotBiomass(sim, min_w = 10, max_w = 1000)
#' }
setGeneric('plotBiomass', function(object, ...)
    standardGeneric('plotBiomass'))

#' Plot the biomass using a \code{MizerSim} object.
#' @rdname plotBiomass
setMethod('plotBiomass', signature(object='MizerSim'),
    function(object, print_it=TRUE, start_time=as.numeric(dimnames(object@n)[[1]][1]), end_time = as.numeric(dimnames(object@n)[[1]][dim(object@n)[1]]), ...){
        b <- getBiomass(object, ...)
        names(dimnames(b))[names(dimnames(b))=="sp"] <- "Species"
        if(start_time >= end_time){
            stop("start_time must be less than end_time")
        }
        b <- b[(as.numeric(dimnames(b)[[1]]) >= start_time) & (as.numeric(dimnames(b)[[1]]) <= end_time),,drop=FALSE]
        bm <- melt(b)
        # Force Species column to be a character (if numbers used - may be
        # interpreted as integer and hence continuous)
        bm$Species <- as.character(bm$Species)
        # Due to log10, need to set a minimum value, seems like a feature in ggplot
        min_value <- 1e-300
        bm <- bm[bm$value >= min_value,]
        p <- ggplot(bm) + geom_line(aes(x=time,y=value, colour=Species, linetype=Species)) + scale_y_continuous(trans="log10", name="Biomass") + scale_x_continuous(name="Time") 
        if (nrow(object@params@species_params)>12){
        p <- ggplot(bm) + geom_line(aes(x=time,y=value, group=Species)) + scale_y_continuous(trans="log10", name="Biomass") + scale_x_continuous(name="Time") 
        }
        if (print_it)
            print(p)
        return(p)
    }
)

#' Plot the total yield of each species through time
#'
#' After running a projection, the total yield of each species across all 
#' fishing gears can be plotted against time. This plot is pretty easy to do by
#' hand. It just gets the biomass using the \code{\link{getYield}} method and
#' plots using the ggplot2 package. You can then fiddle about with colours and
#' linetypes etc. Just look at the source code for details.
#' 
#' @param object An object of class \code{MizerSim}
#' @param print_it Display the plot, or just return the ggplot2 object
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
#' sim <- project(params, effort=1, t_max=20, t_save = 2)
#' plotYield(sim)
#' }
setGeneric('plotYield', function(object, ...)
    standardGeneric('plotYield'))

#' Plot the yield using a \code{MizerSim} object.
#' @rdname plotYield
setMethod('plotYield', signature(object='MizerSim'),
    function(object, print_it = TRUE, ...){
	y <- getYield(object, ...)
	names(dimnames(y))[names(dimnames(y))=="sp"] <- "Species"
	ym <- melt(y)
	p <- ggplot(ym) + geom_line(aes(x=time,y=value, colour=Species, linetype=Species)) + scale_y_continuous(trans="log10", name="Yield") + scale_x_continuous(name="Time") 
    if (nrow(object@params@species_params)>12){
        p <- ggplot(ym) + geom_line(aes(x=time,y=value, group=Species)) + scale_y_continuous(trans="log10", name="Yield") + scale_x_continuous(name="Time") 
    }
    if (print_it)
        print(p)
	return(p)
    }
)

#' Plot the total yield of each species by gear through time
#'
#' After running a projection, the total yield of each species by fishing gear
#' can be plotted against time. This plot is pretty easy to do by hand. It just
#' gets the biomass using the \code{\link{getYieldGear}} method and plots using
#' the ggplot2 package. You can then fiddle about with colours and linetypes
#' etc. Just look at the source code for details.
#' 
#' @param object An object of class \code{MizerSim}
#' @param print_it Display the plot, or just return the ggplot2 object
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
#' sim <- project(params, effort=1, t_max=20, t_save = 2)
#' plotYieldGear(sim)
#' }
setGeneric('plotYieldGear', function(object, ...)
    standardGeneric('plotYieldGear'))

#' Plot the yield of each gear using a \code{MizerSim} object.
#' @rdname plotYieldGear
setMethod('plotYieldGear', signature(object='MizerSim'),
    function(object, print_it=TRUE, ...){
	y <- getYieldGear(object, ...)
	names(dimnames(y))[names(dimnames(y))=="sp"] <- "Species"
	ym <- melt(y)
	p <- ggplot(ym) + geom_line(aes(x=time,y=value, colour=Species, linetype=gear)) + scale_y_continuous(trans="log10", name="Yield") + scale_x_continuous(name="Time") 
    if (nrow(object@params@species_params)>12){
        p <- ggplot(ym) + geom_line(aes(x=time,y=value, group=Species)) + scale_y_continuous(trans="log10", name="Yield") + scale_x_continuous(name="Time") 
    }
    if (print_it)
        print(p)
	return(p)
    }
)

#' Plot the abundance spectra of each species and the background population
#' 
#' After running a projection, the spectra of the abundance of each species and
#' the background population can be plotted. The abundance is averaged over the
#' specified time range (a single value for the time range can be used to plot a
#' single time step). The abundance can be in terms of numbers or biomass,
#' depending on the \code{biomass} argument.
#' 
#' @param object An object of class \code{MizerSim}.
#' @param time_range The time range (either a vector of values, a vector of min
#'   and max time, or a single value) to average the abundances over. Default is
#'   the final time step.
#' @param min_w Minimum weight to be plotted (useful for truncating the
#'   background spectrum). Default value is a hundredth of the minimum size
#'   value of the community.
#' @param biomass A boolean value. Should the biomass spectrum (TRUE) be plotted
#'   or the abundance in numbers (FALSE). Default is TRUE.
#' @param print_it Display the plot, or just return the ggplot2 object
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
#' }
setGeneric('plotSpectra', function(object, ...)
    standardGeneric('plotSpectra'))

#' Plot the abundance spectra using a \code{MizerSim} object.
#' @rdname plotSpectra
setMethod('plotSpectra', signature(object='MizerSim'),
    function(object, time_range = max(as.numeric(dimnames(object@n)$time)), min_w =min(object@params@w)/100, biomass = TRUE, print_it = TRUE, ...){
        time_elements <- get_time_elements(object,time_range)
        spec_n <- apply(object@n[time_elements,,,drop=FALSE],c(2,3), mean)
        background_n <- apply(object@n_pp[time_elements,,drop=FALSE],2,mean)
        y_axis_name = "Abundance"
        if (biomass){
            spec_n <- sweep(spec_n,2,object@params@w,"*")
            background_n <- background_n * object@params@w_full
            y_axis_name = "Biomass"
        }
        # Make data.frame for plot
        plot_dat <- data.frame(value = c(spec_n), Species = dimnames(spec_n)[[1]], w = rep(object@params@w, each=nrow(object@params@species_params)))
        plot_dat <- rbind(plot_dat, data.frame(value = c(background_n), Species = "Background", w = object@params@w_full))
        # lop off 0s in background and apply min_w
        plot_dat <- plot_dat[(plot_dat$value > 0) & (plot_dat$w >= min_w),]
        p <- ggplot(plot_dat) + geom_line(aes(x=w, y = value, colour = Species, linetype=Species)) + scale_x_continuous(name = "Size", trans="log10") + scale_y_continuous(name = y_axis_name, trans="log10")
        if (nrow(object@params@species_params)>12){
            p <- ggplot(plot_dat) + geom_line(aes(x=w, y = value, group = Species)) + scale_x_continuous(name = "Size", trans="log10") + scale_y_continuous(name = y_axis_name, trans="log10")
        }
        if (print_it)
            print(p)
        return(p)
    }
)


#' Plot the feeding level of each species by size
#' 
#' After running a projection, plot the feeding level of each species by size. 
#' The feeding level is averaged over the specified time range (a single value
#' for the time range can be used).
#' 
#' @param object An object of class \code{MizerSim}.
#' @param time_range The time range (either a vector of values, a vector of min
#'   and max time, or a single value) to average the abundances over. Default is
#'   the final time step.
#' @param print_it Display the plot, or just return the ggplot2 object
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
setGeneric('plotFeedingLevel', function(object, ...)
    standardGeneric('plotFeedingLevel'))

#' PLot the feeding level using a \code{MizerSim} object.
#' @rdname plotFeedingLevel
setMethod('plotFeedingLevel', signature(object='MizerSim'),
    function(object, time_range = max(as.numeric(dimnames(object@n)$time)), print_it = TRUE, ...){
        feed_time <- getFeedingLevel(object=object, time_range=time_range, drop=FALSE, ...)
        feed <- apply(feed_time, c(2,3), mean)
        plot_dat <- data.frame(value = c(feed), Species = dimnames(feed)[[1]], w = rep(object@params@w, each=nrow(object@params@species_params)))
        p <- ggplot(plot_dat) + geom_line(aes(x=w, y = value, colour = Species, linetype=Species)) + scale_x_continuous(name = "Size", trans="log10") + scale_y_continuous(name = "Feeding Level", limits=c(0,1))
        if (nrow(object@params@species_params)>12){
            p <- ggplot(plot_dat) + geom_line(aes(x=w, y = value, group = Species)) + scale_x_continuous(name = "Size", trans="log10") + scale_y_continuous(name = "Feeding Level", limits=c(0,1))
        }
        if (print_it)
            print(p)
        return(p)
    }
)

#' Plot M2 of each species by size
#' 
#' After running a projection, plot M2 of each species by size. M2 is averaged
#' over the specified time range (a single value for the time range can be used
#' to plot a single time step).
#' 
#' @param object An object of class \code{MizerSim}
#' @param time_range The time range (either a vector of values, a vector of min
#'   and max time, or a single value) to average the abundances over. Default is
#'   the final time step.
#' @param print_it Display the plot, or just return the ggplot2 object
#' @param ... Other arguments to pass to \code{getM2} method.
#'
#' @return A ggplot2 object
#' @export
#' @seealso \code{\link{getM2}}
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
setGeneric('plotM2', function(object, ...)
    standardGeneric('plotM2'))

#' Plot M2 using a \code{MizerSim} object.
#' @rdname plotM2
setMethod('plotM2', signature(object='MizerSim'),
    function(object, time_range = max(as.numeric(dimnames(object@n)$time)), print_it = TRUE, ...){
	m2_time <- getM2(object, time_range=time_range, drop=FALSE, ...)
	m2 <- apply(m2_time, c(2,3), mean)
	plot_dat <- data.frame(value = c(m2), Species = dimnames(m2)[[1]], w = rep(object@params@w, each=nrow(object@params@species_params)))
	p <- ggplot(plot_dat) + geom_line(aes(x=w, y = value, colour = Species, linetype=Species)) + scale_x_continuous(name = "Size", trans="log10") + scale_y_continuous(name = "M2", limits=c(0,max(plot_dat$value)))
    if (nrow(object@params@species_params)>12){
        p <- ggplot(plot_dat) + geom_line(aes(x=w, y = value, group = Species)) + scale_x_continuous(name = "Size", trans="log10") + scale_y_continuous(name = "M2", limits=c(0,max(plot_dat$value)))
    }
    if (print_it)
        print(p)
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
#' @param object An object of class \code{MizerSim}.
#' @param time_range The time range (either a vector of values, a vector of min
#'   and max time, or a single value) to average the abundances over. Default is
#'   the final time step.
#' @param print_it Display the plot, or just return the ggplot2 object
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
setGeneric('plotFMort', function(object, ...)
    standardGeneric('plotFMort'))

#' Plot total fishing mortality using a \code{MizerSim} object.
#' @rdname plotFMort
setMethod('plotFMort', signature(object='MizerSim'),
    function(object, time_range = max(as.numeric(dimnames(object@n)$time)), print_it = TRUE, ...){
	f_time <- getFMort(object, time_range=time_range, drop=FALSE, ...)
	f <- apply(f_time, c(2,3), mean)
	plot_dat <- data.frame(value = c(f), Species = dimnames(f)[[1]], w = rep(object@params@w, each=nrow(object@params@species_params)))
	p <- ggplot(plot_dat) + geom_line(aes(x=w, y = value, colour = Species, linetype=Species)) + scale_x_continuous(name = "Size", trans="log10") + scale_y_continuous(name = "Total fishing mortality", limits=c(0,max(plot_dat$value)))
    if (nrow(object@params@species_params)>12){
        p <- ggplot(plot_dat) + geom_line(aes(x=w, y = value, group = Species)) + scale_x_continuous(name = "Size", trans="log10") + scale_y_continuous(name = "Total fishing mortality", limits=c(0,max(plot_dat$value)))
    }
    if (print_it)
        print(p)
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
	p3 <- plotBiomass(x,print_it = FALSE,...)
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
	print(p2+ theme(legend.position="right", legend.key.size=unit(0.1,"cm")), vp = vplayout(3,1:2))
    }
)

