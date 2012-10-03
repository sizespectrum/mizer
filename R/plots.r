# Plotting methods for the MizerSim class

# Copyright 2012 Finlay Scott and Julia Blanchard. Distributed under the GPL 2 or later
# Maintainer: Finlay Scott, CEFAS

# Biomass through time
# Pretty easy - user could do it themselves by hand for fine tuning if necessary
 
#' Plot the biomass of each species through time
#'
#' After running a projection, the biomass of each species can plotted against time. The biomass is calculated within user defined limits (see \code{\link{getBiomass}}.
#' This plot is pretty easy to do by hand. It just gets the biomass using the \code{\link{getBiomass}} method and plots uses ggplot2. You can then fiddle about with colours and linetypes etc. Just look at the source code for details.
#' 
#' @param object An object of class \code{MizerSim}
#' @param min_w minimum weight of species to be used in the calculation
#' @param max_w maximum weight of species to be used in the calculation
#' @param min_l minimum length of species to be used in the calculation
#' @param max_l maximum length of species to be used in the calculation
#'
#' @return A ggplot2 object
#' @export
#' @docType methods
#' @rdname plotBiomass-methods
setGeneric('plotBiomass', function(object, ...)
    standardGeneric('plotBiomass'))
#' @rdname plotBiomass-methods
#' @aliases plotBiomass,MizerSim-method
setMethod('plotBiomass', signature(object='MizerSim'),
    function(object, ...){
	b <- getBiomass(object, ...)
	bm <- melt(b)
	p <- ggplot(bm) + geom_line(aes(x=time,y=value, colour=sp, linetype=sp)) + scale_y_continuous(trans="log10", name="Biomass") + scale_x_continuous(name="Time") #+ scale_linetype_discrete(name="Species")
	print(p)
	return(p)
    })
# Do we need a window() method that returns a truncated MizerSim object?
# For the instantaneous plots (e.g. spectra, feeding level, F etc) we need to decide what time we are plotting
# Final time step
# A single time step
# Average of last x time steps
# Average of x time steps
# One argument needed: time_period
# Data calculated as mean of that time period (maybe an option for summary function: mean, geomean etc?)
# If just a single value then mean of that 
# How does user specify?
    # Time (as in actual years, not number of time steps)?
    # Row elements?
    # start and end time, or full range?
# Rules: time period has to be contiguous, i.e. cannot do c(1:5, 8:10)
# Argument can be character, or numeric. But always forced as.character.

# get_time_elements
# internal function to get the array element references of the time dimension for the time based slots of a MizerSim object
# time_range can be character or numeric
# Necessary to include a slot_name argument because the effort and abundance slots have different time dimensions
get_time_elements <- function(sim,time_range,slot_name="n"){
    if (!(slot_name %in% c("n","effort")))
	stop("'slot_name' argument should be 'n' or 'effort'")
    if (!is(sim,"MizerSim"))
	stop("First argument to get_time_elements function must be of class MizerSim")
    time_range <- range(as.numeric(time_range))
    # Check that time range is even in object
    sim_time_range <- range(as.numeric(dimnames(slot(sim,slot_name))$time))
    if ((time_range[1] < sim_time_range[1]) | (time_range[2] > sim_time_range[2]))
	stop("Time range is outside the time range of the modell")
    time_elements <- (as.numeric(dimnames(slot(sim,slot_name))$time) >= time_range[1]) & (as.numeric(dimnames(slot(sim,slot_name))$time) <= time_range[2])
    names(time_elements) <- dimnames(slot(sim,slot_name))$time
    return(time_elements)
}


#' Plot the abundance spectra of each species and the background population
#'
#' After running a projection, the spectra of the abundance of each species and the background population can be plotted.
#' The abundance is averaged over the specified time range. A single value for the time range can be used to explore the spectra at a single point in time.
#' 
#' @param object An object of class \code{MizerSim}
#' @param min_w Minimum weight to be plotted (useful for truncating the background spectrum). Default value is a tenth of the minimum size value of the community. 
#' @param time_range The time range (either a vector of values, a vector of min and max time, or a single value) to average the abundances over. Default is the final time step.
#'
#' @return A ggplot2 object
#' @export
#' @docType methods
#' @rdname plotNSpectra-methods
setGeneric('plotNSpectra', function(object, ...)
    standardGeneric('plotNSpectra'))
#' @rdname plotNSpectra-methods
#' @aliases plotNSpectra,MizerSim-method
setMethod('plotNSpectra', signature(object='MizerSim'),
    function(object, time_range = max(as.numeric(dimnames(object@n)$time)), min_w =min(object@params@w)/10, ...){
	time_elements <- get_time_elements(object,time_range)
	spec_n <- apply(object@n[time_elements,,,drop=FALSE],c(2,3), mean)
	background_n <- apply(object@n_pp[time_elements,,drop=FALSE],2,mean)
	# Need this to be real w not names w - real ugly - can we clean this up?
	dimnames(spec_n)$w <- object@params@w
	plot_dat <- rbind(melt(spec_n), cbind(sp = "Background", w = object@params@w_full, value = background_n))
	plot_dat$w <- as.numeric(plot_dat$w)
	plot_dat$value <- as.numeric(plot_dat$value)
	# lop off 0s in background and apply min_w
	plot_dat <- plot_dat[(plot_dat$value > 0) & (plot_dat$w >= min_w),]
	names(plot_dat)[names(plot_dat)=="sp"] <- "Species"
	p <- ggplot(plot_dat) + geom_line(aes(x=w, y = value, colour = Species, linetype=Species)) + scale_x_continuous(name = "Size", trans="log10") + scale_y_continuous(name = "Abundance", trans="log10")
	print(p)
	return(p)
    })


