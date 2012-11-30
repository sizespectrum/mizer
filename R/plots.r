# Plotting methods for the MizerSim class

# Copyright 2012 Finlay Scott and Julia Blanchard. Distributed under the GPL 2 or later
# Maintainer: Finlay Scott, CEFAS

# Biomass through time
# Pretty easy - user could do it themselves by hand for fine tuning if necessary
 
#' Plot the biomass of each species through time
#'
#' After running a projection, the biomass of each species can plotted against time. The biomass is calculated within user defined limits (see \code{\link{getBiomass}}).
#' This plot is pretty easy to do by hand. It just gets the biomass using the \code{\link{getBiomass}} method and plots using ggplot2. You can then fiddle about with colours and linetypes etc. Just look at the source code for details.
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
#' @examples
#' data(species_params_gears)
#' data(inter)
#' params <- MizerParams(species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=20, t_save = 2)
#' plotBiomass(sim)
#' plotBiomass(sim, min_l = 10, max_l = 25)
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

#' Plot the abundance spectra of each species and the background population
#'
#' After running a projection, the spectra of the abundance of each species and the background population can be plotted.
#' The abundance is averaged over the specified time range (a single value for the time range can be used).
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
#' @examples
#' data(species_params_gears)
#' data(inter)
#' params <- MizerParams(species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=20, t_save = 2)
#' plotNSpectra(sim)
#' plotNSpectra(sim, min_w = 1e-6)
#' plotNSpectra(sim, time_range = 10:20)
setMethod('plotNSpectra', signature(object='MizerSim'),
    function(object, time_range = max(as.numeric(dimnames(object@n)$time)), min_w =min(object@params@w)/10, ...){
	time_elements <- get_time_elements(object,time_range)
	spec_n <- apply(object@n[time_elements,,,drop=FALSE],c(2,3), mean)
	background_n <- apply(object@n_pp[time_elements,,drop=FALSE],2,mean)
#	# Need this to be real w not names w - real ugly - can we clean this up?
#	dimnames(spec_n)$w <- object@params@w
#	plot_dat <- rbind(melt(spec_n), cbind(sp = "Background", w = object@params@w_full, value = background_n))
#	plot_dat$w <- as.numeric(plot_dat$w)
#	plot_dat$value <- as.numeric(plot_dat$value)
	plot_dat <- data.frame(value = c(spec_n), Species = dimnames(spec_n)[[1]], w = rep(object@params@w, each=nrow(object@params@species_params)))
	plot_dat <- rbind(plot_dat, data.frame(value = c(background_n), Species = "Background", w = object@params@w_full))
	# lop off 0s in background and apply min_w
	plot_dat <- plot_dat[(plot_dat$value > 0) & (plot_dat$w >= min_w),]
	p <- ggplot(plot_dat) + geom_line(aes(x=w, y = value, colour = Species, linetype=Species)) + scale_x_continuous(name = "Size", trans="log10") + scale_y_continuous(name = "Abundance", trans="log10")
	print(p)
	return(p)
    })


#' Plot the feeding level of each species by size 
#'
#' After running a projection, plot the feeding level of each species by size.
#' The feeding level is averaged over the specified time range (a single value for the time range can be used).
#' 
#' @param object An object of class \code{MizerSim}
#' @param time_range The time range (either a vector of values, a vector of min and max time, or a single value) to average the abundances over. Default is the final time step.
#'
#' @return A ggplot2 object
#' @export
#' @docType methods
#' @rdname plotFeedingLevel-methods
setGeneric('plotFeedingLevel', function(object, ...)
    standardGeneric('plotFeedingLevel'))
#' @rdname plotFeedingLevel-methods
#' @aliases plotFeedingLevel,MizerSim-method
#' @examples
#' data(species_params_gears)
#' data(inter)
#' params <- MizerParams(species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=20, t_save = 2)
#' #plotFeedingLevel(sim)
#' #plotFeedingLevel(sim, time_range = 10:20)
setMethod('plotFeedingLevel', signature(object='MizerSim'),
    function(object, time_range = max(as.numeric(dimnames(object@n)$time)), ...){
	feed_time <- getFeedingLevel(object, time_range=time_range, .drop=FALSE, ...)
	feed <- apply(feed_time, c(2,3), mean)
	plot_dat <- data.frame(value = c(feed), Species = dimnames(feed)[[1]], w = rep(object@params@w, each=nrow(object@params@species_params)))
	p <- ggplot(plot_dat) + geom_line(aes(x=w, y = value, colour = Species, linetype=Species)) + scale_x_continuous(name = "Size", trans="log10") + scale_y_continuous(name = "Feeding Level", lim=c(0,1))
	print(p)
	return(p)
    })

#' Plot M2 of each species by size 
#'
#' After running a projection, plot M2 of each species by size.
#' M2 is averaged over the specified time range (a single value for the time range can be used).
#' 
#' @param object An object of class \code{MizerSim}
#' @param time_range The time range (either a vector of values, a vector of min and max time, or a single value) to average the abundances over. Default is the final time step.
#'
#' @return A ggplot2 object
#' @export
#' @docType methods
#' @rdname plotM2-methods
setGeneric('plotM2', function(object, ...)
    standardGeneric('plotM2'))
#' @rdname plotM2-methods
#' @aliases plotM2,MizerSim-method
#' @examples
#' data(species_params_gears)
#' data(inter)
#' params <- MizerParams(species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=20, t_save = 2)
#' #plotM2(sim)
#' #plotM2(sim, time_range = 10:20)
setMethod('plotM2', signature(object='MizerSim'),
    function(object, time_range = max(as.numeric(dimnames(object@n)$time)), ...){
	m2_time <- getM2(object, time_range=time_range, .drop=FALSE, ...)
	m2 <- apply(m2_time, c(2,3), mean)
	plot_dat <- data.frame(value = c(m2), Species = dimnames(m2)[[1]], w = rep(object@params@w, each=nrow(object@params@species_params)))
	p <- ggplot(plot_dat) + geom_line(aes(x=w, y = value, colour = Species, linetype=Species)) + scale_x_continuous(name = "Size", trans="log10") + scale_y_continuous(name = "M2", lim=c(0,max(plot_dat$value)))
	print(p)
	return(p)
    })

#' Plot total fishing mortality of each species by size 
#'
#' After running a projection, plot the total fishing mortality of each species by size.
#' The total fishing mortality is averaged over the specified time range (a single value for the time range can be used).
#' 
#' @param object An object of class \code{MizerSim}
#' @param time_range The time range (either a vector of values, a vector of min and max time, or a single value) to average the abundances over. Default is the final time step.
#'
#' @return A ggplot2 object
#' @export
#' @docType methods
#' @rdname plotFMort-methods
setGeneric('plotFMort', function(object, ...)
    standardGeneric('plotFMort'))
#' @rdname plotFMort-methods
#' @aliases plotFMort,MizerSim-method
#' @examples
#' data(species_params_gears)
#' data(inter)
#' params <- MizerParams(species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=20, t_save = 2)
#' #plotFMort(sim)
#' #plotFMort(sim, time_range = 10:20)
setMethod('plotFMort', signature(object='MizerSim'),
    function(object, time_range = max(as.numeric(dimnames(object@n)$time)), ...){
	f_time <- getFMort(object, time_range=time_range, .drop=FALSE, ...)
	f <- apply(f_time, c(2,3), mean)
	plot_dat <- data.frame(value = c(f), Species = dimnames(f)[[1]], w = rep(object@params@w, each=nrow(object@params@species_params)))
	p <- ggplot(plot_dat) + geom_line(aes(x=w, y = value, colour = Species, linetype=Species)) + scale_x_continuous(name = "Size", trans="log10") + scale_y_continuous(name = "Total fishing mortality", lim=c(0,max(plot_dat$value)))
	print(p)
	return(p)
    })

