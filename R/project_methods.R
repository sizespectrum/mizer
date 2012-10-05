# Methods used for projecting for the size based modelling package

# Copyright 2012 Finlay Scott and Julia Blanchard. Distributed under the GPL 2 or later
# Maintainer: Finlay Scott, CEFAS

# Calculate the amount of food exposed to each predator by predator size

#' getPhiPrey method for the size based model
#'
#' Calculates the amount of food exposed to each predator by predator size
#' @param object An \code{MizerParams} object
#' @param n A matrix of species abundance (species x size)
#' @param n_pp A vector of the background abundance by size
#'
#' @return A two dimensional array (predator species x predator size) 
#' @export
#' @docType methods
#' @rdname getPhiPrey-methods
setGeneric('getPhiPrey', function(object, n, n_pp,...)
    standardGeneric('getPhiPrey'))

#' @rdname getPhiPrey-methods
#' @aliases getPhiPrey,MizerParams,matrix,numeric-method
setMethod('getPhiPrey', signature(object='MizerParams', n = 'matrix', n_pp='numeric'),
    function(object, n, n_pp, ...){
	# Check n dims
	if(dim(n)[1] != dim(object@interaction)[1])
	    stop("n does not the right number of species (first dimension)")
	if(dim(n)[2] != length(object@w))
	    stop("n does not have the right number of size groups (second dimension)")
	if(length(n_pp) != length(object@w_full))
	    stop("n_pp does not have the right number of size groups")
	# n_eff_prey is the total prey abundance by size exposed to each predator
	# (prey not broken into species - here we are just working out how much a predator eats - not which species are being eaten - that is in the mortality calculation
	n_eff_prey <- sweep(object@interaction %*% n, 2, object@w * object@dw, "*") 
	# Quick reference to just the fish part of the size spectrum
	idx_sp <- (length(object@w_full) - length(object@w) + 1):length(object@w_full)
	# predKernal is predator x predator size x prey size
	# So multiply 3rd dimension of predKernal by the prey abundance
	# Then sum over 3rd dimension to get total eaten by each predator by predator size
	phi_prey_species <- rowSums(sweep(object@pred_kernel[,,idx_sp],c(1,3),n_eff_prey,"*"),dims=2)
	# Eating the background
	phi_prey_background <- rowSums(sweep(object@pred_kernel,3,object@dw_full*object@w_full*n_pp,"*"),dims=2)
	return(phi_prey_species+phi_prey_background)
})


# Feeding level
# The amount of food consumed by a predator, by each predator size

#' getFeedingLevel method for the size based model
#'
#' Calculates the amount of food consumed by a predator by predator size based on food availability, search volume and maximum intake
#' @param object A MizerParams or MizerSim object
#' @param n A matrix of species abundance (species x size). Only used if \code{object} argument is of type MizerParam.
#' @param n_pp A vector of the background abundance by size. Only used if \code{object} argument is of type MizerParam.
#' @param time_range The time range (either a vector of values, a vector of min and max time, or a single value) to average the abundances over. Default is the whole time range. Only used if \code{object} argument is of type MizerSim.
#' @param .drop should extra dimensions of length 1 in the output be dropped, simplifying the output. Defaults to TRUE  
#'
#' @return A two dimensional array (predator species x predator size) 
#' @export
#' @docType methods
#' @rdname getFeedingLevel-methods
setGeneric('getFeedingLevel', function(object, n, n_pp, ...)
    standardGeneric('getFeedingLevel'))

#' @rdname getFeedingLevel-methods
#' @aliases getFeedingLevel,MizerParams,matrix,numeric-method
setMethod('getFeedingLevel', signature(object='MizerParams', n = 'matrix', n_pp='numeric'),
    function(object, n, n_pp, ...){
	phi_prey <- getPhiPrey(object, n=n, n_pp=n_pp)
	# encountered food = available food * search volume
	encount <- object@search_vol * phi_prey
	# calculate feeding level
	f <- encount/(encount + object@intake_max)
	return(f)
})

#' @rdname getFeedingLevel-methods
#' @aliases getFeedingLevel,MizerSim,missing,missing-method
setMethod('getFeedingLevel', signature(object='MizerSim', n = 'missing', n_pp='missing'),
    function(object, time_range=dimnames(object@n)$time, .drop=TRUE, ...){
	time_elements <- get_time_elements(object,time_range)
	feed_time <- aaply(which(time_elements), 1, function(x){
			feed <- getFeedingLevel(object@params, n=object@n[x,,], n_pp = object@n_pp[x,])
			return(feed)}, .drop=.drop)
	return(feed_time)
})

# Predation rate

#' getPredRate method for the size based model
#'
#' Calculates the predation rate on each prey species by predator and prey size
#' @param object An \code{MizerParams} object
#' @param n A matrix of species abundance (species x size)
#' @param n_pp A vector of the background abundance by size
#'
#' @return A three dimensional array (predator species x predator size x prey size) 
#' @export
#' @docType methods
#' @rdname getPredRate-methods
setGeneric('getPredRate', function(object, n, n_pp,...)
    standardGeneric('getPredRate'))

#' @rdname getPredRate-methods
#' @aliases getPredRate,MizerParams,matrix,numeric-method
setMethod('getPredRate', signature(object='MizerParams', n = 'matrix', n_pp='numeric'),
    function(object, n, n_pp, ...){
	n_total_in_size_bins <- sweep(n, 2, object@dw, '*')
	f <- getFeedingLevel(object, n=n, n_pp=n_pp)
	pred_rate <- sweep(object@pred_kernel,c(1,2),(1-f)*object@search_vol*n_total_in_size_bins,"*")
	return(pred_rate)
})


# getM2
# This uses the predation rate which is also used in M2background
# Too much overlap? Inefficient? Same thing is calculated twice

#' getM2 method for the size based model
#'
#' Calculates the predation mortality on each prey species by prey size
#' @param object An \code{MizerParams} object
#' @param n A matrix of species abundance (species x size)
#' @param n_pp A vector of the background abundance by size
#' @param time_range The time range (either a vector of values, a vector of min and max time, or a single value) to average the abundances over. Default is the whole time range. Only used if \code{object} argument is of type MizerSim.
#' @param .drop should extra dimensions of length 1 in the output be dropped, simplifying the output. Defaults to TRUE  
#'
#' @return A two dimensional array (prey species x prey size) 
#' @export
#' @docType methods
#' @rdname getM2-methods
setGeneric('getM2', function(object, n, n_pp,...)
    standardGeneric('getM2'))

#' @rdname getM2-methods
#' @aliases getM2,MizerParams,matrix,numeric-method
setMethod('getM2', signature(object='MizerParams', n = 'matrix', n_pp='numeric'),
    function(object, n, n_pp, ...){
	pred_rate <- getPredRate(object,n=n,n_pp=n_pp)
	# get the element numbers that are just species
	idx_sp <- (length(object@w_full) - length(object@w) + 1):length(object@w_full)
	# Interaction is predator x prey so need to transpose so it is prey x pred
	# Sum pred_kernel over predator sizes to give total predation rate of each predator on each prey size
	m2 <- t(object@interaction) %*% colSums(aperm(pred_rate, c(2,1,3)),dims=1)[,idx_sp]
	return(m2)
})
#' @rdname getM2-methods
#' @aliases getM2,MizerSim,missing,missing-method
setMethod('getM2', signature(object='MizerSim', n = 'missing', n_pp='missing'),
    function(object, time_range=dimnames(object@n)$time, .drop=TRUE, ...){
	time_elements <- get_time_elements(object,time_range)
	m2_time <- aaply(which(time_elements), 1, function(x){
			m2 <- getM2(object@params, n=object@n[x,,], n_pp = object@n_pp[x,])
			return(m2)}, .drop=.drop)
	return(m2_time)
})


#' getM2Background method for the size based model
#'
#' Calculates the predation mortality on the background spectrum by prey size
#' @param object An \code{MizerParams} object
#' @param n A matrix of species abundance (species x size)
#' @param n_pp A vector of the background abundance by size
#'
#' @return A vector (prey size) 
#' @export
#' @docType methods
#' @rdname getM2Background-methods
setGeneric('getM2Background', function(object, n, n_pp,...)
    standardGeneric('getM2Background'))

#' @rdname getM2Background-methods
#' @aliases getM2Background,MizerParams,matrix,numeric-method
setMethod('getM2Background', signature(object='MizerParams', n = 'matrix', n_pp='numeric'),
    function(object, n, n_pp, ...){
	predRate <- getPredRate(object,n=n,n_pp=n_pp)
	M2background <- colSums(predRate,dims=2)
	return(M2background)
})

# FMort just uses a vector of efforts (length = nGear) so isn't ready for a full effortArray - could be though
# Get the fishing mortality by gear x species x size
# And maybe also by time depending on whether effort is an array

#' Get the fishing mortality by time, gear, species and size
#'
#' Calculates the fishing mortality by gear, species and size at each time step in the \code{effort} argument.
#' Fishing mortality = catchability x selectivity x effort
#'
#' @param object A \code{MizerParams} object or an \code{MizerSim} object
#' @param effort The effort of each fishing gear. Only needed if the object argument is of class \code{MizerParams}. Can be an two dimensional array (time x gear), a vector of length equal to the number of gears, or a single numeric value (each gear has the same effort). If the object argument is of class \code{MizerSim} then the effort slot of the \code{MizerSim} object is used and this argument is not used. 
#'
#' @return An array. If effort argument has a time dimension, output array has four dimensions (time x gear x species x size). If effort argument does not have a time dimension, output array has three dimensions (gear x species x size).
#' @export
#' @docType methods
#' @rdname getFMortGear-methods
setGeneric('getFMortGear', function(object, effort, ...)
    standardGeneric('getFMortGear'))

#' @rdname getFMortGear-methods
#' @aliases getFMortGear,MizerParams,numeric-method
setMethod('getFMortGear', signature(object='MizerParams', effort = 'numeric'),
    function(object, effort, ...){
	no_gear <- dim(object@catchability)[1]
	# If a single value, just repeat it for all gears
	if(length(effort) == 1)
	    effort <- rep(effort, no_gear)
	if (length(effort) != no_gear)
	    stop("Effort must be a single value or a vector as long as the number of gears\n")
	# turn to array and call next method
	effort <- array(effort,dim=c(1,no_gear))
	fmort_gear <- getFMortGear(object,effort)
	return(fmort_gear)
    }
)

#' @rdname getFMortGear-methods
#' @aliases getFMortGear,MizerParams,matrix-method
setMethod('getFMortGear', signature(object='MizerParams', effort = 'matrix'),
    function(object, effort, ...){
	no_gear <- dim(object@catchability)[1]
	if (dim(effort)[2] != no_gear)
	    stop("Effort array must have a single value or a vector as long as the number of gears for each time step\n")
	# F = sel * q * effort
	sel_q <- sweep(object@selectivity, c(1,2), object@catchability, "*")
	# Kinda nasty! ends up with 4D array 
	fmort_gear <- aaply(effort, 1, function(x,sel_q) sweep(sel_q, c(2,3), x, "*"), sel_q=sel_q)
	return(fmort_gear)
    }
)

# Returns the fishing mortality: time * gear * species * size
#' @rdname getFMortGear-methods
#' @aliases getFMortGear,MizerSim,missing-method
setMethod('getFMortGear', signature(object='MizerSim', effort='missing'),
    function(object, ...){
	f_mort_gear <- getFMortGear(object@params, object@effort, ...)
	return(f_mort_gear)
})

# Total fishing mortality from all gears
# species x size and maybe also by time if effort is time based

#' Get the total fishing mortality from all fishing gears by time, species and size
#'
#' Calculates the fishing mortality from all grars by species and size at each time step in the \code{effort} argument.
#' Fishing mortality from each gear = catchability x selectivity x effort
#' Total fishing mortality is just the sum from each gear
#'
#' @param object A \code{MizerParams} object or an \code{MizerSim} object
#' @param effort The effort of each fishing gear. Only needed if the object argument is of class \code{MizerParams}. Can be an two dimensional array (time x gear), a vector of length equal to the number of gears, or a single numeric value (each gear has the same effort). If the object argument is of class \code{MizerSim} then the effort slot of the \code{MizerSim} object is used and this argument is not used. 
#' @param time_range The time range (either a vector of values, a vector of min and max time, or a single value) to average the abundances over. Default is the whole time range. Only used if \code{object} argument is of type MizerSim.
#'
#' @return An array. If effort argument has a time dimension, output array has three dimensions (time x species x size). If effort argument does not have a time dimension, output array has two dimensions (species x size).
#' @export
#' @docType methods
#' @rdname getFMort-methods
setGeneric('getFMort', function(object, effort, ...)
    standardGeneric('getFMort'))

#' @rdname getFMort-methods
#' @aliases getFMort,MizerParams,numeric-method
setMethod('getFMort', signature(object='MizerParams', effort='numeric'),
    function(object, effort, ...){
	fMortGear <- getFMortGear(object, effort, ...)
	fMort <- apply(fMortGear, c(2,3), sum)
	return(fMort)
})

#' @rdname getFMort-methods
#' @aliases getFMort,MizerParams,matrix-method
setMethod('getFMort', signature(object='MizerParams', effort='matrix'),
    function(object, effort, ...){
	fMortGear <- getFMortGear(object, effort, ...)
	fMort <- apply(fMortGear, c(1,3,4), sum)
	return(fMort)
})

#' @rdname getFMort-methods
#' @aliases getFMort,MizerSim,missing-method
setMethod('getFMort', signature(object='MizerSim', effort='missing'),
    function(object, effort, time_range=dimnames(object@effort)$time, .drop=TRUE, ...){
	time_elements <- get_time_elements(object,time_range, slot="effort")
	fMort <- getFMort(object@params, object@effort, ...)
	return(fMort[time_elements,,,drop=.drop])
})


# get total Z
#' getZ method for the size based model
#'
#' Calculates the total mortality on each prey species by prey size
#' @param object An \code{MizerParams} object
#' @param n A matrix of species abundance (species x size)
#' @param n_pp A vector of the background abundance by size
#' @param effort A numeric vector of the effort by gear or a single numeric effort value which is used for all gears
#'
#' @return A two dimensional array (prey species x prey size) 
#' @export
#' @docType methods
#' @rdname getZ-methods
setGeneric('getZ', function(object, n, n_pp, effort,...)
    standardGeneric('getZ'))

#' @rdname getZ-methods
#' @aliases getZ,MizerParams,matrix,numeric,numeric-method
setMethod('getZ', signature(object='MizerParams', n = 'matrix', n_pp = 'numeric', effort='numeric'),
    function(object, n, n_pp, effort){
	m2 <- getM2(object, n=n, n_pp=n_pp)
	f_mort <- getFMort(object, effort = effort)
	z = sweep(m2 + f_mort,1,object@species_params$z0,"+")
	return(z)
})


# Energy after metabolism and movement
#' getEForReproAndGrowth method for the size based model
#'
#' Calculates the energy available by species and size for reproduction growth after metabolism and movement have been accounted for
#' @param object An \code{MizerParams} object
#' @param n A matrix of species abundance (species x size)
#' @param n_pp A vector of the background abundance by size
#'
#' @return A two dimensional array (prey species x prey size) 
#' @export
#' @docType methods
#' @rdname getEForReproAndGrowth-methods
setGeneric('getEForReproAndGrowth', function(object, n, n_pp, ...)
    standardGeneric('getEForReproAndGrowth'))

#' @rdname getEForReproAndGrowth-methods
#' @aliases getEForReproAndGrowth,MizerParams,matrix,numeric-method
setMethod('getEForReproAndGrowth', signature(object='MizerParams', n = 'matrix', n_pp = 'numeric'),
    function(object, n, n_pp){
	f <- getFeedingLevel(object, n=n, n_pp=n_pp)
	# assimilated intake
	e <- sweep(f * object@intake_max,1,object@species_params$alpha,"*")
	# Subtract basal metabolism and activity 
	e <- e - object@std_metab - object@activity
	e[e<0] <- 0 # Do not allow negative growth
	return(e)
})

# Energy left fot reproduction
# assimilated food intake, less metabolism and activity, split between reproduction and growth

#' getESpawning method for the size based model
#'
#' Calculates the energy available by species and size for reproduction after metabolism and movement have been accounted for
#' @param object An \code{MizerParams} object
#' @param n A matrix of species abundance (species x size)
#' @param n_pp A vector of the background abundance by size
#'
#' @return A two dimensional array (prey species x prey size) 
#' @export
#' @docType methods
#' @rdname getESpawning-methods
setGeneric('getESpawning', function(object,n,n_pp, ...)
    standardGeneric('getESpawning'))

#' @rdname getESpawning-methods
#' @aliases getESpawning,MizerParams,matrix,numeric-method
setMethod('getESpawning', signature(object='MizerParams', n = 'matrix', n_pp = 'numeric'),
    function(object, n, n_pp){
	e <- getEForReproAndGrowth(object,n=n,n_pp=n_pp)
	e_spawning <- object@psi * e 
	return(e_spawning)
})

#' getEGrowth method for the size based model
#'
#' Calculates the energy available by species and size for growth after metabolism and movement have been accounted for
#' @param object An \code{MizerParams} object
#' @param n A matrix of species abundance (species x size)
#' @param n_pp A vector of the background abundance by size
#'
#' @return A two dimensional array (prey species x prey size) 
#' @export
#' @docType methods
#' @rdname getEGrowth-methods
setGeneric('getEGrowth', function(object, n, n_pp, ...)
    standardGeneric('getEGrowth'))


#' @rdname getEGrowth-methods
#' @aliases getEGrowth,MizerParams,matrix,numeric-method
setMethod('getEGrowth', signature(object='MizerParams', n = 'matrix', n_pp = 'numeric'),
    function(object, n, n_pp){
	# Assimilated intake less activity and metabolism
	e <- getEForReproAndGrowth(object,n=n,n_pp=n_pp)
	e_spawning <- getESpawning(object,n=n,n_pp=n_pp)
	# energy for growth is intake - energy for growth
	e_growth <- e - e_spawning
	return(e_growth)
})

#' getRDI method for the size based model
#'
#' Calculates the density independent recruitment (flux entering the smallest size class of each species) before density dependence, by species 
#' @param object An \code{MizerParams} object
#' @param n A matrix of species abundance (species x size)
#' @param n_pp A vector of the background abundance by size
#' @param sex_ratio Proportion of the population that is female. Default value is 0.5
#'
#' @return A numeric vector the length of the number of species 
#' @export
#' @docType methods
#' @rdname getRDI-methods
setGeneric('getRDI', function(object, n, n_pp, ...)
    standardGeneric('getRDI'))

#' @rdname getRDI-methods
#' @aliases getRDI,MizerParams,matrix,numeric-method
setMethod('getRDI', signature(object='MizerParams', n = 'matrix', n_pp = 'numeric'),
    function(object, n, n_pp, sex_ratio = 0.5){
	# Should we put this in the class as part of species_params?
	# Index of the smallest size class for each species
	#w0_idx <- as.vector(tapply(object@species_params$w_min,1:length(object@species_params$w_min),function(w_min,wx) max(which(wx<=w_min)),wx=params@w))
	e_spawning <- getESpawning(object, n=n, n_pp=n_pp)
	e_spawning_pop <- (e_spawning*n) %*% object@dw
	rdi <- sex_ratio*(e_spawning_pop * object@species_params$erepro)/object@w[object@species_params$w_min_idx] 
	return(rdi)
})

#' getRDD method for the size based model
#'
#' Calculates the density dependent recruitment (flux entering the smallest size class of each species) for each species. The density independent recruitment after it has been put through the density dependence function. 
#' @param object An \code{MizerParams} object
#' @param n A matrix of species abundance (species x size)
#' @param n_pp A vector of the background abundance by size
#' @param sex_ratio Proportion of the population that is female. Default value is 0.5
#'
#' @return A numeric vector the length of the number of species 
#' @export
#' @docType methods
#' @rdname getRDD-methods
setGeneric('getRDD', function(object, n, n_pp, ...)
    standardGeneric('getRDD'))

#' @rdname getRDD-methods
#' @aliases getRDD,MizerParams,matrix,numeric-method
setMethod('getRDD', signature(object='MizerParams', n = 'matrix', n_pp = 'numeric'),
    function(object, n, n_pp, sex_ratio = 0.5){
	rdi <- getRDI(object, n=n, n_pp=n_pp, sex_ratio = sex_ratio)
	rdd <- object@srr(rdi = rdi, species_params = object@species_params)
	return(rdd)
})


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

