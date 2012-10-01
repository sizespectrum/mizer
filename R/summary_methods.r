# Summary methods for mizer package

# Copyright 2012 Finlay Scott and Julia Blanchard. 
# Distributed under the GPL 2 or later
# Maintainer: Finlay Scott, CEFAS. finlay.scott@cefas.co.uk

# SSB
# Biomass
# N
# Yield

#' Calculate the SSB of species
#'
#' Calculates the spawning stock biomass (SSB) through time of the species in the \code{MizerSim} class.
#' SSB is calculated as the total mass of all mature individuals
#'
#' @param object An object of class \code{MizerSim}
#'
#' @return An array containing the SSB (time x species)
#' @export
#' @docType methods
#' @rdname getSSB-methods
# @examples
setGeneric('getSSB', function(object, ...)
    standardGeneric('getSSB'))
#' @rdname getSSB-methods
#' @aliases getSSB,MizerSim-method
setMethod('getSSB', signature(object='MizerSim'),
    function(object,  ...){
	ssb <- apply(sweep(sweep(object@n, c(2,3), object@params@psi,"*"), 3, object@params@w * object@params@dw, "*"),c(1,2),sum) 
	return(ssb)
    })

#' Calculate the biomass of species within a size range
#'
#' Calculates the total biomass through time of the species in the \code{MizerSim} class within user defined size limits.
#' The default option is to use the whole size range
#' You can specify minimum and maximum weight or lengths for the species. Lengths take precedence over weights (i.e. if both min_l and min_w are supplied, only min_l will be used)
#'
#' @param object An object of class \code{MizerSim}
#' @param min_w minimum weight of species to be used in the calculation
#' @param max_w maximum weight of species to be used in the calculation
#' @param min_l minimum length of species to be used in the calculation
#' @param max_l maximum length of species to be used in the calculation
#'
#' @return An array containing the biomass (time x species)
#' @export
#' @docType methods
#' @rdname getBiomass-methods
#' @examples
#' data(species_params_gears)
#' data(inter)
#' params <- MizerParams(species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=10)
#' n <- getN(sim, min_l = 10)
setGeneric('getBiomass', function(object, ...)
    standardGeneric('getBiomass'))
#' @rdname getBiomass-methods
#' @aliases getBiomass,MizerSim-method
setMethod('getBiomass', signature(object='MizerSim'),
    function(object, ...){
	size_range <- get_size_range_array(object@params,...)
	biomass <- apply(sweep(sweep(object@n,c(2,3),size_range,"*"),3,object@params@w * object@params@dw, "*"),c(1,2),sum)
	return(biomass)
    })

#' Calculate the total abundance in terms of numbers of species within a size range
#'
#' Calculates the total numbers through time of the species in the \code{MizerSim} class within user defined size limits.
#' The default option is to use the whole size range
#' You can specify minimum and maximum weight or lengths for the species. Lengths take precedence over weights (i.e. if both min_l and min_w are supplied, only min_l will be used)
#'
#' @param object An object of class \code{MizerSim}
#' @param min_w minimum weight of species to be used in the calculation
#' @param max_w maximum weight of species to be used in the calculation
#' @param min_l minimum length of species to be used in the calculation
#' @param max_l maximum length of species to be used in the calculation
#'
#' @return An array containing the total numbers (time x species)
#' @export
#' @docType methods
#' @rdname getN-methods
#' @examples
#' data(species_params_gears)
#' data(inter)
#' params <- MizerParams(species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=10)
#' n <- getN(sim, min_l = 10)
setGeneric('getN', function(object, ...)
    standardGeneric('getN'))
#' @rdname getN-methods
#' @aliases getN,MizerSim-method
setMethod('getN', signature(object='MizerSim'),
    function(object, ...){
	size_range <- get_size_range_array(object@params,...)
	n <- apply(sweep(sweep(sim@n,c(2,3),size_range,"*"),3,sim@params@dw, "*"),c(1,2),sum)
	return(n)
    })

# Helper function that returns an array (no_sp x no_w) of boolean values indicating whether that size bin is within
# the size limits specified by the arguments
# If min_l or max_l are supplied they take precendence over the min_w and max_w
# But you can mix min_l and max_w etc
# Not exported
get_size_range_array <- function(params, min_w = min(params@w), max_w = max(params@w), min_l = NULL, max_l = NULL){
    no_sp <- nrow(params@species_params)
    if(is.null(min_w) | is.null(min_w))
	if (any(!c("a","b") %in% names(params@species_params)))
	    stop("species_params slot must have columns 'a' and 'b' for length-weight conversion")
    if(!is.null(min_l))
	min_w <- params@species_params$a * min_l ^ params@species_params$b
    else min_w <- rep(min_w,no_sp)
    if(!is.null(max_l))
	max_w <- params@species_params$a * max_l ^ params@species_params$b
    else max_w <- rep(max_w,no_sp)

    min_n <- aaply(min_w, 1, function(x) params@w >= x)
    max_n <- aaply(max_w, 1, function(x) params@w <= x)
    size_n <- min_n & max_n
    # Add dimnames?
    dimnames(size_n) <- list(sp = params@species_params$species, w = signif(params@w,3)) 
    return(size_n)
}

#' Summary method 
#'
#' Outputs a general summary of the structure and content of the object
#'
#' @export
#' @docType methods
#' @rdname summary-methods
#' @aliases summary,MizerParams-method
#'
#' @examples
#' params <- MizerParams(object=3, species_names = c("cod", "haddock", "whiting"))
#' summary(params)
setMethod("summary", signature(object="MizerParams"),
    function(object, ...){
	#cat("An object of class \"", as.character(class(object)), "\" with:\n", sep="")
	cat("An object of class \"", as.character(class(object)), "\" \n", sep="")
	cat("Community size spectrum:\n")
	cat("\tminimum size:\t", signif(min(object@w)), "\n", sep="")
	cat("\tmaximum size:\t", signif(max(object@w)), "\n", sep="")
	cat("\tno. size bins:\t", length(object@w), "\n", sep="")
	# Length of background? 
	cat("Background size spectrum:\n")
	cat("\tminimum size:\t", signif(min(object@w_full)), "\n", sep="")
	cat("\tmaximum size:\t", signif(max(object@w_full)), "\n", sep="")
	cat("\tno. size bins:\t", length(object@w_full), "\n", sep="")
	# w range - min, max, number of w
	# w background min max
	# no species and names and wInf,  - not all these wMat, beta, sigma
	# no gears, gear names catching what
	cat("Species details:\n")
	cat("\tSpecies\t\tw_inf\n")
	for (i in 1:nrow(object@species_params))
	    cat("\t",as.character(object@species_params$species)[i], "\t\t ",signif(object@species_params$w_inf[i],3), "\n", sep="")
	cat("Fishing gear details:\n")
	cat("\tGear\t\t\tTarget species\n")
	for (i in 1:dim(object@catchability)[1]){
	    cat("\t",dimnames(object@catchability)$gear[i], "\t\t",dimnames(params@catchability)$sp[params@catchability[i,]>0], "\n", sep=" ") 
	}
})

#' @rdname summary-methods
#' @aliases summary,MizerSim-method
#' @examples
#' params <- MizerParams(object=3, species_names = c("cod", "haddock", "whiting"))
#' sim <- project(params, effort=1, t_max=5)
#' summary(sim)
setMethod("summary", signature(object="MizerSim"),
    function(object, ...){
	#cat("An object of class \"", as.character(class(object)), "\" with:\n", sep="")
	cat("An object of class \"", as.character(class(object)), "\" \n", sep="")
	cat("Parameters:\n")
	summary(object@params)
	cat("Simulation parameters:\n")
	# Need to store t_max and dt in a description slot? Or just in simulation time parameters? Like a list?
	cat("\tFinal time step: ", max(as.numeric(dimnames(object@n)$time)), "\n", sep="")
	cat("\tOutput stored every ", as.numeric(dimnames(object@n)$time)[2] - as.numeric(dimnames(object@n)$time)[1], " time units\n", sep="")
})


