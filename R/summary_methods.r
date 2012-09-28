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
#'
#' @param object An object of class \code{MizerSim}
#' @param other params to set up size range (min_l, max_l, min_w, max_w)
#' @param min_w minimum weight of species to be used in the calculation
#' @param max_w maximum weight of species to be used in the calculation
#' @param min_l minimum length of species to be used in the calculation
#' @param max_l maximum length of species to be used in the calculation
#'
#' @return An array containing the biomass (time x species)
#' @export
#' @docType methods
#' @rdname getBiomass-methods
# @examples
setGeneric('getBiomass', function(object, ...)
    standardGeneric('getBiomass'))
#' @rdname getSSB-methods
#' @aliases getSSB,MizerSim-method
setMethod('getBiomass', signature(object='MizerSim'),
    function(object, min_w, max_w, min_l, max_l,...){
	#ssb <- apply(sweep(sweep(object@n, c(2,3), object@params@psi,"*"), 3, object@params@w * object@params@dw, "*"),c(1,2),sum) 
	#return(ssb)
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

