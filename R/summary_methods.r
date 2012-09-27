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
#' @param other params to set up size range
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


