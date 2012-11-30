# Class specification and constructors for the simulation class

# Copyright 2012 Finlay Scott and Julia Blanchard. 
# Distributed under the GPL 2 or later
# Maintainer: Finlay Scott, CEFAS. finlay.scott@cefas.co.uk


# Validity check
valid_MizerSim <- function(object){
    errors <- character()
    validObject(object@params)
    # array dimensions
    if(length(dim(object@n)) != 3){
	msg <- "n slot must have three dimensions"
	errors <- c(errors, msg)
    }
    if(length(dim(object@effort)) != 2){
	msg <- "effort slot must have two dimensions"
	errors <- c(errors, msg)
    }
    if(length(dim(object@n_pp)) != 2){
	msg <- "n_pp slot must have two dimensions"
	errors <- c(errors, msg)
    }
    # Check time dimension is good - size, dim name, and names
    if(!all(c(dim(object@n)[1]-1,dim(object@n_pp)[1]-1) == dim(object@effort)[1])){
	msg <- "First dimension of effort, n and n_pp slots must be the same length (n and n_pp are + 1 longer)"
	errors <- c(errors, msg)
    }
    if(!all(c(names(dimnames(object@n))[1], names(dimnames(object@n_pp))[1], names(dimnames(object@effort))[1]) == "time")){
	msg <- "First dimension of effort, n and n_pp slots must be called 'time'"
	errors <- c(errors, msg)
    }
    if(!all(c(names(dimnames(object@n))[1], names(dimnames(object@n_pp))[1]) == 
	names(dimnames(object@effort))[1])){
	msg <- "First dimension of effort, n and n_pp slots must have the same names"
	errors <- c(errors, msg)
    }
    # species dimension of n
    if(dim(object@n)[2] != dim(object@params@psi)[1]){
	msg <- "Second dimension of n slot must have same length as the species names in the params slot"
	errors <- c(errors, msg)
    }
    if(names(dimnames(object@n))[2] != "sp"){
	msg <- "Second dimension of n slot must be called 'sp'"
	errors <- c(errors, msg)
    }
    if(!all(names(dimnames(object@n))[2] == names(dimnames(object@params@psi))[1])){
	msg <- "Second dimension of n slot must have same species names as in the params slot"
	errors <- c(errors, msg)
    }
    # w dimension of n
    if(dim(object@n)[3] != length(object@params@w)){
	msg <- "Third dimension of n slot must have same length as w in the params slot"
	errors <- c(errors, msg)
    }
    if(names(dimnames(object@n))[3] != "w"){
	msg <- "Third dimension of n slot must be called 'w'"
	errors <- c(errors, msg)
    }
    if(!all(names(dimnames(object@n))[3] == names(dimnames(object@params@psi))[2])){
	msg <- "Third dimension of n slot must have same size names as in the params slot"
	errors <- c(errors, msg)
    }
    # w dimension of n_pp
    if(dim(object@n_pp)[2] != length(object@params@w_full)){
	msg <- "Second dimension of n_pp slot must have same length as w_full in the params slot"
	errors <- c(errors, msg)
    }
    if(names(dimnames(object@n_pp))[2] != "w"){
	msg <- "Second dimension of n_pp slot must be called 'w'"
	errors <- c(errors, msg)
    }
    if(!all(dimnames(object@n_pp)$w == names(object@params@rr_pp))){
	msg <- "Second dimension of n_pp slot must have same size names as rr_pp in the params slot"
	errors <- c(errors, msg)
    }
    # gear dimension of effort
    if(dim(object@effort)[2] != dim(object@params@catchability)[1]){
	msg <- "Second dimension of effort slot must have same number of gears as in the params slot"
	errors <- c(errors, msg)
    }
    if(names(dimnames(object@effort))[2] != "gear"){
	msg <- "Second dimension of effort slot must be called 'gear'"
	errors <- c(errors, msg)
    }
    if(!all(names(dimnames(object@effort))[2] == names(dimnames(object@params@catchability)[1]))){
	msg <- "Second dimension of effort slot must have same gear names as in the params slot"
	errors <- c(errors, msg)
    }
    if (length(errors) == 0) TRUE else errors
}




#' MizerSim
#' 
#' A class that holds the results of a projecting an \code{MizerParams} object through time
#'
#' \section{Slots}{
#'     \describe{
#'         \item{params}{An object of type \code{MizerParams}}
#'         \item{n}{Array to store the projected community populations by time, species and size}
#'         \item{effort}{Array that stores the fishing effort through time by time and gear. SEE DETAILS}
#'         \item{n_pp}{Array to store the projected background population by time and size}
#'     }
#' }
#' @name MizerSim-class
#' @rdname MizerSim-class
#' @docType class
#' @export
setClass("MizerSim",
    representation(
	params = "MizerParams",
	n = "array",
	effort = "array",
	n_pp = "array"
    ),
    prototype = prototype(
	params = new("MizerParams"),
	n = array(NA,dim=c(1,1,1), dimnames = list(time = NULL, sp = NULL, w = NULL)),
	effort = array(NA,dim=c(1,1), dimnames = list(time = NULL, gear = NULL)),
	n_pp = array(NA,dim=c(1,1), dimnames = list(time = NULL, w = NULL))
    ),
    validity = valid_MizerSim
)

setValidity("MizerSim", valid_MizerSim)
remove(valid_MizerSim)


# Constructors

#' Constructors for objects of \code{MizerSim} class
#'
#' A range of constructors for the \code{MizerSim} class. Some are more user friendly than others and some are really only used 
#' internally to help construct \code{MizerSim objects} of the right dimension
#' 
#' @param object a MizerParams object 
#' @param ... other arguments
#'
#' @return An object of type \code{MizerSim}
#' @note Something on the constructor
#' @export
#' @docType methods
#' @rdname MizerSim-methods
#'
setGeneric('MizerSim', function(object, ...)
    standardGeneric('MizerSim'))

#' @rdname MizerSim-methods
#' @aliases MizerSim,MizerParams-method
setMethod('MizerSim', signature(object='MizerParams'),
    function(object, t_max = 100, t_save=1, ...){
	if((t_max %% t_save) != 0)
	    stop("t_max must be divisible by t_save with no remainder")
	no_sp <- nrow(object@species_params)
	species_names <- dimnames(object@psi)$sp
	no_w <- length(object@w)
	w_names <- dimnames(object@psi)$w
	t_dim_n <- 1 + (t_max / t_save) # N is 1 bigger because it holds the initial population
	t_dim_effort <- (t_max / t_save)
	array_n <- array(NA, dim = c(t_dim_n, no_sp, no_w), dimnames = list(time = seq(from = 0, to = t_max, by = t_save), sp = species_names, w = w_names))

	no_gears <- dim(object@selectivity)[1]
	gear_names <- dimnames(object@selectivity)$gear
	array_effort <- array(NA, dim = c(t_dim_effort, no_gears), dimnames = list(time = seq(from = t_save, to = t_max, by = t_save), gear = gear_names))

	no_w_full <- length(object@w_full)
	w_full_names <- names(object@rr_pp)
	array_n_pp <- array(NA, dim = c(t_dim_n, no_w_full), dimnames = list(time=seq(from=0, to = t_max, by = t_save), w = w_full_names))

	sim <- new('MizerSim',
		   n = array_n, 
		   effort = array_effort,
		   n_pp = array_n_pp,
		   params = object)
	return(sim)
    }
)

