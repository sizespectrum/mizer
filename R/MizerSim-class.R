# Class specification and constructors for the simulation class

# Copyright 2012 Finlay Scott and Julia Blanchard.
# Copyright 2018 Gustav Delius and Richard Southwell.
# Development has received funding from the European Commission's Horizon 2020 
# Research and Innovation Programme under Grant Agreement No. 634495 
# for the project MINOUW (http://minouw-project.eu/).
# Distributed under the GPL 3 or later 
# Maintainer: Gustav Delius, University of York, <gustav.delius@york.ac.uk>

# Validity check
valid_MizerSim <- function(object) {
    errors <- character()
    validObject(object@params)
    # array dimensions
    if (length(dim(object@n)) != 3) {
	msg <- "n slot must have three dimensions"
	errors <- c(errors, msg)
    }
    if (length(dim(object@effort)) != 2) {
	msg <- "effort slot must have two dimensions"
	errors <- c(errors, msg)
    }
    if (length(dim(object@n_pp)) != 2) {
	msg <- "n_pp slot must have two dimensions"
	errors <- c(errors, msg)
    }
    # Check time dimension is good - size, dim name, and names
    if (!all(c(dim(object@n)[1], 
               dim(object@n_pp)[1],
               length(object@n_other)) == dim(object@effort)[1])) {
	msg <- "First dimension of effort, n, n_pp and n_other slots must be the same length."
	errors <- c(errors, msg)
    }
    if (!all(c(names(dimnames(object@n))[1],
               names(dimnames(object@n_pp))[1],
               names(dimnames(object@effort))[1]) == "time")) {
	msg <- "First dimension of effort, n and n_pp slots must be called 'time'."
	errors <- c(errors, msg)
    }
    # species dimension of n
    if (dim(object@n)[2] != dim(object@params@psi)[1]) {
	msg <- "Second dimension of n slot must have same length as the species names in the params slot"
	errors <- c(errors, msg)
    }
    if (names(dimnames(object@n))[2] != "sp") {
	msg <- "Second dimension of n slot must be called 'sp'"
	errors <- c(errors, msg)
    }
    if (!all(names(dimnames(object@n))[2] == names(dimnames(object@params@psi))[1])) {
	msg <- "Second dimension of n slot must have same species names as in the params slot"
	errors <- c(errors, msg)
    }
    # w dimension of n
    if (dim(object@n)[3] != length(object@params@w)) {
	msg <- "Third dimension of n slot must have same length as w in the params slot"
	errors <- c(errors, msg)
    }
    if (names(dimnames(object@n))[3] != "w") {
	msg <- "Third dimension of n slot must be called 'w'"
	errors <- c(errors, msg)
    }
    if (!all(names(dimnames(object@n))[3] == names(dimnames(object@params@psi))[2])) {
	msg <- "Third dimension of n slot must have same size names as in the params slot"
	errors <- c(errors, msg)
    }
    # w dimension of n_pp
    if (dim(object@n_pp)[2] != length(object@params@w_full)) {
	msg <- "Second dimension of n_pp slot must have same length as w_full in the params slot"
	errors <- c(errors, msg)
    }
    if (names(dimnames(object@n_pp))[2] != "w") {
	msg <- "Second dimension of n_pp slot must be called 'w'"
	errors <- c(errors, msg)
    }
    if (!all(dimnames(object@n_pp)$w == names(object@params@rr_pp))) {
	msg <- "Second dimension of n_pp slot must have same size names as rr_pp in the params slot"
	errors <- c(errors, msg)
    }
    
    # gear dimension of effort
    if (dim(object@effort)[2] != dim(object@params@catchability)[1]) {
	msg <- "Second dimension of effort slot must have same number of gears as in the params slot"
	errors <- c(errors, msg)
    }
    if (names(dimnames(object@effort))[2] != "gear") {
	msg <- "Second dimension of effort slot must be called 'gear'"
	errors <- c(errors, msg)
    }
    if (!all(names(dimnames(object@effort))[2] == names(dimnames(object@params@catchability)[1]))) {
	msg <- "Second dimension of effort slot must have same gear names as in the params slot"
	errors <- c(errors, msg)
    }
    if (length(errors) == 0) TRUE else errors
}

# Soundtrack: Yob - Quantum Mystic
#### Class definition ####
#' A class to hold the results of a simulation
#' 
#' A class that holds the results of projecting a \linkS4class{MizerParams}
#' object through time using the \code{\link{project}} function.
#' 
#' A new \code{MizerSim} object can be created with the \code{\link{MizerSim}}
#' constructor, but you will never have to do that because the object is
#' created automatically by the \code{\link{project}} function when needed.
#' 
#' The arrays all have named dimensions. The names of the "time" dimension are
#' numeric and denote the time in years. The names of the "sp" dimension are the
#' same as the species name in the order specified in the species_params data
#' frame. The names of the "gear" dimension are the names of the gears, in the
#' same order as specified when setting up the \code{MizerParams} object.
#' 
#' There are several \code{link{summary_functions}} and
#' \code{\link{plotting_functions}} available to explore the contents of a
#' \code{MizerSim} object.
#' 
#' @slot params An object of type \linkS4class{MizerParams}.
#' @slot n Array that stores the projected community population abundances by
#'   time, species and size.
#' @slot n_pp Array that stores the projected plankton abundance by time and
#'   size.
#' @slot n_other List of lists where the outer list is named by time and the
#'   inner lists are named by the names of the other ecosystem components.
#' @slot effort Array that stores the fishing effort by time and gear.
#' 
#' @export
setClass(
    "MizerSim",
    slots = c(
        params = "MizerParams",
        n = "array",
        effort = "array",
        n_pp = "array",
        n_other = "list"
    )
)

setValidity("MizerSim", valid_MizerSim)
remove(valid_MizerSim)


#' Constructor for the \code{MizerSim} class
#' 
#' A constructor for the \code{MizerSim} class. This is used by the
#' \code{project} function to create \code{MizerSim} objects of the right
#' dimensions. It is not necessary for users to use this constructor.
#' 
#' @param params a \linkS4class{MizerParams} object
#' @param t_dimnames Numeric vector that is used for the time dimensions of the
#'   slots. Default = NA.
#' @param t_max The maximum time step of the simulation. Only used if t_dimnames
#'   = NA. Default value = 100.
#' @param t_save How often should the results of the simulation be stored. Only
#'   used if t_dimnames = NA. Default value = 1.
#'   
#' @return An object of type \linkS4class{MizerSim}
#' @export
MizerSim <- function(params, t_dimnames = NA, t_max = 100, t_save = 1) {
    # If the dimnames for the time dimension not passed in, calculate them
    # from t_max and t_save
    if (any(is.na(t_dimnames))) {
        t_dimnames <- seq(from = 0, to = t_max, by = t_save)
    }
    if (!is.numeric(t_dimnames)) {
        stop("The t_dimnames argument must be numeric.")
    }
    if (is.unsorted(t_dimnames)) {
        stop("The t_dimnames argument should be increasing.")
    }
    no_sp <- nrow(params@species_params)
    species_names <- dimnames(params@psi)$sp
    no_w <- length(params@w)
    w_names <- dimnames(params@psi)$w
    no_t <- length(t_dimnames)
    array_n <- array(NA, dim = c(no_t, no_sp, no_w), 
                     dimnames = list(time = t_dimnames, 
                                     sp = species_names, w = w_names))
    
    no_gears <- dim(params@selectivity)[1]
    gear_names <- dimnames(params@selectivity)$gear
    array_effort <- array(NA, dim = c(no_t, no_gears), 
                          dimnames = list(time = t_dimnames, 
                                          gear = gear_names))
    
    no_w_full <- length(params@w_full)
    w_full_names <- names(params@rr_pp)
    array_n_pp <- array(NA, dim = c(no_t, no_w_full), 
                        dimnames = list(time = t_dimnames, 
                                        w = w_full_names))
    
    other_names <- names(params@other_dynamics)
    no_other <- length(other_names)
    other <- rep(list(NA), no_other)
    names(other) <- other_names
    list_n_other <- rep(list(other), no_t)
    names(list_n_other) <- t_dimnames
    
    sim <- new('MizerSim',
               params = params,
               n = array_n,
               n_pp = array_n_pp,
               n_other = list_n_other,
               effort = array_effort)
    return(sim)
}

#' Consumer size spectra at end of simulation
#' 
#' @param sim A MizerSim object
#' @return A matrix (species x size) holding the consumer number densities at
#'   the end of the simulation
#' @export
final_n <- function(sim) {
    assert_that(is(sim, "MizerSim"))
    n <- sim@params@initial_n  # Needed to get the right dimnames
    n[] <- sim@n[dim(sim@n)[1], , ]
    n
}

#' Plankton size spectrum at end of simulation
#' 
#' @param sim A MizerSim object
#' @return A vector holding the plankton number densities at the end of the
#'   simulation for all size classes
#' @export
final_n_pp <- function(sim) {
    assert_that(is(sim, "MizerSim"))
    sim@n_pp[dim(sim@n_pp)[1], ]
}

#' Values of other ecosystem components at end of simulation
#' 
#' @param sim A MizerSim object
#' @return A named list holding the values of other ecosystem components at the
#'   end of the simulation
#' @export
final_n_other <- function(sim) {
    assert_that(is(sim, "MizerSim"))
    sim@n_other[[length(sim@n_other)]]
}

#' Times for which simulation results are available
#' 
#' @param sim A MizerSim object
#' @return A numeric vectors of the times (in years) at which simulation results
#'   have been stored in the MizerSim object.
#' @export
times <- function(sim) {
    as.numeric(names(sim@n_other))
}