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
    if (!all(c(dim(object@n)[[1]], 
               dim(object@n_pp)[[1]],
               dim(object@n_other)[[1]]) == dim(object@effort)[[1]])) {
	msg <- "First dimension of effort, n, n_pp and n_other slots must be the same length."
	errors <- c(errors, msg)
    }
    if (!all(c(names(dimnames(object@n))[[1]],
               names(dimnames(object@n_pp))[[1]],
               names(dimnames(object@n_other))[[1]],
               names(dimnames(object@effort))[[1]]) == "time")) {
	msg <- "First dimension of effort, n, n_pp and n_other slots must be called 'time'."
	errors <- c(errors, msg)
    }
    # species dimension of n
    if (dim(object@n)[[2]] != dim(object@params@psi)[[1]]) {
	msg <- "Second dimension of n slot must have same length as the species names in the params slot"
	errors <- c(errors, msg)
    }
    if (names(dimnames(object@n))[[2]] != "sp") {
	msg <- "Second dimension of n slot must be called 'sp'"
	errors <- c(errors, msg)
    }
    if (!all(names(dimnames(object@n))[[2]] == names(dimnames(object@params@psi))[[1]])) {
	msg <- "Second dimension of n slot must have same species names as in the params slot"
	errors <- c(errors, msg)
    }
    # w dimension of n
    if (dim(object@n)[[3]] != length(object@params@w)) {
	msg <- "Third dimension of n slot must have same length as w in the params slot"
	errors <- c(errors, msg)
    }
    if (names(dimnames(object@n))[[3]] != "w") {
	msg <- "Third dimension of n slot must be called 'w'"
	errors <- c(errors, msg)
    }
    if (!all(names(dimnames(object@n))[[3]] == names(dimnames(object@params@psi))[[2]])) {
	msg <- "Third dimension of n slot must have same size names as in the params slot"
	errors <- c(errors, msg)
    }
    # w dimension of n_pp
    if (dim(object@n_pp)[[2]] != length(object@params@w_full)) {
	msg <- "Second dimension of n_pp slot must have same length as w_full in the params slot"
	errors <- c(errors, msg)
    }
    if (names(dimnames(object@n_pp))[[2]] != "w") {
	msg <- "Second dimension of n_pp slot must be called 'w'"
	errors <- c(errors, msg)
    }
    if (!all(dimnames(object@n_pp)$w == names(object@params@rr_pp))) {
	msg <- "Second dimension of n_pp slot must have same size names as rr_pp in the params slot"
	errors <- c(errors, msg)
    }
    # component dimension of n_other
    if (dim(object@n_other)[[2]] != length(object@params@other_dynamics)) {
        msg <- "Second dimension of n_other slot must have same length as other_dynamics in the params slot"
        errors <- c(errors, msg)
    }
    if (names(dimnames(object@n_other))[[2]] != "component") {
        msg <- "Second dimension of n_other slot must be called 'component'"
        errors <- c(errors, msg)
    }
    if (!all(dimnames(object@n_pp)$component == names(object@params@other_dynamics))) {
        msg <- "Second dimension of n_other slot must have same component names as other_dynamics in the params slot"
        errors <- c(errors, msg)
    }
    
    # gear dimension of effort
    if (dim(object@effort)[[2]] != dim(object@params@catchability)[[1]]) {
	msg <- "Second dimension of effort slot must have same number of gears as in the params slot"
	errors <- c(errors, msg)
    }
    if (names(dimnames(object@effort))[[2]] != "gear") {
	msg <- "Second dimension of effort slot must be called 'gear'"
	errors <- c(errors, msg)
    }
    if (!all(names(dimnames(object@effort))[[2]] == names(dimnames(object@params@catchability)[[1]]))) {
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
#' object through time using [project()].
#' 
#' A new `MizerSim` object can be created with the [MizerSim()]
#' constructor, but you will never have to do that because the object is
#' created automatically by [project()] when needed.
#' 
#' As a user you should never have to access the slots of a MizerSim object
#' directly. Instead there are a range of functions to extract the information.
#' [N()] and [NResource()] return arrays with the saved abundances of
#' the species and the resource population at size respectively. [effort()]
#' returns the fishing effort of each gear through time. 
#' [times()] returns the vector of times at which simulation results
#' were stored and [idxFinalT()] returns the index with which to access
#' specifically the value at the final time in the arrays returned by the other
#' functions. T[params()] returns the `MizerParams` object that was
#' passed in to `project()`. There are also several
#' [summary_functions] and [plotting_functions]
#' available to explore the contents of a `MizerSim` object.
#' 
#' The arrays all have named dimensions. The names of the `time` dimension
#' denote the time in years. The names of the `w` dimension are weights in grams
#' rounded to three significant figures. The names of the `sp` dimension are the
#' same as the species name in the order specified in the species_params data
#' frame. The names of the `gear` dimension are the names of the gears, in the
#' same order as specified when setting up the `MizerParams` object.
#' 
#' Extensions of mizer can use the `n_other` slot to store the abundances of
#' other ecosystem components and these extensions should provide their own
#' functions for accessing that information.
#' 
#' The `MizerSim` class has changed since previous versions of mizer. To use
#' a `MizerSim` object created by a previous version, you need to upgrade it
#' with [upgradeSim()].
#' 
#' @slot params An object of type \linkS4class{MizerParams}.
#' @slot n Three-dimensional array (time x species x size) that stores the 
#'   projected community number densities.
#' @slot n_pp An array (time x size) that stores the projected resource number
#'   densities.
#' @slot n_other A list array (time x component) that stores the projected
#'   values for other ecosystem components.
#' @slot effort An array (time x gear) that stores the fishing effort by time and 
#'   gear.
#' 
#' @export
setClass(
    "MizerSim",
    slots = c(
        params = "MizerParams",
        n = "array",
        effort = "array",
        n_pp = "array",
        n_other = "array"
    )
)

setValidity("MizerSim", valid_MizerSim)
remove(valid_MizerSim)


#' Constructor for the `MizerSim` class
#' 
#' A constructor for the `MizerSim` class. This is used by 
#' [project()] to create `MizerSim` objects of the right
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
    
    component_names <- names(params@other_dynamics)
    no_components <- length(component_names)
    list_n_other <- rep(list(NA), no_t * no_components)
    dim(list_n_other) <- c(no_t, no_components)
    dimnames(list_n_other) <- list(time = t_dimnames,
                                   component = component_names)
    
    sim <- new('MizerSim',
               params = params,
               n = array_n,
               n_pp = array_n_pp,
               n_other = list_n_other,
               effort = array_effort)
    return(sim)
}

#' Time series of size spectra
#' 
#' Fetch the simulation results for the size spectra over time.
#' 
#' @param sim A MizerSim object
#' @return For `N()`: A three-dimensional array (time x species x size) with the
#'   number density of consumers
#' @export
N <- function(sim) {
    sim@n
}

#' @rdname N
#' @return For `NResource()`: An array (time x size) with the number density of resource
#' @export
NResource <- function(sim) {
    sim@n_pp
}


#' Size spectra at end of simulation
#' 
#' @param sim A MizerSim object
#' @return For `finalN()`: An array (species x size) holding the consumer
#'   number densities at the end of the simulation
#' @export
finalN <- function(sim) {
    assert_that(is(sim, "MizerSim"))
    n <- sim@params@initial_n  # Needed to get the right dimnames
    n[] <- sim@n[dim(sim@n)[[1]], , ]
    n
}

#' @rdname finalN
#' @return For `finalNResource()`: A vector holding the resource number densities at
#'   the end of the simulation for all size classes
#' @export
finalNResource <- function(sim) {
    assert_that(is(sim, "MizerSim"))
    sim@n_pp[dim(sim@n_pp)[[1]], ]
}

#' Time index at end of simulation
#' 
#' @param sim A MizerSim object
#' @return An integer giving the index for extracting the
#'   results for the final time step
#' @export
#' @examples
#' \dontrun{
#' sim <- project(NS_params, t_max = 12, t_save = 0.5)
#' idx <- idxFinalT(sim)
#' idx
#' # This coincides with
#' length(times(sim))
#' # and corresponds to the final time
#' times(sim)[idx]
#' # We can use this index to extract the result at the final time
#' identical(N(sim)[idx, , ], finalN(sim))
#' identical(NResource(sim)[idx, ], finalNResource(sim))
#' }
idxFinalT <- function(sim) {
    assert_that(is(sim, "MizerSim"))
    dim(sim@n_pp)[[1]]
}


#' Times for which simulation results are available
#' 
#' @param sim A MizerSim object
#' @return A numeric vectors of the times (in years) at which simulation results
#'   have been stored in the MizerSim object.
#' @export
times <- function(sim) {
    as.numeric(dimnames(sim@n)$t)
}

#' Fishing effort used in simulation
#' 
#' Note that the array returned may not be exactly the same as the `effort`
#' argument that was passed in to `project()`. This is because only the saved
#' effort is stored (the frequency of saving is determined by the argument
#' `t_save`).
#' 
#' @param sim A MizerSim object
#' @return An array (time x gear) that stores the fishing effort by time and 
#'   gear.
#' @export
effort <- function(sim) {
    sim@effort
}

#' Extract the parameter object underlying a simulation
#' 
#' @param sim A MizerSim object
#' @return The MizerParams object that was used to run the simulation
#' @export
params <- function(sim) {
    sim@params
}
