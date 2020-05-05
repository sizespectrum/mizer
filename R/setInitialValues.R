#' Set initial values to final values of a simulation
#' 
#' Takes the final values from a simulation in a MizerSim object and stores them
#' as initial values in a MizerParams object.
#'
#' @param params A [MizerParams()] object
#' @param sim A `MizerSim` object.
#'   
#' @return The `params` object with updated initial values and initial effort, 
#'   taken from the
#'   values at the final time of the simulation in `sim`. Because of the way the
#'   R language works, `setInitialValues()` does not make the changes to the
#'   params object that you pass to it but instead returns a new params object.
#'   So to affect the change you call the function in the form
#'   `params <- setInitialValues(params, sim)`.
#' @export
#' @family functions for setting parameters
#' @examples
#' \dontrun{
#' params <- NS_params
#' sim <- project(params, t_max = 20, effort = 0.5)
#' params <- setInitialValues(params, sim)
#' }
setInitialValues <- function(params, sim) {
    assert_that(is(params, "MizerParams"),
                is(sim, "MizerSim"))
    no_t <- dim(sim@n)[1]
    if (!identical(dim(sim@n)[2:3], dim(params@initial_n))) {
        stop("The consumer size spectrum of the simulation in `sim` has a ",
             "different size from that in `params`.")
    }
    if (!identical(dim(sim@n_pp[no_t, ]), dim(params@initial_n_pp))) {
        stop("The resource size spectrum of the simulation in `sim` has a ",
             "different size from that in `params`.")
    }
    if (!identical(length(sim@n_pp[no_t, ]), length(params@initial_n_pp))) {
        stop("The number of other components in the simulation in `sim` is ",
             "different from that in `params`.")
    }
    if (!identical(length(sim@effort[no_t, ]), length(params@initial_effort))) {
        stop("The number of gears in the simulation in `sim` is ",
             "different from that in `params`.")
    }
    if (!identical(dimnames(sim@effort)[[2]], names(params@initial_effort))) {
        stop("The gears in the simulation in `sim` have different names ",
             "from those in `params`.")
    }
    params@initial_n[] <- sim@n[no_t, , ]
    params@initial_n_pp[] <- sim@n_pp[no_t, ]
    params@initial_n_other <- sim@n_other[no_t, ]
    params@initial_effort <- sim@effort[no_t, ]
    params
}

#' Initial values for fish spectra
#' 
#' Values used as starting values for simulations with `project()`.
#' 
#' @param params A MizerParams object
#' @param value A matrix with dimensions species x size holding the initial
#'   number densities for the fish spectra.
#' @export
`initialN<-` <- function(params, value) {
    if (!is(params, "MizerParams")) {
        stop("You can only assign an initial N to a MizerParams object. ",
             params, " is of class ", class(params), ".")
    }
    assert_that(identical(dim(value), dim(params@initial_n)),
                all(value >= 0))
    if (!is.null(dimnames(value)) &&
        !identical(dimnames(value), dimnames(params@initial_n))) {
        warning("The dimnames do not match. I will ignore them.")
    }
    params@initial_n[] <- value
    params
}

#' @rdname initialN-set
#' @param object An object of class MizerParams or MizerSim
#' @export
initialN <- function(object) {
    if (is(object, "MizerParams")) {
        return(object@initial_n)
    }
    if (is(object, "MizerSim")) {
        return(object@params@initial_n)
    }
}

#' Initial value for resource spectrum
#' 
#' Value used as starting value for simulations with `project()`.
#' 
#' @param params A MizerParams object
#' @param value A vector with the initial number densities for the resource
#'   spectrum
#' @export
`initialNResource<-` <- function(params, value) {
    if (!is(params, "MizerParams")) {
        stop("You can only assign an initial N to a MizerParams object. ",
             params, " is of class ", class(params), ".")
    }
    assert_that(identical(dim(value), dim(params@initial_n_pp)),
                all(value >= 0))
    if (!is.null(dimnames(value)) &&
        !identical(dimnames(value), dimnames(params@initial_n_pp))) {
        warning("The dimnames do not match. I will ignore them.")
    }
    params@initial_n_pp[] <- value
    params
}

#' @rdname initialNResource-set
#' @param object An object of class MizerParams or MizerSim
#' @export
initialNResource <- function(object) {
    if (is(object, "MizerParams")) {
        return(object@initial_n_pp)
    }
    if (is(object, "MizerSim")) {
        return(object@params@initial_n_pp)
    }
}
