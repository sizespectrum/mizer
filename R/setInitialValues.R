#' Set initial values to final values of a simulation
#' 
#' Takes the final values from a simulation in a MizerSim object and stores them
#' as initial values in a MizerParams object.
#'
#' @param params A [MizerParams()] object
#' @param sim A `MizerSim` object.
#'   
#' @return The `params` object with updated initial values, taken from the
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
    assert_that(identical(dim(sim@n[no_t, , ]), dim(params@initial_n)),
                identical(dim(sim@n_pp[no_t, ]), dim(params@initial_n_pp)),
                identical(dim(sim@n_other[no_t, ]), length(params@initial_n_other)))
    params@initial_n[] <- sim@n[no_t, , ]
    params@initial_n_pp[] <- sim@n_pp[no_t, ]
    params@initial_n_other <- sim@n_other[no_t, ]
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
`initial_N<-` <- function(params, value) {
    assert_that(is(params, "MizerParams"),
                identical(dim(value), dim(params@initial_n)),
                all(value >= 0))
    if (!is.null(dimnames(value)) &&
        !identical(dimnames(value), dimnames(params@initial_n))) {
        warning("The dimnames do not match. I will ignore them.")
    }
    params@initial_n[] <- value
    params
}

#' @rdname initial_N-set
#' @export
initial_N <- function(params) {
    params@initial_n
}

#' Initial value for plankton spectrum
#' 
#' Value used as starting value for simulations with `project()`.
#' 
#' @param params A MizerParams object
#' @param value A vector with the initial number densities for the plankton
#'   spectrum
#' @export
`initial_N_pp<-` <- function(params, value) {
    assert_that(is(params, "MizerParams"),
                identical(dim(value), dim(params@initial_n_pp)),
                all(value >= 0))
    if (!is.null(dimnames(value)) &&
        !identical(dimnames(value), dimnames(params@initial_n_pp))) {
        warning("The dimnames do not match. I will ignore them.")
    }
    params@initial_n_pp[] <- value
    params
}

#' @rdname initial_N_pp-set
#' @export
initial_N_pp <- function(params) {
    params@initial_n_pp
}
