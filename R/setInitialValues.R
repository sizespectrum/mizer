#' Set initial values to values from a simulation
#' 
#' This is used to use the results from one simulation as the starting values
#' for another simulation.
#' 
#' The initial abundances (for both species and resource) in the `params`
#' object are set to the abundances in a MizerSim object, averaged over
#' a range of times. Similarly, the initial effort in the `params` object is
#' set to the effort in the MizerSim object, again averaged over that range
#' of times.
#' When no time range is specified, the initial values are taken from the final
#' time step of the simulation.
#' 
#' If the model described by `sim` and `params` has additional components
#' created with [setComponent()] then the values of these components are also
#' averaged and copied to `params`.
#' 
#' The MizerSim object must come from a model with the same set of species and
#' gears and other components and the same size bins as the MizerParams object.
#' Otherwise an error is raised.
#'
#' @param params A `MizerParams` object in which to set the initial values
#' @param sim A `MizerSim` object from which to take the values.
#' @param time_range The time range to average the abundances over. Can be a
#'   vector of values, a vector of min and max time, or a single value. Only the
#'   range of times is relevant, i.e., all times between the smallest and
#'   largest will be selected.  Default is the final time step.
#' @param geometric_mean `r lifecycle::badge("experimental")`
#'   If TRUE then the average of the abundances over the
#'   time range is a geometric mean instead of the default arithmetic mean. This
#'   does not affect the average of the effort or of other components, which is
#'   always arithmetic.
#'   
#' @return The `params` object with updated initial values and initial effort. 
#'   Because of the way the
#'   R language works, `setInitialValues()` does not make the changes to the
#'   params object that you pass to it but instead returns a new params object.
#'   So to affect the change you call the function in the form
#'   `params <- setInitialValues(params, sim)`.
#' @export
#' @family functions for setting parameters
#' @examples
#' \donttest{
#' params <- NS_params
#' sim <- project(params, t_max = 20, effort = 0.5)
#' params <- setInitialValues(params, sim)
#' }
setInitialValues <- function(params, sim, time_range, geometric_mean = FALSE) {
    assert_that(is(params, "MizerParams"),
                is(sim, "MizerSim"),
                is.flag(geometric_mean))
    no_t <- dim(sim@n)[1]
    if (!identical(dim(sim@n)[2:3], dim(params@initial_n))) {
        stop("The consumer size spectrum of the simulation in `sim` has a ",
             "different size from that in `params`.")
    }
    if (!identical(params@species_params$species, 
                   sim@params@species_params$species)) {
        stop("The species in `sim` have different names from those in `params`.")
    }
    if (!identical(length(sim@n_pp[no_t, ]), length(params@initial_n_pp))) {
        stop("The resource size spectrum of the simulation in `sim` has a ",
             "different size from that in `params`.")
    }
    if (!identical(length(sim@n_other[no_t, ]), length(params@initial_n_other))) {
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
    if (missing(time_range)) {
        time_range  <- max(as.numeric(dimnames(sim@n)$time))
    }
    time_elements <- get_time_elements(sim, time_range)
    mean_fn <- mean
    if (geometric_mean) {
        mean_fn <- function(x) {
            exp(mean(log(x)))
        }
    }
    params@initial_n[] <-
        apply(sim@n[time_elements, , , drop = FALSE], c(2, 3), mean_fn)
    params@initial_n_pp[] <-
        apply(sim@n_pp[time_elements, , drop = FALSE], 2, mean_fn)
    # Below we have to work around the fact that "+" does not do
    # componentwise addition of lists
    mizer_add <- function(l1, l2) {
        ifelse(is.list(l1), Map("+", l1, l2), l1 + l2)
    }
    params@initial_n_other[] <-
        apply(sim@n_other[time_elements, , drop = FALSE], 2,
              function(l) Reduce(mizer_add, l) / length(l),
              simplify = FALSE)
    
    params@initial_effort[] <-
        apply(sim@effort[time_elements, , drop = FALSE], 2, mean)
    
    params@time_modified <- lubridate::now()
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
    
    params@time_modified <- lubridate::now()
    params
}

#' @rdname initialN-set
#' @param object An object of class MizerParams or MizerSim
#' @return A matrix with dimensions species x size holding the initial number
#'   densities for the fish spectra.
#' @export
#' @examples 
#' # Doubling abundance of Cod in the initial state of the North Sea model
#' params <- NS_params
#' initialN(params)["Cod", ] <- 2 * initialN(params)["Cod", ]
#' # Calculating the corresponding initial biomass
#' biomass <- initialN(params)["Cod", ] * dw(NS_params) * w(NS_params)
#' # Of course this initial state will no longer be a steady state
#' params <- steady(params)
initialN <- function(object) {
    if (is(object, "MizerParams")) {
        params <- validParams(object)
        return(params@initial_n)
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
#' @examples
#' # Doubling resource abundance in the initial state of the North Sea model
#' params <- NS_params
#' initialNResource(params) <- 2 * initialNResource(params)
#' # Of course this initial state will no longer be a steady state
#' params <- steady(params)
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
    
    params@time_modified <- lubridate::now()
    params
}

#' @rdname initialNResource-set
#' @param object An object of class MizerParams or MizerSim
#' @return A vector with the initial number densities for the resource
#'   spectrum
#' @export
initialNResource <- function(object) {
    if (is(object, "MizerParams")) {
        params <- validParams(object)
        return(params@initial_n_pp)
    }
    if (is(object, "MizerSim")) {
        return(object@params@initial_n_pp)
    }
}
