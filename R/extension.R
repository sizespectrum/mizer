#' Set own rate function to replace mizer rate function
#' 
#' At each time step during a simulation with the [project()] function, mizer
#' needs to calculate the instantaneous values of the various rates. By
#' default it calls the [mizerRates()] function which creates a list with the
#' following components:
#' * `encounter` from [mizerEncounter()]
#' * `feeding_level` from [mizerFeedingLevel()]
#' * `pred_rate` from [mizerPredRate()]
#' * `pred_mort` from [mizerPredMort()]
#' * `fishing_mort` from [mizerFMort()]
#' * `mort` from [mizerMort()]
#' * `plankton_mort` from [mizerPlanktonMort()]
#' * `e` from [mizerEReproAndGrowth()]
#' * `e_repro` from [mizerERepro()]
#' * `e_growth` from [mizerEGrowth()]
#' * `rdi` from [mizerRDI()]
#' * `rdd` from [BevertonHoltRDD()]
#' 
#' You can modify these in two ways.
#' 
#' @param params A `MizerParams` object
#' @param rate Name of the rate for which a new function is to be set.
#' @param fun Name of the function to use to calculate the rate.
#' @md
#' @export
setRateFunction <- function(params, rate = "Rates", fun) {
    assert_that(is(params, "MizerParams"),
                is.string(rate),
                is.string(fun),
                is.function(get(fun)))
    if (!(rate %in% names(params@rates_funcs))) {
        stop("The `rate` argument must be one of ", 
             toString(names(params@rates_funcs)), ".")
    }
    f <- get0(fun, mode = "function")
    if (is.null(f)) {
        stop(fun, " should be a function")
    }
    # TODO: put some code to test that the function has the right kind of
    # arguments
    params@rates_funcs[[rate]] <- fun
    
    validObject(params)
    params
}

#' @rdname setRateFunction
#' @export
getRateFunction <- function(params, rate = "Rates") {
    assert_that(is(params, "MizerParams"),
                is.string(rate))
    validObject(params)
    if (rate == "All") {
        return(params@rates_funcs)
    }
    if (!(rate %in% names(params@rates_funcs))) {
        stop("The `rate` argument must be one of ", 
             toString(names(params@rates_funcs)), ".")
    }
    params@rates_funcs[[rate]]
}

#' Add a dynamical ecosystem component
#' 
#' @param params A MizerParams object
#' @param component Name of the component
#' @param initial_value Initial value of the component
#' @param dynamics_fun Name of function to calculate value at the next time step
#' @param encounter_fun Name of function to calculate contribution to encounter
#'   rate
#' @param mortality_fun Name of function to calculate contribution to the
#'   predation mortality rate.
#' @param component_params Named list of parameters needed by the component
#'   functions.
#' @return For `setComponent`: The updated MizerParams object
#' @export
setComponent <- function(params, component, initial_value,
                         encounter_fun, mortality_fun,  
                         dynamics_fun, component_params) {
    assert_that(is(params, "MizerParams"),
                is.string(component),
                is.string(dynamics_fun),
                is.string(encounter_fun),
                is.string(mortality_fun),
                is.function(get0(dynamics_fun)),
                is.function(get0(encounter_fun)),
                is.function(get0(mortality_fun)),
                is.list(component_params))
    params@other_dynamics[[component]] <- dynamics_fun
    params@other_pred_mort[[component]] <- pred_mort_fun
    params@other_encounter[[component]] <- encounter_fun
    params@other_params[[component]] <- other_params
    initial_n_other(params)[[component]] <- initial_value
}

#' Get information about other ecosystem components
#' 
#' @param params A MizerParams object
#' @param component Name of the component of interest. If missing, a list of
#'   all components will be returned.
#' @return For `getComponent`: A list with the entries `initial_value`, `dynamics_fun`,
#'   `encounter_fun`, `morality_fun`, `component_params`. If `component` is
#'   missing, then a list of lists for all components is returned.
#' @md
#' @rdname setComponent
#' @export
getComponent <- function(params, component) {
    if (missing(component)) {
        lapply(names(params@other_dynamics),
               function(x) getComponent(params, x))
    }
    comp_list <- list(
        initial_value = params@initial_n_other[[component]],
        component_params = params@other_params[[component]],
        dynamics_fun = params@other_dynamics[[component]],
        mortality_fun = params@other_mort[[component]],
        encounter_fun = params@other_encounter[[component]]
    )
}


#' Initial values for other ecosystem components
#' 
#' Values used as starting values for simulations with `project()`.
#' 
#' @param params A MizerParams object
#' @param value A named list with the initial values of other ecosystem components
#' @export
`initial_n_other<-` <- function(params, value) {
    assert_that(is(params, "MizerParams"),
                is.list(value))
    params@initial_n_other <- value
    params
}

#' @rdname initial_n_other-set
#' @export
initial_n_other <- function(params) {
    params@initial_n_other
}


#' Time series of other components
#' 
#' Fetch the simulation results for other components over time.
#' 
#' @param sim A MizerSim object
#' @return A list array (time x component) that stores the projected
#'   values for other ecosystem components.
#' @export
n_other <- function(sim, component) {
    if (missing(component)) {
        return(sim@n_other)
    }
    
}


#' Values of other ecosystem components at end of simulation
#' 
#' @param sim A MizerSim object
#' @return A named list holding the values of other ecosystem components at the
#'   end of the simulation
#' @export
final_n_other <- function(sim) {
    assert_that(is(sim, "MizerSim"))
    sim@n_other[dim(sim@n)[[1]], ]
}