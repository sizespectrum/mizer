#' Set own rate function to replace mizer rate function
#' 
#' If the way mizer calculates a fundamental rate entering the model is
#' not flexible enough for you (for example if you need to introduce time
#' dependence) then you can write your own functions for calculating that
#' rate and use `setRateFunction()` to register it with mizer.
#' 
#' At each time step during a simulation with the [project()] function, mizer
#' needs to calculate the instantaneous values of the various rates. By
#' default it calls the [mizerRates()] function which creates a list with the
#' following components:
#' * `encounter` from [mizerEncounter()]
#' * `feeding_level` from [mizerFeedingLevel()]
#' * `pred_rate` from [mizerPredRate()]
#' * `pred_mort` from [mizerPredMort()]
#' * `f_mort` from [mizerFMort()]
#' * `mort` from [mizerMort()]
#' * `resource_mort` from [mizerResourceMort()]
#' * `e` from [mizerEReproAndGrowth()]
#' * `e_repro` from [mizerERepro()]
#' * `e_growth` from [mizerEGrowth()]
#' * `rdi` from [mizerRDI()]
#' * `rdd` from [BevertonHoltRDD()]
#' 
#' For each of these you can substitute your own function. So for example if
#' you have written your own function for calculating the total mortality
#' rate and have called it `myMort` and have a mizer model stored in a
#' MizerParams object called `params` that you want to run with your new
#' mortality rate, then you would call
#' ```
#' params <- setRateFunction(params, "Mort", "myMort")
#' ```
#' In general if you want to replace a function `mizerSomeRateFunc()` with
#' a function `myVersionOfThis()` you would call
#' ```
#' params <- setRateFunction(params, "SomeRateFunc", "myVersionOfThis")
#' ```
#' 
#' In some extreme cases you may need to swap out the entire `mizerRates()`
#' function for your own function called `myRates()`. That you can do with
#' ```
#' params <- setRateFunction(params, "Rates", "myRates")
#' ```
#' 
#' Your new rate functions may need their own model parameters. These you 
#' can store in `other_params(params)`. For example
#' ```
#' other_params(params)$my_param <- 42
#' ```
#' Note that your own rate functions need to be defined in the global
#' environment or in a package. If they are defined within a function then
#' mizer will not find them.
#' 
#' @param params A MizerParams object
#' @param rate Name of the rate for which a new function is to be set.
#' @param fun Name of the function to use to calculate the rate.
#' @return For `setRateFunction()`: An updated MizerParams object
#' @export
setRateFunction <- function(params, rate, fun) {
    assert_that(is(params, "MizerParams"),
                is.string(rate),
                is.string(fun))
    if (!(rate %in% names(params@rates_funcs))) {
        stop("The `rate` argument must be one of ", 
             toString(names(params@rates_funcs)), ".")
    }
    if (!exists(fun, mode = "function")) {
        stop("There is no function called '", fun, "'.")
    }
    # TODO: put some code to test that the function has the right kind of
    # arguments
    params@rates_funcs[[rate]] <- fun
    
    validObject(params)
    
    params@time_modified <- lubridate::now()
    params
}

#' @rdname setRateFunction
#' @return For `getRateFunction()`: The name of the registered rate function for
#'   the requested `rate`, or the list of all rate functions if called without
#'   `rate` argument.
#' @export
getRateFunction <- function(params, rate) {
    params <- validParams(params)
    if (missing(rate)) {
        return(params@rates_funcs)
    }
    if (!(rate %in% names(params@rates_funcs))) {
        stop("The `rate` argument must be one of ", 
             toString(names(params@rates_funcs)), ".")
    }
    params@rates_funcs[[rate]]
}

#' @rdname setRateFunction
#' @return For `other_params()`: A named list with all the parameters for which
#'   you have set values.
#' @export
other_params <- function(params) {
    assert_that(is(params, "MizerParams"))
    params@other_params$other
}

#' @rdname setRateFunction
#' @param value Values for other parameters
#' @export
`other_params<-` <- function(params, value) {
    assert_that(is(params, "MizerParams"))
    if (!is.list(value) || is.null(names(value))) {
        stop("other_params should be a named list.")
    }
    # We save the value in the $other slot in order to make it impossible for
    # the user to overwrite component parameters by mistake.
    params@other_params$other <- value
    
    params@time_modified <- lubridate::now()
    params
}

#' Add a dynamical ecosystem component
#' 
#' By default, mizer models any number of size-resolved consumer species
#' and a single size-resolved resource spectrum. Your model may require
#' additional components, like for example detritus or carrion or multiple
#' resources or .... This function allows you to set up such components.
#' 
#' The component can be a number, a vector, an array, a list, or any other
#' data structure you like. 
#' 
#' If you set a component with a new name, the new component will be added
#' to the existing components. If you set a component with an existing name,
#' that component will be overwritten. You can remove a component with
#' `removeComponent()`.
#' 
#' @param params A MizerParams object
#' @param component Name of the component
#' @param initial_value Initial value of the component
#' @param dynamics_fun Name of function to calculate value at the next time step
#' @param encounter_fun Name of function to calculate contribution to encounter
#'   rate. Optional.
#' @param mort_fun Name of function to calculate contribution to the
#'   mortality rate. Optional.
#' @param component_params Object holding the parameters needed by the component
#'   functions. This could for example be a named list of parameters. Optional.
#' @return The updated MizerParams object
#' @export
setComponent <- function(params, component, initial_value,
                         dynamics_fun, 
                         encounter_fun, mort_fun,  
                         component_params) {
    assert_that(is(params, "MizerParams"),
                is.string(component),
                is.string(dynamics_fun),
                is.function(get0(dynamics_fun)))
    params@other_dynamics[[component]] <- dynamics_fun
    params@initial_n_other[[component]] <- initial_value
    # TODO: Add checks that the functions have the right arguments and
    # return values
    if (!missing(mort_fun)) {
        if (!is.null(mort_fun) && !is.function(get0(mort_fun))) {
            stop("`mort_fun` needs to be NULL or a function.")
        }
        params@other_mort[[component]] <- mort_fun
    }
    if (!missing(encounter_fun)) {
        if (!is.null(encounter_fun) && !is.function(get0(encounter_fun))) {
            stop("`encounter_fun` needs to be NULL or a function.")
        }
        params@other_encounter[[component]] <- encounter_fun
    }
    if (!missing(component_params)) {
        params@other_params[[component]] <- component_params
    }
    
    params@time_modified <- lubridate::now()
    params
}

#' @rdname setComponent
#' @export
removeComponent <- function(params, component) {
    if (!component %in% names(params@other_dynamics)) {
        stop("There is no component named ", component)
    }
    params@other_dynamics[[component]] <- NULL
    params@other_encounter[[component]] <- NULL
    params@other_mort[[component]] <- NULL
    params@other_params[[component]] <- NULL
    params@initial_n_other[[component]] <- NULL
    
    params@time_modified <- lubridate::now()
    params
}


#' Get information about other ecosystem components
#' 
#' @param params A MizerParams object
#' @param component Name of the component of interest. If missing, a list of
#'   all components will be returned.
#' @return A list with the entries `initial_value`, `dynamics_fun`,
#'   `encounter_fun`, `mort_fun`, `component_params` for the requested
#'   component. If the requested component does not exist, `NULL` is returned.
#'   If no `component` argument is given, then a list of lists for all
#'   components is returned.
#' @export
getComponent <- function(params, component) {
    assert_that(is(params, "MizerParams"))
    if (missing(component)) {
        l <- lapply(names(params@other_dynamics),
                    function(x) getComponent(params, x))
        names(l) <- names(params@other_dynamics)
        return(l)
    }
    if (!component %in% names(params@other_dynamics)) {
        return(NULL)
    }
    list(initial_value = initialNOther(params)[[component]],
         dynamics_fun = params@other_dynamics[[component]],
         mort_fun = params@other_mort[[component]],
         encounter_fun = params@other_encounter[[component]],
         component_params = params@other_params[[component]]
    )
}


#' Initial values for other ecosystem components
#'
#' Values used as starting values for simulations with `project()`.
#'
#' @param params A MizerParams object
#' @param value A named list with the initial values of other ecosystem
#'   components
#' @export
#' @seealso [initialNResource()], [initialN()]
`initialNOther<-` <- function(params, value) {
    assert_that(is(params, "MizerParams"),
                is.list(value))
    components <- names(params@other_dynamics)
    missing <- !(names(value) %in% components)
    if (any(missing)) {
        stop("The following components do not exist: ", names(value)[missing])
    }
    extra <- !(components %in% names(value))
    if (any(extra)) {
        stop("Missing values for components ", components[extra])
    }
    params@initial_n_other <- value
    
    params@time_modified <- lubridate::now()
    params
}

#' @param object An object of class MizerParams or MizerSim
#' @rdname initialNOther-set
#' @return A named list with the initial values of other ecosystem
#'   components
#' @export
initialNOther <- function(object) {
    if (is(object, "MizerParams")) {
        params <- validParams(object)
        return(params@initial_n_other)
    }
    if (is(object, "MizerSim")) {
        return(object@params@initial_n_other)
    } else {
        stop("The argument needs to be a MizerSim or a MizerParams object.")
    }
}

#' Time series of other components
#' 
#' Fetch the simulation results for other components over time.
#' 
#' @param sim A MizerSim object
#' @return A list array (time x component) that stores the projected values for
#'   other ecosystem components.
#' @export
NOther <- function(sim) {
    return(sim@n_other)
}


#' Values of other ecosystem components at end of simulation
#' 
#' @param sim A MizerSim object
#' @return A named list holding the values of other ecosystem components at the
#'   end of the simulation
#' @export
finalNOther <- function(sim) {
    assert_that(is(sim, "MizerSim"))
    n_other <- sim@n_other[dim(sim@n)[[1]], ]
    names(n_other) <- dimnames(sim@n_other)$component
    n_other
}
