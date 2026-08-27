#' Check that a rate function returns the correct output dimensions
#'
#' Called by [setRateFunction()] to verify that a candidate rate function
#' returns an array (or vector/list) of the correct dimensions for the
#' requested rate.
#'
#' @param params A MizerParams object
#' @param rate Name of the rate being replaced, e.g. `"Encounter"`.
#' @param fun Name of the candidate function to validate.
#' @return Invisibly `NULL`. Called for its side-effect of stopping with an
#'   informative error if the output has the wrong shape.
#' @keywords internal
.checkRateFunctionOutput <- function(params, rate, fun) {
    no_sp  <- nrow(params@species_params)
    no_w   <- length(params@w)
    no_w_full <- length(params@w_full)

    n       <- params@initial_n
    n_pp    <- params@initial_n_pp
    n_other <- params@initial_n_other
    t       <- 0
    effort  <- params@initial_effort

    f <- get(fun, mode = "function")

    # For functions that depend on prerequisite rates, compute the current rates
    # first (using whatever functions are already registered).
    if (!(rate %in% c("Encounter", "Rates"))) {
        rates <- tryCatch(
            getRates(params, n = n, n_pp = n_pp, n_other = n_other,
                     effort = effort, t = t),
            error = function(e) {
                warning("Could not validate '", fun, "' because the current ",
                        "model rates could not be computed: ",
                        conditionMessage(e), call. = FALSE)
                NULL
            }
        )
        if (is.null(rates)) return(invisible(NULL))
    }

    # Call the candidate function with appropriate test inputs.
    result <- tryCatch(
        switch(rate,
            Encounter =
                f(params, n = n, n_pp = n_pp, n_other = n_other, t = t),
            FeedingLevel =
                f(params, n = n, n_pp = n_pp, n_other = n_other, t = t,
                  encounter = rates$encounter),
            EReproAndGrowth =
                f(params, n = n, n_pp = n_pp, n_other = n_other, t = t,
                  encounter = rates$encounter,
                  feeding_level = rates$feeding_level),
            ERepro =
                f(params, n = n, n_pp = n_pp, n_other = n_other, t = t,
                  e = rates$e),
            EGrowth =
                f(params, n = n, n_pp = n_pp, n_other = n_other, t = t,
                  e = rates$e, e_repro = rates$e_repro),
            Diffusion =
                f(params, n = n, n_pp = n_pp, n_other = n_other, t = t,
                  feeding_level = rates$feeding_level),
            PredRate =
                f(params, n = n, n_pp = n_pp, n_other = n_other, t = t,
                  feeding_level = rates$feeding_level),
            PredMort =
                f(params, n = n, n_pp = n_pp, n_other = n_other, t = t,
                  pred_rate = rates$pred_rate),
            FMort =
                f(params, n = n, n_pp = n_pp, n_other = n_other, t = t,
                  effort = effort,
                  e_growth = rates$e_growth, pred_mort = rates$pred_mort),
            Mort =
                f(params, n = n, n_pp = n_pp, n_other = n_other, t = t,
                  f_mort = rates$f_mort, pred_mort = rates$pred_mort),
            RDI =
                f(params, n = n, n_pp = n_pp, n_other = n_other, t = t,
                  e_growth = rates$e_growth, mort = rates$mort,
                  e_repro = rates$e_repro, diffusion = rates$diffusion),
            RDD =
                f(rdi = rates$rdi, species_params = params@species_params,
                  params = params, t = t),
            ResourceMort =
                f(params, n = n, n_pp = n_pp, n_other = n_other, t = t,
                  pred_rate = rates$pred_rate),
            Rates =
                f(params, n = n, n_pp = n_pp, n_other = n_other, t = t,
                  effort = effort,
                  rates_fns = lapply(params@rates_funcs, get))
        ),
        error = function(e) {
            stop("The function '", fun, "' failed when called with test inputs: ",
                 conditionMessage(e), call. = FALSE)
        }
    )

    # --- Check return dimensions ---
    if (rate == "Rates") {
        if (!is.list(result)) {
            stop("The function '", fun, "' must return a list, not a ",
                 class(result)[[1L]], ".", call. = FALSE)
        }
        required <- c("encounter", "feeding_level", "pred_rate", "pred_mort",
                      "f_mort", "mort", "resource_mort", "e", "e_repro",
                      "e_growth", "diffusion", "rdi", "rdd")
        missing  <- setdiff(required, names(result))
        if (length(missing) > 0L) {
            stop("The list returned by '", fun,
                 "' is missing the following components: ",
                 toString(missing), ".", call. = FALSE)
        }
    } else if (rate %in% c("RDI", "RDD")) {
        if (length(result) != no_sp) {
            stop("The function '", fun, "' must return a vector of length ",
                 no_sp, " (one value per species), but returned length ",
                 length(result), ".", call. = FALSE)
        }
    } else if (rate == "ResourceMort") {
        if (length(result) != no_w_full) {
            stop("The function '", fun, "' must return a vector of length ",
                 no_w_full,
                 " (one value per size bin in the full size grid), ",
                 "but returned length ", length(result), ".", call. = FALSE)
        }
    } else {
        expected_dim <- if (rate == "PredRate") c(no_sp, no_w_full)
                        else                    c(no_sp, no_w)
        actual_dim   <- dim(result)
        size_desc    <- if (rate == "PredRate")
            paste0(no_sp, " x ", no_w_full, " (species x full size grid)")
        else
            paste0(no_sp, " x ", no_w, " (species x size grid)")
        if (is.null(actual_dim) || length(actual_dim) < 2L ||
            actual_dim[[1L]] != expected_dim[[1L]] ||
            actual_dim[[2L]] != expected_dim[[2L]]) {
            if (!is.null(actual_dim)) {
                stop("The function '", fun,
                     "' must return a 2D array of dimensions ", size_desc,
                     " but returned dimensions ",
                     paste(actual_dim, collapse = " x "), ".", call. = FALSE)
            } else {
                stop("The function '", fun,
                     "' must return a 2D array of dimensions ", size_desc,
                     " but returned a non-array object of class ",
                     class(result)[[1L]], ".", call. = FALSE)
            }
        }
    }

    invisible(NULL)
}

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
#' * `diffusion` from [mizerDiffusion()]
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
#' @section Avoid rates that jump as a function of abundance:
#' Make sure your rate function depends **continuously** on the abundances `n`,
#' `n_pp` and `n_other`. It is tempting to write a rate that switches abruptly
#' on the state of the model — a fishery that closes when a stock falls below a
#' limit reference point, a predator that switches diet when its preferred prey
#' becomes scarce, a mortality that kicks in below a critical condition. Such a
#' rate breaks the assumption underlying every one of mizer's time-stepping
#' methods, which freeze the rates during each density update and so cannot see
#' a threshold being crossed within a step.
#'
#' The symptoms are quiet: mizer issues no warning, but the trajectory keeps
#' changing as you refine `dt`, the Newton solver stalls, and [getStability()]
#' reports a confident but meaningless answer. Choosing the L-stable
#' `method = "tr_bdf2"` in [project()] does not help, because the difficulty
#' lies in the frozen rates rather than in the linear solve.
#'
#' The remedy is to give the switch a finite width, using a linear ramp between
#' two thresholds or a logistic transition, which is usually the more realistic
#' model anyway. Rates built with `max()` or `min()` are continuous but not
#' differentiable; these are much less troublesome, costing some accuracy but
#' not correctness. See the [Discontinuous rate
#' functions](https://sizespectrum.org/mizer/articles/discontinuous_rates.html)
#' article for the full story, the diagnostics and the fix.
#'
#' @param params A MizerParams object
#' @param rate Name of the rate for which a new function is to be set.
#' @param fun Name of the function to use to calculate the rate.
#' @return For `setRateFunction()`: An updated MizerParams object
#' @seealso "Extending mizer":
#'   [guide to extending mizer](https://sizespectrum.org/mizer/articles/guide-extend-mizer.html);
#'   [Discontinuous rate
#'   functions](https://sizespectrum.org/mizer/articles/discontinuous_rates.html)
#' @export
#' @family extension tools
setRateFunction <- function(params, rate, fun) {
    params <- validParams(params)
    assert_that(is.string(rate),
                is.string(fun))
    if (!(rate %in% names(params@rates_funcs))) {
        stop("The `rate` argument must be one of ",
             toString(names(params@rates_funcs)), ".")
    }
    if (!exists(fun, mode = "function")) {
        stop("There is no function called '", fun, "'.")
    }
    .checkRateFunctionOutput(params, rate, fun)
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
#' @return For `other_params()`: The user-defined parameters stored in
#'   `other_params(params)`, or `NULL` if none have been set. This excludes any
#'   component-specific parameters stored via [setComponent()].
#' @export
other_params <- function(params) {
    assert_that(is(params, "MizerParams"))
    params@other_params$other
}

#' @rdname setRateFunction
#' @param value A named list of user-defined parameters to store in
#'   `other_params(params)`.
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

#' Extra contributions to the mortality and encounter rates
#'
#' Besides the rates that mizer calculates itself, a model can carry extra
#' contributions to the mortality rate and to the encounter rate. Each
#' contribution is an R function, registered here by name, that mizer calls at
#' every time step: [getMort()] adds the result of every function listed in
#' `other_mort(params)` to the mortality rate and [getEncounter()] adds the
#' result of every function listed in `other_encounter(params)` to the encounter
#' rate.
#'
#' Use these when the extra contribution depends on the state of the model, as
#' a starvation mortality or a density-dependent mortality does. A contribution
#' that is simply a fixed array is better set with [ext_mort()] or
#' [ext_encounter()], which is what mizer's own external mortality and external
#' encounter use.
#'
#' Assigning `NULL` to an entry removes that contribution:
#' ```
#' other_mort(params)[["starvation"]] <- NULL
#' ```
#'
#' @section How your function is called:
#' Each registered function is called as
#' ```
#' fun(params, n, n_pp, n_other, t, component, ...)
#' ```
#' and has to return an array with the same dimensions as the rate it
#' contributes to, that is species x size. The `component` argument holds the
#' name under which the function is registered, so a single implementation can
#' serve several entries. Always give your function a `...` argument so that it
#' tolerates being passed arguments it does not use.
#'
#' Make sure the function depends *continuously* on the abundances, for the
#' reasons set out in [setRateFunction()].
#'
#' @section Contributions belonging to a component:
#' A component set up with [setComponent()] can contribute to these rates
#' through that function's `mort_fun` and `encounter_fun` arguments. Such an
#' entry belongs to its component: [setComponent()] sets it, [removeComponent()]
#' removes it and [getComponent()] reports it. The accessors here therefore
#' leave those entries alone. `other_mort()` and `other_encounter()` list only
#' the contributions that do not belong to a component, and assigning through
#' them preserves the ones that do. This mirrors [other_params()], which
#' likewise hides the parameters belonging to a component, and it makes it
#' impossible to wipe out a component's contribution by assigning a whole list.
#'
#' @param params A MizerParams object
#' @param value A named list of function names, or `NULL` to remove all
#'   contributions that do not belong to a component. You choose the names, but
#'   they must not clash with the name of a component.
#' @return A named list with the names of the functions contributing to the
#'   rate, excluding any that belong to a component.
#' @export
#' @family extension tools
#' @seealso [setComponent()], [other_params()], [ext_mort()], [ext_encounter()]
#' @examples
#' # An extra mortality that grows with the total biomass in the community.
#' # Like any custom rate function it has to live in the global environment
#' # or in a package, so that mizer can find it by name.
#' crowdingMort <- function(params, n, n_pp, n_other, t, component, ...) {
#'     biomass <- sum(n %*% (params@w * params@dw))
#'     # Same dimensions as the mortality rate: species x size
#'     0 * params@mu_b + 1e-14 * biomass
#' }
#' params <- NS_params
#' other_mort(params)[["crowding"]] <- "crowdingMort"
#' other_mort(params)
#'
#' # The contribution is included in the mortality rate
#' range(getMort(params) - getMort(NS_params))
#'
#' # and can be removed again
#' other_mort(params)[["crowding"]] <- NULL
#' other_mort(params)
other_mort <- function(params) {
    assert_that(is(params, "MizerParams"))
    freeRateContributions(params, "other_mort")
}

#' @rdname other_mort
#' @export
`other_mort<-` <- function(params, value) {
    setRateContributions(params, "other_mort", value,
                         noun = "mortality", arg = "mort_fun")
}

#' @rdname other_mort
#' @export
other_encounter <- function(params) {
    assert_that(is(params, "MizerParams"))
    freeRateContributions(params, "other_encounter")
}

#' @rdname other_mort
#' @export
`other_encounter<-` <- function(params, value) {
    setRateContributions(params, "other_encounter", value,
                         noun = "encounter", arg = "encounter_fun")
}

#' The entries of `other_mort` or `other_encounter` that belong to a component
#'
#' The `other_mort` and `other_encounter` slots are keyed by name, and an entry
#' whose name is that of a component created with [setComponent()] belongs to
#' that component. Those entries are owned by `setComponent()`,
#' `removeComponent()` and `getComponent()`; the rest are owned by
#' [other_mort()] and [other_encounter()], which is why the accessors have to
#' split the list.
#'
#' @param params A MizerParams object.
#' @param slot_name Either `"other_mort"` or `"other_encounter"`.
#' @return The named list of function names, restricted to the entries that do
#'   (`componentRateContributions()`) or do not (`freeRateContributions()`)
#'   belong to a component.
#' @noRd
componentRateContributions <- function(params, slot_name) {
    funs <- slot(params, slot_name)
    funs[names(funs) %in% names(params@other_dynamics)]
}

#' @noRd
freeRateContributions <- function(params, slot_name) {
    funs <- slot(params, slot_name)
    funs[!names(funs) %in% names(params@other_dynamics)]
}

#' Validate and store the rate contributions that do not belong to a component
#'
#' Shared implementation of `other_mort<-()` and `other_encounter<-()`. The
#' entries belonging to a component are preserved unchanged, so that assigning a
#' whole list cannot wipe out what [setComponent()] put there.
#'
#' @param params A MizerParams object.
#' @param slot_name Either `"other_mort"` or `"other_encounter"`.
#' @param value The named list of function names to store, or `NULL`.
#' @param noun The word for the rate, used in error messages.
#' @param arg The name of the corresponding [setComponent()] argument, used in
#'   error messages.
#' @return The updated MizerParams object.
#' @noRd
setRateContributions <- function(params, slot_name, value, noun, arg) {
    assert_that(is(params, "MizerParams"))
    if (is.null(value)) value <- list()
    if (!is.list(value)) {
        stop("`", slot_name, "` should be a named list of function names.")
    }
    if (length(value) > 0 && !rlang::is_named(value)) {
        stop("All entries of `", slot_name, "` need to be named.")
    }
    clash <- intersect(names(value), names(params@other_dynamics))
    if (length(clash) > 0) {
        stop("`", toString(clash), "` ",
             if (length(clash) > 1) "are components" else "is a component",
             ". Set the ", noun, " function of a component with ",
             "setComponent(", arg, " = ).")
    }
    for (i in seq_along(value)) {
        fun <- value[[i]]
        if (!is.string(fun) || !is.function(get0(fun))) {
            stop("The entry `", names(value)[[i]], "` of `", slot_name,
                 "` needs to be the name of a function.")
        }
    }
    slot(params, slot_name) <-
        c(componentRateContributions(params, slot_name), value)

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
#' the `initial_value` and `dynamics_fun` are overwritten, while the optional
#' `encounter_fun`, `mort_fun` and `component_params` are only changed if the
#' corresponding arguments are supplied. You can remove a component with
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
#' @param colour Line colour to use for the component in plots. Defaults to
#'   `"grey"`.
#' @param linetype Line type to use for the component in plots. Defaults to
#'   `"solid"`.
#' @return The updated MizerParams object
#' @seealso "Extending mizer":
#'   [guide to extending mizer](https://sizespectrum.org/mizer/articles/guide-extend-mizer.html)
#' @export
#' @family extension tools
setComponent <- function(params, component, initial_value,
                         dynamics_fun,
                         encounter_fun, mort_fun,
                         component_params,
                         colour = "grey", linetype = "solid") {
    assert_that(is(params, "MizerParams"),
                is.string(component),
                is.string(dynamics_fun),
                is.function(get0(dynamics_fun)))
    # A rate contribution registered with `other_mort()` or `other_encounter()`
    # under this name would be taken over by the component and thereby become
    # invisible to those accessors, so refuse the name instead.
    if (component %in% c(names(freeRateContributions(params, "other_mort")),
                         names(freeRateContributions(params,
                                                     "other_encounter")))) {
        stop("There is already a rate contribution registered under the name '",
             component, "'. Remove it with `other_mort()` or ",
             "`other_encounter()` first, or choose a different component name.")
    }
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
    params <- setColours(params, stats::setNames(list(colour), component))
    params <- setLinetypes(params, stats::setNames(list(linetype), component))

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


#' @param params A MizerParams object
#' @param component Name of the component of interest. If missing, a list of
#'   all components will be returned.
#' @return For `getComponent`: A list with the entries `initial_value`, `dynamics_fun`,
#'   `encounter_fun`, `mort_fun`, `component_params` for the requested
#'   component. If the requested component does not exist, `NULL` is returned.
#'   If no `component` argument is given, then a list of lists for all
#'   components is returned.
#' @export
#' @rdname setComponent
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
#' @family extension tools
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

#' @rdname initialNOther-set
#' @param object An object of class MizerParams or MizerSim
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
#' @return For `NOther`: A list array indexed by time and component that stores the projected
#'   values for other ecosystem components.
#' @export
#' @family extension tools
NOther <- function(sim) {
    return(sim@n_other)
}


#' Values of other ecosystem components at end of simulation
#'
#' @param sim A MizerSim object
#' @return For `finalNOther`: A named list holding the values of other ecosystem components at the
#'   end of the simulation
#' @export
#' @rdname NOther
finalNOther <- function(sim) {
    assert_that(is(sim, "MizerSim"))
    n_other <- sim@n_other[dim(sim@n)[[1]], ]
    names(n_other) <- dimnames(sim@n_other)$component
    n_other
}
