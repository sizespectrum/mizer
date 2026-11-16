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

    f <- resolve_model_callback(fun, "Rate callback")

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
                  rates_fns = projectRateFunctions(params))
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
    if (rate == "RDD" && grepl("::", fun, fixed = TRUE)) {
        resolve_model_callback(fun, "Reproduction rate callback")
    } else if (!exists(fun, mode = "function")) {
        stop("There is no function called '", fun, "'.")
    }
    .checkRateFunctionOutput(params, rate, fun)
    params@rates_funcs[[rate]] <- fun
    params <- refresh_builtin_tuning(params)

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
#'   they must be unique and must not clash with the name of a component.
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
    if (anyDuplicated(names(value)) > 0) {
        stop("All entries of `", slot_name, "` need unique names.")
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
    # Keep component entries in their existing positions and fill the
    # free-standing positions from `value` in its order. This preserves the
    # complete registry when the getter's result is assigned back unchanged.
    funs <- slot(params, slot_name)
    is_component <- names(funs) %in% names(params@other_dynamics)
    updated <- list()
    value_i <- 1L
    for (i in seq_along(funs)) {
        if (is_component[[i]]) {
            updated <- c(updated, funs[i])
        } else if (value_i <= length(value)) {
            updated <- c(updated, value[value_i])
            value_i <- value_i + 1L
        }
    }
    if (value_i <= length(value)) {
        updated <- c(updated, value[value_i:length(value)])
    }
    slot(params, slot_name) <- updated
    if (slot_name == "other_mort" && !isS4(params)) {
        params <- refresh_builtin_mortality(params)
    }

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
#' A new component needs an `initial_value`, supplied separately from its
#' definition. Supply at least one of `dynamics_fun` and `rate_fun`, directly
#' or through `definition`. For an existing component, omitted arguments retain
#' their existing values. Explicit `NULL` removes an optional callback, provided
#' at least one of the two dynamics capabilities remains. Remove a component
#' with `removeComponent()`.
#'
#' A `definition` is a reusable named list containing `dynamics_fun`,
#' `rate_fun`, `encounter_fun`, `mort_fun`, `colour`, `linetype` and/or
#' `integration_order`, `convergence_fun`, `role`, `stationary`,
#' `equation_rates_fun` and/or `tuning`.
#' Explicit arguments
#' override its entries. Supplying a definition replaces the component's
#' previous definition: omitted optional callbacks are removed and omitted
#' styles default to `"grey"` and `"solid"`. Without `definition`, explicitly
#' supplied arguments patch the existing definition instead. Initial state and
#' `component_params` belong to the component instance and are retained when
#' omitted, including when its definition is replaced.
#'
#' The assembled definition is stored in the model. Subsequent changes to the
#' list used as `definition` do not change that model; call `setComponent()`
#' again to apply an update. Callback function names are stored, rather than
#' their implementations, so the functions still need to be available when
#' the model is used.
#'
#' @section Equations and custom steps:
#' A `rate_fun` returns the component's derivative, with the numeric shape and
#' names of its state preserved. It receives `params`, `n`, `n_pp`, `n_other`,
#' `rates`, `t`, `component` and `...`. Mizer applies Euler under first-order
#' projection and Heun under the second-order methods. This is explicit
#' integration, including when the fish use TR-BDF2.
#'
#' A `dynamics_fun` returns the new state and owns the integration. It is the
#' appropriate capability for stiff or non-numeric states. When both callbacks
#' are supplied, the custom step performs projection and the derivative remains
#' available for diagnostics. Callback names may be unqualified or use
#' `package::exported_function`.
#'
#' Custom steps may accept `predicted_state` and `endpoint_rates`. During a
#' second-order corrector these contain the predicted endpoint state (a list
#' with `n`, `n_pp`, `n_other` and `t`) and its rate list. They are `NULL` when
#' no prediction is available. The ordinary state arguments still describe the
#' start of the step. Prediction data do not by themselves establish the
#' accuracy or stability of the custom integration.
#'
#' `integration_order` declares the order of a custom step; omission means
#' unknown, not an assumed second order. Replacing `dynamics_fun` clears an
#' inherited order unless an order is explicitly supplied with the replacement.
#' Mizer chooses the integrator for equations-only components, so they must not
#' declare a custom integration order.
#'
#' @section Steady-state assessment:
#' A `convergence_fun` returns a data frame with `measure`, `value`, `tolerance`
#' and `units` columns. It receives the model, state, complete rate list, time,
#' derivative when available and component name. Accept `...` for unused
#' arguments. Mizer compares absolute measured values to their tolerances; the
#' callback does not return its own pass/fail verdict.
#'
#' Components are dynamic by default, including legacy components using
#' `constant_other`. Declare `role = "prescribed"` explicitly when the entire
#' component is a prescribed input. Temporarily holding a dynamic component
#' fixed during analysis requires conditional assessment scope. Replacing a
#' dynamics callback removes its inherited convergence callback and resets an
#' inherited prescribed role to dynamic unless those declarations are explicitly
#' supplied again with the replacement.
#'
#' @param params A MizerParams object
#' @param component Name of the component
#' @param initial_value Initial value of the component. Required for a new
#'   component; retained when omitted for an existing component. Must not be
#'   `NULL`; use `removeComponent()` to remove a component.
#' @param dynamics_fun Name of function to calculate value at the next time
#'   step. Required unless `rate_fun` is supplied.
#' @param rate_fun Name of function returning the component's derivative.
#'   Required unless `dynamics_fun` is supplied.
#' @param integration_order Optional positive integer declaring the order of a
#'   custom `dynamics_fun`. `NULL` removes the declaration.
#' @param convergence_fun Name of a convergence-measurement callback, or `NULL`
#'   to remove it. A dynamic component without this support makes a whole-model
#'   assessment incomplete.
#' @param role Either `"dynamic"` (the default) or `"prescribed"`.
#' @param prescribed_mask_fun Named structural input-mask callback, or NULL.
#'   For a dynamic numeric component, return logical values with exactly its
#'   state dimensions and names; TRUE marks prescribed inputs. The callback
#'   receives params and optional component/owner, without time or state.
#' @param autonomous TRUE for equations and contributions without explicit time
#'   dependence, FALSE for forced providers, or NULL for undeclared support.
#'   This declaration is independent of stationary tuning and convergence.
#' @param stationary Explicit TRUE/FALSE stationary-tuning declaration.
#' @param equation_rates_fun Callback exposing effective equation quantities.
#' @param tuning Named strategies and their default, or NULL to remove support.
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
#' @param biomass_fun Named biomass measurement callback, or NULL to remove
#'   measurement support. This is independent of convergence and dynamics.
#' @param biomass_filters Character vector of explicitly supported filters:
#'   min_w, max_w, min_l, max_l and use_cutoff. Defaults to none.
#' @param transformations Named transformation declarations. Update hooks may
#'   modify this component's initial state and parameters; additional paths
#'   must be explicit. Use action = 'none' for supported no-action operations.
#' @param definition A named list specifying the component's callbacks and
#'   plotting styles. Supplying a list replaces its previous definition before
#'   explicitly supplied arguments are applied. State and component parameters
#'   are supplied separately.
#'   Updating only state or component parameters does not re-adopt capabilities
#'   after a global provider change; supply the relevant capability explicitly.
#' @return The updated MizerParams object
#' @seealso "Extending mizer":
#'   [guide to extending mizer](https://sizespectrum.org/mizer/articles/guide-extend-mizer.html)
#' @export
#' @family extension tools
setComponent <- function(params, component, initial_value,
                         dynamics_fun,
                         encounter_fun, mort_fun,
                         component_params,
                         colour = "grey", linetype = "solid",
                         definition, rate_fun, integration_order,
                         convergence_fun, role, stationary,
                         equation_rates_fun, tuning, transformations,
                         biomass_fun, biomass_filters, prescribed_mask_fun,
                         autonomous) {
    params <- validParams(params)
    assert_that(is.string(component), nzchar(component))
    if (!missing(initial_value) && is.null(initial_value)) {
        stop("`initial_value` must not be NULL. Use `removeComponent()` ",
             "to remove a component.")
    }
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
    existing <- existing_component_definition(params, component)
    if (is.null(existing) && missing(initial_value)) {
        stop("`initial_value` is required for a new component.")
    }
    overrides <- list()
    if (!missing(dynamics_fun)) overrides["dynamics_fun"] <- list(dynamics_fun)
    if (!missing(rate_fun)) overrides["rate_fun"] <- list(rate_fun)
    if (!missing(encounter_fun)) {
        overrides["encounter_fun"] <- list(encounter_fun)
    }
    if (!missing(mort_fun)) overrides["mort_fun"] <- list(mort_fun)
    if (!missing(colour)) overrides["colour"] <- list(colour)
    if (!missing(linetype)) overrides["linetype"] <- list(linetype)
    if (!missing(integration_order)) {
        overrides["integration_order"] <- list(integration_order)
    }
    if (!missing(convergence_fun)) {
        overrides["convergence_fun"] <- list(convergence_fun)
    }
    if (!missing(role)) overrides["role"] <- list(role)
    if (!missing(stationary)) overrides["stationary"] <- list(stationary)
    if (!missing(equation_rates_fun)) overrides["equation_rates_fun"] <-
        list(equation_rates_fun)
    if (!missing(tuning)) overrides["tuning"] <- list(tuning)
    if (!missing(transformations)) overrides["transformations"] <-
        list(transformations)
    if (!missing(biomass_fun)) overrides["biomass_fun"] <- list(biomass_fun)
    if (!missing(biomass_filters)) overrides["biomass_filters"] <-
        list(biomass_filters)
    if (!missing(prescribed_mask_fun)) overrides["prescribed_mask_fun"] <-
        list(prescribed_mask_fun)
    if (!missing(autonomous)) overrides["autonomous"] <- list(autonomous)
    assembled <- assemble_component_definition(
        existing = existing,
        definition = if (missing(definition)) NULL else definition,
        overrides = overrides,
        replace = !missing(definition)
    )

    params@other_dynamics[[component]] <- if (is.null(assembled$dynamics_fun)) {
        "component_rate_dynamics"
    } else {
        assembled$dynamics_fun
    }
    params@other_mort[[component]] <- assembled$mort_fun
    params@other_encounter[[component]] <- assembled$encounter_fun
    if (!missing(initial_value)) {
        params@initial_n_other[component] <- list(initial_value)
    }
    if (!missing(component_params)) {
        params@other_params[[component]] <- component_params
    }
    if (!is.null(assembled$colour)) {
        params <- setColours(params,
                             stats::setNames(list(assembled$colour), component))
    }
    if (!is.null(assembled$linetype)) {
        params <- setLinetypes(
            params, stats::setNames(list(assembled$linetype), component)
        )
    }
    # The plot setters ignore invalid styles with a warning. Record their
    # accepted values, including absent styles, rather than rejected inputs.
    accepted <- legacy_component_definition(params, component)
    assembled[c("colour", "linetype")] <- accepted[c("colour", "linetype")]
    params@capabilities$components[[component]] <- assembled
    if (is.null(existing) || !missing(definition) ||
            any(c("stationary", "equation_rates_fun", "tuning") %in%
                    names(overrides))) {
        params$capabilities$components[[component]]$tuning_binding <-
            tuning_provider_binding(params)
    }
    if (is.null(existing) || !missing(definition) ||
            !missing(transformations)) {
        params$capabilities$components[[component]]$transformation_binding <-
            transformation_provider_binding(params)
    }
    if (is.null(existing) || !missing(definition) || !missing(role) ||
            !missing(prescribed_mask_fun) ||
            !identical(existing$dynamics_fun, assembled$dynamics_fun) ||
            !identical(existing$rate_fun, assembled$rate_fun)) {
        params$capabilities$components[[component]]$prescribed_binding <-
            state_provider_binding(params, paste0("component:", component))
    }
    if (is.null(existing) || !missing(definition) || !missing(autonomous)) {
        params$capabilities$components[[component]]$temporal_binding <-
            temporal_provider_binding(params, paste0("component:", component))
    }
    providers <- c("dynamics_fun", "rate_fun", "encounter_fun", "mort_fun")
    provider_changed <- is.null(existing) || any(vapply(
        providers,
        function(field) !identical(existing[[field]], assembled[[field]]),
        logical(1)
    ))
    if (provider_changed) params <- refresh_builtin_tuning(params)

    params@time_modified <- lubridate::now()
    params
}

#' @rdname setComponent
#' @export
removeComponent <- function(params, component) {
    params <- validParams(params)
    if (!component %in% names(params@other_dynamics)) {
        stop("There is no component named ", component)
    }
    params@other_dynamics[[component]] <- NULL
    params@other_encounter[[component]] <- NULL
    params@other_mort[[component]] <- NULL
    params@other_params[[component]] <- NULL
    params@initial_n_other[[component]] <- NULL
    params@capabilities$components[[component]] <- NULL
    params <- refresh_builtin_tuning(params)

    params@time_modified <- lubridate::now()
    params
}


#' @param params A MizerParams object
#' @param component Name of the component of interest. If missing, a list of
#'   all components will be returned.
#' @return For `getComponent`: A list with the entries `initial_value`,
#'   `dynamics_fun`, `rate_fun`, `encounter_fun`, `mort_fun`,
#'   `component_params`, `convergence_fun`, `role`, and `definition` for the
#'   requested component.
#'   If the requested component does not exist, `NULL` is returned.
#'   If no `component` argument is given, then a list of lists for all
#'   components is returned.
#' @export
#' @rdname setComponent
getComponent <- function(params, component) {
    params <- validParams(params)
    if (missing(component)) {
        l <- lapply(names(params@other_dynamics),
                    function(x) getComponent(params, x))
        names(l) <- names(params@other_dynamics)
        return(l)
    }
    if (!component %in% names(params@other_dynamics)) {
        return(NULL)
    }
    definition <- existing_component_definition(params, component)
    dynamics_fun <- params@other_dynamics[[component]]
    if (identical(dynamics_fun, "component_rate_dynamics")) {
        dynamics_fun <- NULL
    }
    list(
        initial_value = initialNOther(params)[[component]],
        dynamics_fun = dynamics_fun,
        rate_fun = definition$rate_fun,
        mort_fun = params@other_mort[[component]],
        encounter_fun = params@other_encounter[[component]],
        component_params = params@other_params[[component]],
        convergence_fun = definition$convergence_fun,
        role = definition$role,
        definition = definition
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
