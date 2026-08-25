# Finding a steady state ------------------------------------------------------
#
# Two things can be varied when putting a model onto a steady state, and the two
# exported functions here are named after them rather than after the algorithm
# they use.
#
# `tuneSteadyState()` holds the inputs to the fish dynamics — the reproduction
# rate, the resource and the other components — at the values it was given while
# it solves for the consumer spectra, and then adjusts the parameters that
# generate those inputs (`erepro`/`R_max` and `cc_pp`) so that the held values
# are themselves steady. That is the calibration job.
#
# `findSteadyState()` changes no parameter at all: reproduction, the resource
# and the spectra settle together wherever the model's own dynamics take them.
#
# Each takes a `solver` argument choosing between running the dynamics
# (`project_until_settled()`, in `steady.R`) and solving the steady-state
# equation directly (`newton_steady_state()`, in `steadyNewton.R`).


#' Tune a model so that the state it is in becomes a steady state
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Solves for the consumer size spectra while holding the reproduction rate
#' (RDD), the resource and any other components at the values stored in
#' `params`, and then adjusts the parameters that generate those held values so
#' that they are steady too. This is the function to use while setting up and
#' calibrating a model: holding the inputs to the fish dynamics fixed is what
#' makes the search reliable.
#'
#' @details
#' Concretely, three things are held fixed during the search and two parameters
#' are re-derived afterwards:
#'
#' * The reproduction rate is pinned at `getRDD(params)`. Afterwards, if the
#'   model uses Beverton-Holt reproduction, [setBevertonHolt()] restores the
#'   density dependence so that it reproduces exactly that rate at the new
#'   spectra. Use `preserve` to say which of the reproduction parameters should
#'   be kept as it was.
#' * The resource abundance is held at `initialNResource(params)`. Afterwards
#'   the resource capacity `cc_pp` is recomputed so that this abundance is a
#'   steady state of the resource under the new spectra. If the capacity had
#'   been set by hand (frozen), it is rebalanced and thereby unfrozen.
#' * Other components are held constant throughout and are not adjusted.
#'
#' So the model you get back is at a fixed point of the *full* dynamics, with
#' reproduction and the resource free, and [getStability()] can be applied to it
#' directly. Contrast [findSteadyState()], which changes no parameter and
#' instead lets the reproduction rate and the resource move to wherever the
#' parameters you already have put them.
#'
#' Holding those inputs fixed is what makes the search reliable, but it does not
#' make the result certain: the state that is stored is only as close to a fixed
#' point as the solver's tolerance allowed, and with `solver = "project"` the
#' run may instead have stopped on a limit cycle or on a species going extinct.
#' Check the result rather than assuming it; see the section below.
#'
#' # Choosing a solver
#'
#' `solver = "project"` (the default) runs the dynamics until they settle, via
#' [projectUntilSettled()], using [distanceMaxRelRDI()] as its distance
#' function. It needs no extra packages and works with any resource dynamics.
#'
#' `solver = "newton"` solves the steady-state equation directly with a
#' Newton-type root finder from the `nleqslv` package. It converges even when
#' the steady state is dynamically **unstable**, where the time-stepping solver
#' cannot, and it discovers the support of the steady state automatically. It
#' starts from the spectra in `initialN(params)`, so a reasonable initial guess
#' still matters — for example the spectra from a nearby stable
#' parameterisation, or the (diverging) output of `solver = "project"`.
#'
#' Because the resource is held fixed either way, `solver = "newton"` here does
#' not need the `steady_<resource_dynamics>()` companion that
#' [findSteadyState()] uses when the resource is one of the unknowns.
#'
#' @template section_check_steady
#'
#' @param params A \linkS4class{MizerParams} object
#' @param solver The solver to use: `"project"` to run the dynamics until they
#'   settle, `"newton"` to solve the steady-state equation directly. See
#'   *Choosing a solver*.
#' @param effort The fishing effort to use throughout. By default the initial
#'   effort stored in `params`.
#' @param preserve `r lifecycle::badge("experimental")`
#'   Specifies whether the `reproduction_level` should be preserved (default)
#'   or the maximum reproduction rate `R_max` or the reproductive
#'   efficiency `erepro`. See [setBevertonHolt()] for an explanation
#'   of the `reproduction_level`.
#' @param info_level Controls the amount of information messages that are shown.
#'   Higher levels lead to more messages, `info_level = 0` gives silence. The
#'   default is taken from the `mizer_info_level` option, see
#'   [default_info_level()].
#' @param ... Arguments for the chosen solver.
#'
#'   With `solver = "project"`: `t_max`, `t_check`, `dt`,
#'   `distance_tol`, `residual_tol`, `amplitude_tol`, `amp_rel_tol`,
#'   `extinction_threshold`, `progress_bar` and `method`, all as described in
#'   [projectUntilSettled()]. There is no `t_save`, because no trajectory is
#'   returned.
#'   Note that `distance_tol` here defaults to `0.1 * dt` and measures the
#'   largest relative change in egg production, because the distance function is
#'   [distanceMaxRelRDI()]. `residual_tol` is judged on the model as the search
#'   sees it, with reproduction, the resource and the other components pinned;
#'   the residual reported in the result is measured again on the model that is
#'   actually returned.
#'
#'   With `solver = "newton"`: `solver_tol` (default `1e-6`), a tolerance on
#'   the scaled steady-state residual passed to [nleqslv::nleqslv()]. It was
#'   called `residual_tol` before mizer 3.3, a name that now belongs to the
#'   biomass drift criterion above; `maxit`
#'   (default `200`); `jacobian`, either `"update"` (default, the Jacobian is
#'   computed once and then updated cheaply each iteration — `nleqslv`'s
#'   `"Broyden"`) or `"recompute"` (a numerical Jacobian at every iteration —
#'   `nleqslv`'s `"Newton"`); `global`, the globalisation strategy (default
#'   `"dbldog"`, a robust double-dogleg trust region); and `verbose` to trace
#'   the iterations.
#'
#' @return A `MizerParams` object with the initial state replaced by the steady
#'   state found, and with `erepro`/`R_max` and `cc_pp` adjusted as described
#'   above. It carries a `"convergence"` attribute describing the solution
#'   found; see [projectUntilSettled()]. Check it: convergence is not
#'   guaranteed.
#' @seealso [findSteadyState()], [projectUntilSettled()],
#'   [steadySingleSpecies()], [isSteady()], [getSteadyResidual()],
#'   [getStability()]
#' @export
#' @examples
#' \donttest{
#' params <- newTraitParams()
#' species_params(params)$gamma[5] <- 3000
#' params <- tuneSteadyState(params)
#' plotSpectra(params)
#' }
tuneSteadyState <- function(params, solver = c("project", "newton"),
                            effort = params@initial_effort,
                            preserve = c("reproduction_level", "erepro",
                                         "R_max"),
                            info_level = default_info_level(), ...) {
    UseMethod("tuneSteadyState")
}

#' @export
tuneSteadyState.MizerParams <- function(params,
                                        solver = c("project", "newton"),
                                        effort = params@initial_effort,
                                        preserve = c("reproduction_level",
                                                     "erepro", "R_max"),
                                        info_level = default_info_level(),
                                        ...) {
    solver <- match.arg(solver)
    with_info_level(info_level = info_level, {
        params <- validParams(params)
        effort <- validEffortVector(effort, params = params)
        params@initial_effort <- effort
        warn_other_components_fixed(params, paste(
            "The model returned is therefore at a fixed point of the",
            "consumer-resource dynamics; check the `residual` entry of its",
            "`\"convergence\"` attribute, which does cover the components,",
            "before treating it as a steady state of the whole model."))
        if (solver == "project") {
            tune_steady_project(params, effort = effort, preserve = preserve,
                                info_level = info_level, ...)
        } else {
            tune_steady_newton(params, effort = effort, preserve = preserve,
                               ...)
        }
    })
}


#' Find the steady state of a model
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Puts the model onto a steady state of its own dynamics, changing no
#' parameter: the reproduction rate, the resource and the consumer spectra all
#' settle together at whatever the parameters you already have imply.
#'
#' @details
#' This is the counterpart of [tuneSteadyState()], which instead holds the
#' reproduction rate and the resource at the values you supplied and adjusts
#' `erepro`/`R_max` and `cc_pp` to make those values steady. Use this function
#' when the parameters are the thing you want to keep — when asking what state a
#' given model settles into, for instance under a changed fishing effort — and
#' [tuneSteadyState()] while calibrating.
#'
#' Nothing being held fixed means the search has more ways to end badly. With
#' `solver = "project"` the run can settle on a limit cycle or drive a species
#' extinct rather than reach a fixed point, and with either solver reproduction
#' can collapse. Check the `"convergence"` attribute rather than assuming a
#' fixed point was reached; see the section below.
#'
#' # Choosing a solver
#'
#' `solver = "project"` (the default) runs the dynamics until they settle, via
#' [projectUntilSettled()], and takes the final state. It is exactly that
#' function with the trajectory thrown away; call [projectUntilSettled()]
#' instead if you want to watch the approach.
#'
#' `solver = "newton"` solves the steady-state equation directly with a
#' Newton-type root finder from the `nleqslv` package, so it converges even when
#' the steady state is dynamically **unstable**, where the time-stepping solver
#' diverges away from it. This is the natural entry point for a stability
#' analysis with [getStability()].
#'
#' The Newton solver treats the resource densities as unknowns alongside the
#' fish and appends a resource steady-state equation to the system, so the
#' resource density and the feeding levels it implies are self-consistent even
#' where consumers are satiated. It obtains that equation from a function named
#' `steady_<resource_dynamics>()`, which returns the resource equilibrium
#' implied by the rates at the current state. Mizer supplies companions for
#' [resource_semichemostat()] and [resource_logistic()]. A custom resource
#' dynamics function can support the Newton solver by supplying its own
#' companion with the same naming and argument convention.
#'
#' Resource densities are solved for in log space and must therefore remain
#' positive. For logistic resource dynamics this solver finds the positive
#' branch only; use `solver = "project"` when a resource size class can be
#' depleted to zero.
#'
#' The Newton iteration also needs the residual \eqn{F(N)} to be continuous. A
#' custom rate function registered with [setRateFunction()] that jumps as a
#' function of the abundances makes \eqn{F} discontinuous, and where the
#' equilibrium lies on the switching threshold there is no root at all, because
#' neither branch is in equilibrium there. The solver then stalls (`nleqslv`
#' termination code 3) and returns an iterate pinned to the threshold. See
#' [Discontinuous rate
#' functions](https://sizespectrum.org/mizer/articles/discontinuous_rates.html).
#'
#' It also respects the active transport scheme: if the experimental
#' second-order scheme is enabled (see [second_order_w()]) it solves the
#' steady-state equation of that scheme. With the van Leer reconstruction the
#' residual is only Lipschitz, so the iteration converges to a fixed point of the
#' dynamics but not to machine precision. The unlimited `"centred"`
#' reconstruction admits an undamped odd-even mode at a steady state with no
#' physical diffusion, giving an ill-conditioned steady-state Jacobian for which
#' the solver is not expected to converge.
#'
#' @template section_check_steady
#'
#' @inheritParams tuneSteadyState
#' @param ... Arguments for the chosen solver.
#'
#'   With `solver = "project"`: `distance_func`, `t_max`, `t_check`, `dt`,
#'   `distance_tol`, `residual_tol`, `amplitude_tol`, `amp_rel_tol`,
#'   `extinction_threshold`, `progress_bar` and `method`, all as described in
#'   [projectUntilSettled()]. There is no `t_save`, because no trajectory is
#'   returned.
#'
#'   With `solver = "newton"`: `extinction_floor` (default `1e-6`), the relative
#'   abundance below which a species counts as extinct, plus `solver_tol`,
#'   `maxit`, `jacobian`, `global` and `verbose` as described in
#'   [tuneSteadyState()].
#'
#' @return A `MizerParams` object with the initial state replaced by the steady
#'   state found and no parameter changed. It carries a `"convergence"`
#'   attribute describing the solution found; see [projectUntilSettled()].
#' @seealso [tuneSteadyState()], [projectUntilSettled()], [isSteady()],
#'   [getSteadyResidual()], [getStability()]
#' @export
#' @examples
#' \donttest{
#' params <- findSteadyState(NS_params, solver = "newton")
#' plotSpectra(params)
#' }
findSteadyState <- function(params, solver = c("project", "newton"),
                            effort = params@initial_effort,
                            info_level = default_info_level(), ...) {
    UseMethod("findSteadyState")
}

#' @export
findSteadyState.MizerParams <- function(params,
                                        solver = c("project", "newton"),
                                        effort = params@initial_effort,
                                        info_level = default_info_level(),
                                        ...) {
    solver <- match.arg(solver)
    with_info_level(info_level = info_level, {
        params <- validParams(params)
        effort <- validEffortVector(effort, params = params)
        params@initial_effort <- effort
        if (solver == "project") {
            # Nothing is held fixed here: the projection advances the other
            # components along with everything else, so any dynamics are fine.
            project_until_settled(params, effort = effort,
                                  info_level = info_level, ...,
                                  return_sim = FALSE)
        } else {
            warn_other_components_fixed(params, paste(
                "Use `solver = \"project\"`, which advances them like",
                "everything else, if they need to settle too."))
            find_steady_newton(params, effort = effort, ...)
        }
    })
}


#### Solver branches ####

#' Tune to a steady state by running the dynamics
#'
#' The `solver = "project"` branch of [tuneSteadyState()], and the whole of the
#' superseded [steady()]. The two differ only in that `steady()` can ask for the
#' `MizerSim` of the run, which is why `return_sim` lives here rather than on
#' [tuneSteadyState()].
#'
#' The defaults of the projection arguments are the ones `steady()` shipped
#' with, which are not those of [projectUntilSettled()]: `distance_tol` is
#' `0.1 * dt`
#' rather than `0.1 * t_check` because the distance function is
#' [distanceMaxRelRDI()] rather than [distanceSSLogN()], and `amp_rel_tol` is
#' `0.01` rather than `0.1`. `t_save` matters only to `steady(return_sim =
#' TRUE)`, and defaults to `t_check` there so that the trajectory keeps the
#' spacing that function has always given it.
#'
#' @inheritParams tuneSteadyState
#' @inheritParams projectUntilSettled
#' @param return_sim Whether to return the `MizerSim` of the run instead of the
#'   tuned `MizerParams`.
#' @param require_steady Whether the biomass drift must also be within
#'   `residual_tol` before the run may stop; `FALSE` restores the stopping rule
#'   that [steady()] shipped with. See `project_until_settled()`.
#' @return A `MizerParams`, or a `MizerSim` if `return_sim = TRUE`, carrying the
#'   `"convergence"` attribute.
#' @noRd
tune_steady_project <- function(params, effort, preserve,
                                t_max = 100, t_check = 15 * dt, dt = 0.1,
                                t_save = t_check, distance_tol = 0.1 * dt,
                                residual_tol = steady_residual_tol(),
                                amplitude_tol = 0.01, amp_rel_tol = 0.01,
                                extinction_threshold = 1e-6,
                                return_sim = FALSE, require_steady = TRUE,
                                progress_bar = TRUE,
                                info_level = default_info_level(),
                                method = c("euler", "predictor_corrector",
                                           "tr_bdf2")) {
    method <- normalise_project_method(method)

    old_reproduction_level <- NULL
    old_R_max <- NULL
    old_erepro <- NULL
    if (params@rates_funcs$RDD == "BevertonHoltRDD") {
        preserve <- match.arg(preserve, c("reproduction_level", "erepro",
                                          "R_max"))
        old_reproduction_level <- reproduction_level(params)
        old_R_max <- params@species_params$R_max
        old_erepro <- params@species_params$erepro
    }

    # Force the reproduction to stay at the current level
    params@species_params$constant_reproduction <- getRDD(params)
    old_rdd_fun <- params@rates_funcs$RDD
    params@rates_funcs$RDD <- "constantRDD"

    # Force resource to stay at current level
    old_resource_dynamics <- params@resource_dynamics
    params@resource_dynamics <- "resource_constant"

    # Force other components to stay at current level
    old_other_dynamics <- params@other_dynamics
    for (res in names(params@other_dynamics)) {
        params@other_dynamics[[res]] <- "constant_other"
    }

    object <- project_until_settled(params,
                                    effort = effort,
                                    distance_func = distanceMaxRelRDI,
                                    t_check = t_check,
                                    t_max = t_max,
                                    dt = dt,
                                    t_save = t_save,
                                    distance_tol = distance_tol,
                                    residual_tol = residual_tol,
                                    amplitude_tol = amplitude_tol,
                                    amp_rel_tol = amp_rel_tol,
                                    extinction_threshold = extinction_threshold,
                                    return_sim = return_sim,
                                    require_steady = require_steady,
                                    progress_bar = progress_bar,
                                    info_level = info_level,
                                    method = method)
    # Capture the convergence diagnostic before the setter functions below
    # return fresh objects that drop attributes; it is re-attached at the end.
    conv <- attr(object, "convergence")
    if (return_sim) {
        params <- object@params
    } else {
        params <- object
    }
    # Restore original RDD and dynamics
    params@rates_funcs$RDD <- old_rdd_fun
    params@other_dynamics <- old_other_dynamics
    params@species_params$constant_reproduction <- NULL

    params <- rebalance_resource(params, old_resource_dynamics)
    params <- restore_reproduction(params, preserve, old_reproduction_level,
                                   old_R_max, old_erepro)

    # The residual recorded during the run was measured on the model it was
    # handed, which had reproduction, the resource and the other components
    # pinned. Restoring the real dynamics above changes the model, so the
    # residual is measured again on what is actually being returned, and the
    # verdict on the state is derived from that measurement.
    #
    # This is also why a cycle detected during the search does not make the
    # returned model a limit-cycle attractor: the search was run with the inputs
    # to the fish dynamics pinned, so what oscillated was a constrained system,
    # and the model handed back is not that system. `termination` still records
    # that the search ended on a detected cycle.
    conv$residual <- tryCatch(steady_biomass_drift(params),
                              error = function(e) NA_real_)
    conv$attractor <- steady_attractor(conv$residual, residual_tol)

    if (return_sim) {
        object@params <- params
        attr(object, "convergence") <- conv
        return(object)
    }
    params@time_modified <- lubridate::now()
    attr(params, "convergence") <- conv
    params
}

#' Tune to a steady state with the Newton solver
#'
#' The `solver = "newton"` branch of [tuneSteadyState()]. The reproduction rate
#' is held constant and the resource is left out of the system altogether, so
#' the unknowns are the consumer densities alone; the capacity is rebalanced
#' afterwards, exactly as in the projecting branch.
#'
#' @inheritParams tuneSteadyState
#' @return A `MizerParams` carrying the `"convergence"` attribute.
#' @noRd
tune_steady_newton <- function(params, effort, preserve,
                               verbose = FALSE, solver_tol = 1e-6,
                               residual_tol = steady_residual_tol(),
                               maxit = 200,
                               jacobian = c("update", "recompute"),
                               global = "dbldog") {
    old_reproduction_level <- NULL
    old_R_max <- NULL
    old_erepro <- NULL
    if (params@rates_funcs$RDD == "BevertonHoltRDD") {
        preserve <- match.arg(preserve, c("reproduction_level", "erepro",
                                          "R_max"))
        old_reproduction_level <- reproduction_level(params)
        old_R_max <- params@species_params$R_max
        old_erepro <- params@species_params$erepro
    }
    old_resource_dynamics <- params@resource_dynamics

    sol <- newton_steady_state(params, effort = effort,
                               rdd_const = getRDD(params),
                               resource = "fixed",
                               verbose = verbose,
                               solver_tol = solver_tol, maxit = maxit,
                               jacobian = jacobian, global = global)

    params@initial_n[] <- sol$n
    params@initial_n_pp[] <- sol$n_pp

    params <- rebalance_resource(params, old_resource_dynamics)
    params <- restore_reproduction(params, preserve, old_reproduction_level,
                                   old_R_max, old_erepro)

    finish_newton(params, sol, effort, residual_tol)
}

#' Find a steady state with the Newton solver
#'
#' The `solver = "newton"` branch of [findSteadyState()]. Reproduction runs
#' through the model's own RDI and RDD functions and the resource densities are
#' unknowns of the system, so no parameter is touched.
#'
#' @inheritParams findSteadyState
#' @param extinction_floor The relative abundance floor below which a species is
#'   considered extinct.
#' @return A `MizerParams` carrying the `"convergence"` attribute.
#' @noRd
find_steady_newton <- function(params, effort, extinction_floor = 1e-6,
                               verbose = FALSE, solver_tol = 1e-6,
                               residual_tol = steady_residual_tol(),
                               maxit = 200,
                               jacobian = c("update", "recompute"),
                               global = "dbldog") {
    sol <- newton_steady_state(params, effort = effort, rdd_const = NULL,
                               resource = "solve",
                               extinction_floor = extinction_floor,
                               verbose = verbose,
                               solver_tol = solver_tol, maxit = maxit,
                               jacobian = jacobian, global = global)

    params@initial_n[] <- sol$n
    params@initial_n_pp[] <- sol$n_pp

    finish_newton(params, sol, effort, residual_tol)
}


#### Shared pieces ####

#' Report that the analysis covers the fish and the resource only
#'
#' `tuneSteadyState()`, the Newton solver and the stability analyses all treat
#' the components registered with [setComponent()] as fixed inputs: the search
#' holds them at their stored values and the Jacobian has no rows for them.
#' Where those components have dynamics of their own, the result is a statement
#' about a subsystem, and this is what says so.
#'
#' It is a report rather than a refusal, for the same reason that
#' `warn_if_not_steady()` is: the analysis is still the right one whenever the
#' components are slaved to the fish or move on a far slower timescale, and the
#' user is the one who knows which. Refusing would also lock every model built
#' with [setComponent()] out of these functions. What the user must not be
#' allowed to do is *assume* the components were covered, which is why this is a
#' warning rather than a quiet note, and why the returned `residual` — which
#' does cover them — is reported alongside.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param context A sentence naming what is being claimed, appended to the
#'   report.
#' @return `TRUE` if a report was raised, `FALSE` otherwise, invisibly.
#' @noRd
warn_other_components_fixed <- function(params, context) {
    dynamic <- names(params@other_dynamics)[
        !vapply(params@other_dynamics,
                function(f) identical(f, "constant_other"), logical(1))]
    if (length(dynamic) == 0) return(invisible(FALSE))
    signal_info("other_components", paste0(
        "The component", if (length(dynamic) > 1) "s " else " ",
        paste0("`", dynamic, "`", collapse = ", "),
        if (length(dynamic) > 1) " have" else " has",
        " dynamics of their own, and mizer's steady-state and stability ",
        "machinery covers the consumers and the resource only: ",
        if (length(dynamic) > 1) "they are" else "it is",
        " held at the stored value throughout. ", context),
        level = 1, severity = "warning", unhandled = "show")
    invisible(TRUE)
}

#' Rebalance the resource so that its preserved abundance is steady
#'
#' Restores the original resource dynamics and derives the capacity from the
#' (preserved) rate, so that the abundance the search held fixed is a genuine
#' steady state of the resource under the new spectra. A manually set (frozen)
#' capacity would otherwise block this, leaving the resource off its fixed
#' point, so the frozen mark is cleared first — exactly as earlier versions of
#' mizer did when they rebalanced the resource in `steady()`.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param resource_dynamics The resource dynamics to restore.
#' @return The `MizerParams` object with the capacity rebalanced.
#' @noRd
rebalance_resource <- function(params, resource_dynamics) {
    # Balancing is what makes the preserved resource abundance a steady state of
    # the restored dynamics. `setResource()` can only do it where a
    # `balance_<dynamics>()` function exists, and otherwise says nothing, so a
    # custom resource would silently come back off its own fixed point.
    if (!is.function(get0(paste0("balance_", resource_dynamics)))) {
        signal_info("cc_pp", paste0(
            "There is no `balance_", resource_dynamics, "()` function, so the ",
            "resource capacity could not be rebalanced and the preserved ",
            "resource abundance need not be a steady state of `",
            resource_dynamics, "()`. The reported residual measures how far ",
            "off it is."),
            level = 1, severity = "warning", unhandled = "show")
    }
    comment(params@cc_pp) <- NULL
    setResource(params, resource_dynamics = resource_dynamics)
}

#' Restore density-dependent reproduction after a constant-RDD search
#'
#' The search held the reproduction rate constant; this puts the Beverton-Holt
#' density dependence back so that it reproduces exactly that rate at the new
#' spectra, keeping whichever of the reproduction parameters `preserve` names.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param preserve One of `"reproduction_level"`, `"R_max"` or `"erepro"`.
#' @param old_reproduction_level,old_R_max,old_erepro The values recorded before
#'   the search.
#' @return The `MizerParams` object with reproduction restored.
#' @noRd
restore_reproduction <- function(params, preserve, old_reproduction_level,
                                 old_R_max, old_erepro) {
    if (params@rates_funcs$RDD != "BevertonHoltRDD") {
        return(params)
    }
    if (preserve == "reproduction_level") {
        reproduction_level(params) <- old_reproduction_level
    } else if (preserve == "R_max") {
        params <- setBevertonHolt(params, R_max = old_R_max)
    } else {
        params <- setBevertonHolt(params, erepro = old_erepro)
    }
    params
}

#' Report on a Newton solve and attach its convergence diagnostic
#'
#' Gives the Newton solvers the same `"convergence"` attribute that
#' `project_until_settled()` attaches, so that callers reading it — `scanModel()`
#' among them — do not have to know which solver produced the object. The
#' fields that only mean something for a run over time (`distance`, `years`,
#' `period`, `amplitude`) are `NA`.
#'
#' @param params A \linkS4class{MizerParams} object holding the solution.
#' @param sol The list returned by `newton_steady_state()`.
#' @param effort The fishing effort the solve used.
#' @param residual_tol The biomass drift below which the solution counts as a
#'   fixed point.
#' @return The `MizerParams` object with the attribute attached.
#' @noRd
finish_newton <- function(params, sol, effort,
                          residual_tol = steady_residual_tol()) {
    # A species that had no density to begin with is not news of an extinction:
    # nothing died during the solve, it was already absent and was held at zero
    # rather than solved for. Reporting it as an extinction would send the user
    # looking for a collapse that did not happen.
    if (!is.null(sol$absent) && any(sol$absent)) {
        signal_info("convergence", paste0(
            "The following species were already absent and were held at zero ",
            "rather than solved for: ",
            paste(names(sol$absent)[sol$absent], collapse = ", "), "."),
            level = 1, severity = "warning", unhandled = "show")
    }
    if (any(sol$extinct)) {
        warning("The following species went extinct and were set to zero: ",
                paste(names(sol$extinct)[sol$extinct], collapse = ", "))
    }
    if (sol$termcd > 2) {
        warning("The Newton solver did not converge (nleqslv termination code ",
                sol$termcd, ": ", sol$message,
                "). Returning the best iterate found.", call. = FALSE)
    }

    params@time_modified <- lubridate::now()

    # Report how close to a fixed point the solver actually got. This is the
    # honest measure of the result: with the van Leer reconstruction the
    # residual is only Lipschitz and cannot reach machine precision, so the
    # number is not always as small as `solver_tol` suggests.
    residual <- tryCatch(steady_biomass_drift(params, effort = effort),
                         error = function(e) NA_real_)
    if (is.finite(residual)) {
        signal_info("convergence", paste0(
            "The biomasses of the solution change at up to ",
            signif(residual, 2), " per year."),
            level = 3, unhandled = "show")
    }

    # `converged` is about the root finder and `attractor` about the state it
    # returned. They can disagree in both directions: a solve that stopped on
    # its own criterion can still sit at a state that drifts (the van Leer
    # residual is only Lipschitz), and a solve that gave up can leave a
    # perfectly good fixed point behind.
    converged <- sol$termcd <= 2
    termination <- if (any(sol$extinct)) {
        "extinction"
    } else if (converged) {
        "solver_converged"
    } else {
        "solver_failed"
    }
    attr(params, "convergence") <- convergence_result(
        termination = termination,
        converged = converged,
        attractor = steady_attractor(residual, residual_tol),
        residual = residual,
        extinct = if (any(sol$extinct)) {
            names(sol$extinct)[sol$extinct]
        } else {
            character(0)
        }
    )
    params
}
