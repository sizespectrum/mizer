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
#' Solves for a model's steady state using declared tuning strategies. The
#' built-in defaults hold final recruitment and the resource abundance fixed,
#' then adjust reproduction parameters and resource capacity. Extensions declare
#' which states their strategies hold or evolve, what parameters may change,
#' and how the held targets are balanced. The restored candidate is checked
#' against every effective equation rate before its final steady assessment.
#'
#' It is not the function for asking where a model with given parameters
#' settles. Because the reproduction rate is held at its current value, it
#' cannot reach a steady state in which a species' reproduction rate is very
#' different, such as an alternative stable state or a state that exists
#' because predation changes the reproduction of the prey. Because it
#' re-derives `erepro`/`R_max` and `cc_pp`, two states found with it do not
#' share parameters. Use [findSteadyState()] for those questions.
#'
#' @details
#' Each provider supplies named strategies and a default. Use `strategies` to
#' choose alternatives, and `joint` for declared strategies owning coupled
#' processes. The selected strategies must cover the whole tuning operation
#' without overlapping parameter ownership. Missing support, unknown stationary
#' behaviour and unsupported solver combinations stop before the search.
#'
#' Recruitment targets come from the complete effective rate pipeline,
#' including named replacements and extension wrappers. A held target replaces
#' the final recruitment entering the equations, so wrappers are not applied a
#' second time to a total that already includes their contribution.
#'
#' Balancers run once after the constrained search. They must preserve the
#' candidate state and all effective target equation rates, including growth,
#' diffusion, total mortality and final recruitment. Mizer independently checks
#' the released model; it does not repair changed rates by solving the consumers
#' again or by choosing an undeclared numerical balancing method.
#'
#' The `"convergence"` attribute separates `search_assessment`, measured with
#' the selected state and rate constraints, from `assessment`, measured after
#' restoring and balancing the model. Whole-model assessment is the default.
#' Select `scope = "conditional"` with explicit `fixed` owners for a subsystem;
#' those inputs must also be held by the selected tuning strategies.
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
#' function. It needs no extra packages. Every required dynamic state must
#' provide convergence support.
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
#' not need the resource to be a semichemostat, unlike in [findSteadyState()]
#' where the resource is one of the unknowns. The Newton solver does not solve
#' extension component states; dynamic components require explicit conditional
#' exclusions with `fixed`, or use the projection solver.
#'
#' @templateVar rederived TRUE
#' @template section_check_steady
#'
#' @param params A \linkS4class{MizerParams} object
#' @param solver The solver to use: `"project"` to run the dynamics until they
#'   settle, `"newton"` to solve the steady-state equation directly. See
#'   *Choosing a solver*.
#' @param effort The fishing effort to use throughout. By default the initial
#'   effort stored in `params`.
#' @param preserve `r lifecycle::badge("experimental")`
#'   An explicit reproduction-preservation choice, such as `"erepro"`,
#'   `"R_max"` or `"reproduction_level"` for Beverton-Holt reproduction.
#'   Omitting this argument selects the provider's declared default strategy.
#'   Other laws support their own controls; unsupported choices error. Do not
#'   combine this selector with an explicit reproduction entry in `strategies`.
#' @param strategies Optional named character vector selecting a strategy for
#'   each specified owner, such as `c(resource = "evolve")` or
#'   `c("component:pool" = "hold")`. Omitted owners use their declared defaults.
#' @param joint Names of joint strategies registered with [setTuningStrategy()].
#'   `NULL` selects declared joint defaults; `character()` selects none. A joint
#'   strategy replaces the local defaults of its owners. Conflicting explicit
#'   local selections or overlapping joint ownership error.
#' @param info_level Controls the amount of information messages that are shown.
#'   Higher levels lead to more messages, `info_level = 0` gives silence. The
#'   default is taken from the `mizer_info_level` option, see
#'   [default_info_level()].
#' @param ... Arguments for the chosen solver.
#'
#'   With `solver = "project"`: `t_max`, `t_check`, `dt`,
#'   `distance_tol`, `residual_tol`, `amplitude_tol`, `amp_rel_tol`,
#'   `extinction_threshold`, `progress_bar`, `method`, `subcycle`,
#'   `max_courant` and `max_substeps`, all as described in
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
#'   the per-capita rate of change passed to [nleqslv::nleqslv()]. It was called
#'   `residual_tol` before mizer 3.3, a name now used by the assessment
#'   criterion above; `maxit`
#'   (default `200`); `jacobian`, either `"update"` (default, the Jacobian is
#'   computed once and then updated cheaply each iteration — `nleqslv`'s
#'   `"Broyden"`) or `"recompute"` (a numerical Jacobian at every iteration —
#'   `nleqslv`'s `"Newton"`); `global`, the globalisation strategy (default
#'   `"dbldog"`, a robust double-dogleg trust region); and `verbose` to trace
#'   the iterations.
#'
#' @inheritParams getSteadyAssessment
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
                            info_level = default_info_level(),
                            scope = c("whole", "conditional"),
                            fixed = character(),
                            criterion = c("state", "biomass"),
                            tolerances = NULL, atol = 1e-12,
                            strategies = NULL, joint = NULL, ...) {
    UseMethod("tuneSteadyState")
}

#' @export
tuneSteadyState.MizerParams <- function(params,
                                        solver = c("project", "newton"),
                                        effort = params@initial_effort,
                                        preserve = c("reproduction_level",
                                                     "erepro", "R_max"),
                                        info_level = default_info_level(),
                                        scope = c("whole", "conditional"),
                                        fixed = character(),
                                        criterion = c("state", "biomass"),
                                        tolerances = NULL, atol = 1e-12,
                                        strategies = NULL, joint = NULL, ...) {
    if (missing(preserve)) preserve <- NULL
    solver <- match.arg(solver)
    with_info_level(info_level = info_level, {
        params <- validParams(params)
        effort <- validEffortVector(effort, params = params)
        params@initial_effort <- effort
        if (solver == "project") {
            tune_steady_project(params, effort = effort, preserve = preserve,
                                info_level = info_level, scope = scope,
                                fixed = fixed, criterion = criterion,
                                tolerances = tolerances, atol = atol,
                                strategies = strategies, joint = joint, ...)
        } else {
            tune_steady_newton(params, effort = effort, preserve = preserve,
                               scope = scope, fixed = fixed,
                               criterion = criterion, tolerances = tolerances,
                               atol = atol, strategies = strategies,
                               joint = joint, ...)
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
#' fish and appends the resource steady-state equation to the system, so the
#' resource density and the feeding levels it implies are self-consistent even
#' where consumers are satiated. That equation is written for the default
#' semichemostat resource dynamics, so `solver = "newton"` stops with an error
#' for any other dynamic `resource_dynamics`; use `solver = "project"` there,
#' or explicitly fix the resource in a conditional analysis. Dynamic extension
#' components likewise require explicit conditional exclusions for Newton.
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
#'   `extinction_threshold`, `progress_bar`, `method`, `subcycle`,
#'   `max_courant` and `max_substeps`, all as described in
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
                            info_level = default_info_level(),
                            scope = c("whole", "conditional"),
                            fixed = character(),
                            criterion = c("state", "biomass"),
                            tolerances = NULL, atol = 1e-12, ...) {
    UseMethod("findSteadyState")
}

#' @export
findSteadyState.MizerParams <- function(params,
                                        solver = c("project", "newton"),
                                        effort = params@initial_effort,
                                        info_level = default_info_level(),
                                        scope = c("whole", "conditional"),
                                        fixed = character(),
                                        criterion = c("state", "biomass"),
                                        tolerances = NULL, atol = 1e-12, ...) {
    solver <- match.arg(solver)
    with_info_level(info_level = info_level, {
        params <- validParams(params)
        effort <- validEffortVector(effort, params = params)
        params@initial_effort <- effort
        if (solver == "project") {
            # Projection advances every state except explicitly fixed inputs.
            project_until_settled(params, effort = effort,
                                  info_level = info_level, scope = scope,
                                  fixed = fixed, criterion = criterion,
                                  tolerances = tolerances, atol = atol, ...,
                                  return_sim = FALSE)
        } else {
            find_steady_newton(params, effort = effort, scope = scope,
                               fixed = fixed, criterion = criterion,
                               tolerances = tolerances, atol = atol, ...)
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
                                           "tr_bdf2"),
                                subcycle = FALSE, max_courant = 1,
                                max_substeps = 10000,
                                scope = c("whole", "conditional"),
                                fixed = character(),
                                criterion = c("state", "biomass"),
                                tolerances = NULL, atol = 1e-12,
                                strategies = NULL, joint = NULL) {
    method <- normalise_project_method(method)
    prepare_solver_assessment(params, residual_tol, scope, fixed, criterion,
                              tolerances, atol, require_steady = require_steady)
    preserve <- tuning_preserve_selection(preserve)
    plan <- prepare_tuning(params, strategies = strategies, joint = joint,
                           solver = "project", effort = effort,
                           preserve = preserve)
    runtime <- prepare_tuning_runtime(params, plan, fixed, "project")
    search_tolerances <- tuning_search_tolerances(tolerances,
                                                  runtime$hold_states)
    object <- project_until_settled(
        runtime$params, effort = effort, distance_func = distanceMaxRelRDI,
        t_check = t_check, t_max = t_max, dt = dt, t_save = t_save,
        distance_tol = distance_tol, residual_tol = residual_tol,
        amplitude_tol = amplitude_tol, amp_rel_tol = amp_rel_tol,
        extinction_threshold = extinction_threshold, return_sim = return_sim,
        require_steady = require_steady, progress_bar = progress_bar,
        info_level = info_level, method = method, subcycle = subcycle,
        max_courant = max_courant, max_substeps = max_substeps,
        scope = runtime$scope, fixed = runtime$hold_states,
        criterion = criterion, tolerances = search_tolerances, atol = atol,
        .rates_fns = runtime$rates_fns,
        .rate_constraints = runtime$rate_constraints
    )
    conv <- attr(object, "convergence")
    candidate <- if (return_sim) object@params else object
    candidate@resource_dynamics <- params@resource_dynamics
    candidate@other_dynamics <- params@other_dynamics
    rates <- tuning_target_rates(candidate, runtime$rates_fns, effort,
                                 t = conv$years)
    targets <- capture_tuning_targets(candidate, plan, rates = rates,
                                      t = conv$years)
    candidate <- finish_tuning(candidate, plan, targets)

    conv$search_assessment <- conv$assessment
    assessment_plan <- prepare_solver_assessment(
        candidate, residual_tol, scope, fixed, criterion, tolerances, atol,
        require_steady = require_steady
    )
    conv$assessment <- steady_assessment(candidate, prepared = assessment_plan,
                                         effort = effort, t = conv$years)
    conv$scope <- conv$assessment$scope
    conv$residual <- attr(conv$assessment, "drift_report")$drift
    conv$attractor <- steady_attractor(conv$assessment)
    candidate@time_modified <- lubridate::now()
    if (return_sim) {
        object@params <- candidate
        attr(object, "convergence") <- conv
        return(object)
    }
    attr(candidate, "convergence") <- conv
    candidate
}

# Preserve the legacy wrapper's unevaluated default without imposing a
# Beverton-Holt control on another provider's declared strategy.
tuning_preserve_selection <- function(preserve) {
    if (identical(preserve, c("reproduction_level", "erepro", "R_max"))) {
        NULL
    } else {
        preserve
    }
}

prepare_tuning_runtime <- function(params, plan, fixed, solver) {
    hold <- plan$hold_states
    evolve <- plan$evolve_states
    if ("consumers" %in% hold) {
        stop("Tuning solvers cannot hold consumers fixed.", call. = FALSE)
    }
    if (length(setdiff(fixed, hold))) {
        stop("Conditional fixed inputs must be held by the selected tuning ",
             "strategies: ", paste(setdiff(fixed, hold), collapse = ", "),
             ". Select a compatible strategy.", call. = FALSE)
    }
    if (solver == "newton" && any(startsWith(evolve, "component:"))) {
        stop("The Newton solver cannot evolve component states. Select a ",
             "strategy that holds them or use `solver = 'project'`.",
             call. = FALSE)
    }
    if (solver == "newton" && "resource" %in% evolve &&
            sub("^mizer::", "", params@resource_dynamics) !=
                "resource_semichemostat") {
        stop("The Newton solver can evolve only a semichemostat resource. ",
             "Select a strategy that holds it or use `solver = 'project'`.",
             call. = FALSE)
    }
    rates_fns <- projectRateFunctions(params)
    constraints <- character()
    if (isTRUE(plan$hold_recruitment)) {
        rates_fns <- pin_recruitment_rates(rates_fns, plan$recruitment)
        constraints <- "consumers/recruitment"
    }
    if ("resource" %in% hold) params@resource_dynamics <- "resource_constant"
    components <- hold[startsWith(hold, "component:")]
    for (owner in components) {
        component <- substring(owner, nchar("component:") + 1L)
        params@other_dynamics[[component]] <- "constant_other"
    }
    list(params = params, rates_fns = rates_fns, hold_states = hold,
         rate_constraints = constraints,
         scope = if (length(hold) || length(constraints)) {
             "conditional"
         } else {
             "whole"
         })
}

tuning_search_tolerances <- function(tolerances, hold) {
    if (!length(tolerances) || !length(hold)) return(tolerances)
    excluded <- vapply(names(tolerances), function(id) {
        any(startsWith(id, paste0(hold, "/")))
    }, logical(1))
    tolerances <- tolerances[!excluded]
    if (!length(tolerances)) NULL else tolerances
}

tuning_target_rates <- function(params, rates_fns, effort, t = 0) {
    rates_fns$Rates(params, n = params@initial_n,
                    n_pp = params@initial_n_pp,
                    n_other = params@initial_n_other, t = t, effort = effort,
                    rates_fns = rates_fns)
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
                               global = "dbldog",
                               scope = c("whole", "conditional"),
                               fixed = character(),
                               criterion = c("state", "biomass"),
                               tolerances = NULL, atol = 1e-12,
                               strategies = NULL, joint = NULL) {
    prepare_solver_assessment(params, residual_tol, scope, fixed, criterion,
                              tolerances, atol)
    require_newton_equations(params)
    preserve <- tuning_preserve_selection(preserve)
    plan <- prepare_tuning(params, strategies = strategies, joint = joint,
                           solver = "newton", effort = effort,
                           preserve = preserve)
    runtime <- prepare_tuning_runtime(params, plan, fixed, "newton")
    resource <- if ("resource" %in% plan$evolve_states) "solve" else "fixed"
    recruitment <- if (isTRUE(plan$hold_recruitment)) plan$recruitment else NULL
    sol <- newton_steady_state(runtime$params, effort = effort,
                               rdd_const = recruitment, resource = resource,
                               verbose = verbose, solver_tol = solver_tol,
                               maxit = maxit, jacobian = jacobian,
                               global = global, rates_fns = runtime$rates_fns)
    candidate <- runtime$params
    candidate@initial_n[] <- sol$n
    candidate@initial_n_pp[] <- sol$n_pp
    search_plan <- prepare_solver_assessment(
        candidate, residual_tol, runtime$scope, runtime$hold_states,
        criterion, tuning_search_tolerances(tolerances, runtime$hold_states),
        atol, rate_constraints = runtime$rate_constraints
    )
    search_assessment <- steady_assessment(
        candidate, prepared = search_plan, effort = effort,
        rates_fns = runtime$rates_fns
    )
    candidate@resource_dynamics <- params@resource_dynamics
    candidate@other_dynamics <- params@other_dynamics
    rates <- tuning_target_rates(candidate, runtime$rates_fns, effort)
    targets <- capture_tuning_targets(candidate, plan, rates = rates)
    candidate <- finish_tuning(candidate, plan, targets)
    candidate <- finish_newton(candidate, sol, effort, residual_tol,
                               scope = scope, fixed = fixed,
                               criterion = criterion, tolerances = tolerances,
                               atol = atol)
    conv <- attr(candidate, "convergence")
    conv$search_assessment <- search_assessment
    attr(candidate, "convergence") <- conv
    candidate
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
                               global = "dbldog",
                               scope = c("whole", "conditional"),
                               fixed = character(),
                               criterion = c("state", "biomass"),
                               tolerances = NULL, atol = 1e-12) {
    assessment_plan <- prepare_solver_assessment(
        params, residual_tol, scope, fixed, criterion, tolerances, atol,
        newton = TRUE
    )
    resource_mode <- if ("resource" %in%
                             c(assessment_plan$fixed,
                               assessment_plan$prescribed)) {
        "fixed"
    } else {
        "solve"
    }
    sol <- newton_steady_state(params, effort = effort, rdd_const = NULL,
                               resource = resource_mode,
                               extinction_floor = extinction_floor,
                               verbose = verbose,
                               solver_tol = solver_tol, maxit = maxit,
                               jacobian = jacobian, global = global)

    params@initial_n[] <- sol$n
    params@initial_n_pp[] <- sol$n_pp

    finish_newton(params, sol, effort, residual_tol, scope = scope,
                  fixed = fixed, criterion = criterion,
                  tolerances = tolerances, atol = atol)
}


#### Shared pieces ####

#' Report that the analysis covers the fish and the resource only
#'
#' Used by stability analyses whose Jacobian contains consumers and resource
#' states only. Dynamic components are fixed inputs to that analysis; their
#' separate steady-state assessment does not add rows to the Jacobian.
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
        " dynamics of their own, and this stability analysis covers the ",
        "consumers and the resource only: ",
        if (length(dynamic) > 1) "they are" else "it is",
        " held at the stored value throughout and ",
        if (length(dynamic) > 1) "are" else "is",
        " not included in the biomass drift that mizer reports. See ",
        "`attr(getSteadyResidual(params), \"other\")` for ",
        if (length(dynamic) > 1) "their" else "its",
        " rate of change. ", context),
        level = 1, severity = "warning", unhandled = "show")
    invisible(TRUE)
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
#' @param residual_tol Relative tolerance for the built-in assessment measures.
#' @return The `MizerParams` object with the attribute attached.
#' @noRd
finish_newton <- function(params, sol, effort,
                          residual_tol = steady_residual_tol(),
                          scope = c("whole", "conditional"),
                          fixed = character(),
                          criterion = c("state", "biomass"),
                          tolerances = NULL, atol = 1e-12) {
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
    assessment_plan <- prepare_solver_assessment(
        params, residual_tol, scope, fixed, criterion, tolerances, atol
    )
    assessment <- steady_assessment(params, prepared = assessment_plan,
                                    effort = effort)
    residual <- attr(assessment, "drift_report")$drift
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
        attractor = steady_attractor(assessment),
        residual = residual,
        assessment = assessment,
        extinct = if (any(sol$extinct)) {
            names(sol$extinct)[sol$extinct]
        } else {
            character(0)
        }
    )
    params
}
