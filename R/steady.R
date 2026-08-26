#' Measure distance between current and previous state in terms of RDI
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' This function can be used in [projectUntilSettled()] to decide when sufficient
#' convergence to steady state has been achieved.
#'
#' @param params MizerParams
#' @param current A named list with entries `n`, `n_pp` and `n_other`
#'   describing the current state
#' @param previous A named list with entries `n`, `n_pp` and `n_other`
#'   describing the previous state
#' @return The largest absolute relative change in rdi:
#'   `max(abs((current_rdi - previous_rdi) / previous_rdi))`. If any entry of
#'   `previous_rdi` is zero, the result can be infinite.
#' @family distance functions
#' @concept helper
#' @export
distanceMaxRelRDI <- function(params, current, previous) {
    UseMethod("distanceMaxRelRDI")
}
#' @export
distanceMaxRelRDI.MizerParams <- function(params, current, previous) {
    current_rdi <- getRDI(params, n = current$n, n_pp = current$n_pp,
                          n_other = current$n_other)
    previous_rdi <- getRDI(params, n = previous$n, n_pp = previous$n_pp,
                           n_other = previous$n_other)
    d <- max(abs((current_rdi - previous_rdi) / previous_rdi))
    d[is.nan(d)] <- Inf
    d
}

#' Measure distance between current and previous state in terms of fish abundances
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Calculates the sum squared difference between log(N) in current and previous
#' state. This function can be used in [projectUntilSettled()] to decide when
#' sufficient convergence to steady state has been achieved.
#'
#' Only the size classes that hold fish are measured. A class whose density is
#' zero in either state has no log to take, and one holding a negligible share
#' of its species' biomass has a log that never stops moving: above a size where
#' growth stops the density decays exponentially, so `log(n)` falls by the same
#' amount between every pair of states, for ever. Left in the sum, a single such
#' class holding \eqn{10^{-92}} g of fish can hold the distance above any
#' tolerance indefinitely and stop [projectUntilSettled()] from ever converging,
#' while the biomass drift correctly reports a fixed point. Excluding it changes
#' nothing for a model that does not have one, because the classes that hold the
#' fish are exactly the classes that are kept.
#'
#' @param params MizerParams
#' @param current A named list with entries `n`, `n_pp` and `n_other`
#'   describing the current state
#' @param previous A named list with entries `n`, `n_pp` and `n_other`
#'   describing the previous state
#' @param biomass_share_cutoff `r lifecycle::badge("experimental")`
#'   A finite number between 0 and 1, inclusive, giving the share of a species'
#'   biomass that a size class must hold in the current state to be measured.
#'   `0` measures every class with a nonzero density in both states, which is
#'   what this function did before mizer 3.3.
#' @param ... Unused. Accepted because [projectUntilSettled()] forwards its own
#'   `...` to both the rate functions and the distance function, so an argument
#'   meant for one arrives at the other.
#' @return The sum of squares of the difference in the logs of the fish
#'   abundances `n`, over the size classes that hold fish in both states:
#'   `sum((log(current$n) - log(previous$n))^2)`
#' @family distance functions
#' @concept helper
#' @export
distanceSSLogN <- function(params, current, previous,
                           biomass_share_cutoff = steady_share_cutoff(), ...) {
    UseMethod("distanceSSLogN")
}
#' @export
distanceSSLogN.MizerParams <- function(params, current, previous,
                                       biomass_share_cutoff =
                                           steady_share_cutoff(), ...) {
    # The cutoff is applied to the current state, the one the run has reached
    # and the one the verdict is about. A class filling up from nothing joins
    # the sum as soon as it holds anything; a dying one leaves it as soon as it
    # does not.
    sel <- current$n > 0 & previous$n > 0 &
        !negligible_cells(params, current$n, biomass_share_cutoff)
    sum((log(current$n[sel]) - log(previous$n[sel]))^2)
}

#' Project the dynamics until they settle
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Run the full dynamics, as in [project()], but stop once the run has settled:
#' either the change has slowed down sufficiently, in the sense that the
#' distance between states `t_check` years apart is less than `distance_tol` and
#' the state has stopped drifting, or the run has been recognised as being on a
#' limit cycle. You determine how the distance
#' is calculated.
#'
#' Nothing is held fixed, so the run can only ever end up on an attractor of the
#' dynamics, and that need not be a fixed point: besides a limit cycle it may
#' stop on a species going extinct, or simply at `t_max`. "Settled" is therefore
#' the most this function claims; the state it leaves behind is not necessarily
#' a steady state. See the section below on how to check, and use
#' [findSteadyState()] if what you want is the steady state itself rather than
#' the trajectory leading to it.
#'
#' @details
#' # How the run is organised
#'
#' The dynamics are advanced with time step `dt` exactly as in [project()].
#' Every `t_check` years the function pauses to decide whether to stop, so
#' `t_check` sets how often the stopping criteria are evaluated and also the
#' interval over which change is measured. You should not normally need to set it: it defaults
#' to `15 * dt`, which is an odd multiple of the time step for the reason given
#' below. The run ends at the latest after `t_max` years.
#'
#' Independently of that, the state is stored in the returned `MizerSim` every
#' `t_save` years, exactly as in [project()], and a cheap scalar summary of the
#' state (the biomass of each species) is recorded after every time step. That
#' finely resolved series is what the limit-cycle detection works on, so that a
#' cycle can be found and its period measured even when that period bears no
#' simple relation to `t_check`. The three intervals are independent of each
#' other; `t_check` and `t_save` need only be multiples of `dt`.
#'
#' At each check the following tests are made, in this order.
#'
#' ## 1. Extinction
#'
#' If the reproduction rate (RDD) of any species has fallen below
#' `extinction_threshold` times its value at the start of the run, or has become
#' `NA`, that species is deemed to be on its way to extinction. A warning naming
#' the affected species is issued and the run stops with
#' `type = "extinction"`. Because the criterion is relative to the initial
#' reproduction, a species that starts with zero reproduction is flagged
#' immediately, whereas in [tuneSteadyState()], where reproduction is held constant, a
#' healthy species is never flagged.
#'
#' ## 2. Limit cycle
#'
#' The recorded biomass series is examined to see whether the run has settled
#' onto a limit cycle. If it has, the run stops with `type = "cycle"` and the
#' period and amplitude of the cycle are reported.
#'
#' It is made at every check, whether or not the state looks converged
#' by the measure below. A cycle whose period divides `t_check` puts the two
#' states that the distance function compares at the same phase, so it would
#' otherwise be reported as a fixed point of zero width. The detection works on
#' the biomass series sampled at every time step instead, which is blind to
#' `t_check`.
#'
#' ## 3. Convergence to a fixed point
#'
#' Two things have to hold for the run to stop on a fixed point.
#'
#' First, `distance_func` is called with the state at the previous check and the
#' state at the current one, i.e. with two states
#' `t_check` years apart, and the number it returns must be less than
#' `distance_tol`. What
#' "distance" means is entirely up to that function: the default
#' [distanceSSLogN()] uses the sum of squared changes in log abundance, while
#' [tuneSteadyState()] instead passes [distanceMaxRelRDI()], which uses the
#' largest relative change in egg production.
#'
#' Second, the state actually reached must be a fixed point: the largest
#' relative rate of biomass change there, as measured by [getSteadyResidual()],
#' must be at most `residual_tol`. The distance criterion on its own says only
#' that the state stopped moving on the scale of the distance function, which is
#' a different question — a distance function can be insensitive to the very
#' motion that is left. When the distance criterion is met but the drift is not,
#' the run carries on rather than declaring a fixed point.
#'
#' When both hold the run stops with `termination = "residual_tolerance"`. That
#' is deliberately not called `"steady"`: `residual_tol` is a working tolerance
#' rather than a proof, and the `residual` entry of the `"convergence"`
#' attribute reports the drift that was actually reached.
#'
#' Even so, `t_check` should be an *odd* multiple of `dt`, which is why it
#' defaults to `15 * dt`. A period-2 cycle (period `2 * dt`), the most common
#' numerical oscillation, is otherwise sampled at the same phase at every check,
#' and its amplitude can sit below `amplitude_tol` where the cycle detection
#' deliberately ignores it.
#'
#' If none of the three checks fires before `t_max` is reached, the run stops
#' with `type = "not_converged"`. In every case the outcome is recorded in the
#' `"convergence"` attribute of the returned object, described under *Value*
#' below.
#'
#' # How a limit cycle is detected
#'
#' The detection uses the community-total biomass, on a log scale and with its
#' mean removed, as a scalar signal, sampled after every time step. At least 20
#' steps are needed before any cycle can be reported.
#'
#' 1. **Candidate period.** The autocorrelation function of the signal is
#'    computed up to a lag of half the length of the series, and the first local
#'    maximum with an autocorrelation above `0.5` is taken as the candidate
#'    period. If there is no such peak, or the peak is at a lag of one sample,
#'    no cycle is reported.
#' 2. **Enough history.** The series must cover at least three full candidate
#'    periods. Otherwise the check is deferred to a later block, when more
#'    history has accumulated.
#' 3. **Amplitude.** For each of the last three period-long windows, the
#'    amplitude is measured as the largest over species of the relative
#'    peak-to-trough biomass range `(max - min) / mean`. The amplitude in the
#'    most recent window must exceed `amplitude_tol`; a smaller oscillation is
#'    considered negligible and the state is left to be treated as a fixed
#'    point.
#' 4. **Settled.** The amplitudes of the three successive windows must agree
#'    with each other to within `amp_rel_tol`, and the most recent amplitude must not be
#'    smaller than the oldest by more than `amp_rel_tol`.
#'
#' The last condition is what distinguishes a genuine limit cycle from a slowly
#' decaying spiral towards a stable fixed point: the spiral loses amplitude from
#' one period to the next, the cycle does not. The distinction is necessarily
#' imperfect when the decay is extremely slow, because over any finite run such
#' a spiral is indistinguishable from a cycle. If you need a definitive answer,
#' use [getStability()] on the fixed point found by
#' [findSteadyState()], which
#' works out the eigenvalues of the linearised dynamics instead of watching
#' a trajectory.
#'
#' The reported `period` is a multiple of `dt`, so it is only resolved to that
#' accuracy; reduce `dt` if you need the period more precisely.
#'
#' @template section_check_steady
#' @param params A \linkS4class{MizerParams} object
#' @param effort The fishing effort to be used throughout the simulation.
#'   This is validated by [validEffortVector()] and can therefore be `NULL`, a
#'   single numeric value used for all gears, an unnamed numeric vector with one
#'   entry per gear, or a named numeric vector for some or all gears.
#' @param distance_func A function that will be called at every check with both
#'   the previous and the new state and that should return a number
#'   that in some sense measures the distance between the states. By default
#'   this uses the function [distanceSSLogN()] that you can use as a model for your
#'   own distance function.
#' @param t_check The interval in years at which the run pauses to check whether
#'   it has settled, and hence also the interval over which `distance_func`
#'   measures change. Must be a positive multiple of `dt`. The default
#'   `15 * dt` is an odd multiple of the time step, which is what lets a
#'   period-2 cycle be seen; you should rarely need to change it.
#' @param t_max The maximum number of years to run the simulation. Default is 100.
#' @param dt The time step to use in `project()`.
#' @param t_save The interval in years at which the state is stored in the
#'   returned `MizerSim`, as in [project()]. Must be a positive multiple of
#'   `dt`, but need bear no relation to `t_check`. Default is 1. The state the
#'   run settles on is always the final time point, even when the run stops
#'   between two saves.
#' @param distance_tol The run stops when the number returned by `distance_func`
#'   for two states `t_check` years apart drops below `distance_tol`, provided
#'   the drift criterion below is also met. Its meaning therefore depends on the
#'   distance function you supply. It was called `tol` before mizer 3.3.
#' @param residual_tol `r lifecycle::badge("experimental")`
#'   The largest relative rate of biomass change, in 1/year, that a state may
#'   have and still be called a fixed point. This is the criterion of
#'   [isSteady()] and is measured with [getSteadyResidual()], so it has the same
#'   meaning for every model and every distance function. Default `0.05`.
#'
#'   It is a backstop against a distance function that has gone quiet while the
#'   model is still moving, not the main line of defence against a cycle: an
#'   oscillation of relative amplitude \eqn{A} and period \eqn{T} drifts at up
#'   to \eqn{2\pi A/T} per year, which for a small `amplitude_tol` can be below
#'   any sensible `residual_tol`. Cycles that small are caught by the cycle
#'   detection above, which is why that runs unconditionally.
#' @param amplitude_tol `r lifecycle::badge("experimental")`
#'   The minimum relative biomass amplitude for a persistent oscillation to be
#'   reported as a limit cycle rather than treated as an (effectively steady)
#'   fixed point. This is a fraction of mean biomass and is kept separate from
#'   `distance_tol` (which measures convergence to a fixed point on a different
#'   scale).
#'   Default `0.01`.
#' @param amp_rel_tol `r lifecycle::badge("experimental")`
#'   Maximum relative change of amplitude between successive periods for the
#'   cycle to count as settled. Default `0.01`.
#' @param extinction_threshold `r lifecycle::badge("experimental")`
#'   A species is treated as going extinct, stopping the run, once its
#'   reproduction rate (RDD) falls below this fraction of its value at the start
#'   of the run. For example the default `1e-6` treats a species as extinct once
#'   its reproduction has collapsed to a millionth of its initial level. Because
#'   it is relative to the initial reproduction, a species that starts with zero
#'   reproduction is flagged immediately, and (in [tuneSteadyState()], where
#'   reproduction is held constant) a healthy species is never flagged.
#' @param progress_bar A shiny progress object to implement a progress bar in a
#'   shiny app. Default FALSE.
#' @param info_level Controls the amount of information messages that are shown.
#'   Higher levels lead to more messages, `info_level = 0` gives silence. The
#'   default is taken from the `mizer_info_level` option, see
#'   [default_info_level()].
#' @param method The numerical method to use for the consumer density update.
#'   See [project()].
#' @param ... Further arguments will be passed on to your distance function.
#'
#' @return A `MizerSim` object containing the states saved every `t_save` years,
#'   with the state the run settled on as its final time point. That last
#'   interval is shorter than the others when the run stops between two saves.
#'   Use [finalParams()] to
#'   extract that final state as a `MizerParams` object, or call
#'   [findSteadyState()] instead, which returns it directly.
#'
#'   The returned object carries an attribute `"convergence"` describing the
#'   solution the run settled on, a named list with entries. The first three
#'   answer three different questions and should not be read as one:
#'   \describe{
#'     \item{`termination`}{Why the run stopped: `"residual_tolerance"` (both
#'       convergence criteria were met), `"distance_tolerance"` (the distance
#'       function was satisfied but the state is still drifting — reachable with
#'       a loose `residual_tol`, and from the superseded [steady()], which stops
#'       on the distance criterion alone), `"cycle_detected"` (a limit cycle),
#'       `"time_limit"` (still changing at `t_max`) or `"extinction"` (a species
#'       died out). The steady-state finders can also return
#'       `"solver_converged"` and `"solver_failed"` from the Newton solver.}
#'     \item{`converged`}{Logical, `TRUE` when the run stopped on a criterion of
#'       its own rather than running out of time or losing a species. This is a
#'       statement about the numerics, not about the state.}
#'     \item{`attractor`}{What the state that was reached actually is:
#'       `"fixed_point"` when the biomass drift is within `residual_tol`,
#'       `"limit_cycle"` when a cycle was detected, and `NA` when it is neither
#'       — a run stopped in mid-flight, or a species on its way out. This is the
#'       entry to test before treating a result as a steady state.}
#'     \item{`distance`}{The final value returned by `distance_func`.}
#'     \item{`residual`}{The largest absolute relative rate of biomass change,
#'       in 1/year, at the state that was reached. For each consumer species
#'       this is \eqn{(dB_i/dt) / B_i}, a biomass-weighted aggregate over its
#'       size classes; the resource is treated the same way and the other
#'       components contribute their own relative rates. It is **not** the
#'       largest cellwise value of [getSteadyResidual()], which is dominated by
#'       fast-turnover size classes holding almost no mass and should not be
#'       reduced to its maximum. This is the quantity `residual_tol` is a
#'       tolerance on, so the two cannot mean different things. Unlike
#'       `distance`, which compares two states `t_check` apart on whatever scale
#'       the distance function uses, this measures how far the state actually is
#'       from being a fixed point.}
#'     \item{`years`}{The number of years simulated. `NA` for a direct solve.}
#'     \item{`period`}{For a limit cycle, its period in years; otherwise `NA`.}
#'     \item{`amplitude`}{For a limit cycle, the largest per-species relative
#'       peak-to-trough biomass amplitude; otherwise `NA`.}
#'     \item{`extinct`}{Character vector naming any species that went extinct
#'       during the run, or `character(0)` if none.}
#'   }
#' @seealso [findSteadyState()], [tuneSteadyState()], [isSteady()],
#'   [getSteadyResidual()], [distanceSSLogN()], [distanceMaxRelRDI()],
#'   [getStability()]
#' @export
projectUntilSettled <- function(params,
                                effort = params@initial_effort,
                                distance_func = distanceSSLogN,
                                t_check = 15 * dt,
                                t_max = 100,
                                dt = 0.1,
                                t_save = 1,
                                distance_tol = 0.1 * t_check,
                                residual_tol = steady_residual_tol(),
                                amplitude_tol = 0.01,
                                amp_rel_tol = 0.1,
                                extinction_threshold = 1e-6,
                                progress_bar = TRUE,
                                info_level = default_info_level(),
                                method = c("euler", "predictor_corrector",
                                           "tr_bdf2"), ...) {
    UseMethod("projectUntilSettled")
}

#' @export
projectUntilSettled.MizerParams <- function(params,
                                effort = params@initial_effort,
                                distance_func = distanceSSLogN,
                                t_check = 15 * dt,
                                t_max = 100,
                                dt = 0.1,
                                t_save = 1,
                                distance_tol = 0.1 * t_check,
                                residual_tol = steady_residual_tol(),
                                amplitude_tol = 0.01,
                                amp_rel_tol = 0.1,
                                extinction_threshold = 1e-6,
                                progress_bar = TRUE,
                                info_level = default_info_level(),
                                method = c("euler", "predictor_corrector",
                                           "tr_bdf2"), ...) {
    project_until_settled(params = params, effort = effort,
                          distance_func = distance_func,
                          t_check = t_check, t_max = t_max, dt = dt,
                          t_save = t_save, distance_tol = distance_tol,
                          residual_tol = residual_tol,
                          amplitude_tol = amplitude_tol,
                          amp_rel_tol = amp_rel_tol,
                          extinction_threshold = extinction_threshold,
                          progress_bar = progress_bar, info_level = info_level,
                          method = method, ..., return_sim = TRUE)
}

#' The `"convergence"` attribute, in one place
#'
#' Every steady-state search — the projecting one here and the Newton solves in
#' `steadyState.R` — attaches a description of what it settled on, and callers
#' such as [scanModel()] read it without knowing which solver produced it. The
#' three verdicts it carries answer three different questions and are
#' deliberately not collapsed into one:
#'
#' * `termination` — *why the search stopped*. A fact about the run.
#' * `converged` — *whether the solver met its own criterion*. A fact about the
#'   numerics.
#' * `attractor` — *what the state that was reached is*. A fact about the model,
#'   and the only one of the three that may be used to claim a fixed point.
#'
#' The old `type`/`settled` pair conflated them, which is how a limit cycle
#' could be reported as a converged fixed point.
#'
#' @param termination One of `"residual_tolerance"`, `"cycle_detected"`,
#'   `"time_limit"`, `"extinction"`, `"solver_converged"` or `"solver_failed"`.
#' @param converged Whether the solver met its own convergence criterion.
#' @param attractor `"fixed_point"`, `"limit_cycle"` or `NA`.
#' @param distance The final value of the distance function, or `NA`.
#' @param residual The biomass drift at the state reached, in 1/year.
#' @param years The number of years simulated, or `NA` for a direct solve.
#' @param period,amplitude The period and relative amplitude of a limit cycle.
#' @param extinct Character vector naming species that went extinct.
#' @return A named list.
#' @noRd
convergence_result <- function(termination, converged, attractor,
                               distance = NA_real_, residual = NA_real_,
                               years = NA_real_, period = NA_real_,
                               amplitude = NA_real_,
                               extinct = character(0)) {
    list(termination = termination,
         converged = converged,
         attractor = attractor,
         distance = distance,
         residual = residual,
         years = years,
         period = period,
         amplitude = amplitude,
         extinct = extinct)
}

#' Whether a state counts as a fixed point
#'
#' The single place that turns a measured biomass drift into the claim that a
#' state is a fixed point, so that the claim means the same thing whichever
#' solver produced the state and whatever criterion that solver stopped on.
#'
#' @param residual The biomass drift, from `steady_biomass_drift()`.
#' @param residual_tol The tolerance to judge it against.
#' @return `"fixed_point"` or `NA_character_`.
#' @noRd
steady_attractor <- function(residual, residual_tol) {
    if (is.finite(residual) && residual <= residual_tol) {
        "fixed_point"
    } else {
        NA_character_
    }
}

#' The engine behind [projectUntilSettled()] and the steady-state finders
#'
#' Holds the whole block-by-block run. It is kept separate from the exported
#' [projectUntilSettled()] for one reason: with `return_sim = FALSE` it never
#' allocates the `MizerSim` arrays at all, and [tuneSteadyState()],
#' [findSteadyState()] and [scanModel()] all run it in that mode, the last of
#' them in a loop. The exported function is the `return_sim = TRUE` face of it.
#'
#' @inheritParams projectUntilSettled
#' @param return_sim Whether to build and return a `MizerSim` of the run rather
#'   than a `MizerParams` holding only the final state.
#' @param require_steady Whether the biomass drift has to be within
#'   `residual_tol` before the run may stop on the distance criterion. `TRUE`
#'   for everything mizer exports; the superseded [steady()] and
#'   [projectToSteady()] pass `FALSE`, because released code was written against
#'   a stopping rule that used the distance function alone and a run that
#'   suddenly takes ten times as long is a bad way to learn about a new
#'   criterion. It changes when the run stops and nothing else: the drift is
#'   measured either way, and `attractor` is derived from it either way, so a
#'   wrapper that stops early still reports honestly what it stopped on.
#' @return A `MizerSim` or a `MizerParams`, in either case carrying the
#'   `"convergence"` attribute.
#' @noRd
project_until_settled <- function(params,
                                effort = params@initial_effort,
                                distance_func = distanceSSLogN,
                                t_check = 15 * dt,
                                t_max = 100,
                                dt = 0.1,
                                t_save = 1,
                                distance_tol = 0.1 * t_check,
                                residual_tol = steady_residual_tol(),
                                amplitude_tol = 0.01,
                                amp_rel_tol = 0.1,
                                extinction_threshold = 1e-6,
                                progress_bar = TRUE,
                                info_level = default_info_level(),
                                method = c("euler", "predictor_corrector",
                                           "tr_bdf2"), ...,
                                  return_sim = FALSE,
                                  require_steady = TRUE) {
    with_info_level(info_level = info_level, {
    # `tol` used to be the name of `distance_tol`. It would otherwise be passed
    # on to the distance function and silently ignored, which is the one way
    # this rename could cost someone a wrong answer rather than an error.
    if ("tol" %in% names(list(...))) {
        stop("`tol` has been renamed to `distance_tol`, because the run now ",
             "also has to meet a tolerance on the biomass drift, ",
             "`residual_tol`. See `?projectUntilSettled`.", call. = FALSE)
    }
    params <- validParams(params)
    method <- normalise_project_method(method)
    effort <- validEffortVector(effort, params = params)
    params@initial_effort <- effort
    assert_that(t_max >= t_check,
                distance_tol > 0,
                is.number(residual_tol), residual_tol > 0,
                amplitude_tol > 0,
                extinction_threshold >= 0)
    if ((t_check < dt) ||
        !isTRUE(all.equal((t_check - round(t_check / dt) * dt), 0))) {
        stop("t_check must be a positive multiple of dt")
    }
    # A run that stops before the first save would otherwise return a sim
    # holding only its starting state, so the save interval is capped at the
    # length of the run, exactly as in project().
    if (t_max < t_save) {
        t_save <- t_max
    }
    if ((t_save < dt) || !isTRUE(all.equal((t_save - round(t_save / dt) * dt), 0))) {
        stop("t_save must be a positive multiple of dt")
    }
    # `t_check` and `t_save` are independent grids: each has to land on the time
    # steps, but neither has to be a multiple of the other. The run advances one
    # `dt` at a time and tests both, so they need no common divisor beyond `dt`.
    steps_total <- round(t_max / dt)
    save_steps  <- round(t_save / dt)
    check_steps <- round(t_check / dt)
    t_dimnames <- seq(0, t_max, by = t_save)

    if (is(progress_bar, "Progress")) {
        # We have been passed a shiny progress object
        progress_bar$set(message = "Finding steady state", value = 0)
        proginc <- 1 / ceiling(t_max / t_check)
    }

    if (return_sim) {
        # create MizerSim object. The run stops at a `t_check` boundary, which
        # need not be on the `t_save` grid, so one extra slot is kept for the
        # state the run settled on; both the slot and any unused ones are
        # trimmed at the end.
        sim <- MizerSim(params, t_dimnames = c(t_dimnames, t_max + t_save))
        sim@sim_params <- list(method = method, dt = dt)
        sim@n[1, , ] <- params@initial_n
        sim@n_pp[1, ] <- params@initial_n_pp
        sim@n_other[1, ] <- params@initial_n_other
        sim@effort[1, ] <- params@initial_effort
        # Index of the last row written, so that the settled state can be
        # appended after the final save without overwriting it.
        sim_idx <- 1L
    }

    # get functions
    resource_dynamics_fn <- get(params@resource_dynamics)
    other_dynamics_fns <- lapply(params@other_dynamics, get)
    rates_fns <- projectRateFunctions(params)
    r <- rates_fns$Rates(
        params, n = params@initial_n,
        n_pp = params@initial_n_pp,
        n_other = params@initial_n_other,
        t = 0,
        effort = effort, rates_fns = rates_fns, ...)

    previous <- list(n = params@initial_n,
                     n_pp = params@initial_n_pp,
                     n_other = params@initial_n_other,
                     rates = r)

    # Reference reproduction rate for the relative extinction test below.
    rdd_start <- r$rdd

    # Record a cheap scalar summary (per-species biomass) after every time step,
    # so that a limit cycle can be detected and characterised at the finest
    # resolution the run has, whatever the period and whatever `t_check` is. It
    # costs one matrix product per step against a whole rate evaluation, and it
    # is why the cycle detection needs no tuning knob of its own.
    wdw <- params@w * params@dw
    bio_series <- matrix(NA_real_, nrow = steps_total + 1,
                         ncol = nrow(params@species_params))
    bio_series[1, ] <- as.numeric(previous$n %*% wdw)

    cycle <- NULL
    success <- FALSE
    extinct <- FALSE
    distance <- NA_real_
    residual <- NA_real_
    current <- previous
    # The run advances one time step at a time and tests two independent grids:
    # the state is stored every `save_steps` steps and the stopping criteria are
    # evaluated every `check_steps` steps. Stepping singly rather than in blocks
    # is what frees the two grids from having to be multiples of each other; it
    # gives bit-identical results, because each step picks up where the last one
    # left off, and costs a few percent. The time passed is the absolute time
    # since the start of the run, exactly as project() passes it, so that a rate
    # or component function that depends on `t` sees the same monotonically
    # increasing time here as it would there.
    for (k in seq_len(steps_total)) {
        current <- project_simple(params, n = current$n, n_pp = current$n_pp,
                                  n_other = current$n_other,
                                  t = (k - 1) * dt,
                                  dt = dt, steps = 1,
                                  effort = params@initial_effort,
                                  resource_dynamics_fn = resource_dynamics_fn,
                                  other_dynamics_fns = other_dynamics_fns,
                                  rates_fns = rates_fns,
                                  method = method)
        bio_series[k + 1, ] <- as.numeric(current$n %*% wdw)

        if (return_sim && k %% save_steps == 0) {
            # Store result
            sim_idx <- sim_idx + 1L
            sim@n[sim_idx, , ] <- current$n
            sim@n_pp[sim_idx, ] <- current$n_pp
            sim@n_other[sim_idx, ] <-
                unserialize(serialize(current$n_other, NULL))
            sim@effort[sim_idx, ] <- params@initial_effort
        }

        if (k %% check_steps != 0) next

        # advance shiny progress bar
        if (is(progress_bar, "Progress")) {
            progress_bar$inc(amount = proginc)
        }

        # A species whose reproduction has collapsed to a tiny fraction of its
        # starting value is going extinct, so stop.
        extinct <- is.na(current$rates$rdd) |
            current$rates$rdd <= extinction_threshold * rdd_start
        if (any(extinct)) {
            warning(paste(params@species_params$species[extinct], collapse = ", "),
                    " are going extinct.")
            success <- FALSE
            distance <- NA
            break
        }

        distance <- distance_func(params,
                                  current = current,
                                  previous = previous, ...)

        # Check for a limit cycle whether or not the distance criterion is met.
        # A cycle whose period divides `t_check` is sampled at the same phase by
        # every call to the distance function, so on that evidence alone it is
        # indistinguishable from a fixed point; the detection below works on the
        # biomass series sampled at every time step and is blind to `t_check`.
        cycle <- detect_limit_cycle(bio_series[seq_len(k + 1), , drop = FALSE],
                                    dt, amplitude_tol,
                                    amp_rel_tol = amp_rel_tol)
        if (!is.null(cycle)) {
            break
        }

        if (distance < distance_tol) {
            # The distance function has gone quiet. Whether that means a fixed
            # point has been reached is a different question, asked of the state
            # itself rather than of the distance function's scale. It is asked
            # here whatever `require_steady` says, because the answer is
            # reported either way; only whether it can hold the run open
            # depends on that argument.
            residual <- tryCatch(
                steady_biomass_drift(params, effort = effort,
                                     n = current$n, n_pp = current$n_pp,
                                     n_other = current$n_other),
                error = function(e) NA_real_)
            # An unmeasurable drift is no reason to keep running: the distance
            # criterion the user asked for has been met, and `residual` records
            # that the check could not be made.
            if (!require_steady || !is.finite(residual) ||
                residual <= residual_tol) {
                success <- TRUE
                break
            }
        }
        previous <- current
    }

    years <- k * dt

    params@initial_n[] <- current$n
    params@initial_n_pp[] <- current$n_pp
    params@initial_n_other[] <- current$n_other

    # How far the state that was reached actually is from a fixed point. The
    # distance function above only compares two states `t_check` apart, which
    # is a proxy; this is the thing itself. It is recorded even for a cycle or a
    # non-converged run, where it says how far off the run stopped. On a
    # successful run it has already been measured at this very state, as the
    # second half of the convergence criterion, so it is not measured twice.
    if (!success) {
        residual <- tryCatch(steady_biomass_drift(params, effort = effort),
                             error = function(e) NA_real_)
    }
    residual_txt <- if (is.finite(residual)) {
        paste0(" The biomasses change at up to ", signif(residual, 2),
               " per year.")
    } else {
        ""
    }

    if (!is.null(cycle)) {
        termination <- "cycle_detected"
        converged <- TRUE
        attractor <- "limit_cycle"
        signal_info("convergence", paste0(
            "Settled onto a limit cycle of period ",
            signif(cycle$period, 3), " years (relative amplitude ",
            signif(cycle$amplitude, 2), ") after ", years, " years."),
            unhandled = "show")
    } else if (success) {
        # The drift was measured at this very state, so the verdict is not
        # re-derived; it still goes through the same function, so that "fixed
        # point" cannot come to mean two things. Where the run stopped on the
        # distance criterion alone the termination says so, rather than naming a
        # tolerance that was not part of the test.
        attractor <- steady_attractor(residual, residual_tol)
        converged <- TRUE
        termination <- if (identical(attractor, "fixed_point")) {
            "residual_tolerance"
        } else {
            "distance_tolerance"
        }
        # The distance function being satisfied only says that the state stopped
        # moving on that function's own scale, so the message says that and no
        # more. The residual is reported every time, not only when it is large:
        # it is the evidence for whether this is a fixed point, and a number the
        # user only ever sees when something has gone wrong is a number they
        # never learn to look at. Only the instruction to act on it is
        # conditional. This stays an "info": the user asked for convergence at
        # their `distance_tol` and got it, so nothing they asked for failed to
        # happen.
        # The drift is only unmeasurable when `steady_biomass_drift()` itself
        # failed, in which case the run stopped on the distance criterion alone
        # and the user should know that half the test did not happen.
        caveat <- if (!is.finite(residual)) {
            paste0(" The biomass drift could not be measured, so only the ",
                   "distance criterion was checked.")
        } else if (residual > steady_residual_tol()) {
            paste0(" Reduce the tolerance on the distance function to ",
                   "converge further.")
        } else {
            ""
        }
        signal_info("convergence",
                    paste0("Reached the convergence tolerance after ", years,
                           " years.", residual_txt, caveat),
                    unhandled = "show")
    } else {
        termination <- if (any(extinct)) "extinction" else "time_limit"
        converged <- FALSE
        # Running out of time says nothing about where the run had got to. If
        # the state it stopped at is not drifting, it is a fixed point, however
        # unhappy the distance function was about it.
        attractor <- if (any(extinct)) {
            NA_character_
        } else {
            steady_attractor(residual, residual_tol)
        }
        # A run that has met the distance criterion but is still drifting is
        # the case the two-part test exists for. Saying so is the difference
        # between a puzzling long run and a diagnosis.
        why <- if (!any(extinct) && is.finite(distance) &&
                   distance < distance_tol && is.finite(residual) &&
                   residual > residual_tol) {
            paste0("The distance function returned ", signif(distance, 3),
                   ", which is below the distance tolerance, but the ",
                   "biomasses are still changing at up to ",
                   signif(residual, 2), " per year, which is above the ",
                   "residual tolerance, so this state is not a fixed point.")
        } else if (identical(attractor, "fixed_point")) {
            # The opposite case, and the one that reads as a contradiction
            # until it is spelled out: the distance function is still
            # unsatisfied while the state it is unsatisfied about has stopped
            # drifting. `converged` and `attractor` answer different questions,
            # so say which is which rather than leaving the two entries of the
            # convergence attribute to be compared.
            paste0("The distance function returned ", signif(distance, 3),
                   ", which is above the distance tolerance, but the state ",
                   "reached is a fixed point: the biomasses change at only ",
                   signif(residual, 2), " per year. The distance function is ",
                   "measuring motion that the biomasses do not see, so ",
                   "`attractor` is \"fixed_point\" even though `converged` is ",
                   "FALSE. Loosen `distance_tol`, or use a distance function ",
                   "that ignores what this one is still seeing.")
        } else {
            paste0("Value returned by the distance function was: ", distance)
        }
        signal_info("convergence", paste0(
            "Simulation run did not converge after ", years, " years. ", why),
            unhandled = "show")
    }

    convergence <- convergence_result(
        termination = termination,
        converged = converged,
        attractor = attractor,
        distance = distance,
        residual = residual,
        years = years,
        period = if (!is.null(cycle)) cycle$period else NA_real_,
        amplitude = if (!is.null(cycle)) cycle$amplitude else NA_real_,
        extinct = if (any(extinct)) {
            params@species_params$species[extinct]
        } else {
            character(0)
        }
    )

    if (return_sim) {
        sim@params <- params
        times <- t_dimnames[seq_len(sim_idx)]
        if (k %% save_steps != 0) {
            # The run stopped between two saves, so the state it settled on is
            # appended as the final time point. project() does the same when
            # `t_max` is not a multiple of `t_save`: the last interval is
            # shorter than the others rather than the final state being lost.
            sim_idx <- sim_idx + 1L
            sim@n[sim_idx, , ] <- current$n
            sim@n_pp[sim_idx, ] <- current$n_pp
            sim@n_other[sim_idx, ] <-
                unserialize(serialize(current$n_other, NULL))
            sim@effort[sim_idx, ] <- params@initial_effort
            times <- c(times, years)
        }
        sel <- seq_len(sim_idx)
        sim@n <- sim@n[sel, , , drop = FALSE]
        sim@n_pp <- sim@n_pp[sel, , drop = FALSE]
        sim@n_other <- sim@n_other[sel, , drop = FALSE]
        sim@effort <- sim@effort[sel, , drop = FALSE]
        dimnames(sim@n)$time <- times
        dimnames(sim@n_pp)$time <- times
        dimnames(sim@n_other)$time <- times
        dimnames(sim@effort)$time <- times
        attr(sim, "convergence") <- convergence
        return(sim)
    } else {
        params@time_modified <- lubridate::now()
        attr(params, "convergence") <- convergence
        return(params)
    }
    })
}

#' Detect a limit cycle from a fine-resolution biomass time series
#'
#' Used by [projectUntilSettled()] to decide whether a run that is not settling to a
#' fixed point has instead converged onto a limit cycle. Works on the per-species
#' biomass sampled at every time step, so it is independent of
#' whether the cycle period is commensurate with the check interval `t_check`.
#'
#' The community-total log-biomass is used as a scalar signal. Its
#' autocorrelation gives a candidate period (the first autocorrelation peak). The
#' oscillation is accepted as a settled limit cycle only if it has already
#' persisted for three periods: the relative amplitude of the three successive
#' period-windows must agree within `amp_rel_tol`, must exceed `amplitude_tol`,
#' and must show no net decay across the three periods. The no-decay condition is what
#' distinguishes a genuine limit cycle from a slowly-decaying spiral toward a
#' stable fixed point (which has a spectral radius just below 1). Discrimination
#' is necessarily imperfect when the spectral radius is extremely close to 1,
#' because such a spiral is indistinguishable from a cycle over any finite run.
#'
#' @param bio Numeric matrix of per-species biomass, one row per sample.
#' @param t_sample The sampling interval, so period `= lag * t_sample`.
#' @param amplitude_tol Minimum relative amplitude for an oscillation to count as
#'   a cycle.
#' @param acf_threshold Minimum autocorrelation at the candidate period.
#' @param amp_rel_tol Maximum relative change of amplitude between successive
#'   periods for the cycle to count as settled.
#' @return `NULL` if no settled cycle is detected, otherwise a list with the
#'   `period` (in the same time units as `t_sample`) and the relative
#'   `amplitude`.
#' @noRd
detect_limit_cycle <- function(bio, t_sample, amplitude_tol,
                               acf_threshold = 0.5,
                               amp_rel_tol = 0.1) {
    n <- nrow(bio)
    if (n < 20) return(NULL)
    
    window_length <- max(20, ceiling(n / 2))
    recent_bio <- bio[(n - window_length + 1):n, , drop = FALSE]
    
    s <- log(rowSums(recent_bio))
    s <- s - mean(s)
    if (all(abs(s) < .Machine$double.eps)) return(NULL)
    
    lag_max <- floor(window_length / 2)
    ac <- stats::acf(s, lag.max = lag_max, plot = FALSE, demean = TRUE)$acf[, 1, 1]
    w <- find_first_acf_peak(ac, acf_threshold)
    # Need three full periods of history to confirm a settled cycle.
    if (is.na(w) || w < 2 || n < 3 * w) return(NULL)
    amp_old <- amp_window(bio[(n - 3 * w + 1):(n - 2 * w), , drop = FALSE])
    amp_mid <- amp_window(bio[(n - 2 * w + 1):(n - w), , drop = FALSE])
    amp_new <- amp_window(bio[(n - w + 1):n, , drop = FALSE])
    if (amp_new <= amplitude_tol) return(NULL)
    # Successive periods must have matching amplitude ...
    if (abs(amp_new - amp_mid) / amp_new > amp_rel_tol) return(NULL)
    if (abs(amp_mid - amp_old) / amp_mid > amp_rel_tol) return(NULL)
    # ... and there must be no net decay (which would signal a decaying spiral).
    if (amp_new < amp_old * (1 - amp_rel_tol)) return(NULL)
    list(period = w * t_sample, amplitude = amp_new)
}

#' Largest per-species relative peak-to-trough amplitude in a window
#' @param bio Numeric matrix of per-species biomass over the window.
#' @return The maximum over species of `(max - min) / mean`.
#' @noRd
amp_window <- function(bio) {
    rng <- apply(bio, 2, function(x) {
        m <- mean(x)
        if (!is.finite(m) || m <= 0) return(0)
        (max(x) - min(x)) / m
    })
    max(rng)
}

#' First local maximum of an autocorrelation vector
#'
#' `ac[1]` is the lag-0 autocorrelation; a peak at vector position `k`
#' corresponds to lag `k - 1`.
#' @param ac Autocorrelation values indexed from lag 0.
#' @param threshold Minimum autocorrelation for a peak to count.
#' @return The lag of the first qualifying local maximum, or `NA`.
#' @noRd
find_first_acf_peak <- function(ac, threshold) {
    n <- length(ac)
    if (n < 3) return(NA_integer_)
    for (k in 2:(n - 1)) {
        if (ac[k] > ac[k - 1] && ac[k] >= ac[k + 1] && ac[k] > threshold) {
            return(k - 1L)
        }
    }
    NA_integer_
}



#' Helper function to keep other components constant
#'
#' @param params MizerParams object
#' @param n_other Abundances of other components
#' @param component Name of the component that is being updated
#' @param ... Unused
#' @return The current value of the component
#' @export
#' @concept helper
constant_other <- function(params, n_other, component, ...) {
    n_other[[component]]
}

#' Helper function to assure validity of species argument
#'
#' If the species argument contains invalid species, then these are
#' ignored but a warning is issued.
#'
#' @param object A MizerSim or MizerParams object from which the species
#'   should be selected.
#' @param species The species to be selected. Optional. By default all target
#'   species are selected. A vector of species names, or a
#'   numeric vector with the species indices, or a logical vector indicating for
#'   each species whether it is to be selected (TRUE) or not.
#' @param return.logical Whether the return value should be a logical vector.
#'   Default FALSE.
#' @param error_on_empty Whether to throw an error if there are zero valid
#'   species. Default FALSE.
#'
#' @return A vector of species names, in the same order as specified in the
#'   'species' argument. If 'return.logical = TRUE' then a logical vector is
#'   returned instead, with length equal to the number of species, with
#'   TRUE entry for each selected species.
#' @export
#' @concept helper
valid_species_arg <- function(object, species = NULL, return.logical = FALSE,
                              error_on_empty = FALSE) {
    if (is(object, "MizerSim")) {
        params <- object@params
    } else if (is(object, "MizerParams")) {
        params <- object
    } else {
        stop("The first argument must be a MizerSim or MizerParams object.")
    }
    assert_that(is.flag(return.logical),
                is.flag(error_on_empty))
    all_species <- dimnames(params@initial_n)$sp
    no_sp <- nrow(params@species_params)
    # Set species if missing to list of all non-background species
    if (is.null(species)) {
        species <- dimnames(params@initial_n)$sp[!params@species_params$is_background]
        if (length(species) == 0) {  # There are no non-background species.
            if (error_on_empty) {
                stop("No species have been selected.")
            }
            if (return.logical) {
                return(rep(FALSE, no_sp))
            } else {
                vector("character")
            }
        }
    }
    if (is.logical(species)) {
        if (length(species) != no_sp) {
            stop("The boolean `species` argument has the wrong length.")
        }
        if (!any(species) && error_on_empty) {
            stop("No species have been selected.")
        }
        if (return.logical) {
            return(species)
        }
        return(all_species[species])
    }
    if (is.numeric(species)) {
        if (!all(species %in% (1:no_sp))) {
            warning("A numeric 'species' argument should only contain the ",
                    "integers 1 to ", no_sp, ".")
        }
        species.logical <- 1:no_sp %in% species
        if (!any(species.logical) && error_on_empty) {
            stop("No species have been selected.")
        }
        if (return.logical) {
            return(species.logical)
        }
        return(all_species[species[species %in% (1:no_sp)]])
    }
    invalid <- setdiff(species, all_species)
    if (length(invalid) > 0) {
        warning("The following species do not exist: ",
                toString(invalid), ".")
    }
    species <- intersect(species, all_species)
    if (length(species) == 0 && error_on_empty) {
        stop("No species have been selected.")
    }
    if (return.logical) {
        return(all_species %in% species)
    }
    species
}

#' Helper function to assure validity of gears argument
#'
#' If the gears argument contains invalid gears, then these are
#' ignored but a warning is issued.
#'
#' @param object A MizerSim or MizerParams object from which the gears
#'   should be selected.
#' @param gears The gears to be selected. Optional. By default all gears
#'   are selected. A vector of gear names.
#' @param error_on_empty Whether to throw an error if there are zero valid
#'   gears. Default FALSE.
#'
#' @return A vector of gear names in the same order as supplied in `gears`,
#'   with invalid names removed. If `gears` is `NULL`, all gears are returned
#'   in the order stored in the model.
#' @export
#' @concept helper
valid_gears_arg <- function(object, gears = NULL,
                            error_on_empty = FALSE) {
    if (is(object, "MizerSim")) {
        params <- object@params
    } else if (is(object, "MizerParams")) {
        params <- object
    } else {
        stop("The first argument must be a MizerSim or MizerParams object.")
    }
    assert_that(is.flag(error_on_empty))
    all_gears <- unique(params@gear_params$gear)
    no_gear <- length(all_gears)
    # Set gears if missing to list of all gears
    if (is.null(gears)) {
        gears <- all_gears
    }
    invalid <- setdiff(gears, all_gears)
    if (length(invalid) > 0) {
        warning("The following gears do not exist: ",
                toString(invalid), ".")
    }
    gears <- intersect(gears, all_gears)
    if (length(gears) == 0 && error_on_empty) {
        stop("No gears have been selected.")
    }
    gears
}
