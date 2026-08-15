#' Measure distance between current and previous state in terms of RDI
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' This function can be used in [projectToSteady()] to decide when sufficient
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
#' state. This function can be used in [projectToSteady()] to decide when
#' sufficient convergence to steady state has been achieved.
#'
#' @param params MizerParams
#' @param current A named list with entries `n`, `n_pp` and `n_other`
#'   describing the current state
#' @param previous A named list with entries `n`, `n_pp` and `n_other`
#'   describing the previous state
#' @return The sum of squares of the difference in the logs of the (nonzero)
#'   fish abundances `n`, ignoring entries where either state has zero
#'   abundance:
#'   `sum((log(current$n) - log(previous$n))^2)`
#' @family distance functions
#' @concept helper
#' @export
distanceSSLogN <- function(params, current, previous) {
    UseMethod("distanceSSLogN")
}
#' @export
distanceSSLogN.MizerParams <- function(params, current, previous) {
    sel <- current$n > 0 & previous$n > 0
    sum((log(current$n[sel]) - log(previous$n[sel]))^2)
}

#' Project to steady state
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Run the full dynamics, as in [project()], but stop once the change has slowed
#' down sufficiently, in the sense that the distance between states at
#' successive time steps is less than `tol`. You determine how the distance is
#' calculated.
#'
#' @details
#' # How the run is organised
#'
#' The simulation is not run in one go but is broken into blocks of `t_per`
#' years. Within a block the dynamics are advanced with time step `dt` exactly
#' as in [project()]. At the end of each block the function decides whether to
#' stop, so `t_per` sets how often the stopping criteria are evaluated and also
#' the interval over which change is measured. The run ends at the latest after
#' `t_max` years, i.e. after `floor(t_max / t_per)` blocks.
#'
#' Independently of the blocks, a cheap scalar summary of the state (the biomass
#' of each species) is recorded every `t_save` years. This finely resolved
#' series is what the limit-cycle detection works on, so that a cycle can be
#' found and its period measured even when that period bears no simple relation
#' to `t_per`.
#'
#' At the end of each block the following checks are made, in this order.
#'
#' ## 1. Extinction
#'
#' If the reproduction rate (RDD) of any species has fallen below
#' `extinction_threshold` times its value at the start of the run, or has become
#' `NA`, that species is deemed to be on its way to extinction. A warning naming
#' the affected species is issued and the run stops with
#' `type = "extinction"`. Because the criterion is relative to the initial
#' reproduction, a species that starts with zero reproduction is flagged
#' immediately, whereas in [steady()], where reproduction is held constant, a
#' healthy species is never flagged.
#'
#' ## 2. Convergence to a fixed point
#'
#' `distance_func` is called with the state at the end of the previous block and
#' the state at the end of the current block, i.e. with two states `t_per` years
#' apart. If the number it returns is less than `tol`, the run stops with
#' `type = "steady"`. What "distance" means is entirely up to that function:
#' the default [distanceSSLogN()] uses the sum of squared changes in log
#' abundance, while [steady()] instead passes [distanceMaxRelRDI()], which uses
#' the largest relative change in egg production.
#'
#' Note that this test only compares two states one `t_per` apart; it cannot by
#' itself tell a fixed point from a cycle whose period divides `t_per`. This is
#' why `t_per` should be chosen as an *odd* multiple of `dt`: a period-2 cycle
#' (period `2 * dt`), which is the most common numerical oscillation, would
#' otherwise be sampled at the same phase in every block and would look
#' perfectly converged.
#'
#' ## 3. Limit cycle
#'
#' If the state is still changing more than `tol`, the recorded biomass series
#' is examined to see whether the run has instead settled onto a limit cycle. If
#' it has, the run stops with `type = "cycle"` and the period and amplitude of
#' the cycle are reported.
#'
#' If none of the three checks fires before `t_max` is reached, the run stops
#' with `type = "not_converged"`. In every case the outcome is recorded in the
#' `"convergence"` attribute of the returned object, described under *Value*
#' below.
#'
#' # How a limit cycle is detected
#'
#' The detection uses the community-total biomass, on a log scale and with its
#' mean removed, as a scalar signal, sampled every `t_save` years. At least 20
#' samples are needed before any cycle can be reported.
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
#'    with each other to within 10%, and the most recent amplitude must not be
#'    smaller than the oldest by more than 10%.
#'
#' The last condition is what distinguishes a genuine limit cycle from a slowly
#' decaying spiral towards a stable fixed point: the spiral loses amplitude from
#' one period to the next, the cycle does not. The distinction is necessarily
#' imperfect when the decay is extremely slow, because over any finite run such
#' a spiral is indistinguishable from a cycle. If you need a definitive answer,
#' use [getStability()] on the fixed point found by [steadyNewton()], which
#' works out the spectral radius of the linearised dynamics instead of watching
#' a trajectory.
#'
#' The reported `period` is a multiple of `t_save`, so it is only resolved to
#' that accuracy; reduce `t_save` if you need the period more precisely.
#'
#' @inheritParams steady
#' @param effort The fishing effort to be used throughout the simulation.
#'   This is validated by [validEffortVector()] and can therefore be `NULL`, a
#'   single numeric value used for all gears, an unnamed numeric vector with one
#'   entry per gear, or a named numeric vector for some or all gears.
#' @param distance_func A function that will be called after every `t_per` years
#'   with both the previous and the new state and that should return a number
#'   that in some sense measures the distance between the states. By default
#'   this uses the function [distanceSSLogN()] that you can use as a model for your
#'   own distance function.
#' @param tol The run stops when the number returned by `distance_func` for two
#'   states `t_per` years apart drops below `tol`. Its meaning therefore depends
#'   on the distance function you supply.
#' @param info_level Controls the amount of information messages that are shown.
#'   Higher levels lead to more messages, `info_level = 0` gives silence. The
#'   default is taken from the `mizer_info_level` option, see
#'   [default_info_level()].
#' @param method The numerical method to use for the consumer density update.
#'   See [project()].
#' @param ... Further arguments will be passed on to your distance function.
#'
#' @return If `return_sim = FALSE`, a `MizerParams` object with the initial
#'   state replaced by the final state found by the steady-state search. If
#'   `return_sim = TRUE`, a `MizerSim` object containing the intermediate states
#'   saved every `t_per` years.
#'
#'   In either case the returned object carries an attribute `"convergence"`
#'   describing the solution the run settled on, a named list with entries:
#'   \describe{
#'     \item{`type`}{One of `"steady"` (a stable fixed point), `"cycle"` (a
#'       limit cycle), `"not_converged"` (still changing at `t_max`) or
#'       `"extinction"` (a species died out).}
#'     \item{`converged`}{Logical, `TRUE` for `"steady"` and `"cycle"`.}
#'     \item{`distance`}{The final value returned by `distance_func`.}
#'     \item{`residual`}{The largest per-capita rate of change, in 1/year, at the
#'       state that was reached, as returned by [getSteadyResidual()]. Unlike
#'       `distance`, which compares two states `t_per` apart on whatever scale
#'       the distance function uses, this measures how far the state actually is
#'       from being a fixed point.}
#'     \item{`years`}{The number of years simulated.}
#'     \item{`period`}{For a limit cycle, its period in years; otherwise `NA`.}
#'     \item{`amplitude`}{For a limit cycle, the largest per-species relative
#'       peak-to-trough biomass amplitude; otherwise `NA`.}
#'   }
#'   This mirrors how [steadyNewton()] attaches an `"stability"` attribute.
#' @seealso [distanceSSLogN()], [distanceMaxRelRDI()], [steadyNewton()],
#'   [getStability()]
#' @export
projectToSteady <- function(params,
                            effort = params@initial_effort,
                            distance_func = distanceSSLogN,
                            t_per = 1.5,
                            t_max = 100,
                            dt = 0.1,
                            t_save = dt,
                            tol = 0.1 * t_per,
                            amplitude_tol = 0.01,
                            extinction_threshold = 1e-6,
                            return_sim = FALSE,
                            progress_bar = TRUE,
                            info_level = default_info_level(),
                            method = c("euler", "predictor_corrector", "tr_bdf2"), ...) {
    UseMethod("projectToSteady")
}
#' @export
projectToSteady.MizerParams <- function(params,
                            effort = params@initial_effort,
                            distance_func = distanceSSLogN,
                            t_per = 1.5,
                            t_max = 100,
                            dt = 0.1,
                            t_save = dt,
                            tol = 0.1 * t_per,
                            amplitude_tol = 0.01,
                            extinction_threshold = 1e-6,
                            return_sim = FALSE,
                            progress_bar = TRUE,
                            info_level = default_info_level(),
                            method = c("euler", "predictor_corrector", "tr_bdf2"), ...) {
    with_info_level(info_level = info_level, {
    params <- validParams(params)
    method <- normalise_project_method(method)
    effort <- validEffortVector(effort, params = params)
    params@initial_effort <- effort
    assert_that(t_max >= t_per,
                tol > 0,
                amplitude_tol > 0,
                extinction_threshold >= 0)
    if ((t_per < dt) || !isTRUE(all.equal((t_per - round(t_per / dt) * dt), 0))) {
        stop("t_per must be a positive multiple of dt")
    }
    if ((t_save < dt) || !isTRUE(all.equal((t_save - round(t_save / dt) * dt), 0))) {
        stop("t_save must be a positive multiple of dt")
    }
    if ((t_per < t_save) ||
        !isTRUE(all.equal((t_per - round(t_per / t_save) * t_save), 0))) {
        stop("t_per must be a positive multiple of t_save")
    }
    t_dimnames <-  seq(0, t_max, by = t_per)

    if (is(progress_bar, "Progress")) {
        # We have been passed a shiny progress object
        progress_bar$set(message = "Finding steady state", value = 0)
        proginc <- 1 / ceiling(t_max/t_per)
    }

    if (return_sim) {
        # create MizerSim object
        sim <- MizerSim(params, t_dimnames =  t_dimnames)
        sim@sim_params <- list(method = method, dt = dt)
        sim@n[1, , ] <- params@initial_n
        sim@n_pp[1, ] <- params@initial_n_pp
        sim@n_other[1, ] <- params@initial_n_other
        sim@effort[1, ] <- params@initial_effort
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

    # Record a cheap scalar summary (per-species biomass) at the fine `t_save`
    # resolution so that a limit cycle can be detected and characterised even
    # when its period is incommensurate with `t_per`.
    wdw <- params@w * params@dw
    steps_per_save  <- round(t_save / dt)
    saves_per_block <- round(t_per / t_save)
    max_saves <- (length(t_dimnames) - 1) * saves_per_block + 1
    bio_series <- matrix(NA_real_, nrow = max_saves,
                         ncol = nrow(params@species_params))
    bio_series[1, ] <- as.numeric(previous$n %*% wdw)
    save_idx <- 1L

    cycle <- NULL
    success <- FALSE
    extinct <- FALSE
    distance <- NA_real_
    for (i in 2:length(t_dimnames)) {
        # advance shiny progress bar
        if (is(progress_bar, "Progress")) {
            progress_bar$inc(amount = proginc)
        }
        # Step the block in `t_save`-sized sub-steps, recording the biomass
        # summary after each. Passing the within-block time keeps the result
        # bit-identical to a single project_simple() call over the whole block.
        current <- previous
        for (s in seq_len(saves_per_block)) {
            current <- project_simple(params, n = current$n, n_pp = current$n_pp,
                                      n_other = current$n_other,
                                      t = (s - 1) * t_save,
                                      dt = dt, steps = steps_per_save,
                                      effort = params@initial_effort,
                                      resource_dynamics_fn = resource_dynamics_fn,
                                      other_dynamics_fns = other_dynamics_fns,
                                      rates_fns = rates_fns,
                                      method = method)
            save_idx <- save_idx + 1L
            bio_series[save_idx, ] <- as.numeric(current$n %*% wdw)
        }
        if (return_sim) {
            # Store result
            sim@n[i, , ] <- current$n
            sim@n_pp[i, ] <- current$n_pp
            sim@n_other[i, ] <- unserialize(serialize(current$n_other, NULL))
            sim@effort[i, ] <- params@initial_effort
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
        success <- distance < tol
        if (success == TRUE) {
            break
        }
        # Not settling to a fixed point: check whether we are on a limit cycle.
        cycle <- detect_limit_cycle(bio_series[seq_len(save_idx), , drop = FALSE],
                                    t_save, amplitude_tol)
        if (!is.null(cycle)) {
            break
        }
        previous <- current
    }

    years <- (i - 1) * t_per
    converged <- FALSE

    params@initial_n[] <- current$n
    params@initial_n_pp[] <- current$n_pp
    params@initial_n_other[] <- current$n_other

    # How far the state that was reached actually is from a fixed point. The
    # distance function above only compares two states `t_per` apart, which is a
    # proxy; this is the thing itself. It is recorded even for a cycle or a
    # non-converged run, where it says how far off the run stopped.
    residual <- tryCatch(steady_biomass_drift(params, effort = effort),
                         error = function(e) NA_real_)
    residual_txt <- if (is.finite(residual)) {
        paste0(" A biomass is still changing at up to ",
               signif(residual, 2), " per year.")
    } else {
        ""
    }

    if (!is.null(cycle)) {
        type <- "cycle"
        converged <- TRUE
        signal_info("convergence", paste0(
            "Converged to a limit cycle of period ",
            signif(cycle$period, 3), " years (relative amplitude ",
            signif(cycle$amplitude, 2), ") after ", years, " years."),
            unhandled = "show")
    } else if (success) {
        type <- "steady"
        converged <- TRUE
        # The distance function being satisfied only says that the state stopped
        # moving on that function's own scale. Report the residual alongside it,
        # and flag the case where the two disagree — a run declared converged
        # whose biomasses are still visibly moving, which would otherwise be
        # invisible. This stays an "info": the user asked for convergence at
        # their `tol` and got it, so nothing they asked for failed to happen.
        caveat <- if (is.finite(residual) && residual > steady_residual_tol()) {
            paste0(residual_txt, " Reduce `tol` to converge further.")
        } else {
            ""
        }
        signal_info("convergence",
                    paste0("Convergence was achieved in ", years, " years.",
                           caveat),
                    unhandled = "show")
    } else {
        type <- if (any(extinct)) "extinction" else "not_converged"
        signal_info("convergence", paste0(
            "Simulation run did not converge after ", years,
            " years. Value returned by the distance function was: ", distance),
            unhandled = "show")
    }

    convergence <- list(
        type = type,
        converged = converged,
        distance = distance,
        residual = residual,
        years = years,
        period = if (type == "cycle") cycle$period else NA_real_,
        amplitude = if (type == "cycle") cycle$amplitude else NA_real_
    )

    if (return_sim) {
        sim@params <- params
        sel <- 1:i
        sim@n <- sim@n[sel, , , drop = FALSE]
        sim@n_pp <- sim@n_pp[sel, , drop = FALSE]
        sim@n_other <- sim@n_other[sel, , drop = FALSE]
        sim@effort <- sim@effort[sel, , drop = FALSE]
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
#' Used by [projectToSteady()] to decide whether a run that is not settling to a
#' fixed point has instead converged onto a limit cycle. Works on the per-species
#' biomass sampled at the fine `t_save` resolution, so it is independent of
#' whether the cycle period is commensurate with the block length `t_per`.
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
#' @param bio Numeric matrix of per-species biomass, one row per saved time step.
#' @param t_save The sampling interval, so period `= lag * t_save`.
#' @param amplitude_tol Minimum relative amplitude for an oscillation to count as
#'   a cycle.
#' @param acf_threshold Minimum autocorrelation at the candidate period.
#' @param amp_rel_tol Maximum relative change of amplitude between successive
#'   periods for the cycle to count as settled.
#' @return `NULL` if no settled cycle is detected, otherwise a list with the
#'   `period` (in the same time units as `t_save`) and the relative `amplitude`.
#' @noRd
detect_limit_cycle <- function(bio, t_save, amplitude_tol,
                               acf_threshold = 0.5,
                               amp_rel_tol = 0.1) {
    n <- nrow(bio)
    if (n < 20) return(NULL)
    s <- log(rowSums(bio))
    s <- s - mean(s)
    if (all(abs(s) < .Machine$double.eps)) return(NULL)
    lag_max <- floor(n / 2)
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
    list(period = w * t_save, amplitude = amp_new)
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

#' Set initial values to a steady state for the model
#'
#' The steady state is found by running the dynamics while keeping reproduction,
#' resource and other components constant until the size spectra no longer
#' change much (or until time `t_max` is reached, if earlier).
#'
#' If the model use Beverton-Holt reproduction then the reproduction parameters
#' are set to values that give the level of reproduction observed in that
#' steady state. The `preserve` argument can be used to specify which of the
#' reproduction parameters should be preserved.
#'
#' The resource is rebalanced so that the returned state is a steady state of
#' the resource as well as of the consumers. The resource abundance is
#' preserved while the capacity is recomputed to keep it steady. If the resource
#' capacity had been set manually (frozen), it is rebalanced and thereby
#' unfrozen.
#'
#' @param params A \linkS4class{MizerParams} object
#' @param t_max The maximum number of years to run the simulation. Default is 100.
#' @param t_per The simulation is broken up into shorter runs of `t_per` years,
#'   after each of which we check for convergence. Default value is 1.5. This
#'   should be chosen as an odd multiple of the timestep `dt` in order to be
#'   able to detect period 2 cycles.
#' @param dt The time step to use in `project()`.
#' @param t_save The interval at which a cheap per-species biomass summary is
#'   recorded for limit-cycle detection. Must be a positive multiple of `dt` and
#'   a divisor of `t_per`. Smaller values resolve the cycle period more finely at
#'   a small extra cost. Default is `dt`.
#' @param tol The simulation stops when the relative change in the egg
#'   production RDI over `t_per` years is less than `tol` for every species.
#' @param amplitude_tol `r lifecycle::badge("experimental")`
#'   The minimum relative biomass amplitude for a persistent oscillation to be
#'   reported as a limit cycle rather than treated as an (effectively steady)
#'   fixed point. This is a fraction of mean biomass and is kept separate from
#'   `tol` (which measures convergence to a fixed point on a different scale).
#'   Default `0.01`.
#' @param extinction_threshold `r lifecycle::badge("experimental")`
#'   A species is treated as going extinct, stopping the run, once its
#'   reproduction rate (RDD) falls below this fraction of its value at the start
#'   of the run. For example the default `1e-6` treats a species as extinct once
#'   its reproduction has collapsed to a millionth of its initial level. Because
#'   it is relative to the initial reproduction, a species that starts with zero
#'   reproduction is flagged immediately, and (in [steady()], where reproduction
#'   is held constant) a healthy species is never flagged.
#' @param return_sim If TRUE, the function returns the MizerSim object holding
#'   the result of the simulation run, saved at intervals of `t_per`. If FALSE (default) the function returns
#'   a MizerParams object with the "initial" slots set to the steady state.
#' @param preserve `r lifecycle::badge("experimental")`
#'   Specifies whether the `reproduction_level` should be preserved (default)
#'   or the maximum reproduction rate `R_max` or the reproductive
#'   efficiency `erepro`. See [setBevertonHolt()] for an explanation
#'   of the `reproduction_level`.
#' @param progress_bar A shiny progress object to implement a progress bar in a
#'   shiny app. Default FALSE.
#' @param info_level Controls the amount of information messages that are shown.
#'   Higher levels lead to more messages, `info_level = 0` gives silence. The
#'   default is taken from the `mizer_info_level` option, see
#'   [default_info_level()].
#' @param method The numerical method to use for the consumer density update.
#'   See [project()].
#' @return If `return_sim = FALSE`, a `MizerParams` object with the initial
#'   state replaced by the steady state. If `return_sim = TRUE`, a `MizerSim`
#'   object containing the intermediate states saved every `t_per` years. The
#'   returned object carries a `"convergence"` attribute describing the solution
#'   found (steady state, limit cycle, or non-convergence); see
#'   [projectToSteady()].
#' @export
#' @examples
#' \donttest{
#' params <- newTraitParams()
#' species_params(params)$gamma[5] <- 3000
#' params <- steady(params)
#' plotSpectra(params)
#' }
steady <- function(params, t_max = 100, t_per = 1.5, dt = 0.1, t_save = dt,
                   tol = 0.1 * dt, amplitude_tol = 0.01,
                   extinction_threshold = 1e-6, return_sim = FALSE,
                   preserve = c("reproduction_level", "erepro", "R_max"),
                   progress_bar = TRUE,
                   info_level = default_info_level(),
                   method = c("euler", "predictor_corrector", "tr_bdf2")) {
    UseMethod("steady")
}

#' @export
steady.MizerParams <- function(params, t_max = 100, t_per = 1.5, dt = 0.1,
                   t_save = dt,
                   tol = 0.1 * dt, amplitude_tol = 0.01,
                   extinction_threshold = 1e-6, return_sim = FALSE,
                   preserve = c("reproduction_level", "erepro", "R_max"),
                   progress_bar = TRUE,
                   info_level = default_info_level(),
                   method = c("euler", "predictor_corrector", "tr_bdf2")) {
    with_info_level(info_level = info_level, {
    method <- normalise_project_method(method)

    if (params@rates_funcs$RDD == "BevertonHoltRDD") {
        preserve <- match.arg(preserve)
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

    object <- projectToSteady(params,
                              distance_func = distanceMaxRelRDI,
                              t_per = t_per,
                              t_max = t_max,
                              dt = dt,
                              t_save = t_save,
                              tol = tol,
                              amplitude_tol = amplitude_tol,
                              extinction_threshold = extinction_threshold,
                              return_sim = return_sim,
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

    # Restore the original resource dynamics and rebalance so that the
    # converged state is a genuine steady state of the resource as well as of
    # the consumers. Balancing derives the resource capacity from the (preserved)
    # rate; a manually set (frozen) capacity would otherwise block this, leaving
    # the resource off its fixed point. We therefore clear the frozen mark on the
    # capacity first, so that it is rebalanced from the rate exactly as earlier
    # versions of mizer did when they rebalanced the resource in `steady()`.
    comment(params@cc_pp) <- NULL
    params <- setResource(params, resource_dynamics = old_resource_dynamics)

    if (params@rates_funcs$RDD == "BevertonHoltRDD") {
        if (preserve == "reproduction_level") {
            reproduction_level(params) <- old_reproduction_level
        } else if (preserve == "R_max") {
            params <- setBevertonHolt(params,
                                      R_max = old_R_max)
        } else {
            params <- setBevertonHolt(params, erepro = old_erepro)
        }
    }

    # The residual recorded by projectToSteady() was measured on the model it
    # was handed, which had reproduction, the resource and the other components
    # pinned. Restoring the real dynamics above changes the model, so the
    # residual is measured again on what is actually being returned.
    conv$residual <- tryCatch(steady_biomass_drift(params),
                              error = function(e) NA_real_)

    if (return_sim) {
        object@params <- params
        attr(object, "convergence") <- conv
        return(object)
    } else {
        params@time_modified <- lubridate::now()
        attr(params, "convergence") <- conv
        return(params)
    }
    })
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
