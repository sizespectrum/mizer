# Steady-state residual diagnostic -------------------------------------------
#
# `steady()` and `steadyNewton()` put a model *onto* its steady state. This file
# provides the diagnostic that asks whether a model *is* on it, so that the
# instruction "re-run steady() after any match.../calibrate... step" can be
# checked rather than remembered.

#' Default tolerance for calling a model steady
#'
#' The relative rate of biomass change, in units of 1/year, below which mizer
#' treats a model as being at its steady state. It has to sit above what a
#' converged [steady()] actually leaves behind and below the drift that the
#' mistake it is meant to catch produces. Measured on the North Sea model with
#' `steady_biomass_drift()`:
#'
#' | State | Drift (1/year) |
#' |---|---|
#' | `steadyNewton()` | 4e-7 |
#' | `steady(tol = 1e-4)` | 2e-5 |
#' | `steady()` at its default `tol` | 3e-3 |
#' | as shipped in `NS_params` | 1e-2 |
#' | `newTraitParams()`, freshly built | 0.14 |
#' | after `matchBiomasses()` | 2 |
#' | after `matchGrowth()` | 4 |
#'
#' `0.05` therefore sits more than an order of magnitude above a default
#' `steady()` run and one to three orders below a model knocked off its steady
#' state. Note that [steady()] stops on the relative change in egg
#' production rather than on this drift, so how close it gets is governed by its
#' own `tol` argument; tighten that if you need to settle further.
#'
#' Everything that judges steadiness — the `summary()` line, the [project()]
#' check, the guards in [getStability()] and [getLimitCycleSim()] — goes through
#' this function, so the tolerance has a single definition.
#'
#' @return The tolerance, a single number.
#' @noRd
steady_residual_tol <- function() {
    0.05
}

#' How far a model is from its steady state
#'
#' `r lifecycle::badge("experimental")`
#' Returns the rate at which the abundances would change if the model were
#' projected forward from its current initial state, relative to those
#' abundances. At a steady state this is zero, so it answers the question that
#' every calibration workflow otherwise has to remember to ask: *is this model
#' still at its steady state?*
#'
#' The value is a **per-capita rate of change, in units of 1/year**:
#' \deqn{R_i(w) = \frac{1}{N_i(w)}\frac{dN_i(w)}{dt}.}
#' A value of `1e-8` means nothing is moving. A value of `0.05` means that size
#' class would change by about 5% over the first year of a projection, and
#' `-0.05` that it would shrink by about that much. The sign is therefore the
#' direction the model would drift.
#'
#' For the consumers this is exact, not a finite-difference approximation: the
#' backward-Euler transport coefficients used by [project()] satisfy
#' \eqn{A N - S = -dt\,dN/dt} identically, so evaluating them at `dt = 1` gives
#' the instantaneous rate with no time-discretisation error. The resource and
#' other components have arbitrary user-supplied dynamics functions, so their
#' rates are obtained by taking one short step of length `dt`, accurate to
#' \eqn{O(dt)}.
#'
#' Everything is evaluated at the model's own stored state — `initialN()`,
#' `initialNResource()`, `initialNOther()` — using the model's own reproduction
#' function and its own `resource_dynamics`. Nothing is substituted or held
#' fixed. The number therefore answers exactly "if I called [project()] now,
#' would anything move?", which is why it works for every model rather than only
#' for the semichemostat resource that [steadyNewton()] requires.
#'
#' ## Reading the result
#'
#' The returned array is an [ArraySpeciesBySize] object, so it prints,
#' summarises and plots itself:
#'
#' ```r
#' res <- getSteadyResidual(params)
#' summary(res)                  # per-species minimum, mean and maximum
#' plot(res)                     # which species, and at which sizes
#' ```
#'
#' The plot is the diagnostic one: a model that is off steady state is usually
#' off in one species, or one part of the size range, and the plot says which.
#'
#' Size classes with no fish in them carry no information about steadiness — the
#' relative rate of change of a zero density is undefined — so they are returned
#' as `NA`. Use `na.rm = TRUE` in any summary, as the examples above do.
#'
#' ## Do not reduce this to its maximum
#'
#' `max(abs(res))` is a tempting single-number verdict and a misleading one. The
#' per-capita rate of a single size class is dominated by the fastest-relaxing
#' cells, and near the egg size those turn over in hours: a model settled for
#' every practical purpose can carry a cell rate of \eqn{10^4}/year there while
#' nothing observable moves. Under the second-order scheme (see
#' [second_order_w()]) this is severe enough to reverse the ordering between a
#' converged model and one that has just been knocked off its steady state.
#'
#' What mizer's own checks — the `summary()` line, and
#' `project(check_steady = TRUE)` — judge instead is the relative rate of change
#' of each species' *biomass*, which weights each size class by the mass it
#' holds, and is the drift the user would actually see in [plotBiomass()]. Use
#' this array to find out *where* a model is unsteady, and those checks to find
#' out *whether* it is.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param effort The fishing effort at which to evaluate the residual. By
#'   default the initial effort stored in `params`, which is the effort the
#'   model's steady state belongs to.
#' @param dt The step length used for the resource and other components, whose
#'   dynamics functions are only available as one-step maps. Smaller is more
#'   accurate. Not used for the consumers, whose rate is exact.
#' @return An [ArraySpeciesBySize] object (species x size) of per-capita rates
#'   of change in 1/year, `NA` where the density is zero. It carries two
#'   further attributes:
#'   \describe{
#'     \item{`resource`}{The per-capita rate of change of the resource, a
#'       numeric vector over `w_full`, `NA` where the resource density is zero.}
#'     \item{`other`}{A named list with one entry per other component, holding
#'       its per-capita rate of change, or `NA` for a component whose state is
#'       not numeric.}
#'   }
#' @seealso [isSteady()], [steady()], [steadyNewton()], [getStability()]
#' @export
#' @family summary functions
#' @concept summary_function
#' @examples
#' summary(getSteadyResidual(NS_params))
#' \donttest{
#' # Matching biomasses moves the model off its steady state, and the plot
#' # shows which species and which sizes have moved.
#' params <- NS_params
#' species_params(params)$biomass_observed <-
#'     c(0.8, 61, 12, 35, 1.6, 20, 10, 7.6, 135, 60, 30, 78)
#' species_params(params)$biomass_cutoff <- 10
#' params <- calibrateBiomass(params)
#' params <- matchBiomasses(params)
#' plot(getSteadyResidual(params))
#' }
getSteadyResidual <- function(params, effort = params@initial_effort,
                              dt = 1e-4) {
    rates <- steady_rates(params, effort = effort, dt = dt)

    n <- rates$n
    residual <- rates$dNdt / n
    residual[n == 0] <- NA_real_
    resource <- rates$dn_pp_dt / rates$n_pp
    resource[rates$n_pp == 0] <- NA_real_

    residual <- ArraySpeciesBySize(residual,
                                   value_name = "Steady-state residual",
                                   units = "1/year", params = rates$params)
    attr(residual, "resource") <- resource
    attr(residual, "other") <- rates$other
    residual
}

#' Check whether a model is at steady state
#'
#' `r lifecycle::badge("experimental")`
#' Returns `TRUE` if the model is at its steady state (within a specified
#' tolerance), `FALSE` otherwise.
#'
#' Steadiness is judged by computing the relative rate of change of biomass
#' across all consumer species, resource, and other components (see
#' [getSteadyResidual()]). If the largest biomass drift is less than or equal to
#' `tol`, the model is considered to be at steady state.
#'
#' @param params A \linkS4class{MizerParams} object or an extension thereof.
#' @param tol Tolerance for the relative rate of biomass change in 1/year.
#'   Defaults to `0.05` (5% change per year).
#' @param effort The fishing effort at which to evaluate steadiness. By default
#'   the initial effort stored in `params`.
#' @param ... Additional arguments passed to methods.
#' @return `TRUE` if the model's biomass drift is within `tol`, `FALSE`
#'   otherwise.
#' @seealso [getSteadyResidual()], [steady()], [steadyNewton()]
#' @export
#' @examples
#' isSteady(NS_params)
#'
#' \donttest{
#' # Moving a species abundance off its steady state makes isSteady() FALSE
#' params <- NS_params
#' initialN(params)[1, ] <- initialN(params)[1, ] * 2
#' isSteady(params)
#' }
isSteady <- function(params, tol = 0.05, effort = params@initial_effort, ...) {
    UseMethod("isSteady")
}

#' @rdname isSteady
#' @usage NULL
#' @export
isSteady.MizerParams <- function(params, tol = 0.05,
                                 effort = params@initial_effort, ...) {
    params <- validParams(params)
    drift <- tryCatch(steady_biomass_drift(params, effort = effort, ...),
                      error = function(e) NA_real_)
    is.finite(drift) && drift <= tol
}

#' The absolute rates of change of every state variable
#'
#' The shared computation behind [getSteadyResidual()] and
#' `steady_biomass_drift()`. It returns the rates of change themselves rather
#' than the relative ones, because the absolute rate is well defined where the
#' density is zero — a size class with no fish in it can still be filling up —
#' whereas the relative rate there is `0/0`. The two callers then divide, or
#' integrate, as each needs.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param effort The fishing effort to evaluate at.
#' @param dt The step length for the components whose dynamics are only
#'   available as a one-step map.
#' @return A list with the validated `params`, the state (`n`, `n_pp`), the
#'   consumer and resource rates of change (`dNdt`, `dn_pp_dt`) and the list
#'   `other` of relative rates of change of the other components.
#' @noRd
steady_rates <- function(params, effort = params@initial_effort, dt = 1e-4) {
    params <- validParams(params)
    effort <- validEffortVector(effort, params = params)
    assert_that(is.number(dt), dt > 0)

    n <- params@initial_n
    n_pp <- params@initial_n_pp
    n_other <- params@initial_n_other

    # Consumers. `rdd = NULL` asks for the model's own reproduction function to
    # be used, so that the residual reflects the dynamics the user would get.
    # `consumer_residual()` returns A N - S, which is -dt dN/dt at dt = 1, so
    # the sign is flipped to give the rate of change itself.
    dNdt <- -consumer_residual(params, n = n, n_pp = n_pp, n_other = n_other,
                               effort = effort, rdd = NULL)

    # Resource. Its dynamics function is a one-step map, so difference it.
    r <- mizer_rates_subset(params, n = n, n_pp = n_pp, n_other = n_other,
                            t = 0, effort = effort,
                            rates_fns = projectRateFunctions(params),
                            targets = "ResourceMort")
    resource_dynamics_fn <- get(params@resource_dynamics)
    n_pp_new <- resource_dynamics_fn(params, n = n, n_pp = n_pp,
                                     n_other = n_other, rates = r,
                                     t = 0, dt = dt,
                                     resource_rate = params@rr_pp,
                                     resource_capacity = params@cc_pp)
    dn_pp_dt <- (n_pp_new - n_pp) / dt

    # Other components, likewise. Their state can be any object at all, so only
    # a numeric one can be differenced, and there is no size grid to integrate
    # over, so these stay relative rates.
    other <- list()
    for (component in names(params@other_dynamics)) {
        current <- n_other[[component]]
        if (!is.numeric(current)) {
            other[[component]] <- NA_real_
            next
        }
        new <- get(params@other_dynamics[[component]])(
            params, n = n, n_pp = n_pp, n_other = n_other, rates = r,
            t = 0, dt = dt, component = component)
        rate <- (new - current) / (dt * current)
        rate[current == 0] <- NA_real_
        other[[component]] <- rate
    }

    list(params = params, n = n, n_pp = n_pp,
         dNdt = dNdt, dn_pp_dt = dn_pp_dt, other = other)
}

#' How fast the biomasses in a model are drifting
#'
#' The scalar that every steady-state check in mizer is phrased in terms of: the
#' largest relative rate of change of any species' biomass, of the resource
#' biomass, or of any other component. Having one function compute it means
#' `summary()`, [project()], [steady()] and the guards in [getStability()] all
#' report the same number for the same model.
#'
#' ## Why biomass rather than the largest cell
#'
#' The obvious scalar is `max(abs(getSteadyResidual(params)))`, and it is the
#' wrong one. The per-capita rate of change of a *single size class* is
#' dominated by the fastest-relaxing cells, whose residence time near the egg
#' size is hours; a model settled for every practical purpose can carry a cell
#' rate of \eqn{10^4}/year there while nothing observable moves. Under the
#' second-order scheme this is severe enough to invert the ordering: a converged
#' `steady()` run scores *worse* on the cell maximum than a model that has just
#' been knocked off its steady state by `matchGrowth()`.
#'
#' Weighting by biomass removes that: fast cells holding no mass contribute
#' nothing, and the number that comes out is the one the user would actually see
#' drift in [plotBiomass()]. Measured on the North Sea model, it separates the
#' settled models (\eqn{10^{-7}} to \eqn{10^{-3}}) from a model that has been
#' knocked off its steady state (2 to 4) by three orders of magnitude.
#'
#' The integrals go through [sizeIntegral()], so they use whichever quadrature
#' scheme the model is on.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param ... Passed on to `steady_rates()`.
#' @return A single number, in units of 1/year. `NA` if the model has no biomass
#'   at all.
#' @noRd
steady_biomass_drift <- function(params, ...) {
    rates <- steady_rates(params, ...)
    params <- rates$params

    # Consumers: (dB_i/dt) / B_i, with both integrals taken the same way.
    biomass <- as.numeric(sizeIntegral(params, weighting = params@w, n = rates$n))
    dBdt <- as.numeric(sizeIntegral(params, weighting = params@w, n = rates$dNdt))
    consumers <- dBdt[biomass > 0] / biomass[biomass > 0]

    # Resource: the same, over the full grid it lives on.
    wdw <- params@w_full * params@dw_full
    resource_biomass <- sum(rates$n_pp * wdw)
    resource <- if (resource_biomass > 0) {
        sum(rates$dn_pp_dt * wdw) / resource_biomass
    } else {
        numeric(0)
    }

    values <- c(consumers, resource,
                unlist(lapply(rates$other, as.numeric)))
    values <- values[is.finite(values)]
    if (length(values) == 0) return(NA_real_)
    max(abs(values))
}

#' Report that a function has moved the model off its steady state
#'
#' The `match...()` functions rescale abundances or growth parameters per
#' species, which is not a symmetry of the model, so whatever steady state the
#' model was on it is no longer on. This is the report that says so, turning the
#' instruction the documentation used to have to give into something the package
#' says at the moment it becomes true.
#'
#' Note that the `calibrate...()` functions and [scaleModel()] do *not* need
#' this: they apply one overall scaling factor to the whole model, which is an
#' exact symmetry, and leave the residual unchanged to the last digit.
#'
#' `severity = "info"` because mizer did exactly what it was asked to do and is
#' reporting a consequence, which is the rule in
#' `.claude/skills/info-signals.md`; and `level = 3` because in the calibration
#' loop this is expected and the user is about to run [steady()] anyway.
#'
#' @param fname The name of the calling function, without brackets.
#' @return `NULL` invisibly.
#' @noRd
signal_off_steady <- function(fname) {
    signal_info("steady", paste0(
        "`", fname, "()` has rescaled the model and so moved it off its ",
        "steady state. Run `steady()` to settle it again. You can check with ",
        "`getSteadyResidual()`."),
        level = 3, unhandled = "show")
}

#' Report a model that is not at its steady state
#'
#' The shared body of the steady-state checks. Computes the residual, and if it
#' is above tolerance raises the report through the info-signal mechanism, so
#' `info_level` governs it like every other thing mizer says.
#'
#' `severity = "warning"` because the user is being told that something they
#' probably assumed is not true, and `level = 1` because it survives
#' `info_level = 1`; both follow the rule that a report about an assumption that
#' does not hold has to be hard to miss.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param context A sentence naming what the caller was about to do with the
#'   model, appended to the report.
#' @param tol The tolerance to judge against.
#' @return `TRUE` if a report was raised, `FALSE` otherwise, invisibly.
#' @noRd
warn_if_not_steady <- function(params, context,
                               tol = steady_residual_tol()) {
    # `info_level = 0` means silence, and there is no point paying for the rate
    # evaluation only to say nothing about it.
    if (isTRUE(default_info_level() == 0)) return(invisible(FALSE))
    drift <- tryCatch(steady_biomass_drift(params),
                      error = function(e) NA_real_)
    if (!is.finite(drift) || drift <= tol) return(invisible(FALSE))
    signal_info("steady", paste0(
        "This model is not at its steady state: a biomass is changing at ",
        "up to ", signif(drift, 2), " per year. ", context,
        " Run `steady(params)` if that was not intended, or check ",
        "`getSteadyResidual(params)` to see which species are moving."),
        level = 1, severity = "warning", unhandled = "show")
    invisible(TRUE)
}
