# Steady-state residual diagnostic -------------------------------------------
#
# `tuneSteadyState()` and `findSteadyState()` put a model *onto* its steady
# state. This file
# provides the diagnostic that asks whether a model *is* on it, so that the
# instruction "re-run steady() after any match.../calibrate... step" can be
# checked rather than remembered.

#' Default tolerance for calling a model steady
#'
#' The relative rate of biomass change, in units of 1/year, below which mizer
#' treats a model as being at its steady state. It has to sit above what a
#' converged [tuneSteadyState()] actually leaves behind and below the drift that
#' the
#' mistake it is meant to catch produces. Measured on the North Sea model with
#' `steady_biomass_drift()`:
#'
#' | State | Drift (1/year) |
#' |---|---|
#' | `findSteadyState(solver = "newton")` | 4e-7 |
#' | `steady(tol = 1e-4)` | 2e-5 |
#' | `steady()` at its default `tol` | 3e-3 |
#' | as shipped in `NS_params` | 1e-2 |
#' | `newTraitParams()`, freshly built | 0.14 |
#' | after `matchBiomasses()` | 2 |
#' | after `matchGrowth()` | 4 |
#'
#' `0.05` therefore sits more than an order of magnitude above a default
#' `tuneSteadyState()` run and one to three orders below a model knocked off its
#' steady
#' state. Note that [tuneSteadyState()] stops on the relative change in egg
#' production rather than on this drift, so how close it gets is governed by its
#' own `tol` argument; tighten that if you need to settle further.
#'
#' Everything that judges steadiness — the `summary()` line, the [project()]
#' check, the guards in [getStability()] and [getOscillationModeSim()] — goes through
#' this function, so the tolerance has a single definition.
#'
#' @return The tolerance, a single number.
#' @noRd
steady_residual_tol <- function() {
    0.05
}

#' Below what share of a species' biomass a size class holds nothing
#'
#' A size class can carry a positive density and still hold no fish. Above a
#' size where the growth rate vanishes the density decays exponentially: both
#' `dN/dt` and `N` shrink together, so the per-capita rate stays equal to minus
#' the mortality rate for ever, however small the density becomes. Such a class
#' is not evidence that a model is still moving; it is the arithmetic of a
#' number on its way to zero.
#'
#' This is the share of a species' biomass below which a size class is treated
#' as holding nothing. Relevance is measured as a share of biomass rather than
#' of density, because the density in a size spectrum falls by fifteen orders or
#' more between the egg size and `w_max` for entirely healthy reasons, so no
#' threshold on density can separate a dying trace from real large fish. The bin
#' biomass, \eqn{N_i(w) w \Delta w} against the species total, can. Measured
#' shares:
#'
#' | Size class | Share of the species' biomass |
#' |---|---|
#' | one holding an appreciable number of fish | above 1e-6 |
#' | Saithe at 7 kg in a Datta North Sea reconstruction | 1e-11 |
#' | the inaccessible 17-42 kg Saithe tail there | 1e-137 |
#'
#' `1e-8` therefore sits two orders below anything that holds fish and three or
#' more above the traces it is meant to remove. Erring on the large side is the
#' safe direction: whether a state is a fixed point is decided by
#' `steady_biomass_drift()`, which integrates over *every* size class, cut off
#' or not, so a generous cutoff here cannot make a drifting model look settled.
#' All it can do is stop a class with no mass in it from blocking a convergence
#' test. The residual diagnostic needs no cutoff at all: [getSteadyResidual()]
#' weights each class by the biomass it holds, which gives a trace the weight it
#' deserves without anyone having to name a threshold.
#'
#' @return The cutoff, a single number.
#' @noRd
steady_share_cutoff <- function() {
    1e-8
}

#' The quadrature weight that turns a density into a bin biomass
#'
#' The bin integral \eqn{\int_{bin} w\,dw} that [sizeIntegral()] applies when
#' asked for a biomass, in the one form that the two functions in this file need
#' cellwise rather than summed. Formed with [bin_average_weight()], so it
#' follows whichever quadrature scheme the model is on.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @return A numeric vector over `params@w`.
#' @noRd
biomass_bin_weight <- function(params) {
    bin_average_weight(params@w, params) * params@dw
}

#' Which size classes hold a negligible share of their species' biomass
#'
#' The definition behind the cutoff in [distanceSSLogN()], and so behind which
#' size classes [projectUntilSettled()] measures convergence on.
#'
#' The bin biomass is formed with `biomass_bin_weight()`, so the shares are
#' taken with whichever quadrature scheme the model is on and sum to one under
#' both.
#'
#' A species with no biomass at all has no size class holding a share of it, so
#' all of its classes count as negligible and nothing is divided by zero.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param n Consumer densities (species x size).
#' @param cutoff A finite number between 0 and 1, inclusive, giving the share of
#'   a species' biomass below which a size class is treated as holding nothing.
#'   `0` selects nothing, which is how the callers behaved before this cutoff
#'   existed.
#' @return A logical matrix (species x size), `TRUE` where the size class is
#'   negligible.
#' @noRd
negligible_cells <- function(params, n, cutoff = steady_share_cutoff()) {
    if (!is.number(cutoff) || !is.finite(cutoff) ||
            cutoff < 0 || cutoff > 1) {
        stop("`biomass_share_cutoff` must be a finite number between 0 and ",
             "1 (inclusive).", call. = FALSE)
    }
    if (cutoff == 0) {
        return(matrix(FALSE, nrow = nrow(n), ncol = ncol(n)))
    }
    wdw <- biomass_bin_weight(params)
    bin_biomass <- n * rep(wdw, each = nrow(n))
    total_biomass <- rowSums(bin_biomass)
    negligible <- bin_biomass < cutoff * total_biomass
    # A species with no biomass has no size class holding a share of it. Handle
    # that case explicitly rather than changing the documented strict cutoff
    # into an inclusive one for every other species.
    negligible[total_biomass == 0, ] <- TRUE
    negligible
}

#' How far a model is from its steady state
#'
#' `r lifecycle::badge("experimental")`
#' Returns the rate at which the abundances would change if the model were
#' projected forward from its current initial state, resolved by species and by
#' size. At a steady state this is zero, so it answers the question that every
#' calibration workflow otherwise has to remember to ask: *is this model still
#' at its steady state?*
#'
#' There are two ways to say how fast a size class is changing, and `measure`
#' selects between them. Both are in units of 1/year.
#'
#' ## `measure = "biomass"`, the default
#'
#' How much of its species' biomass each size class is adding or removing per
#' year:
#' \deqn{C_i(w) = \frac{1}{B_i}\,\frac{dN_i(w)}{dt}\,w\,\Delta w,
#'       \qquad B_i = \int N_i(w)\,w\,dw.}
#' The bin weight \eqn{w\,\Delta w} is the one [sizeIntegral()] uses for a
#' biomass, so it follows whichever quadrature scheme the model is on (see
#' [second_order_w()]), and the values **add up over sizes to the relative rate
#' of change of the species' biomass**:
#'
#' ```r
#' rowSums(getSteadyResidual(params))    # (dB_i/dt) / B_i, one per species
#' ```
#'
#' That total is the number [isSteady()], the `summary()` line of a
#' \linkS4class{MizerParams} object and `project(check_steady = TRUE)` all judge
#' a model by, and it is the drift that would actually show up in
#' [plotBiomass()]. This array therefore says *where* a model is unsteady in the
#' same currency that mizer uses to decide *whether* it is, and a size class can
#' only be conspicuous in it if it is moving enough biomass to matter.
#'
#' That last property is why this is the default. A size class near the egg size
#' turns over in hours, and one above the size where growth stops decays
#' exponentially towards zero for ever; both carry enormous *per-capita* rates
#' while holding no mass at all. Weighting by biomass gives them the weight they
#' deserve, which is none, with no need for a threshold below which a class is
#' declared to hold nothing.
#'
#' ## `measure = "per_capita"`
#'
#' The rate of change of each size class relative to its own density,
#' \deqn{R_i(w) = \frac{1}{N_i(w)}\frac{dN_i(w)}{dt}.}
#' A value of `0.05` means that size class would grow by about 5% over the first
#' year of a projection, and `-0.05` that it would shrink by about that much.
#' This is the scale-free reading. It shows a size class whose growth and
#' mortality are out of balance even when the class holds a millionth of its
#' species' biomass, which is a real statement about the structure of a model,
#' but not one about whether anything observable is moving.
#'
#' Do not reduce this measure to `max(abs(...))`. Its extremes belong to the
#' fastest-relaxing cells, which are exactly the ones holding no mass: on
#' `NS_params` the largest per-capita rate is 1.2/year, in a cell holding 2e-8
#' of its species' biomass, while the biomass drift is 0.014/year and the model
#' counts as settled. Under the second-order scheme this is severe enough to
#' reverse the ordering between a converged model and one that has just been
#' knocked off its steady state.
#'
#' ## How the rates are obtained
#'
#' For the consumers `dN/dt` is exact, not a finite-difference approximation:
#' the backward-Euler transport coefficients used by [project()] satisfy
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
#' for the semichemostat resource that `findSteadyState(solver = "newton")`
#' requires.
#'
#' ## Reading the result
#'
#' The returned array is an [ArraySpeciesBySize] object, so it prints,
#' summarises and plots itself:
#'
#' ```r
#' res <- getSteadyResidual(params)
#' rowSums(res)                  # the biomass drift of each species
#' plot(res)                     # which species, and at which sizes
#' ```
#'
#' The plot is the diagnostic one: a model that is off steady state is usually
#' off in one species, or one part of the size range, and the plot says which.
#' Mizer's size grid is logarithmic, so the bin widths are proportional to `w`
#' and the default measure keeps the shape of a density per unit of log size
#' when drawn against a logarithmic size axis: equal areas under the curve are
#' equal contributions to the drift.
#'
#' With `measure = "per_capita"` a size class with no fish in it has no relative
#' rate of change — it is `0/0` — and is returned as `NA`, so pass
#' `na.rm = TRUE` to any summary. The default measure needs no such exception:
#' `dN/dt` is perfectly well defined in a class with no fish in it, which can be
#' filling up, and its contribution to the biomass drift is reported like any
#' other.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param effort The fishing effort at which to evaluate the residual. By
#'   default the initial effort stored in `params`, which is the effort the
#'   model's steady state belongs to.
#' @param dt The step length used for the resource and other components, whose
#'   dynamics functions are only available as one-step maps. Smaller is more
#'   accurate. Not used for the consumers, whose rate is exact.
#' @param measure `r lifecycle::badge("experimental")`
#'   Which rate of change to report. `"biomass"` (the default) gives the
#'   contribution of each size class to the relative rate of change of its
#'   species' biomass, which sums over sizes to the drift that [isSteady()]
#'   judges a model by. `"per_capita"` gives the rate of change of each size
#'   class relative to its own density. See the sections above.
#' @return An [ArraySpeciesBySize] object (species x size) of rates in 1/year:
#'   with `measure = "biomass"`, contributions to each species' relative rate of
#'   biomass change, which sum over sizes to that rate and are `NA` only for a
#'   species with no biomass at all; with `measure = "per_capita"`, per-capita
#'   rates of change, `NA` where the size class holds no fish. It carries two
#'   further attributes:
#'   \describe{
#'     \item{`resource`}{The same measure for the resource, a numeric vector
#'       over `w_full`.}
#'     \item{`other`}{A named list with one entry per other component, holding
#'       its per-capita rate of change, or `NA` for a component whose state is
#'       not numeric.}
#'   }
#' @seealso [isSteady()], [tuneSteadyState()], [findSteadyState()],
#'   [getStability()]
#' @export
#' @family summary functions
#' @concept summary_function
#' @examples
#' # The relative rate of change of each species' biomass, in 1/year
#' rowSums(getSteadyResidual(NS_params))
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
                              dt = 1e-4,
                              measure = c("biomass", "per_capita")) {
    measure <- match.arg(measure)
    rates <- steady_rates(params, effort = effort, dt = dt)
    params <- rates$params
    n <- rates$n

    if (measure == "per_capita") {
        residual <- rates$dNdt / n
        # A class with no density has no rate of change relative to it.
        residual[n == 0] <- NA_real_
        resource <- rates$dn_pp_dt / rates$n_pp
        resource[rates$n_pp == 0] <- NA_real_
        value_name <- "Per-capita rate of change"
        # A ratio of two bin averages belongs to no smaller a size than the bin
        # it was formed in, but it is not itself a bin integral, so it keeps the
        # point representation it has always had.
        representation <- "point"
    } else {
        # The bin biomasses, formed with the quadrature the model is on, so that
        # the row sums agree with `getBiomass()` and with the biomass drift the
        # steady-state checks are phrased in.
        wdw <- biomass_bin_weight(params)
        biomass <- rowSums(n * rep(wdw, each = nrow(n)))
        residual <- rates$dNdt * rep(wdw, each = nrow(n)) / biomass
        # A species with no biomass has no relative rate of change of it. This
        # is the one place the default measure can be undefined, and it is a
        # whole species rather than a size class.
        residual[biomass == 0, ] <- NA_real_
        wdw_full <- params@w_full * params@dw_full
        resource_biomass <- sum(rates$n_pp * wdw_full)
        resource <- if (resource_biomass > 0) {
            rates$dn_pp_dt * wdw_full / resource_biomass
        } else {
            rep(NA_real_, length(wdw_full))
        }
        value_name <- "Biomass drift contribution"
        # Each value is an integral over its bin, so it belongs at the centre of
        # the bin rather than at its left edge.
        representation <- "average"
    }

    residual <- ArraySpeciesBySize(residual, value_name = value_name,
                                   units = "1/year", params = params,
                                   representation = representation)
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
#' @seealso [getSteadyResidual()], [tuneSteadyState()], [findSteadyState()]
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
#' The state is taken from the model unless it is supplied explicitly. The
#' explicit form is what the convergence checks in `project_until_settled()`
#' use: they need the drift at the state the run has just reached, and writing
#' that state into `params` only to read it back would both copy the object and
#' invalidate its validation fingerprint on every block.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param effort The fishing effort to evaluate at.
#' @param dt The step length for the components whose dynamics are only
#'   available as a one-step map.
#' @param n,n_pp,n_other The state to evaluate at. By default the state stored
#'   in `params`.
#' @return A list with the validated `params`, the state (`n`, `n_pp`), the
#'   consumer and resource rates of change (`dNdt`, `dn_pp_dt`) and the list
#'   `other` of relative rates of change of the other components.
#' @noRd
steady_rates <- function(params, effort = params@initial_effort, dt = 1e-4,
                         n = params@initial_n,
                         n_pp = params@initial_n_pp,
                         n_other = params@initial_n_other) {
    params <- validParams(params)
    effort <- validEffortVector(effort, params = params)
    assert_that(is.number(dt), dt > 0)

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
#' `summary()`, [project()], [tuneSteadyState()] and the guards in
#' [getStability()] all
#' report the same number for the same model.
#'
#' ## Why biomass rather than the largest cell
#'
#' The obvious scalar is
#' `max(abs(getSteadyResidual(params, measure = "per_capita")))`, and it is the
#' wrong one. The per-capita rate of change of a *single size class* is
#' dominated by the fastest-relaxing cells, whose residence time near the egg
#' size is hours; a model settled for every practical purpose can carry a cell
#' rate of \eqn{10^4}/year there while nothing observable moves. Under the
#' second-order scheme this is severe enough to invert the ordering: a converged
#' `tuneSteadyState()` run scores *worse* on the cell maximum than a model that
#' has just
#' been knocked off its steady state by `matchGrowth()`.
#'
#' Weighting by biomass removes that: fast cells holding no mass contribute
#' nothing, and the number that comes out is the one the user would actually see
#' drift in [plotBiomass()]. Measured on the North Sea model, it separates the
#' settled models (\eqn{10^{-7}} to \eqn{10^{-3}}) from a model that has been
#' knocked off its steady state (2 to 4) by three orders of magnitude.
#'
#' This is the same weighting that [getSteadyResidual()] applies by default, and
#' the same number: `max(abs(rowSums(getSteadyResidual(params))))` is the
#' consumer part of what this function returns. The difference is only that this
#' function also takes in the resource and the other components, and reduces the
#' whole lot to the single number that the tolerances are stated against.
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
#' loop this is expected and the user is about to run [tuneSteadyState()] anyway.
#'
#' @param fname The name of the calling function, without brackets.
#' @return `NULL` invisibly.
#' @noRd
signal_off_steady <- function(fname) {
    signal_info("steady", paste0(
        "`", fname, "()` has rescaled the model and so moved it off its ",
        "steady state. Run `tuneSteadyState()` to settle it again. You can ",
        "check with ",
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
        " Run `tuneSteadyState(params)` if that was not intended, or check ",
        "`getSteadyResidual(params)` to see which species are moving."),
        level = 1, severity = "warning", unhandled = "show")
    invisible(TRUE)
}
