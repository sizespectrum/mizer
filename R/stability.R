#' Relative finite-difference step scale
#'
#' Returns the scale that [getStability()] multiplies by `h` to get the
#' finite-difference step for each state variable: the variable's own magnitude
#' where that is nonzero, and otherwise the local scale of the spectrum as
#' supplied by `positive_initial_guess()` (log-interpolated from the nonzero
#' neighbours). Cells whose whole row is zero get no local scale from
#' `positive_initial_guess()`; they fall back to `.Machine$double.eps`, which is
#' sufficient to resolve linear invasion dynamics at a zero baseline.
#'
#' A step floored at an *absolute* `.Machine$double.eps` instead would be
#' swamped by the rounding error of the order-`N` outputs it perturbs, and the
#' column of the Jacobian would come out as noise (in practice exactly zero).
#'
#' @param x The state variables (a vector of densities).
#' @param local_scale A strictly positive local scale for each entry of `x`, or
#'   zero where none could be determined.
#' @return A strictly positive vector of step scales, the same length as `x`.
#' @noRd
fd_step_scale <- function(x, local_scale) {
    local_scale[!is.finite(local_scale) | local_scale <= 0] <-
        .Machine$double.eps
    pmax(abs(x), local_scale)
}

# Stability analysis ----------------------------------------------------------

# Assemble explicit coverage before constructing a Jacobian. A declaration that
# removes coordinates does not by itself establish time-independent forcing.
stability_coverage <- function(params, states, discrete = FALSE) {
    required <- names(states)[vapply(states, function(spec) {
        !isTRUE(spec$fixed) || isTRUE(spec$temporal_required)
    }, logical(1))]
    reasons <- new.env(parent = emptyenv())
    reasons$table <- data.frame(owner = character(), capability = character(),
                                reason = character(), stringsAsFactors = FALSE)
    add <- function(owner, capability, reason) {
        reasons$table <- rbind(reasons$table, data.frame(
            owner = owner, capability = capability, reason = reason,
            stringsAsFactors = FALSE
        ))
    }
    autonomous <- TRUE
    for (owner in required) {
        spec <- states[[owner]]
        if (!isTRUE(spec$temporal_current) || !isTRUE(spec$autonomous)) {
            add(owner, "autonomy", spec$temporal_reason %||%
                    "Autonomous dynamics have not been declared.")
            autonomous <- if (identical(spec$autonomous, FALSE) ||
                                  isFALSE(autonomous)) FALSE else NA
        }
        if (isTRUE(spec$fixed)) next
        if (!isTRUE(spec$mask_current)) {
            add(owner, "state", spec$mask_reason %||%
                    "The state declaration is not current.")
        }
        if (!any(spec$dynamic_mask)) next
        if (startsWith(owner, "component:")) {
            add(owner, "jacobian", paste0(
                "Dynamic component coordinates are not supported by this ",
                "Jacobian. Use explicit conditional scope to hold '", owner,
                "' fixed."
            ))
        } else if (!builtin_consumer_projection(params)) {
            add(owner, if (discrete) "step" else "equation",
                "The current projection driver has no supported operator.")
        } else if (!discrete && owner == "resource" &&
                       (!isTRUE(spec$equation_current) ||
                            is.null(spec$rate_fun))) {
            add(owner, "equation", spec$equation_reason %||%
                    "No current resource equation is declared.")
        } else if (!discrete && owner == "resource") {
            reason <- resource_equation_domain(params, spec, spec$dynamic_mask)
            if (!is.null(reason)) add(owner, "equation", reason)
        }
    }
    unavailable <- reasons$table
    missing <- unique(unavailable$owner)
    list(required = required, covered = setdiff(required, missing),
         missing = missing, unavailable = unavailable,
         complete = !length(missing), autonomous = autonomous)
}

#' Shared scaffolding for the stability analyses
#'
#' The masks describe structural state support, not the Newton solver's
#' logarithmic active set. In particular an absent species can invade.
#' @noRd
stability_context <- function(params, effort = params@initial_effort,
                              fn_name = "getStability",
                              map_name = "the rates of change",
                              scope = c("whole", "conditional"),
                              fixed = character(), t = 0,
                              discrete = FALSE) {
    params <- validParams(params)
    effort <- validEffortVector(effort, params = params)
    params@initial_effort <- effort
    if (!is.numeric(t) || is.complex(t) || length(t) != 1L || !is.finite(t)) {
        stop("`t` must be a finite number.", call. = FALSE)
    }
    prepared <- prepare_steady_assessment(params, scope = scope, fixed = fixed)
    states <- prepared$state_capabilities
    coverage <- stability_coverage(params, states, discrete = discrete)
    assessment <- steady_assessment(params, prepared = prepared,
                                    effort = effort, t = t)
    active <- list(mask = states$consumers$dynamic_mask,
                   w_top = support_top_idx(params))
    active_idx <- which(active$mask)
    resource_idx <- which(states$resource$dynamic_mask)
    n_fish_active <- length(active_idx)
    n_npp <- length(resource_idx)
    n_ss <- params@initial_n
    npp_ss <- params@initial_n_pp
    n_vec <- n_ss[active_idx]
    npp_vec <- as.numeric(npp_ss)[resource_idx]
    n_scale <- fd_step_scale(n_vec, positive_initial_guess(
        n_ss, params@w_min_idx, active$w_top
    )[active_idx])
    npp_local <- positive_initial_guess(
        matrix(as.numeric(npp_ss), nrow = 1),
        w_min_idx = 1L, w_top = length(npp_ss)
    )
    npp_scale <- fd_step_scale(npp_vec, as.numeric(npp_local)[resource_idx])
    x0 <- c(n_vec, npp_vec)
    if (any(!is.finite(x0)) || any(x0 < 0)) {
        stop("Stability analysis requires finite nonnegative dynamic states.",
             call. = FALSE)
    }
    unpack <- function(x) {
        n <- n_ss
        n[active_idx] <- x[seq_len(n_fish_active)]
        n_pp <- npp_ss
        n_pp[resource_idx] <- x[n_fish_active + seq_len(n_npp)]
        list(n = n, n_pp = n_pp)
    }
    check <- function(x, where) {
        if (!is.numeric(x) || length(x) != length(x0) ||
                !all(is.finite(x))) {
            stop(fn_name, "(): ", map_name, " returned non-finite values ",
                 "or an invalid state shape ", where, ".", call. = FALSE)
        }
    }
    describe <- function(k) {
        if (k <= n_fish_active) {
            rc <- arrayInd(active_idx[k], dim(n_ss))
            paste0("when perturbing species ",
                   params@species_params$species[rc[1]], " at w = ",
                   signif(params@w[rc[2]], 3))
        } else {
            paste0("when perturbing the resource at w = ",
                   signif(params@w_full[resource_idx[k - n_fish_active]], 3))
        }
    }
    prescribed_masks <- lapply(states, `[[`, "prescribed_mask")
    prescribed <- names(states)[vapply(states, function(spec) {
        !isTRUE(spec$fixed) && !any(spec$dynamic_mask) &&
            isTRUE(spec$mask_current)
    }, logical(1))]
    list(params = params, effort = effort, n_other = params@initial_n_other,
         active = active, active_idx = active_idx, resource_idx = resource_idx,
         n_fish_active = n_fish_active, n_npp = n_npp,
         rates_fns = projectRateFunctions(params),
         flux_limiter = flux_limiter_scheme(params),
         n_ss = n_ss, npp_ss = npp_ss, x0 = x0,
         scale = c(n_scale, npp_scale), unpack = unpack, check = check,
         describe = describe, states = states, coverage = coverage,
         scope = prepared$scope, fixed = prepared$fixed, t = t,
         prescribed = prescribed, prescribed_masks = prescribed_masks,
         assessment = assessment)
}

#' Finite-difference Jacobian of a map on the packed state
#'
#' Differences `f` column by column at the steady state, with a multiplicative
#' step `h * scale[k]`.
#'
#' Every state at which `f` is evaluated satisfies \eqn{N \ge 0}: where a
#' centred step would push a cell negative — which can only happen for a cell at
#' (or below) the floor described in `stability_context()` — the column is
#' differenced forwards from the unperturbed state instead. Such columns are
#' second order in `h`, using two forward perturbations.
#'
#' @param ctx The context from `stability_context()`.
#' @param f The map, a function of the packed state vector returning a vector of
#'   the same length.
#' @param h Relative step size.
#' @return The Jacobian matrix of `f` at the steady state.
#' @noRd
stability_jacobian <- function(ctx, f, h) {
    x <- ctx$x0
    n <- length(x)
    jacobian <- matrix(0, nrow = n, ncol = n)

    # Baseline for the one-sided columns, shared by all of them. It has to be
    # evaluated rather than assumed: the state is a fixed point only to the
    # tolerance of whatever produced it, and that residual would otherwise be
    # divided by eps into every one-sided column.
    base <- f(x)
    ctx$check(base, "at the unperturbed steady state")

    for (k in seq_len(n)) {
        eps <- h * ctx$scale[k]
        x_plus <- x
        x_plus[k] <- x[k] + eps
        if (x[k] - eps < 0) {
            # A centred step would leave the physical cone N >= 0. Difference
            # forwards instead: at the boundary of the cone the one-sided
            # derivative is the right object anyway, because the dynamics never
            # visit the states a centred step would sample. This keeps every
            # rate function evaluation inside N >= 0, so extensions need not be
            # defined at negative abundances.
            x_twice <- x
            x_twice[k] <- x[k] + 2 * eps
            jacobian[, k] <- (-3 * base + 4 * f(x_plus) - f(x_twice)) /
                (2 * eps)
        } else {
            x_minus <- x
            x_minus[k] <- x[k] - eps
            jacobian[, k] <- (f(x_plus) - f(x_minus)) / (2 * eps)
        }
        ctx$check(jacobian[, k], ctx$describe(k))
    }
    jacobian
}

#' The continuous equations on the packed state
#'
#' Resource derivatives are exact built-in equations or explicitly declared
#' callbacks. A resource step is never probed to infer a continuous equation.
#' @noRd
stability_rhs <- function(ctx) {
    params <- ctx$params
    function(x) {
        st <- ctx$unpack(x)
        r <- ctx$rates_fns$Rates(params, n = st$n, n_pp = st$n_pp,
                                 n_other = ctx$n_other, t = ctx$t,
                                 effort = ctx$effort,
                                 rates_fns = ctx$rates_fns)
        dNdt <- if (ctx$n_fish_active) {
            -consumer_residual(params, n = st$n, n_pp = st$n_pp,
                               n_other = ctx$n_other, effort = ctx$effort,
                               rdd = r$rdd, rates_fns = ctx$rates_fns,
                               flux_limiter = ctx$flux_limiter, rates = r)
        } else {
            ctx$n_ss * 0
        }
        dnpp <- if (ctx$n_npp) {
            evaluate_resource_equation(
                params, ctx$states$resource, n = st$n, n_pp = st$n_pp,
                n_other = ctx$n_other, rates = r, t = ctx$t
            )
        } else {
            ctx$npp_ss * 0
        }
        c(dNdt[ctx$active_idx], as.numeric(dnpp)[ctx$resource_idx])
    }
}

#' Mizer's numerical one-step map on the packed state
#'
#' Returns the function that [getDiscreteStability()] differentiates: one step
#' of length `dt` of the dynamics as [project()] takes it with
#' `method = "euler"`, using the same `project_n_loop()` C++ Thomas solver, the
#' same transport coefficients and the model's own `resource_dynamics`.
#'
#' It has to be that map and not one that merely agrees with it to first order
#' in `dt`, because the whole purpose of the function is to say what mizer's
#' solver does at a particular `dt`. An earlier version stepped the resource
#' with backward Euler so that the fast resource modes could not land on a
#' discrete eigenvalue of exactly zero. There is nothing wrong with such an
#' eigenvalue: it is the correct statement that the mode relaxes inside a single
#' step, and it can never be the largest in modulus, so it never touches the
#' spectral radius.
#'
#' @param ctx The context from `stability_context()`.
#' @param dt The step length.
#' @return A function of the packed state vector.
#' @noRd
stability_step <- function(ctx, dt) {
    params  <- ctx$params
    w_top   <- support_top_idx(params)
    resource_dyn <- if (ctx$n_npp) {
        resolve_model_callback(params@resource_dynamics, "resource dynamics")
    }

    function(x) {
        st <- ctx$unpack(x)
        r <- ctx$rates_fns$Rates(params, n = st$n, n_pp = st$n_pp,
                                 n_other = ctx$n_other, t = ctx$t,
                                 effort = ctx$effort,
                                 rates_fns = ctx$rates_fns)
        coefs <- get_transport_coefs(params, n = st$n, g = r$e_growth,
                                     mu = r$mort, dt = dt,
                                     recruitment_flux = r$rdd,
                                     d = r$diffusion,
                                     flux_limiter = ctx$flux_limiter)
        n_projected <- project_n_loop(st$n, coefs$a, coefs$b, coefs$c, coefs$S,
                                      params@w_min_idx)
        args <- list(params = params, n = st$n, n_pp = st$n_pp,
                     n_other = ctx$n_other, rates = r, t = ctx$t, dt = dt,
                     resource_rate = params@rr_pp,
                     resource_capacity = params@cc_pp)
        npp_out <- if (ctx$n_npp) {
            invoke_model_step(resource_dyn, args)
        } else {
            st$n_pp
        }
        c(zero_above_support(n_projected, w_top)[ctx$active$mask],
          as.numeric(npp_out)[ctx$resource_idx])
    }
}

#' Reshape the leading eigenvectors back into abundance space
#'
#' @param evecs The eigenvector matrix, columns already sorted.
#' @param ctx The context from `stability_context()`.
#' @param cols The columns of `evecs` to reshape. Defaults to the first two.
#' @return A list with `$fish`, a complex `(n_species, n_sizes, length(cols))`
#'   array, and `$resource`, a complex `(n_w_full, length(cols))` matrix.
#' @noRd
stability_eigenvectors <- function(evecs, ctx,
                                   cols = seq_len(min(2L, ncol(evecs)))) {
    n_ss  <- ctx$n_ss
    no_sp <- nrow(n_ss)
    no_w  <- ncol(n_ss)
    n_out <- length(cols)

    fish <- array(0 + 0i, dim = c(no_sp, no_w, n_out),
                  dimnames = c(dimnames(n_ss), list(NULL)))
    resource <- matrix(0 + 0i, nrow = length(ctx$npp_ss), ncol = n_out)

    # Fish abundances and resource densities differ by many orders of
    # magnitude, so a normalisation by absolute modulus would be set entirely
    # by the resource block and leave the fish block at ~1e-13. Scale instead
    # by the largest perturbation *relative to the local steady state*, which
    # is the quantity the two blocks have in common. `ctx$x0` is the packed
    # steady state, in the same row order as the eigenvectors.
    x_ss <- ctx$x0
    pos  <- x_ss > 0

    for (k in seq_len(n_out)) {
        # One scalar for the whole state vector, so that the relative
        # amplitude and phase between the fish block and the resource block
        # are those of the mode. Normalising the two blocks separately would
        # rescale them by different factors and destroy exactly the
        # information that makes the resource component worth returning.
        v <- evecs[, cols[k]]
        max_rel <- if (any(pos)) max(Mod(v[pos]) / x_ss[pos]) else 0
        if (max_rel > 0) v <- v / max_rel

        fish_mode <- matrix(0 + 0i, nrow = no_sp, ncol = no_w)
        fish_mode[ctx$active_idx] <- v[seq_len(ctx$n_fish_active)]
        fish[, , k] <- fish_mode
        resource[ctx$resource_idx, k] <-
            v[ctx$n_fish_active + seq_len(ctx$n_npp)]
    }
    rownames(resource) <- names(ctx$npp_ss)
    list(fish = fish, resource = resource)
}

#' Drop the trailing mode dimension of a `stability_eigenvectors()` result
#'
#' @param ev A `stability_eigenvectors()` list holding a single mode.
#' @return A list with `$fish`, a complex `(n_species, n_sizes)` matrix, and
#'   `$resource`, a complex vector of length `n_w_full`.
#' @noRd
single_eigenvector <- function(ev) {
    fish <- ev$fish[, , 1, drop = FALSE]
    dim(fish) <- dim(ev$fish)[1:2]
    dimnames(fish) <- dimnames(ev$fish)[1:2]
    list(fish = fish, resource = ev$resource[, 1])
}

# A local linearisation at a non-equilibrium remains useful diagnostics, but
# cannot establish the stability of a fixed point.
stability_verdict <- function(ctx, value, tolerance, discrete = FALSE) {
    if (!ctx$coverage$complete) return("inconclusive")
    if (!length(ctx$x0)) return("not_applicable")
    if (!isTRUE(ctx$assessment$verdict)) return("inconclusive")
    boundary <- if (discrete) 1 else 0
    if (value < boundary - tolerance) return("stable")
    if (value > boundary + tolerance) return("unstable")
    "inconclusive"
}

stability_model_key <- function(params) {
    digest::digest(params, algo = "xxhash64")
}

stability_result <- function(ctx, tolerance, value, discrete = FALSE) {
    status <- stability_verdict(ctx, value, tolerance, discrete)
    coverage <- ctx$coverage
    autonomous <- coverage$autonomous
    coverage$autonomous <- NULL
    list(status = status,
         stable = switch(status, stable = TRUE, unstable = FALSE, NA),
         scope = ctx$scope, fixed = ctx$fixed, prescribed = ctx$prescribed,
         prescribed_masks = ctx$prescribed_masks, coverage = coverage,
         autonomous = autonomous, time = ctx$t, tolerance = tolerance,
         assessment = ctx$assessment)
}

#' Analyse the dynamic stability of a mizer steady state
#'
#' `r lifecycle::badge("experimental")`
#' Computes the eigenvalues of the linearised dynamics at the steady state
#' stored in `params@initial_n`. These eigenvalues determine whether the steady
#' state is dynamically stable and, where the spectrum contains a complex pair,
#' the period at which the model oscillates.
#'
#' ## Mathematical background
#'
#' Mizer discretises the size axis but not time: on the size grid the model is a
#' system of ordinary differential equations
#' \deqn{\frac{dN}{dt} = F(N, n_{pp}),}
#' where \eqn{F} collects the divergence of the growth flux, the mortality sink
#' and the reproductive influx at the egg size, assembled with the spatial
#' scheme configured via [second_order_w()]. `getStability()` differentiates
#' \eqn{F} directly, by centred finite differences in each state variable, and
#' returns the eigenvalues \eqn{\lambda_i} of the resulting Jacobian
#' \eqn{J = \partial F/\partial N}. The steady state is **stable** when all of
#' them satisfy \eqn{\text{Re}(\lambda_i) < 0} and **unstable** when at least
#' one exceeds 0.
#'
#' No time step enters this calculation. The eigenvalues are a property of the
#' model, not of any solver: they describe the continuous-time dynamics of the
#' semi-discretised model, and are what a simulation with a small enough time
#' step converges to. The stability of the numerical step itself is a separate
#' question, answered by [getDiscreteStability()].
#'
#' A complex-conjugate pair \eqn{\lambda = \sigma \pm i\omega} is an
#' oscillatory mode: a perturbation along it rings with period
#' \deqn{T = \frac{2\pi}{|\omega|} \text{ years,}}
#' growing or decaying as \eqn{e^{\sigma t}}. The pair with the largest
#' \eqn{\sigma} is returned as `leading_oscillatory_eigenvalue`, with its
#' period and eigenvector.
#'
#' That is a statement about the mode, not about a bifurcation. A **Hopf
#' bifurcation** is the event of such a pair *crossing* the imaginary axis, and
#' a single spectrum cannot show a crossing: the leading oscillatory mode of a
#' comfortably stable model can sit far to the left, ringing only as a transient
#' on the way back to the fixed point. Establishing a Hopf bifurcation means
#' watching \eqn{\sigma} pass through zero as a parameter is varied, which is
#' what [scanModel()] is for. Only then is \eqn{T} the period of an emerging
#' limit cycle; otherwise it is the period of a damped (or growing) oscillation.
#'
#' ## What is in the Jacobian
#'
#' The default scope is the whole model. Consumers and resource are supported
#' dynamic coordinates; consumers include absent-species invasion modes. Only
#' structural support limits, explicit prescribed masks and an explicit
#' conditional `fixed` selection remove coordinates. Zero rates, zero capacity
#' and neutral eigenvalues never remove a dynamic coordinate.
#'
#' Continuous resource dynamics require a built-in equation or a declared
#' `rate_fun`; a custom step alone is insufficient. Dynamic component Jacobian
#' coordinates are not yet supported. Missing state, equation or autonomous
#' dynamics declarations return an inconclusive report without eigenvalues.
#' To analyse a subsystem, request `scope = "conditional"` and name the owners
#' to hold fixed. The result records this scope and the excluded coordinates.
#'
#' A stability verdict additionally requires the scoped steady assessment to
#' pass. Otherwise the eigenvalues, when available, are local diagnostics only.
#' Eigenvalues within `stability_tol` of the imaginary axis are inconclusive.
#'
#' Reproduction is a state-dependent rate like any other: the reproduction
#' function stored in `params@rates_funcs$RDD` is evaluated at each perturbed
#' state, so the feedback from the spectra back onto the influx of eggs is part
#' of the Jacobian, exactly as it is part of [project()]. There is no option to
#' pin the reproduction rate at its value at the fixed point. A model in which
#' reproduction really is constant expresses that as a model: with
#' `rates_funcs$RDD = "constantRDD"` the derivative of the reproduction rate is
#' zero and the pinned Jacobian is what the analysis returns.
#'
#' This is why the stability of a steady state depends on the reproduction
#' parameters even though the steady state itself does not.
#' [setBevertonHolt()] moves along a family of `erepro`/`R_max` pairs that all
#' leave the same fixed point, but they do not all leave the same dynamics: at a
#' [reproduction_level()] near 1 the reproduction rate barely responds to the
#' energy invested in it, approaching the constant-reproduction case, while at a
#' level near 0 it follows that energy proportionally. The two ends can differ
#' in their verdict, so the analysis has to read the model rather than take an
#' argument.
#'
#' ## Numerical details
#'
#' The Jacobian is computed numerically using a multiplicative (relative)
#' finite-difference step \eqn{h \cdot N^*}. Where a cell sits at exactly zero
#' and so has no scale of its own, the step is floored at the local scale of the
#' spectrum, interpolated from the nonzero neighbours, so that the cell still
#' gets a resolved derivative rather than a column of rounding error.
#'
#' Every state at which the rates are evaluated satisfies \eqn{N \ge 0}: where
#' a centred step would push a cell negative — which can only happen for a cell
#' at (or below) the floor described above — the column is differenced forwards
#' from the unperturbed state instead. At the boundary of the physical cone the
#' one-sided derivative is the appropriate object anyway, since the dynamics
#' never visit the states a centred step would sample. A rate function
#' registered with [setRateFunction()] therefore never has to be defined at
#' negative abundances. Two forward perturbations give a second-order
#' one-sided derivative at the boundary.
#'
#' @param params A \linkS4class{MizerParams} object whose `initial_n` holds the
#'   steady state to analyse. Typically the output of [findSteadyState()].
#' @param effort The fishing effort to use. By default the initial effort
#'   stored in `params`.
#' @param h Relative step size for centred finite differences. Default `1e-4`.
#'   The result should not depend on this choice. If it does, the dynamics are
#'   not smooth at the state being analysed — see the section below.
#' @param scope `"whole"` (the default) or `"conditional"`.
#' @param fixed Canonical owner IDs held fixed for an explicitly conditional
#'   analysis: `"consumers"`, `"resource"`, or `"component:name"`.
#' @param t Time at which to evaluate the state and rates, in years.
#' @param stability_tol Nonnegative margin used to distinguish strict stability
#'   or instability from a marginal spectrum. Default `1e-8`.
#' @return A named list with the following components:
#'   \describe{
#'     \item{`eigenvalues`}{Complex vector of the continuous-time eigenvalues
#'       \eqn{\lambda_i}, sorted by decreasing real part.}
#'     \item{`max_real_part`}{The largest real part of the eigenvalues:
#'       \eqn{\max_i \text{Re}(\lambda_i)}. A value above `stability_tol`
#'       establishes instability only for complete, autonomous coverage and a
#'       passing scoped steady assessment.}
#'     \item{`stable`}{`TRUE` for strictly stable, `FALSE` for unstable and
#'       `NA` for inconclusive or inapplicable analyses.}
#'     \item{`dominant_period`}{The period (in years) of the dominant
#'       eigenvalue: `2*pi / abs(Im(lambda_1))`. `Inf` for a real
#'       dominant eigenvalue (monotone dynamics).}
#'     \item{`oscillation_period`}{Period (in years) of the oscillatory mode
#'       below, \eqn{2\pi/|\omega|}; `NULL` when no complex eigenvalue exists.
#'       It is the period at which the model rings, and only at a Hopf
#'       bifurcation — where the real part is zero — the period of a limit
#'       cycle.}
#'     \item{`leading_oscillatory_eigenvalue`}{The complex eigenvalue with the
#'       largest real part, or `NULL` when there is none. Its real part is the
#'       rate at which that oscillation grows, so a strongly negative one means
#'       the ringing is a transient, not a cycle the model settles onto.}
#'     \item{`leading_oscillatory_eigenvector`}{Its eigenvector, with `$fish`, a
#'       complex `(n_species, n_sizes)` matrix, and `$resource`, a complex
#'       vector of length `n_w_full`. This is the mode [getOscillationModeSim()]
#'       draws, and it is not in general one of `leading_eigenvectors`: the
#'       dominant mode of the system can be real while the dominant
#'       *oscillatory* mode is well down the spectrum.}
#'     \item{`n_active`}{Dimension of the Jacobian: the number of active fish
#'       consumer and resource coordinates; `NA` if coverage is incomplete.}
#'     \item{`leading_eigenvectors`}{The eigenvectors of the two eigenvalues
#'       with the largest real part, reshaped back into the state space: a list
#'       with `$fish`, a complex array of shape `(n_species, n_sizes, 2)` with
#'       the same species and size dimnames as `params@initial_n`, and
#'       `$resource`, a complex matrix of shape `(n_w_full, 2)`.
#'       Each eigenvector is normalised by a single scalar covering both
#'       blocks, so that the relative amplitude and phase between fish and
#'       resource are those of the mode. The scalar is chosen so that the
#'       largest perturbation *relative to the steady state*,
#'       \eqn{|v_i| / x^*_i}, is 1 somewhere in the state: an absolute
#'       normalisation would be set entirely by the resource, whose densities
#'       dwarf the fish abundances. `Mod(fish[, , 1]) / initialN(params)` is
#'       therefore the relative amplitude pattern, peaking at 1 in whichever
#'       cell swings hardest. The real and imaginary parts of eigenvector 1
#'       span the two-dimensional oscillation plane of the dominant mode.}
#'     \item{`params`}{The validated `params` object the analysis was made at.}
#'     \item{`status`}{`"stable"`, `"unstable"`, `"inconclusive"`, or
#'       `"not_applicable"` when a complete scope has no dynamic coordinates.}
#'     \item{`scope`, `fixed`, `prescribed`, `prescribed_masks`}{The requested
#'       scope and the explicitly excluded states.}
#'     \item{`coverage`}{Required, covered and missing owners, an `unavailable`
#'       reason table, and a logical `complete` flag.}
#'     \item{`autonomous`, `time`, `tolerance`, `assessment`}{Temporal coverage,
#'       evaluation time, stability margin and the scoped steady assessment.}
#'   }
#' @section Requires smooth dynamics:
#' The finite-difference Jacobian is only meaningful if the rates of change are
#' differentiable at \eqn{N^*}. A custom rate function registered with
#' [setRateFunction()] that jumps as a function of the abundances breaks this in
#' two ways. If the state sits on the switching threshold, some perturbations
#' straddle it and pick up the jump, and the reported eigenvalues then vary
#' wildly with `h`. If the state is near but not on the threshold, no
#' perturbation crosses it, and the function silently returns the stability of
#' the single branch the state happens to lie on. Those local eigenvalues do
#' not describe trajectories that cross the discontinuity. A stability verdict
#' additionally requires a passing scoped steady assessment.
#'
#' Re-running with a different `h` is the cheapest check: if the answer moves,
#' do not trust it. See [Discontinuous rate
#' functions](https://sizespectrum.org/mizer/articles/discontinuous_rates.html).
#'
#' @seealso [findSteadyState()], [getDiscreteStability()],
#'   [getOscillationModeSim()]
#' @export
getStability <- function(params,
                         effort = params@initial_effort,
                         h = 1e-4,
                         scope = c("whole", "conditional"),
                         fixed = character(), t = 0,
                         stability_tol = 1e-8) {
    assert_that(is.number(h), is.finite(h), h > 0,
                is.number(stability_tol), is.finite(stability_tol),
                stability_tol >= 0)
    ctx <- stability_context(params, effort = effort, scope = scope,
                             fixed = fixed, t = t)
    evals <- NULL
    evecs <- matrix(complex(), 0, 0)
    max_real_part <- NA_real_
    dominant_period <- NA_real_
    oscillation_period <- NULL
    leading_oscillatory_eigenvalue <- NULL
    osc_vector <- NULL
    vectors <- NULL
    n_active <- NA_integer_
    if (ctx$coverage$complete) {
        n_active <- length(ctx$x0)
        if (n_active) {
            eig <- eigen(stability_jacobian(ctx, stability_rhs(ctx), h))
            ord <- order(Re(eig$values), decreasing = TRUE)
            evals <- as.complex(eig$values[ord])
            evecs <- eig$vectors[, ord, drop = FALSE]
            max_real_part <- Re(evals[1])
            omega <- abs(Im(evals[1]))
            dominant_period <- if (omega < 1e-10) Inf else 2 * pi / omega
            is_complex <- abs(Im(evals)) > 1e-8
            if (any(is_complex)) {
                osc_idx <- which(is_complex)[which.max(Re(evals[is_complex]))]
                leading_oscillatory_eigenvalue <- evals[osc_idx]
                oscillation_period <- 2 * pi / abs(Im(evals[osc_idx]))
                osc_vector <- single_eigenvector(
                    stability_eigenvectors(evecs, ctx, cols = osc_idx)
                )
            }
        } else {
            evals <- complex()
        }
        vectors <- stability_eigenvectors(evecs, ctx)
    }
    result <- c(
        list(eigenvalues = evals, max_real_part = max_real_part,
             dominant_period = dominant_period,
             oscillation_period = oscillation_period,
             leading_oscillatory_eigenvalue = leading_oscillatory_eigenvalue,
             leading_oscillatory_eigenvector = osc_vector,
             n_active = n_active, leading_eigenvectors = vectors,
             params = ctx$params),
        stability_result(ctx, stability_tol, max_real_part)
    )
    attr(result, "stability_key") <- stability_model_key(ctx$params)
    result
}

#' Analyse the stability of mizer's numerical time step
#'
#' `r lifecycle::badge("experimental")`
#' Computes the eigenvalues \eqn{\mu_i} of the linearised one-step-ahead map at
#' the steady state stored in `params@initial_n`, for a given step size `dt`.
#' These describe how mizer's numerical scheme, rather than the model, behaves
#' near the steady state: the map does not amplify perturbations when the
#' spectral radius \eqn{\max_i|\mu_i|} is less than 1.
#'
#' This is the numerical counterpart of [getStability()], which analyses the
#' model itself and involves no time step at all. Use `getStability()` to ask
#' whether the steady state of the *model* is stable, and this function to ask
#' what mizer's solver does at a particular `dt`. The two can disagree, and that
#' disagreement is the point: the implicit transport solve damps oscillations
#' artificially, so a physically unstable steady state can have a spectral
#' radius below 1 at a large `dt`, and the simulation then sits at a state the
#' model does not actually hold.
#'
#' ## The map that is linearised
#'
#' One step is what [project()] takes with `method = "euler"`: the rates are
#' evaluated at the state at the start of the step, and the resulting transport
#' problem is solved implicitly,
#' \deqn{A(N^t, n_{pp}^t)\,N^{t+1} = S(N^t, n_{pp}^t),}
#' with the same `project_n_loop()` C++ Thomas solver and the same spatial
#' scheme ([second_order_w()]) as the regular dynamics.
#'
#' Because the rates are evaluated at the *input* state, the step is not fully
#' implicit, and the discrete eigenvalues therefore cannot be converted into
#' continuous-time eigenvalues by any exact algebraic relation. That conversion
#' is what [getStability()] avoids by differentiating the rates of change
#' themselves.
#'
#' The resource is advanced by the model's own `resource_dynamics` function, the
#' one [project()] calls. Nothing is substituted for it: the map that is
#' differentiated here reproduces a single `project(method = "euler")` step
#' exactly, which is what makes the spectral radius a statement about mizer's
#' solver rather than about a nearby scheme.
#'
#' With excluded coordinates, this is the projected map
#' \eqn{P_{dynamic} \Phi_{dt}(x_{dynamic}, x_{fixed})}: the full Euler step
#' is evaluated at the fixed inputs, then only dynamic output coordinates are
#' retained. Excluded inputs are restored before the next step. This need not
#' equal discretising the constrained continuous equations: an implicit solve
#' can couple excluded and dynamic rows within a step. A provider declaring
#' structural prescribed inputs remains responsible for honouring that claim
#' in its own equations and step map.
#'
#' @inheritParams getStability
#' @param dt The time step size of the one-step map. Default `1`.
#' @return A named list with the following components:
#'   \describe{
#'     \item{`discrete_eigenvalues`}{Complex vector of the eigenvalues
#'       \eqn{\mu_i} of the one-step map, sorted by decreasing modulus.}
#'     \item{`spectral_radius`}{\eqn{\max_i|\mu_i|}. Less than 1 means the
#'       numerical scheme is stable at this `dt`.}
#'     \item{`stable`}{`TRUE` below `1 - stability_tol`, `FALSE` above
#'       `1 + stability_tol`, and otherwise `NA`. Requires complete coverage
#'       and a passing scoped steady assessment.}
#'     \item{`dt`}{The step size the map was evaluated at.}
#'     \item{`n_active`}{Dimension of the Jacobian.}
#'     \item{`leading_eigenvectors`}{The eigenvectors of the two
#'       largest-modulus eigenvalues, in the same shape as for
#'       [getStability()].}
#'     \item{`params`}{The validated `params` object the analysis was made at.}
#'   }
#'   Scope, coverage and assessment fields are the same as [getStability()].
#'   Missing support gives `NULL` eigenvalues and an inconclusive verdict.
#' @inheritSection getStability Requires smooth dynamics
#' @seealso [getStability()], [findSteadyState()]
#' @export
getDiscreteStability <- function(params,
                                 effort = params@initial_effort,
                                 h = 1e-4, dt = 1,
                                 scope = c("whole", "conditional"),
                                 fixed = character(), t = 0,
                                 stability_tol = 1e-8) {
    assert_that(is.number(h), is.finite(h), h > 0,
                is.number(dt), is.finite(dt), dt > 0,
                is.number(stability_tol), is.finite(stability_tol),
                stability_tol >= 0)
    ctx <- stability_context(params, effort = effort, scope = scope,
                             fixed = fixed, t = t, discrete = TRUE,
                             fn_name = "getDiscreteStability",
                             map_name = "the one-step map")
    evals <- NULL
    spectral_radius <- NA_real_
    vectors <- NULL
    n_active <- NA_integer_
    if (ctx$coverage$complete) {
        n_active <- length(ctx$x0)
        if (n_active) {
            eig <- eigen(stability_jacobian(ctx, stability_step(ctx, dt), h))
            ord <- order(Mod(eig$values), decreasing = TRUE)
            evals <- as.complex(eig$values[ord])
            evecs <- eig$vectors[, ord, drop = FALSE]
            spectral_radius <- Mod(evals[1])
        } else {
            evals <- complex()
            evecs <- matrix(complex(), 0, 0)
        }
        vectors <- stability_eigenvectors(evecs, ctx)
    }
    result <- c(
        list(discrete_eigenvalues = evals, spectral_radius = spectral_radius,
             dt = dt, n_active = n_active, leading_eigenvectors = vectors,
             params = ctx$params),
        stability_result(ctx, stability_tol, spectral_radius, discrete = TRUE)
    )
    attr(result, "stability_key") <- stability_model_key(ctx$params)
    result
}
