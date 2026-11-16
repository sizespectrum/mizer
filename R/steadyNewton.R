# Direct steady-state solver -------------------------------------------------
#
# `solver = "project"` finds a steady state by running the dynamics until they
# stop changing. That only works for steady states that are dynamically stable:
# at an unstable steady state the time-stepping diverges away from the fixed
# point.
#
# `solver = "newton"` instead solves the discrete steady-state equation
# `F(N) = 0` directly with a Newton-type root finder (via the `nleqslv`
# package). Because it solves the algebraic equation rather than following the
# time evolution, it converges to the steady state irrespective of its dynamic
# stability. See the "Steady-State Solution" section of the
# `vignette("numerical_details")` for the equation that is being solved.
#
# This file holds `newton_steady_state()`. Stability analysis lives in
# `stability.R`; the exported solver entry points, `tuneSteadyState()` and
# `findSteadyState()`, live in `steadyState.R`.

#' Solve the steady-state equation with a Newton-type root finder
#'
#' The engine behind `solver = "newton"` in [tuneSteadyState()] and
#' [findSteadyState()]. It solves the steady-state equation `F(N) = 0` with a
#' globalised Newton iteration from the `nleqslv` package instead of running the
#' dynamics to convergence, so it converges to the steady state even when that
#' steady state is dynamically **unstable**, a case in which the time-stepping
#' solver fails because it diverges away from the fixed point.
#'
#' Two things vary between the two callers, and they are the two arguments that
#' say what is held fixed:
#'
#' * `rdd_const` is the per-species reproduction rate to hold constant, or
#'   `NULL` to let the model's own RDI and RDD functions set it at each iterate.
#'   [tuneSteadyState()] passes the rate of the state it was given;
#'   [findSteadyState()] passes `NULL`.
#' * `resource = "solve"` adds the resource densities to the unknowns and
#'   appends the semichemostat steady-state equation to the residual, so the
#'   resource density and the feeding levels it implies are self-consistent even
#'   where consumers are satiated. This requires the default semichemostat
#'   resource dynamics and is what [findSteadyState()] uses. `resource = "fixed"`
#'   leaves the resource out of the system entirely, holding it at the density
#'   in `initial_n_pp`; the caller is then responsible for rebalancing the
#'   capacity afterwards, which is what [tuneSteadyState()] does.
#'
#' Restoring density-dependent reproduction after a `rdd_const` solve is the
#' caller's job too, so that both solvers share one implementation of it.
#'
#' The consumer densities are solved for in log space, which both keeps them
#' positive and conditions the otherwise badly-scaled system. The unknowns are
#' the size classes in the full potential grid support (running from the egg
#' size up to the grid truncation limit `support_top_idx()`), regardless of
#' whether they carry non-zero density in the supplied initial spectra. A
#' smoothed log-abundance penalty floor is applied to all size classes, which
#' automatically handles zero-density and tail classes smoothly and prevents
#' singular Jacobians near that floor. Log-density steps are bounded, and if
#' the residual tolerance is not met the solver retries from a frozen-rate
#' steady spectrum, within the remaining iteration budget. This changes only
#' the starting guess: the active support and nonlinear residual stay the same.
#' After convergence, size classes that remain at or near
#' the floor are set to zero. This allows the solver to automatically discover
#' the support of the steady state. Newton's method converges from any starting
#' point in the *root's* basin of attraction (which is unrelated to the dynamic
#' stability of the steady state), so a reasonable initial guess should still be
#' supplied in `initialN(params)`.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param resource `"solve"` to make the resource densities unknowns of the
#'   system, `"fixed"` to leave the resource out of it entirely. Its `initial_n` is used as
#'   the starting guess for the iteration.
#' @param effort The fishing effort vector, already validated.
#' @param rdd_const Per-species reproduction rate to hold constant, or `NULL`
#'   to let the reproduction dynamics run.
#' @param resource Either `"solve"` (the resource densities are unknowns of the
#'   system) or `"fixed"` (the resource is held at `initial_n_pp`).
#' @param extinction_floor The relative abundance floor below which a species is
#'   considered extinct. Only used when `rdd_const` is `NULL`.
#' @param verbose Whether to trace the solver iterations to the console.
#' @param solver_tol Tolerance on the residual, passed to
#'   [nleqslv::nleqslv()] as both `ftol` and `xtol`.
#' @param maxit Maximum total number of iterations for [nleqslv::nleqslv()],
#'   including a retry if needed.
#' @param jacobian Either `"update"` (the Jacobian is computed once and then
#'   updated cheaply on each iteration) or `"recompute"` (a numerical Jacobian
#'   is computed afresh at every iteration).
#' @param global The globalisation strategy passed to [nleqslv::nleqslv()].
#' @return A list with the new `initial_n` and `initial_n_pp`, the logical
#'   vector of species that went extinct, and the `nleqslv` termination code.
#' @noRd
newton_steady_state <- function(params, effort, rdd_const,
                                resource = c("solve", "fixed"),
                                extinction_floor = 1e-6,
                                verbose = FALSE,
                                solver_tol = 1e-6, maxit = 200,
                                jacobian = c("update", "recompute"),
                                global = "dbldog",
                                rates_fns = projectRateFunctions(params)) {
    require_newton_equations(params)
    resource <- match.arg(resource)
    jacobian <- match.arg(jacobian)
    # nleqslv names these after the Jacobian's origin rather than after what it
    # does; ours says what differs.
    nleqslv_method <- c(update = "Broyden", recompute = "Newton")[[jacobian]]
    if (!requireNamespace("nleqslv", quietly = TRUE)) {
        stop("The Newton solver requires the 'nleqslv' package. ",
             "Install it with install.packages('nleqslv').")
    }
    if (resource == "solve" &&
            sub("^mizer::", "", params@resource_dynamics) !=
                "resource_semichemostat") {
        stop("`solver = \"newton\"` can only solve for the resource with the ",
             "default semichemostat resource dynamics ",
             "(resource_dynamics = 'resource_semichemostat'). Use ",
             "`tuneSteadyState()`, which holds the resource fixed, or ",
             "`solver = \"project\"`.")
    }

    n_other <- params@initial_n_other

    active <- steady_active_set(params, resource = resource)
    if (active$n_fish_active == 0) {
        stop("There is no fish density to solve for: every species is absent ",
             "from `initialN(params)`. The Newton solver holds an absent ",
             "species at zero, so there is nothing left to determine.",
             call. = FALSE)
    }

    # Save initial state for relative floor checks
    N0_initial <- params@initial_n

    # The log-space solve needs a strictly positive, well-scaled start. The
    # second-order schemes can leave isolated zeros inside the support (a
    # negativity-floor artefact), which would make the 1/N-scaled residual
    # overflow. Fill those by log-interpolation from the nonzero neighbours.
    N0 <- positive_initial_guess(params@initial_n, params@w_min_idx,
                                 active$w_top)

    residual_fn <- steady_state_residual(params, rdd_const, n_other, effort,
                                         active, extinction_floor = extinction_floor,
                                         n0 = N0, rates_fns = rates_fns)

    x0 <- log(N0[active$mask])
    if (resource == "solve") {
        x0_n_pp <- as.numeric(params@initial_n_pp[active$mask_pp])
        x0_n_pp[x0_n_pp <= 0] <- params@cc_pp[active$mask_pp][x0_n_pp <= 0]
        x0 <- c(x0, log(x0_n_pp))
    }

    # An almost empty tail with negligible inflow has a nearly constant
    # per-capita residual. Its Newton step can otherwise span millions of log
    # units before the structural floor gives that column a useful derivative.
    control <- list(maxit = maxit, ftol = solver_tol, xtol = solver_tol,
                    stepmax = 10,
                    trace = if (isTRUE(verbose)) 1 else 0)

    sol <- nleqslv::nleqslv(x0, residual_fn, method = nleqslv_method,
                            global = global, control = control)

    residual_size <- function(sol) {
        if (all(is.finite(sol$fvec))) max(abs(sol$fvec)) else Inf
    }
    error <- residual_size(sol)
    remaining <- maxit - sol$iter
    # Small-step termination (code 2) need not satisfy the residual tolerance.
    # Keep an already successful starting guess, and never replace a result
    # with a worse one merely because the retry has a different stopping code.
    if (error > solver_tol && remaining > 0) {
        seed <- tryCatch(newton_frozen_guess(params, active, N0, effort,
                                             rates_fns = rates_fns),
                         mizer_steady_nonconvergence = function(e) NULL)
        if (!is.null(seed)) {
            retry_x0 <- x0
            retry_x0[seq_len(active$n_fish_active)] <- log(seed[active$mask])
            control$maxit <- remaining
            retry <- nleqslv::nleqslv(retry_x0, residual_fn,
                                      method = nleqslv_method,
                                      global = global, control = control)
            if (residual_size(retry) < error) {
                sol <- retry
                error <- residual_size(sol)
            }
        }
    }

    unpacked <- active$unpack(sol$x)
    N <- unpacked$N
    n_pp <- unpacked$n_pp
    is_extinct <- rep(FALSE, nrow(N))
    names(is_extinct) <- rownames(N)

    # Check for extinctions if using the relative floor
    if (is.null(rdd_const) && !is.null(extinction_floor) && extinction_floor > 0) {
        extinct_threshold <- extinction_floor * 1.01
        for (i in seq_len(nrow(N))) {
            lo <- params@w_min_idx[i]
            # A species is considered extinct if its recruitment (the egg class)
            # has dropped to the relative abundance floor. Checking any() would
            # produce false positives when the tail naturally shortens.
            if (N0_initial[i, lo] > 0 &&
                N[i, lo] / N0_initial[i, lo] <= extinct_threshold) {
                is_extinct[i] <- TRUE
                N[i, ] <- 0
            }
        }
    }

    # Zero out any tail densities that ended up at or near the structural penalty floor
    for (i in seq_len(nrow(N))) {
        if (!is_extinct[i]) {
            floor_val <- max(N0_initial[i, ]) * 1.1 * 1e-15
            N[i, N[i, ] < floor_val] <- 0
        }
    }

    list(n = N, n_pp = n_pp, extinct = is_extinct, absent = active$absent,
         termcd = sol$termcd, message = sol$message, residual = error)
}

# Rate and verification callbacks cannot replace Newton's residual operator.
require_newton_equations <- function(params) {
    if (!builtin_consumer_projection(params)) {
        stop("The Newton solver requires the built-in consumer transport ",
             "equations. Custom project() or project_simple() dynamics need ",
             "a supported equation operator.", call. = FALSE)
    }
    invisible(TRUE)
}

#' A frozen-rate starting spectrum for a stalled Newton solve
#'
#' Retains the model's transport scheme, diffusion and computational support.
#' In particular, zero growth does not remove a cell from the active set:
#' diffusion or growth at a later iterate may still populate it. Numerically
#' empty cells start just below the existing structural floor, where the
#' log-abundance penalty has a nonzero derivative, instead of extrapolating a
#' substantial density into a tail with almost no inflow.
#'
#' @param params The model at the original starting state.
#' @param active The active-set list from `steady_active_set()`.
#' @param n0 The positive original starting spectrum.
#' @param effort The fishing effort vector.
#' @return A density matrix, strictly positive on the active set.
#' @noRd
newton_frozen_guess <- function(params, active, n0, effort,
                                rates_fns = projectRateFunctions(params)) {
    rates <- rates_fns$Rates(params, n = params@initial_n,
                             n_pp = params@initial_n_pp,
                             n_other = params@initial_n_other,
                             t = 0, effort = effort, rates_fns = rates_fns)
    eggs <- n0[cbind(seq_len(nrow(n0)), params@w_min_idx)]
    seed <- get_steady_state_n(params, rates$e_growth, rates$mort,
                               rates$diffusion, eggs)
    floor <- apply(n0, 1, max) * 1e-15
    seed <- pmax(seed, floor * exp(-0.5))
    seed[!active$mask] <- 0
    seed
}

#' The set of size classes solved for by the Newton solver
#'
#' Builds the logical mask of the (species x size) density matrix that the
#' direct solver treats as unknowns. For each species the unknowns run from the
#' egg size `w_min_idx` up to the grid truncation `support_top_idx()`, which is
#' the whole of the size range the dynamics can put density into.
#'
#' The top is the truncation rather than the highest class that currently
#' carries density, because the tail of the solution need not lie under the tail
#' of the starting guess. The growth rate cannot locate the top either: the main
#' reason fish grow beyond `w_repro_max` is diffusion, and the diffusion rate
#' only grows with `w`, never vanishing. Above the truncation the dynamics hold
#' the density at exactly zero, so those classes can never be unknowns; that is
#' what keeps the band of `log(0)` unknowns that would make the Jacobian
#' singular out of the system.
#'
#' Isolated interior zeros (negativity-floor artefacts of the second-order
#' schemes) below the top are kept in the mask and repaired by
#' `positive_initial_guess()`.
#'
#' A species whose row is zero throughout is left out of the mask altogether. It
#' has no density anywhere, so there is nothing to interpolate a starting guess
#' from and `log()` of its row would put `-Inf` into the starting vector. Zero is
#' in any case a steady state of its own equation, so the species is simply held
#' there; recolonisation is not something this solver can discover.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param resource `"solve"` to make the resource densities unknowns of the
#'   system, `"fixed"` to leave them out of it.
#' @return A list with the logical matrix `mask`, the logical vector `absent`
#'   naming the species that were left out, and a function `unpack(x)` that maps
#'   a vector of log-densities (in `mask` order) to the full density matrix.
#' @noRd
steady_active_set <- function(params, resource = c("solve", "fixed")) {
    resource <- match.arg(resource)
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)

    # The grid truncation caps the support: the dynamics hold the density at zero
    # above this, so it can never be an unknown (see support_top_idx()).
    grid_top <- support_top_idx(params)
    n <- params@initial_n
    w_top <- integer(no_sp)
    absent <- rep(FALSE, no_sp)
    names(absent) <- rownames(n)
    mask <- matrix(FALSE, nrow = no_sp, ncol = no_w)
    for (i in seq_len(no_sp)) {
        lo <- params@w_min_idx[i]
        # Always solve up to the grid truncation limit
        w_top[i] <- grid_top[i]
        # A species that is absent to begin with stays absent: it contributes no
        # unknowns, and its row of the solution is left at zero.
        if (all(n[i, ] <= 0)) {
            absent[i] <- TRUE
            next
        }
        mask[i, lo:w_top[i]] <- TRUE
    }
    
    # With the resource held fixed there are no resource unknowns, so the mask
    # is empty and `unpack()` hands back the supplied density untouched.
    if (resource == "solve") {
        mask_pp <- as.logical(params@cc_pp > 0)
    } else {
        mask_pp <- rep(FALSE, length(params@cc_pp))
    }
    n_fish_active <- sum(mask)

    dn <- dimnames(params@initial_n)
    unpack <- function(x) {
        N <- matrix(0, nrow = no_sp, ncol = no_w, dimnames = dn)
        N[mask] <- exp(x[1:n_fish_active])
        n_pp <- params@initial_n_pp
        if (any(mask_pp)) {
            n_pp[mask_pp] <- exp(x[(n_fish_active + 1):length(x)])
        }
        list(N = N, n_pp = n_pp)
    }
    list(mask = mask, w_top = w_top, mask_pp = mask_pp, absent = absent,
         n_fish_active = n_fish_active, unpack = unpack)
}

#' Strictly positive starting guess for the log-space solve
#'
#' The Newton solve uses log-densities as unknowns and scales the residual by
#' `1/N`, so it needs a starting guess that is strictly positive and reasonably
#' scaled on every active size class. The second-order schemes can leave isolated
#' zeros inside the support (a negativity-floor artefact at a reconstructed
#' over/undershoot); a plain `log()` of such a guess is `-Inf`, and flooring it
#' to a tiny constant instead makes the `1/N`-scaled residual overflow. We repair
#' the guess by interpolating `log(N)` linearly across the gaps from the nonzero
#' neighbours (geometric interpolation on the logarithmic size grid), which is
#' both strictly positive and well scaled. Edges with no nonzero neighbour on one
#' side are extrapolated flat.
#'
#' @param N The current density matrix (species x size).
#' @param w_min_idx Per-species egg-size index.
#' @param w_top Per-species support top from `steady_active_set()`.
#' @return A density matrix that is strictly positive on the active range.
#' @noRd
positive_initial_guess <- function(N, w_min_idx, w_top) {
    for (i in seq_len(nrow(N))) {
        rng <- w_min_idx[i]:w_top[i]
        v <- N[i, rng]
        pos <- v > 0
        if (!any(pos)) next
        floor_val <- max(v) * 1e-20
        if (!all(pos)) {
            idx <- seq_along(v)
            v[!pos] <- exp(stats::approx(idx[pos], log(v[pos]),
                                         xout = idx[!pos], rule = 2)$y)
        }
        v[is.na(v) | v <= 0] <- floor_val
        N[i, rng] <- v
    }
    N
}

#' Unscaled residual of the discrete steady-state equation
#'
#' Evaluates the growth, mortality and diffusion rates at the given state,
#' assembles the steady-state tridiagonal coefficients with
#' `get_transport_coefs()` at `dt = 1` and returns
#' \deqn{F_j = a_j N_{j-1} + b_j N_j + c_j N_{j+1} - S_j.}
#' Since `S = N + recruitment source`, the `+N` cancels the backward-Euler
#' `N^t` term, leaving exactly the steady-state equation
#' \eqn{\tilde A N_{j-1} + \tilde B N_j + \tilde C N_{j+1} - \tilde S = 0}.
#'
#' Because the backward-Euler coefficients are `A = I - dt L`, this residual is
#' \eqn{-dt\,dN/dt} exactly, so at `dt = 1` its negative is the instantaneous
#' rate of change of the density. That is what makes it usable both as the
#' function whose root the Newton solver seeks and as the diagnostic
#' [getSteadyResidual()] reports.
#'
#' This is the single implementation of the transport residual; both callers go
#' through it.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param n Consumer densities (species x size) to evaluate the residual at.
#' @param n_pp Resource density to evaluate the rates at.
#' @param n_other Abundances of other components.
#' @param effort The fishing effort vector.
#' @param rdd Per-species reproduction rate to use. If `NULL`, it is calculated
#'   from the model's own RDI and RDD functions at the given state.
#' @param rates_fns The rate functions, from `projectRateFunctions()`. Passed in
#'   so that a caller evaluating the residual many times need only look them up
#'   once.
#' @param flux_limiter The flux limiter scheme, from `flux_limiter_scheme()`.
#' @param rates The rates at this state, either a full rate evaluation or the
#'   `mizer_rates_subset()` with targets `EGrowth`, `Mort` and `Diffusion`. Pass
#'   them in when the caller has computed them already; otherwise they are
#'   computed here.
#' @return The unscaled residual matrix (species x size).
#' @noRd
consumer_residual <- function(params, n, n_pp, n_other, effort, rdd = NULL,
                              rates_fns = projectRateFunctions(params),
                              flux_limiter = flux_limiter_scheme(params),
                              rates = NULL) {
    no_w <- length(params@w)

    r <- rates
    if (is.null(r)) {
        r <- rates_fns$Rates(params, n = n, n_pp = n_pp, n_other = n_other,
                             t = 0, effort = effort, rates_fns = rates_fns)
    }

    if (is.null(rdd)) {
        rdd <- if (!is.null(r$rdd)) {
            r$rdd
        } else {
            state_rdd(params, n = n, n_pp = n_pp, n_other = n_other,
                      rates = r, effort = effort, rates_fns = rates_fns)
        }
    }

    coefs <- get_transport_coefs(params, n = n, g = r$e_growth,
                                 mu = r$mort, dt = 1,
                                 recruitment_flux = rdd,
                                 d = r$diffusion,
                                 flux_limiter = flux_limiter)

    Nm <- cbind(0, n[, -no_w, drop = FALSE])   # N_{j-1}
    Np <- cbind(n[, -1, drop = FALSE], 0)      # N_{j+1}
    coefs$a * Nm + coefs$b * n + coefs$c * Np - coefs$S
}

#' Reproduction rate at a given state
#'
#' Runs the model's own RDI and RDD functions at the state `(n, n_pp, n_other)`,
#' going through the extension dispatch when the model carries extensions. This
#' is the reproduction rate that the dynamics would use at that state, as
#' opposed to a rate held fixed at its steady-state value. It is shared by
#' `consumer_residual()` and the stability analyses, which have to agree on it.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param n Consumer densities (species x size).
#' @param n_pp Resource density.
#' @param n_other Abundances of other components.
#' @param rates The rates at this state, at least `e_repro`, `e_growth`, `mort`
#'   and `diffusion`, as returned by `mizer_rates_subset()`.
#' @return The per-species density-dependent reproduction rate.
#' @noRd
state_rdd <- function(params, n, n_pp, n_other, rates,
                      effort = params@initial_effort,
                      rates_fns = projectRateFunctions(params)) {
    if (!is.null(rates$rdd)) return(rates$rdd)
    rates_fns$Rates(params, n = n, n_pp = n_pp, n_other = n_other,
                    effort = effort, t = 0, rates_fns = rates_fns)$rdd
}

#' Residual of the discrete steady-state equation
#'
#' Returns a closure `f(x)` suitable for [nleqslv::nleqslv()]. The argument `x`
#' is the vector of log-densities of the active size classes (see
#' `steady_active_set()`). The closure rebuilds the full density matrix,
#' substitutes the analytic resource steady state, and calls
#' `consumer_residual()` for the steady-state residual
#' \eqn{F_j = a_j N_{j-1} + b_j N_j + c_j N_{j+1} - S_j}, to which it adds the
#' scaling, masking and penalty floor that the root finder needs.
#' The residual is divided by `N`, turning it into a per-capita rate of change
#' that is dimensionless and O(1) across the many orders of magnitude spanned by
#' the densities — the natural scaling to pair with the log-space unknowns.
#'
#' The flux-limiter weight is recomputed fresh on every evaluation, so the van
#' Leer / centred nonlinearity is part of the root and the solution matches what
#' [project()] converges to on the same `params`.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param rdd_const Per-species reproduction rate held constant during the solve.
#' @param n_other Abundances of other components (held constant).
#' @param effort The fishing effort vector.
#' @param active The active-set list from `steady_active_set()`.
#' @return A function of the packed log-density vector returning the packed
#'   scaled residual.
#' @noRd
steady_state_residual <- function(params, rdd_const, n_other, effort, active,
                                  extinction_floor = 1e-6, n0 = NULL,
                                  rates_fns = projectRateFunctions(params)) {
    no_w <- length(params@w)
    mask <- active$mask
    mask_pp <- active$mask_pp
    n_fish_active <- active$n_fish_active
    flux_limiter <- flux_limiter_scheme(params)
    force(rates_fns)
    if (is.null(n0)) {
        x0_initial <- params@initial_n
    } else {
        x0_initial <- n0
    }

    # A floor is always active to handle zero-abundance tail classes smoothly
    # and prevent singular Jacobians.
    x0 <- log(x0_initial[mask])
    support_floor <- matrix(0, nrow = nrow(x0_initial), ncol = no_w)
    for (i in seq_len(nrow(x0_initial))) {
        support_floor[i, ] <- log(max(x0_initial[i, ])) + log(1e-15)
    }
    x_floor <- support_floor[mask]

    if (is.null(rdd_const) && !is.null(extinction_floor) && extinction_floor > 0) {
        ext_floor <- x0 + log(extinction_floor)
        x_floor <- pmax(x_floor, ext_floor)
    }

    function(x) {
        # Let nleqslv reject an invalid trial and backtrack before it reaches
        # rate functions or the limiter. Do not catch rate-function errors:
        # those may describe a genuine error in a user's custom dynamics.
        if (any(!is.finite(x))) {
            return(rep(Inf, length(x)))
        }
        unpacked <- active$unpack(x)
        N <- unpacked$N
        n_pp <- unpacked$n_pp
        if (any(!is.finite(N)) || any(N[mask] <= 0) ||
                any(!is.finite(n_pp)) || any(n_pp[mask_pp] <= 0)) {
            return(rep(Inf, length(x)))
        }

        rates <- rates_fns$Rates(params, n = N, n_pp = n_pp,
                                 n_other = n_other, t = 0, effort = effort,
                                 rates_fns = rates_fns)
        recruitment <- if (is.null(rdd_const)) rates$rdd else rdd_const
        res <- consumer_residual(params, n = N, n_pp = n_pp, n_other = n_other,
                                 effort = effort, rdd = recruitment,
                                 rates_fns = rates_fns, rates = rates,
                                 flux_limiter = flux_limiter)

        # Scale to a per-capita rate of change (N > 0 on the active set).
        if (is.null(rdd_const)) {
            r_ss <- (res / x0_initial)[mask]
        } else {
            r_ss <- (res / N)[mask]
        }

        if (!is.null(x_floor)) {
            y <- x[1:n_fish_active] - x_floor
            delta <- 0.001
            penalty <- 0.5 * (y - sqrt(y^2 + delta^2))
            r_ss <- r_ss + penalty
        }
        
        # Resource residual, only when the resource is part of the system.
        if (!any(mask_pp)) {
            return(r_ss)
        }
        resource_mort <- as.numeric(rates$resource_mort)
        r_pp_ss <- (params@rr_pp * params@cc_pp / n_pp -
                        (params@rr_pp + resource_mort))[mask_pp]

        c(r_ss, as.numeric(r_pp_ss))
    }
}
