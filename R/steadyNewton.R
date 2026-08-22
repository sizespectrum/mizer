# Direct steady-state solver -------------------------------------------------
#
# `steady()` finds a steady state by running the dynamics until they stop
# changing. That only works for steady states that are dynamically stable: at an
# unstable steady state the time-stepping diverges away from the fixed point.
#
# `steadyNewton()` instead solves the discrete steady-state equation
# `F(N) = 0` directly with a Newton-type root finder (via the `nleqslv`
# package). Because it solves the algebraic equation rather than following the
# time evolution, it converges to the steady state irrespective of its dynamic
# stability. See the "Steady-State Solution" section of the
# `vignette("numerical_details")` for the equation that is being solved.

#' Find a steady state by directly solving the steady-state equation
#'
#' `r lifecycle::badge("experimental")`
#' This is an alternative to [steady()] that finds the steady state by solving
#' the steady-state equation `F(N) = 0` with a Newton-type root finder instead
#' of running the dynamics to convergence. The advantage is that it converges
#' to the steady state even when that steady state is dynamically **unstable**,
#' a case in which [steady()] fails because the time-stepping diverges away from
#' the fixed point.
#'
#' By default, or when `reproduction = "fixed"`, the function holds the
#' reproduction rate (RDD) constant while solving for the consumer spectra,
#' substitutes the analytic steady state of the resource, and keeps any other
#' components constant. After the spectra have been found it restores
#' density-dependent Beverton-Holt reproduction with [setBevertonHolt()],
#' honouring the `preserve` argument exactly as [steady()] does.
#' If `reproduction = "dynamic"`, the reproduction dynamics are run dynamically,
#' meaning the reproduction rate varies during the solve and the reproduction
#' parameters are not adjusted.
#'
#' The consumer densities are solved for in log space, which both keeps them
#' positive and conditions the otherwise badly-scaled system. The unknowns are
#' the size classes in the full potential grid support (running from the egg
#' size up to the grid truncation limit `support_top_idx()`), regardless of
#' whether they carry non-zero density in the supplied initial spectra. A
#' smoothed log-abundance penalty floor is applied to all size classes, which
#' automatically handles zero-density and tail classes smoothly and prevents
#' singular Jacobians. After convergence, size classes that remain at or near
#' the floor are set to zero. This allows the solver to automatically discover
#' the support of the steady state. The nonlinear system is solved with a
#' globalised Newton iteration from the `nleqslv` package, starting from the
#' current `initial_n`. Newton's method converges from any starting point in the
#' *root's* basin of attraction (which is unrelated to the dynamic stability of
#' the steady state), so a reasonable initial guess should still be supplied in
#' `initialN(params)` — for example the spectra from a nearby stable
#' parameterisation, or the (diverging) output of [steady()].
#'
#' The solver respects the active transport scheme: if the experimental
#' second-order scheme is enabled (see [second_order_w()]) it solves the
#' steady-state equation of that scheme. With the van Leer reconstruction the
#' residual is only Lipschitz, so the iteration converges to a fixed point of the
#' dynamics but not to machine precision. The unlimited `"centred"`
#' reconstruction admits an undamped odd-even mode at a steady state with no
#' physical diffusion, giving an ill-conditioned steady-state Jacobian for which
#' the solver is not expected to converge.
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
#' Only the default semichemostat resource dynamics
#' (`resource_dynamics = "resource_semichemostat"`) are currently supported,
#' because the solver substitutes the analytic resource steady state. For other
#' resource dynamics the function stops with an error.
#'
#' @param params A \linkS4class{MizerParams} object. Its `initial_n` is used as
#'   the starting guess for the iteration.
#' @param effort The fishing effort. By default the initial effort stored in
#'   `params`.
#' @param preserve `r lifecycle::badge("experimental")`
#'   Specifies whether the `reproduction_level` should be preserved (default)
#'   or the maximum reproduction rate `R_max` or the reproductive efficiency
#'   `erepro`. See [setBevertonHolt()] for an explanation of the
#'   `reproduction_level`. This argument is ignored when `reproduction = "dynamic"`.
#' @param reproduction `r lifecycle::badge("experimental")`
#'   If `"fixed"`, the reproduction rate (RDD) is held constant at the initial
#'   value. If `"dynamic"`, the reproduction dynamics are run dynamically and the
#'   reproduction parameters are not adjusted. Default is `"fixed"`.
#' @param extinction_floor `r lifecycle::badge("experimental")`
#'   The relative abundance floor below which a species is considered extinct.
#'   Only used when `reproduction = "dynamic"`. Default is 1e-6.
#' @param verbose If `TRUE` then the solver iterations will be traced and
#'   printed to the console. Default is `FALSE`.
#' @param tol Convergence tolerance passed to [nleqslv::nleqslv()] (both the
#'   function-value tolerance `ftol` and the step tolerance `xtol`).
#' @param maxit Maximum number of iterations for [nleqslv::nleqslv()].
#' @param method The [nleqslv::nleqslv()] method, either `"Newton"` (with a
#'   numerical Jacobian calculated at each iteration) or `"Broyden"` (which
#'   calculates the full Jacobian only once and then only updates it on each
#' @param global The globalisation strategy passed to [nleqslv::nleqslv()].
#'   The default `"dbldog"` (double dogleg) is a robust trust-region method.
#' @param info_level Controls the amount of information messages and warnings
#'   that are shown. Higher levels lead to more messages, `info_level = 0`
#'   gives silence. The default is taken from the `mizer_info_level` option,
#'   see [default_info_level()].
#' @param ... Unused.
#' @return A \linkS4class{MizerParams} object with the initial state set to the
#'   steady state.
#' @seealso [steady()], [steadySingleSpecies()], [getStability()],
#'   [isSteady()], [getSteadyResidual()]
#' @export
#' @examples
#' \donttest{
#' params <- steadyNewton(NS_params)
#' plotSpectra(params)
#' }
steadyNewton <- function(params, ...) {
    UseMethod("steadyNewton")
}

#' @rdname steadyNewton
#' @export
steadyNewton.MizerParams <- function(params,
                                     effort = params@initial_effort,
                                     preserve = c("reproduction_level",
                                                  "erepro", "R_max"),
                                     reproduction = c("fixed", "dynamic"),
                                     extinction_floor = 1e-6,
                                     verbose = FALSE,
                                     tol = 1e-6, maxit = 200,
                                     method = c("Broyden", "Newton"),
                                     global = "dbldog",
                                     info_level = default_info_level(), ...) {
    with_info_level(info_level = info_level, {
    reproduction <- match.arg(reproduction)
    if (reproduction == "fixed") {
        preserve <- match.arg(preserve)
    }
    method <- match.arg(method)
    if (!requireNamespace("nleqslv", quietly = TRUE)) {
        stop("steadyNewton() requires the 'nleqslv' package. ",
             "Install it with install.packages('nleqslv').")
    }
    params <- validParams(params)
    if (params@resource_dynamics != "resource_semichemostat") {
        stop("steadyNewton() currently only supports the default semichemostat ",
             "resource dynamics (resource_dynamics = 'resource_semichemostat'). ",
             "Use steady() for other resource dynamics.")
    }
    effort <- validEffortVector(effort, params = params)
    params@initial_effort <- effort

    if (params@rates_funcs$RDD == "BevertonHoltRDD" && reproduction == "fixed") {
        old_reproduction_level <- reproduction_level(params)
        old_R_max <- params@species_params$R_max
        old_erepro <- params@species_params$erepro
    }

    if (reproduction == "dynamic") {
        rdd_const <- NULL
    } else {
        rdd_const <- getRDD(params)
    }
    n_other <- params@initial_n_other

    active <- steady_active_set(params)

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
                                         N0 = N0)

    x0_fish <- log(N0[active$mask])
    
    x0_n_pp <- as.numeric(params@initial_n_pp[active$mask_pp])
    x0_n_pp[x0_n_pp <= 0] <- params@cc_pp[active$mask_pp][x0_n_pp <= 0]
    x0_pp <- log(x0_n_pp)
    
    x0 <- c(x0_fish, x0_pp)

    if (isTRUE(verbose)) {
        control <- list(maxit = maxit, ftol = tol, xtol = tol, trace = 1)
    } else {
        control <- list(maxit = maxit, ftol = tol, xtol = tol, trace = 0)
    }

    sol <- nleqslv::nleqslv(x0, residual_fn, method = method,
                            global = global, control = control)

    if (sol$termcd > 2) {
        warning("steadyNewton() did not converge (nleqslv termination code ",
                sol$termcd, ": ", sol$message,
                "). Returning the best iterate found.", call. = FALSE)
    }

    unpacked <- active$unpack(sol$x)
    N <- unpacked$N
    n_pp <- unpacked$n_pp
    is_extinct <- rep(FALSE, nrow(N))
    names(is_extinct) <- rownames(N)

    # Check for extinctions if using the relative floor
    if (reproduction == "dynamic" && !is.null(extinction_floor) && extinction_floor > 0) {
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
        if (any(is_extinct)) {
            warning("The following species went extinct and were set to zero: ",
                    paste(names(is_extinct)[is_extinct], collapse = ", "))
        }
    }

    # Zero out any tail densities that ended up at or near the structural penalty floor
    for (i in seq_len(nrow(N))) {
        if (!is_extinct[i]) {
            floor_val <- max(N0_initial[i, ]) * 1.1 * 1e-15
            N[i, N[i, ] < floor_val] <- 0
        }
    }

    n_pp <- resource_steady_semichemostat(params, N, n_other)

    params@initial_n[] <- N
    params@initial_n_pp[] <- n_pp

    # Restore density-dependent reproduction, just as steady() does.
    if (params@rates_funcs$RDD == "BevertonHoltRDD" && reproduction == "fixed") {
        if (preserve == "reproduction_level") {
            reproduction_level(params) <- old_reproduction_level
        } else if (preserve == "R_max") {
            params <- setBevertonHolt(params, R_max = old_R_max)
        } else if (preserve == "erepro") {
            params <- setBevertonHolt(params, erepro = old_erepro)
        }
    }

    params@time_modified <- lubridate::now()

    # Report how close to a fixed point the solver actually got. This is the
    # honest measure of the result: with the van Leer reconstruction the
    # residual is only Lipschitz and cannot reach machine precision, and the
    # analytic resource substitution used during the solve is not quite
    # self-consistent, so the number is not always as small as `tol` suggests.
    residual <- tryCatch(steady_biomass_drift(params, effort = effort),
                         error = function(e) NA_real_)
    if (is.finite(residual)) {
        signal_info("convergence", paste0(
            "The biomasses of the solution change at up to ",
            signif(residual, 2), " per year."),
            level = 3, unhandled = "show")
    }

    params
    })
}

#' The set of size classes solved for by steadyNewton()
#'
#' Builds the logical mask of the (species x size) density matrix that the
#' direct solver treats as unknowns. For each species the unknowns run from the
#' egg size `w_min_idx` up to the highest class that actually carries density in
#' the current `initial_n`, capped at the grid truncation `support_top_idx()`.
#'
#' The support is read off the abundances rather than from `w_max` on purpose.
#' Users routinely set `w_max` far larger than necessary (so the grid need not
#' change when a parameter change produces larger fish), which would otherwise
#' put a whole band of structurally-zero classes — `log(0)` unknowns that make
#' the Jacobian singular — into the system. The growth rate alone cannot locate
#' the top either: the main reason fish grow beyond `w_repro_max` is diffusion,
#' and the diffusion rate only grows with `w`, never vanishing. The abundance is
#' therefore the only reliable indicator of where the (possibly diffusion-fed)
#' tail has died away. The distinction is clean: without diffusion the classes
#' above where growth stops receive no inflow and are *exactly* zero, while with
#' diffusion the tail is genuinely non-zero up to where it decays away (or to the
#' grid truncation).
#'
#' Isolated interior zeros (negativity-floor artefacts of the second-order
#' schemes) below the top are kept in the mask and repaired by
#' `positive_initial_guess()`; only a trailing run of zeros is dropped.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @return A list with the logical matrix `mask` and a function `unpack(x)` that
#'   maps a vector of log-densities (in `mask` order) to the full density
#'   matrix.
#' @noRd
steady_active_set <- function(params) {
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)

    # The grid truncation caps the support: the dynamics hold the density at zero
    # above this, so it can never be an unknown (see support_top_idx()).
    grid_top <- support_top_idx(params)
    n <- params@initial_n
    w_top <- integer(no_sp)
    mask <- matrix(FALSE, nrow = no_sp, ncol = no_w)
    for (i in seq_len(no_sp)) {
        lo <- params@w_min_idx[i]
        # Always solve up to the grid truncation limit
        w_top[i] <- grid_top[i]
        mask[i, lo:w_top[i]] <- TRUE
    }
    
    mask_pp <- as.logical(params@cc_pp > 0)
    n_fish_active <- sum(mask)

    dn <- dimnames(params@initial_n)
    unpack <- function(x) {
        N <- matrix(0, nrow = no_sp, ncol = no_w, dimnames = dn)
        N[mask] <- exp(x[1:n_fish_active])
        n_pp <- params@initial_n_pp
        n_pp[mask_pp] <- exp(x[(n_fish_active + 1):length(x)])
        list(N = N, n_pp = n_pp)
    }
    list(mask = mask, w_top = w_top, mask_pp = mask_pp, 
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

#' Relative finite-difference step scale
#'
#' Returns the scale that [getStability()] multiplies by `h` to get the
#' finite-difference step for each state variable: the variable's own magnitude
#' where that is nonzero, and otherwise the local scale of the spectrum as
#' supplied by `positive_initial_guess()` (log-interpolated from the nonzero
#' neighbours). Cells whose whole row is zero get no local scale from
#' `positive_initial_guess()`; they fall back to `.Machine$double.eps`, which is
#' harmless because the one-step map is exactly linear about a zero baseline.
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

#' Analytic semichemostat resource steady state
#'
#' Returns the resource number density that is in equilibrium with the consumer
#' densities `N`. For the semichemostat resource the equilibrium is
#' `n_pp* = rr_pp * cc_pp / (rr_pp + mu_R(N))`, where the resource predation
#' mortality `mu_R` depends only on the consumer densities (not on the resource
#' density itself), so the substitution is exact. This is the `n_steady` of
#' [resource_semichemostat()].
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param N Consumer densities (species x size).
#' @param n_other Abundances of other components.
#' @return The steady-state resource number density vector.
#' @noRd
resource_steady_semichemostat <- function(params, N, n_other) {
    n_pp <- params@initial_n_pp
    for (i in 1:8) {
        mu_R <- as.numeric(getResourceMort(params, n = N,
                                           n_pp = n_pp,
                                           n_other = n_other, t = 0))
        mur <- params@rr_pp + mu_R
        n_pp_new <- params@rr_pp * params@cc_pp / mur
        # Where both rate and mortality vanish the steady state is undetermined;
        # keep the initial value.
        sel <- !is.finite(n_pp_new)
        n_pp_new[sel] <- params@initial_n_pp[sel]
        n_pp <- n_pp_new
    }
    n_pp
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
#' function whose root [steadyNewton()] seeks and as the diagnostic
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
#' @param rates The rates at this state, as returned by `mizer_rates_subset()`
#'   with targets `EGrowth`, `Mort` and `Diffusion`. Pass them in when the
#'   caller has computed them already; otherwise they are computed here.
#' @return The unscaled residual matrix (species x size).
#' @noRd
consumer_residual <- function(params, n, n_pp, n_other, effort, rdd = NULL,
                              rates_fns = projectRateFunctions(params),
                              flux_limiter = flux_limiter_scheme(params),
                              rates = NULL) {
    no_w <- length(params@w)

    r <- rates
    if (is.null(r)) {
        r <- mizer_rates_subset(params, n = n, n_pp = n_pp, n_other = n_other,
                                t = 0, effort = effort, rates_fns = rates_fns,
                                targets = c("EGrowth", "Mort", "Diffusion"))
    }

    if (is.null(rdd)) {
        rdd <- state_rdd(params, n = n, n_pp = n_pp, n_other = n_other,
                         rates = r)
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
state_rdd <- function(params, n, n_pp, n_other, rates) {
    if (usesExtensionDispatch(params)) {
        rdi <- projectRDI(params, n = n, n_pp = n_pp, n_other = n_other, t = 0,
                          e_repro = rates$e_repro, e_growth = rates$e_growth,
                          mort = rates$mort, diffusion = rates$diffusion)
        projectRDD(params, rdi = rdi, species_params = params@species_params,
                   t = 0)
    } else {
        f_rdi <- get(params@rates_funcs$RDI)
        rdi <- f_rdi(params, n = n, n_pp = n_pp, n_other = n_other, t = 0,
                     e_repro = rates$e_repro, e_growth = rates$e_growth,
                     mort = rates$mort, diffusion = rates$diffusion)
        f_rdd <- get(params@rates_funcs$RDD)
        f_rdd(rdi = rdi, species_params = params@species_params,
              params = params, t = 0)
    }
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
                                  extinction_floor = 1e-6, N0 = NULL) {
    no_w <- length(params@w)
    mask <- active$mask
    mask_pp <- active$mask_pp
    n_fish_active <- active$n_fish_active
    flux_limiter <- flux_limiter_scheme(params)
    rates_fns <- projectRateFunctions(params)
    
    if (is.null(N0)) {
        x0_initial <- params@initial_n
    } else {
        x0_initial <- N0
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
        unpacked <- active$unpack(x)
        N <- unpacked$N
        n_pp <- unpacked$n_pp

        res <- consumer_residual(params, n = N, n_pp = n_pp, n_other = n_other,
                                 effort = effort, rdd = rdd_const,
                                 rates_fns = rates_fns,
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
        
        # Resource residual
        mu_R <- as.numeric(getResourceMort(params, n = N,
                                           n_pp = n_pp,
                                           n_other = n_other, t = 0))
        r_pp_ss <- (params@rr_pp * params@cc_pp / n_pp - (params@rr_pp + mu_R))[mask_pp]

        c(r_ss, as.numeric(r_pp_ss))
    }
}

# Stability analysis ----------------------------------------------------------

#' Shared scaffolding for the stability analyses
#'
#' [getStability()] and [getDiscreteStability()] differ only in the map they
#' linearise: the continuous rates of change in the one case, mizer's numerical
#' one-step map in the other. Everything around that map is the same, and this
#' function assembles it: the validated `params`, the active set, the packed
#' state vector with its finite-difference step scales, and the closures that
#' unpack a packed state, name a cell and check a result for non-finite values.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param reproduction Whether reproduction is held fixed or run dynamically.
#' @param effort The fishing effort to use.
#' @param fn_name The name of the calling function, for error messages.
#' @param map_name What the linearised map is called, for error messages.
#' @return A list holding the validated `params`, the state and the closures
#'   described above.
#' @noRd
stability_context <- function(params,
                              reproduction = c("fixed", "dynamic"),
                              effort = params@initial_effort,
                              fn_name = "getStability",
                              map_name = "the rates of change") {
    reproduction <- match.arg(reproduction)
    params <- validParams(params)
    effort <- validEffortVector(effort, params = params)
    params@initial_effort <- effort

    # The whole calculation linearises the dynamics *at* the stored state. If
    # that state is not a fixed point, the eigenvalues describe the
    # neighbourhood of a point the model is not sitting at, and the verdict on
    # stability means nothing.
    warn_if_not_steady(params, paste(
        "The eigenvalues below describe the dynamics near a fixed point, so on",
        "a state that is not one they are not meaningful. Use `steadyNewton()`",
        "first."))

    rdd_const <- if (reproduction == "dynamic") NULL else getRDD(params)
    n_other      <- params@initial_n_other
    active       <- steady_active_set(params)
    active_idx   <- which(active$mask)
    flux_limiter <- flux_limiter_scheme(params)
    rates_fns    <- projectRateFunctions(params)

    N_ss   <- params@initial_n
    npp_ss <- params@initial_n_pp
    N_vec  <- N_ss[active$mask]
    n_fish_active <- length(N_vec)

    # Finite-difference step sizes. The step is relative, `h * N`, so a cell at
    # exactly zero would get a zero step: the difference quotient then drops
    # the cell from the Jacobian and puts a spurious zero eigenvalue in place of
    # its true decay rate. Such cells do occur: an isolated negativity-floor
    # artefact of the second-order schemes, or a tail class that steadyNewton()
    # zeroed at its structural floor. We therefore floor the step at the local
    # scale of the spectrum, log-interpolated across the gaps from the nonzero
    # neighbours exactly as the Newton solve does. A row that is zero throughout
    # (an extinct species) has no local scale either, but there the dynamics are
    # exactly linear about a zero baseline, so any step recovers the derivative
    # and the absolute floor is kept.
    N_scale <- fd_step_scale(N_vec,
                             positive_initial_guess(N_ss, params@w_min_idx,
                                                    active$w_top)[active$mask])

    npp_vec <- as.numeric(npp_ss)
    n_npp   <- length(npp_vec)
    # Step scale for the resource, floored the same way as for the fish.
    # The resource carries structural zeros above `w_pp_cutoff`, where the
    # capacity vanishes; log-interpolation with flat extrapolation gives
    # them the scale of the last nonzero class.
    npp_local <- positive_initial_guess(matrix(npp_vec, nrow = 1),
                                        w_min_idx = 1L, w_top = n_npp)
    npp_scale <- fd_step_scale(npp_vec, as.numeric(npp_local))
    x0    <- c(N_vec, npp_vec)
    scale <- c(N_scale, npp_scale)

    # Packed state -> the arrays the rate functions want.
    unpack <- function(x) {
        N <- N_ss
        N[active_idx] <- x[seq_len(n_fish_active)]
        n_pp <- npp_ss
        n_pp[] <- x[n_fish_active + seq_len(n_npp)]
        list(N = N, n_pp = n_pp)
    }

    # A custom rate function can be non-finite at a state that the differencing
    # visits. Fail loudly, naming the cell, rather than letting a NaN travel on
    # into eigen() and come back as a meaningless spectrum.
    check <- function(x, where) {
        if (!all(is.finite(x))) {
            stop(fn_name, "(): ", map_name, " returned non-finite values ",
                 where, ". Every state it evaluates satisfies N >= 0, so this ",
                 "points to a rate function that is not finite there.",
                 call. = FALSE)
        }
    }

    describe <- function(k) {
        if (k <= n_fish_active) {
            rc <- arrayInd(active_idx[k], dim(N_ss))
            paste0("when perturbing species ",
                   params@species_params$species[rc[1]], " at w = ",
                   signif(params@w[rc[2]], 3))
        } else {
            paste0("when perturbing the resource at w = ",
                   signif(params@w_full[k - n_fish_active], 3))
        }
    }

    list(params = params, effort = effort, rdd_const = rdd_const,
         n_other = n_other, active = active, active_idx = active_idx,
         n_fish_active = n_fish_active, n_npp = n_npp,
         rates_fns = rates_fns, flux_limiter = flux_limiter,
         N_ss = N_ss, npp_ss = npp_ss, x0 = x0, scale = scale,
         unpack = unpack, check = check, describe = describe)
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
#' first order in `h` rather than second.
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
    L <- matrix(0, nrow = n, ncol = n)

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
            L[, k] <- (f(x_plus) - base) / eps
        } else {
            x_minus <- x
            x_minus[k] <- x[k] - eps
            L[, k] <- (f(x_plus) - f(x_minus)) / (2 * eps)
        }
        ctx$check(L[, k], ctx$describe(k))
    }
    L
}

#' The continuous rates of change as a map on the packed state
#'
#' Returns the function that [getStability()] differentiates: the instantaneous
#' rate of change of every state variable. For the consumers this is
#' `-consumer_residual()`, which assembles exactly the semi-discretised
#' \eqn{dN/dt} with the spatial scheme configured in `params`. For the resource
#' it is the analytic semichemostat derivative; any other resource dynamics is
#' only available as a one-step map, so there it is differenced over a short
#' step, as `steady_rates()` does.
#'
#' @param ctx The context from `stability_context()`.
#' @param dt_resource Step length for differencing a resource dynamics function
#'   that is not the semichemostat.
#' @return A function of the packed state vector.
#' @noRd
stability_rhs <- function(ctx, dt_resource = 1e-4) {
    params  <- ctx$params
    targets <- c("EGrowth", "Mort", "Diffusion", "ResourceMort")
    semichemostat <- params@resource_dynamics == "resource_semichemostat"
    resource_dyn  <- if (semichemostat) NULL else get(params@resource_dynamics)

    function(x) {
        st <- ctx$unpack(x)
        r <- mizer_rates_subset(params, n = st$N, n_pp = st$n_pp,
                                n_other = ctx$n_other, t = 0,
                                effort = ctx$effort,
                                rates_fns = ctx$rates_fns, targets = targets)
        dNdt <- -consumer_residual(params, n = st$N, n_pp = st$n_pp,
                                   n_other = ctx$n_other, effort = ctx$effort,
                                   rdd = ctx$rdd_const,
                                   rates_fns = ctx$rates_fns,
                                   flux_limiter = ctx$flux_limiter,
                                   rates = r)
        if (semichemostat) {
            # dn_pp/dt = rr (cc - n_pp) - mu_R n_pp. Where both the
            # replenishment rate and the mortality vanish this is zero, so
            # such a class simply stays put, as resource_semichemostat()
            # also arranges.
            dnpp <- params@rr_pp * params@cc_pp -
                (params@rr_pp + r$resource_mort) * st$n_pp
        } else {
            npp_new <- resource_dyn(params, n = st$N, n_pp = st$n_pp,
                                    n_other = ctx$n_other, rates = r,
                                    t = 0, dt = dt_resource,
                                    resource_rate = params@rr_pp,
                                    resource_capacity = params@cc_pp)
            dnpp <- (npp_new - st$n_pp) / dt_resource
        }
        c(dNdt[ctx$active$mask], as.numeric(dnpp))
    }
}

#' Mizer's numerical one-step map on the packed state
#'
#' Returns the function that [getDiscreteStability()] differentiates: one step
#' of length `dt` of the dynamics as [project()] takes it with
#' `method = "euler"`, using the same `project_n_loop()` C++ Thomas solver and
#' the same transport coefficients.
#'
#' The resource is stepped with backward Euler rather than with the exact
#' semichemostat solution that [resource_semichemostat()] uses, so that the fast
#' resource modes cannot land on a discrete eigenvalue of exactly zero. The two
#' agree to first order in `dt` for the slow modes that decide stability, and
#' both damp the fast modes out within a step.
#'
#' @param ctx The context from `stability_context()`.
#' @param dt The step length.
#' @return A function of the packed state vector.
#' @noRd
stability_step <- function(ctx, dt) {
    params  <- ctx$params
    w_top   <- support_top_idx(params)
    targets <- c("EGrowth", "Mort", "Diffusion", "ResourceMort")
    semichemostat <- params@resource_dynamics == "resource_semichemostat"
    resource_dyn  <- if (semichemostat) NULL else get(params@resource_dynamics)

    function(x) {
        st <- ctx$unpack(x)
        r <- mizer_rates_subset(params, n = st$N, n_pp = st$n_pp,
                                n_other = ctx$n_other, t = 0,
                                effort = ctx$effort,
                                rates_fns = ctx$rates_fns, targets = targets)
        rdd <- ctx$rdd_const
        if (is.null(rdd)) {
            rdd <- state_rdd(params, n = st$N, n_pp = st$n_pp,
                             n_other = ctx$n_other, rates = r)
        }
        coefs <- get_transport_coefs(params, n = st$N, g = r$e_growth,
                                     mu = r$mort, dt = dt,
                                     recruitment_flux = rdd,
                                     d = r$diffusion,
                                     flux_limiter = ctx$flux_limiter)
        N_out <- project_n_loop(st$N, coefs$a, coefs$b, coefs$c, coefs$S,
                                params@w_min_idx)
        if (semichemostat) {
            mur <- params@rr_pp + r$resource_mort
            npp_out <- (st$n_pp + dt * params@rr_pp * params@cc_pp) /
                (1 + dt * mur)
        } else {
            npp_out <- resource_dyn(params, n = st$N, n_pp = st$n_pp,
                                    n_other = ctx$n_other, rates = r,
                                    t = 0, dt = dt,
                                    resource_rate = params@rr_pp,
                                    resource_capacity = params@cc_pp)
        }
        c(zero_above_support(N_out, w_top)[ctx$active$mask],
          as.numeric(npp_out))
    }
}

#' Reshape the leading eigenvectors back into abundance space
#'
#' @param evecs The eigenvector matrix, columns already sorted.
#' @param ctx The context from `stability_context()`.
#' @return The `leading_eigenvectors` component of the stability list.
#' @noRd
stability_eigenvectors <- function(evecs, ctx) {
    n_leading <- min(2L, ncol(evecs))
    N_ss  <- ctx$N_ss
    no_sp <- nrow(N_ss)
    no_w  <- ncol(N_ss)

    fish <- array(0 + 0i, dim = c(no_sp, no_w, n_leading),
                  dimnames = c(dimnames(N_ss), list(NULL)))
    for (k in seq_len(n_leading)) {
        v <- evecs[seq_len(ctx$n_fish_active), k]
        # Normalise to maximum modulus 1 for comparability.
        max_mod <- max(Mod(v))
        if (max_mod > 0) v <- v / max_mod
        M <- matrix(0 + 0i, nrow = no_sp, ncol = no_w)
        M[ctx$active_idx] <- v
        fish[, , k] <- M
    }
    resource <- matrix(0 + 0i, nrow = ctx$n_npp, ncol = n_leading)
    for (k in seq_len(n_leading)) {
        v <- evecs[ctx$n_fish_active + seq_len(ctx$n_npp), k]
        max_mod <- max(Mod(v))
        if (max_mod > 0) v <- v / max_mod
        resource[, k] <- v
    }
    rownames(resource) <- names(ctx$npp_ss)
    list(fish = fish, resource = resource)
}

#' Analyse the dynamic stability of a mizer steady state
#'
#' `r lifecycle::badge("experimental")`
#' Computes the eigenvalues of the linearised dynamics at the steady state
#' stored in `params@initial_n`. These eigenvalues determine whether the steady
#' state is dynamically stable and, when a Hopf bifurcation is approached, the
#' period of the emergent limit cycle.
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
#' A **Hopf bifurcation** occurs when a complex-conjugate pair of eigenvalues
#' crosses the imaginary axis, giving a limit-cycle period
#' \deqn{T = \frac{2\pi}{|\text{Im}(\lambda)|} \text{ years.}}
#'
#' The resource is a state variable of the system like any other: fish and
#' resource cells are perturbed independently, giving the full coupled Jacobian.
#' Its eigenvalues include both the slow fish modes and a cluster of fast
#' resource-relaxation modes, at \eqn{\lambda \approx -(r_{pp} + \mu_R)}. Any
#' resource dynamics function is supported: the semichemostat derivative is
#' written down analytically, and anything else is differenced over a short step.
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
#' negative abundances. Such columns are first order in `h` rather than second,
#' so they respond slightly more to a change of `h` than the rest.
#'
#' @param params A \linkS4class{MizerParams} object whose `initial_n` holds the
#'   steady state to analyse. Typically the output of [steadyNewton()].
#' @param reproduction Whether the reproduction rate is held fixed (`"fixed"`,
#'   default) or run dynamically (`"dynamic"`) during the evaluation.
#'   Must match the choice used when the steady state was computed.
#' @param effort The fishing effort to use. By default the initial effort
#'   stored in `params`.
#' @param h Relative step size for centred finite differences. Default `1e-4`.
#'   The result should not depend on this choice. If it does, the dynamics are
#'   not smooth at the state being analysed — see the section below.
#' @return A named list with the following components:
#'   \describe{
#'     \item{`eigenvalues`}{Complex vector of the continuous-time eigenvalues
#'       \eqn{\lambda_i}, sorted by decreasing real part.}
#'     \item{`max_real_part`}{The largest real part of the eigenvalues:
#'       \eqn{\max_i \text{Re}(\lambda_i)}. Greater than 0 means unstable.}
#'     \item{`stable`}{Logical: `TRUE` when `max_real_part < 0`.}
#'     \item{`dominant_period`}{The period (in years) of the dominant
#'       eigenvalue: `2*pi / abs(Im(lambda_1))`. `Inf` for a real
#'       dominant eigenvalue (monotone dynamics).}
#'     \item{`hopf_period`}{Period (in years) of the complex eigenvalue
#'       with the largest real part; `NULL` when no complex eigenvalue
#'       exists.  This is the expected limit-cycle period near a Hopf
#'       bifurcation.}
#'     \item{`n_active`}{Dimension of the Jacobian: the number of active fish
#'       cells plus all resource cells.}
#'     \item{`leading_eigenvectors`}{The eigenvectors of the two eigenvalues
#'       with the largest real part, reshaped back into the fish abundance
#'       space: a list with `$fish`, a complex array of shape
#'       `(n_species, n_sizes, 2)` with the same species and size dimnames as
#'       `params@initial_n`, and `$resource`, a complex matrix of shape
#'       `(n_w_full, 2)` for the resource component.
#'       Each eigenvector is normalised so that its maximum modulus equals 1.
#'       The real and imaginary parts of eigenvector 1 span the two-dimensional
#'       oscillation plane of the dominant mode; `Mod()` gives the amplitude
#'       pattern across species and sizes.}
#'     \item{`params`}{The validated `params` object the analysis was made at.}
#'   }
#' @section Requires smooth dynamics:
#' The finite-difference Jacobian is only meaningful if the rates of change are
#' differentiable at \eqn{N^*}. A custom rate function registered with
#' [setRateFunction()] that jumps as a function of the abundances breaks this in
#' two ways. If the state sits on the switching threshold, some perturbations
#' straddle it and pick up the jump, and the reported eigenvalues then vary
#' wildly with `h`. If the state is near but not on the threshold, no
#' perturbation crosses it, and the function silently returns the stability of
#' the single branch the state happens to lie on — which can read as `stable`
#' for a model whose simulations never settle.
#'
#' Re-running with a different `h` is the cheapest check: if the answer moves,
#' do not trust it. See [Discontinuous rate
#' functions](https://sizespectrum.org/mizer/articles/discontinuous_rates.html).
#'
#' @seealso [steadyNewton()], [getDiscreteStability()], [getLimitCycleSim()]
#' @export
getStability <- function(params,
                         reproduction = c("fixed", "dynamic"),
                         effort = params@initial_effort,
                         h = 1e-4) {
    ctx <- stability_context(params, reproduction = reproduction,
                             effort = effort,
                             fn_name = "getStability",
                             map_name = "the rates of change")
    J <- stability_jacobian(ctx, stability_rhs(ctx), h)

    eig <- eigen(J)
    # Sort by decreasing real part (most unstable first).
    ord   <- order(Re(eig$values), decreasing = TRUE)
    evals <- as.complex(eig$values[ord])
    evecs <- eig$vectors[, ord, drop = FALSE]

    max_real_part <- Re(evals[1])
    omega1 <- Im(evals[1])
    dominant_period <- if (abs(omega1) < 1e-10) Inf else 2 * pi / abs(omega1)

    is_complex <- abs(Im(evals)) > 1e-8
    if (any(is_complex)) {
        # The complex eigenvalue closest to (or furthest across) the imaginary
        # axis is the one whose oscillation the dynamics show.
        lam_hopf <- evals[is_complex][which.max(Re(evals[is_complex]))]
        hopf_period <- 2 * pi / abs(Im(lam_hopf))
    } else {
        hopf_period <- NULL
    }

    list(
        eigenvalues          = evals,
        max_real_part        = max_real_part,
        stable               = max_real_part < 0,
        dominant_period      = dominant_period,
        hopf_period          = hopf_period,
        n_active             = length(ctx$x0),
        leading_eigenvectors = stability_eigenvectors(evecs, ctx),
        params               = ctx$params
    )
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
#' The resource is stepped with backward Euler,
#' \deqn{n_{pp}^{t+1} = \frac{n_{pp}^t + dt\,r_{pp}c_{pp}}%
#'                           {1 + dt\,(r_{pp} + \mu_R^t)},}
#' rather than with the exact semichemostat solution that
#' [resource_semichemostat()] uses in `project()`. The two agree to first order
#' in `dt` for the slow modes that decide stability, while backward Euler keeps
#' the fast resource modes away from a discrete eigenvalue of exactly zero.
#'
#' @inheritParams getStability
#' @param dt The time step size of the one-step map. Default `1`.
#' @return A named list with the following components:
#'   \describe{
#'     \item{`discrete_eigenvalues`}{Complex vector of the eigenvalues
#'       \eqn{\mu_i} of the one-step map, sorted by decreasing modulus.}
#'     \item{`spectral_radius`}{\eqn{\max_i|\mu_i|}. Less than 1 means the
#'       numerical scheme is stable at this `dt`.}
#'     \item{`stable`}{Logical: `TRUE` when `spectral_radius < 1`.}
#'     \item{`dt`}{The step size the map was evaluated at.}
#'     \item{`n_active`}{Dimension of the Jacobian.}
#'     \item{`leading_eigenvectors`}{The eigenvectors of the two
#'       largest-modulus eigenvalues, in the same shape as for
#'       [getStability()].}
#'     \item{`params`}{The validated `params` object the analysis was made at.}
#'   }
#' @inheritSection getStability Requires smooth dynamics
#' @seealso [getStability()], [steadyNewton()]
#' @export
getDiscreteStability <- function(params,
                                 reproduction = c("fixed", "dynamic"),
                                 effort = params@initial_effort,
                                 h = 1e-4,
                                 dt = 1) {
    assert_that(is.number(dt), dt > 0)
    ctx <- stability_context(params, reproduction = reproduction,
                             effort = effort,
                             fn_name = "getDiscreteStability",
                             map_name = "the one-step map")
    L <- stability_jacobian(ctx, stability_step(ctx, dt), h)

    eig <- eigen(L)
    ord   <- order(Mod(eig$values), decreasing = TRUE)
    evals <- as.complex(eig$values[ord])
    evecs <- eig$vectors[, ord, drop = FALSE]
    spectral_radius <- Mod(evals[1])

    list(
        discrete_eigenvalues = evals,
        spectral_radius      = spectral_radius,
        stable               = spectral_radius < 1,
        dt                   = dt,
        n_active             = length(ctx$x0),
        leading_eigenvectors = stability_eigenvectors(evecs, ctx),
        params               = ctx$params
    )
}
