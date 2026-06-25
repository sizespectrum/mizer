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
#' Like [steady()], the function holds the reproduction rate (RDD) constant
#' while solving for the consumer spectra, substitutes the analytic steady state
#' of the resource, and keeps any other components constant. After the spectra
#' have been found it restores density-dependent Beverton-Holt reproduction with
#' [setBevertonHolt()], honouring the `preserve` argument exactly as [steady()]
#' does.
#'
#' The consumer densities are solved for in log space, which both keeps them
#' positive and conditions the otherwise badly-scaled system. The unknowns are
#' the size classes that carry non-zero density, taken from the current
#' `initial_n` rather than from `w_max`: the support runs from the egg size up to
#' the highest class with non-zero density in the supplied spectra (capped at the
#' grid). This matters because users routinely set `w_max` much larger than the
#' largest fish, and the structurally-zero classes below such a `w_max` must not
#' enter the log-space system. The nonlinear system is solved with a globalised
#' Newton iteration from the `nleqslv` package, starting from the current
#' `initial_n`. Newton's method converges from any starting point in the *root's*
#' basin of attraction (which is unrelated to the dynamic stability of the steady
#' state), so a reasonable initial guess should be supplied in `initialN(params)`
#' — for example the spectra from a nearby stable parameterisation, or the
#' (diverging) output of [steady()]. The guess also defines the support, so it
#' should be non-zero across the size range that should carry density.
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
#'   `reproduction_level`.
#' @param tol Convergence tolerance passed to [nleqslv::nleqslv()] (both the
#'   function-value tolerance `ftol` and the step tolerance `xtol`).
#' @param maxit Maximum number of iterations for [nleqslv::nleqslv()].
#' @param method The [nleqslv::nleqslv()] method, either `"Newton"` (default, a
#'   full Newton step with a numerical Jacobian — fastest and most accurate when
#'   the starting guess is reasonable) or `"Broyden"`.
#' @param global The globalisation strategy passed to [nleqslv::nleqslv()].
#'   The default `"dbldog"` (double dogleg) is a robust trust-region method.
#' @param ... Unused.
#' @return A \linkS4class{MizerParams} object with the initial state set to the
#'   steady state.
#' @seealso [steady()], [steadySingleSpecies()]
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
                                     tol = 1e-10, maxit = 100,
                                     method = c("Newton", "Broyden"),
                                     global = "dbldog", ...) {
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

    if (params@rates_funcs$RDD == "BevertonHoltRDD") {
        preserve <- match.arg(preserve)
        old_reproduction_level <- getReproductionLevel(params)
        old_R_max <- params@species_params$R_max
        old_erepro <- params@species_params$erepro
    }

    # Hold reproduction at the current level for the duration of the solve,
    # exactly as steady() does with constant_reproduction.
    rdd_const <- getRDD(params)
    n_other <- params@initial_n_other

    active <- steady_active_set(params)
    residual_fn <- steady_state_residual(params, rdd_const, n_other, effort,
                                         active)

    # The log-space solve needs a strictly positive, well-scaled start. The
    # second-order schemes can leave isolated zeros inside the support (a
    # negativity-floor artefact), which would make the 1/N-scaled residual
    # overflow. Fill those by log-interpolation from the nonzero neighbours.
    N0 <- positive_initial_guess(params@initial_n, params@w_min_idx,
                                 active$w_top)
    x0 <- log(N0[active$mask])
    sol <- nleqslv::nleqslv(x0, residual_fn, method = method,
                            global = global,
                            control = list(maxit = maxit, ftol = tol,
                                           xtol = tol))
    # termcd 1 (function convergence) and 2 (x convergence) are successes.
    if (sol$termcd > 2) {
        warning("steadyNewton() did not converge (nleqslv termination code ",
                sol$termcd, ": ", sol$message,
                "). Returning the best iterate found.", call. = FALSE)
    }

    N <- active$unpack(sol$x)
    n_pp <- resource_steady_semichemostat(params, N, n_other)

    params@initial_n[] <- N
    params@initial_n_pp[] <- n_pp

    # Restore density-dependent reproduction, just as steady() does.
    if (params@rates_funcs$RDD == "BevertonHoltRDD") {
        if (preserve == "reproduction_level") {
            params <- setBevertonHolt(params,
                                      reproduction_level = old_reproduction_level)
        } else if (preserve == "R_max") {
            params <- setBevertonHolt(params, R_max = old_R_max)
        } else {
            params <- setBevertonHolt(params, erepro = old_erepro)
        }
    }

    params@time_modified <- lubridate::now()
    params
}

#' The set of size classes solved for by steadyNewton()
#'
#' Builds the logical mask of the (species x size) density matrix that the
#' direct solver treats as unknowns. For each species the unknowns run from the
#' egg size `w_min_idx` up to the highest class that actually carries density in
#' the current `initial_n`, capped at the grid truncation [support_top_idx()].
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
#' [positive_initial_guess()]; only a trailing run of zeros is dropped.
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
        pos <- which(n[i, lo:grid_top[i]] > 0)
        # Highest class carrying density (fall back to the egg class if the
        # species is empty), capped at the grid truncation.
        w_top[i] <- if (length(pos)) lo + max(pos) - 1L else lo
        mask[i, lo:w_top[i]] <- TRUE
    }

    dn <- dimnames(params@initial_n)
    unpack <- function(x) {
        N <- matrix(0, nrow = no_sp, ncol = no_w, dimnames = dn)
        N[mask] <- exp(x)
        N
    }
    list(mask = mask, w_top = w_top, unpack = unpack)
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
#' @param w_top Per-species support top from [steady_active_set()].
#' @return A density matrix that is strictly positive on the active range.
#' @noRd
positive_initial_guess <- function(N, w_min_idx, w_top) {
    for (i in seq_len(nrow(N))) {
        rng <- w_min_idx[i]:w_top[i]
        v <- N[i, rng]
        pos <- v > 0
        if (all(pos) || !any(pos)) next
        idx <- seq_along(v)
        v[!pos] <- exp(stats::approx(idx[pos], log(v[pos]),
                                     xout = idx[!pos], rule = 2)$y)
        N[i, rng] <- v
    }
    N
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
    # Resource mortality is predation by consumers and does not depend on the
    # resource density, so the value of n_pp passed here is irrelevant.
    mu_R <- as.numeric(getResourceMort(params, n = N,
                                       n_pp = params@initial_n_pp,
                                       n_other = n_other, t = 0))
    mur <- params@rr_pp + mu_R
    n_pp <- params@rr_pp * params@cc_pp / mur
    # Where both rate and mortality vanish the steady state is undetermined;
    # keep the current value, exactly as resource_semichemostat() does.
    sel <- !is.finite(n_pp)
    n_pp[sel] <- params@initial_n_pp[sel]
    n_pp
}

#' Residual of the discrete steady-state equation
#'
#' Returns a closure `f(x)` suitable for [nleqslv::nleqslv()]. The argument `x`
#' is the vector of log-densities of the active size classes (see
#' [steady_active_set()]). The closure rebuilds the full density matrix,
#' substitutes the analytic resource steady state, evaluates the growth,
#' mortality and diffusion rates, assembles the steady-state tridiagonal
#' coefficients with [get_transport_coefs()] at `dt = 1`, and returns the
#' steady-state residual
#' \deqn{F_j = a_j N_{j-1} + b_j N_j + c_j N_{j+1} - S_j}
#' (since `S = N + recruitment source`, the `+N` cancels the backward-Euler
#' `N^t` term, leaving exactly the steady-state equation
#' \eqn{\tilde A N_{j-1} + \tilde B N_j + \tilde C N_{j+1} - \tilde S = 0}).
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
#' @param active The active-set list from [steady_active_set()].
#' @return A function of the packed log-density vector returning the packed
#'   scaled residual.
#' @noRd
steady_state_residual <- function(params, rdd_const, n_other, effort, active) {
    no_w <- length(params@w)
    mask <- active$mask
    flux_limiter <- flux_limiter_scheme(params)
    rates_fns <- projectRateFunctions(params)

    function(x) {
        N <- active$unpack(x)
        n_pp <- resource_steady_semichemostat(params, N, n_other)

        r <- mizer_rates_subset(params, n = N, n_pp = n_pp, n_other = n_other,
                                t = 0, effort = effort, rates_fns = rates_fns,
                                targets = c("EGrowth", "Mort", "Diffusion"))

        coefs <- get_transport_coefs(params, n = N, g = r$e_growth,
                                     mu = r$mort, dt = 1,
                                     recruitment_flux = rdd_const,
                                     d = r$diffusion,
                                     flux_limiter = flux_limiter)

        Nm <- cbind(0, N[, -no_w, drop = FALSE])   # N_{j-1}
        Np <- cbind(N[, -1, drop = FALSE], 0)      # N_{j+1}
        res <- coefs$a * Nm + coefs$b * N + coefs$c * Np - coefs$S

        # Scale to a per-capita rate of change (N > 0 on the active set).
        (res / N)[mask]
    }
}
