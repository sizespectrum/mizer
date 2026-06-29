#' Translate the second_order_w flux entry into the internal scheme name
#'
#' The advective-flux scheme is stored directly in the `flux` entry of the
#' `second_order_w` slot as `"upwind"`, `"van_leer"`, or `"centred"`. The
#' internal transport routines use `"none"` for the first-order upwind scheme,
#' so this helper maps `"upwind"` to `"none"`.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @return `"none"`, `"van_leer"` or `"centred"`.
#' @noRd
flux_limiter_scheme <- function(params) {
    flux <- params@second_order_w[["flux"]]
    if (flux == "upwind") "none" else flux
}

#' Uniform log-size grid spacing
#'
#' The mizer size grid is geometric, so `x = log(w)` is uniformly spaced. This
#' helper returns that spacing `h = log(w_{j+1}/w_j)`, the natural step of the
#' second-order finite-volume scheme.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @return The scalar log-size spacing `h`.
#' @noRd
log_dx <- function(params) {
    log(params@w[2] / params@w[1])
}

#' Highest size class retained at the upper boundary (the support top)
#'
#' The size grid is truncated at each species' maximum size `w_max`, and the
#' density is held at zero above it. `w_max` is purely a computational boundary
#' (the top of the size grid); it is **not** the size at which somatic growth
#' stops. Growth slows around the separate parameter `w_repro_max` (where a
#' typical mature individual invests all energy in reproduction), so whether the
#' growth rate is actually zero near `w_max` depends entirely on how large the
#' user chose `w_max`. This cut-off therefore makes no assumption about the rates;
#' it is appropriate when `w_max` is chosen large enough that the density there is
#' negligible.
#'
#' The top retained class depends on the scheme, because the two schemes place the
#' advective inflow differently. With `w_max_idx = sum(w <= w_max)`:
#'
#' * the first-order upwind scheme feeds class `j` from the growth of the class
#'   below, `g(w_{j-1})`, so the class just above the one containing `w_max` can
#'   still carry advected density; the support top is `w_max_idx + 1`. This
#'   reproduces the long-standing (untruncated) mizer behaviour exactly;
#' * the second-order scheme reconstructs the inflow at a class's own lower face,
#'   `g(w_j)`, so its support ends one class lower, at `w_max_idx`.
#'
#' The result is capped at the top of the grid. Because it is fixed by `w_max` and
#' the scheme rather than by the instantaneous growth field, it is well defined
#' even when a frozen, food-limited or otherwise degenerate growth field is
#' supplied (e.g. in the second-order time-steppers' frozen-rate solves).
#'
#' The dynamics impose this as the upper boundary condition (see
#' `get_transport_coefs()`) so that abundance is held at zero above `w_max` even
#' when diffusion would otherwise carry density past it, and [steadyNewton()]
#' solves on exactly this support.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @return A per-species integer vector of the top retained size-class index.
#' @noRd
support_top_idx <- function(params) {
    w <- params@w
    no_w <- length(w)
    w_max <- params@species_params$w_max
    w_max_idx <- vapply(w_max, function(wm) sum(w <= wm), integer(1))
    # The grid is truncated at w_max. First-order upwind feeds a class from the
    # growth of the class below, so its support reaches one bin past w_max (this
    # matches untruncated mizer); the second-order scheme feeds a class from the
    # growth at its own lower face, so it stops at w_max.
    offset <- if (flux_limiter_scheme(params) == "none") 1L else 0L
    pmin(w_max_idx + offset, no_w)
}

#' Sever the coupling to the size classes above the support top
#'
#' Sets `c = 0` at the support-top class `w_top` (see `support_top_idx()`), which
#' decouples the active spectrum from the (held-at-zero) classes above it in the
#' tridiagonal solve. Cutting `c_{top}` zeroes the back-substitution coefficient
#' there, so the active solution up to `w_top` never reads the inactive tail. The
#' tail itself is held at zero after the solve by `zero_above_support()`.
#'
#' For the default no-diffusion scheme this is a no-op, because `c` already
#' vanishes at the top class (nothing diffuses or grows into the class above);
#' with diffusion it stops the diffusive flux above `w_max` from re-entering the
#' active spectrum.
#'
#' @param coefs The list of tridiagonal coefficients `a`, `b`, `c`, `S`.
#' @param w_top Per-species support-top index from `support_top_idx()`.
#' @return The coefficient list with the upper boundary condition applied.
#' @noRd
apply_upper_cutoff <- function(coefs, w_top) {
    top_idx <- cbind(seq_len(nrow(coefs$c)), w_top)
    coefs$c[top_idx] <- 0
    coefs
}

#' Hold abundance at zero above the support top
#'
#' Zeros every size class above the support top `w_top` (see `support_top_idx()`)
#' after a density update. This is the partner of `apply_upper_cutoff()`: the
#' coefficient surgery decouples the active spectrum from these classes during the
#' solve, and this enforces the upper boundary condition `N = 0` above `w_max` on
#' the result. For the default no-diffusion scheme these classes are already zero
#' (nothing flows into them), so this is a no-op; with diffusion it removes the
#' density that would otherwise leak above `w_max`.
#'
#' @param n The updated density matrix (species x size).
#' @param w_top Per-species support-top index from `support_top_idx()`.
#' @return `n` with all classes above the support top set to zero.
#' @noRd
zero_above_support <- function(n, w_top) {
    no_w <- ncol(n)
    w_idx_mat <- matrix(seq_len(no_w), nrow = nrow(n), ncol = no_w, byrow = TRUE)
    n[w_idx_mat > w_top] <- 0
    n
}

#' Helper function to calculate the transport coefficients
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param n An array (species x size) with the number density at the current time step.
#' @param g The growth rate.
#' @param mu The mortality rate.
#' @param dt Time step.
#' @param recruitment_flux Per-species influx at the egg boundary.
#' @param d The diffusion rate.
#' @param flux_limiter Advective-flux scheme: `"none"` (first-order upwind),
#'   `"van_leer"` or `"centred"` (the second-order log-size scheme).
#'
#' This calculates the coefficients A, B, C and S for the linear system
#' A_j * N_{j-1} + B_j * N_j + C_j * N_{j+1} = S_j
#' For details see the [Numerical Details](https://sizespectrum.org/mizer/articles/numerical_details.html#discretised-equation) vignette.
#'
#' @return A list with the coefficients A, B, C and S.
#' @noRd
get_transport_coefs <- function(params, n, g, mu, dt, recruitment_flux, d,
                                flux_limiter = "none") {
    coefs <- if (flux_limiter == "none") {
        get_transport_coefs_upwind(params, n, g, mu, dt, recruitment_flux, d)
    } else {
        get_transport_coefs_logfv(params, n, g, mu, dt, recruitment_flux, d,
                                  flux_limiter)
    }
    # Upper boundary: hold abundance at zero above the maximum size w_max.
    apply_upper_cutoff(coefs, support_top_idx(params))
}

#' Second-order finite-volume transport coefficients
#'
#' Builds the tridiagonal coefficients for the second-order scheme described in
#' the "Reducing the spatial error" section of the
#' [Numerical Details](https://sizespectrum.org/mizer/articles/numerical_details.html)
#' vignette. The finite-volume cells are the size bins \eqn{[w_j, w_{j+1}]}, so
#' the flux divergence keeps the bin width \eqn{\Delta w_j} as its divisor and the
#' bin boundaries are the nodes \eqn{w_j}, exactly as in the first-order scheme.
#' Second order comes from placing each quantity where the finite-volume update
#' needs it:
#' \itemize{
#'   \item the growth velocity is the point value at the bin boundary, \eqn{g_j};
#'   \item the density at the boundary is reconstructed from the two flanking bin
#'     averages, \eqn{N_{j-1} + \tfrac12\chi_j(N_j - N_{j-1})}, which is the
#'     centred value \eqn{\tfrac12(N_{j-1}+N_j)} when \eqn{\chi_j = 1} and pure
#'     upwind when \eqn{\chi_j = 0};
#'   \item the diffusion coefficient is the **bin average** \eqn{d_j} supplied in
#'     `d` ([getDiffusion()] returns the bin average when `bin_average` is on,
#'     gated just like the mortality), so the diffusive flux differences the
#'     bin-averaged products \eqn{d_jN_j} between bin centres.
#' }
#' The advective flux is therefore
#' \eqn{J_j = g_j\,[N_{j-1} + \tfrac12\chi_j(N_j - N_{j-1})]}. The sink uses the
#' supplied mortality `mu`, which is the bin average when `bin_average` is on;
#' both flags together give a fully second-order step. Because \eqn{\chi} is
#' frozen the operator stays tridiagonal with the same \eqn{\Delta w} divisor.
#' @noRd
get_transport_coefs_logfv <- function(params, n, g, mu, dt, recruitment_flux,
                                      d, flux_limiter) {
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)
    w <- params@w
    h <- log_dx(params)
    beta <- w[2] / w[1]
    M <- function(v) matrix(v, nrow = no_sp, ncol = no_w, byrow = TRUE)

    # `d` is the diffusion coefficient d_j as supplied by getDiffusion(): the
    # bin average when bin_average is on (the value the second-order diffusive
    # flux needs, co-located with the bin average N_j) and the point value
    # otherwise. The point-vs-average choice is made there, gated on
    # bin_average, exactly as for the mortality `mu`.
    dR <- cbind(d[, -1, drop = FALSE], d[, no_w, drop = FALSE])  # d_{j+1}
    dL <- cbind(d[, 1, drop = FALSE], d[, -no_w, drop = FALSE])  # d_{j-1}

    # Growth velocity at the faces: g_j at the lower face, g_{j+1} at the upper
    # face (clamped at the outflow top, where the limiter is off anyway).
    gR <- cbind(g[, -1, drop = FALSE], g[, no_w, drop = FALSE])     # g_{j+1}

    # Frozen reconstruction weight chi[, j] at the face below cell j. chi == 1
    # gives the centred reconstruction, chi == 0 first-order upwind.
    chi <- flux_limiter_chi(params, n, g, flux_limiter)
    chiR <- cbind(chi[, -1, drop = FALSE], matrix(0, no_sp, 1))     # chi_{j+1}

    # Per-node 1 / (2 h w) factors for the diffusive flux, with the upper face
    # using w_{j+1} (extrapolated past the top).
    wR <- c(w[-1], w[no_w] * beta)
    f_lo <- M(1 / (2 * h * w))      # belongs to the lower face j
    f_hi <- M(1 / (2 * h * wR))     # belongs to the upper face j+1
    over_dw <- M(dt / params@dw)

    # Interior coefficients (faces at the nodes, divisor Delta w_j).
    a <- -over_dw * (g * (1 - 0.5 * chi) + dL * f_lo)
    b <- 1 + dt * mu +
        over_dw * (gR * (1 - 0.5 * chiR) - 0.5 * chi * g + d * (f_lo + f_hi))
    c <- over_dw * (0.5 * chiR * gR - dR * f_hi)
    dimnames(a) <- dimnames(b) <- dimnames(c) <-
        list(params@species_params$species, NULL)
    c[, no_w] <- 0                  # no cell above the top: pure outflow

    S <- n
    dimnames(S) <- dimnames(a)

    # Egg boundary: the flux from below is the recruitment flux R_dd, so there is
    # no transport from j_start-1 (a = 0) and the lower-face advective and
    # diffusive terms drop out of b, leaving only the upper face j_start+1.
    j_start <- params@w_min_idx
    idxs <- cbind(seq_len(no_sp), j_start)
    b[idxs] <- 1 + dt * mu[idxs] + (dt / params@dw[j_start]) *
        (gR[idxs] * (1 - 0.5 * chiR[idxs]) +
             d[idxs] / (2 * h * wR[j_start]))
    a[idxs] <- 0
    S[idxs] <- S[idxs] + recruitment_flux * dt / params@dw[j_start]

    # Zero out everything below the recruitment size for each species.
    w_idx_mat <- matrix(seq_len(no_w), no_sp, no_w, byrow = TRUE)
    mask_below <- w_idx_mat < j_start
    if (any(mask_below)) {
        a[mask_below] <- 0
        b[mask_below] <- 0
        c[mask_below] <- 0
        S[mask_below] <- 0
    }

    list(a = a, b = b, c = c, S = S)
}

#' First-order upwind transport coefficients (legacy default scheme)
#'
#' The original mizer discretisation: first-order upwind advective flux and a
#' centred diffusive flux, both placed at the lower node `w_j` of bin
#' `[w_j, w_{j+1}]`. This is what `flux_limiter = "none"` selects and is kept
#' bit-for-bit so that the default model reproduces earlier mizer versions.
#' @noRd
get_transport_coefs_upwind <- function(params, n, g, mu, dt, recruitment_flux,
                                       d) {
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)
    dw <- params@dw
    dt_dw <- matrix(dt / dw, nrow = no_sp, ncol = no_w, byrow = TRUE)

    a <- matrix(0, nrow = no_sp, ncol = no_w,
                dimnames = list(params@species_params$species, NULL))
    b <- matrix(0, nrow = no_sp, ncol = no_w, dimnames = dimnames(a))
    c <- matrix(0, nrow = no_sp, ncol = no_w, dimnames = dimnames(a))
    S <- matrix(0, nrow = no_sp, ncol = no_w, dimnames = dimnames(a))

    idx <- 2:no_w
    idx_minus_1 <- idx - 1
    term_diff_minus_1 <- 0.5 * d[, idx_minus_1] /
        matrix(dw[idx_minus_1], nrow = no_sp, ncol = length(idx), byrow = TRUE)
    a[, idx] <- -dt_dw[, idx] * (g[, idx_minus_1] + term_diff_minus_1)

    idx_j <- 1:(no_w - 1)
    term_diff_plus_1 <- 0.5 * d[, idx_j + 1] /
        matrix(dw[idx_j], nrow = no_sp, ncol = length(idx_j), byrow = TRUE)
    c[, idx_j] <- -dt_dw[, idx_j] * term_diff_plus_1

    term_diff_j <- 0.5 * d[, idx] /
        matrix(dw[idx], nrow = no_sp, ncol = length(idx), byrow = TRUE)
    term_diff_j_minus_1 <- 0.5 * d[, idx] /
        matrix(dw[idx_minus_1], nrow = no_sp, ncol = length(idx), byrow = TRUE)
    b[, idx] <- 1 + dt * mu[, idx] +
        dt_dw[, idx] * (g[, idx] + term_diff_j + term_diff_j_minus_1)

    j_start <- params@w_min_idx
    idxs <- cbind(1:no_sp, j_start)
    S[] <- n
    S[idxs] <- S[idxs] + recruitment_flux * dt / params@dw[j_start]
    b[idxs] <- 1 + dt * mu[idxs] +
        dt_dw[idxs] * (g[idxs] + 0.5 * d[idxs] / params@dw[j_start])
    a[idxs] <- 0

    w_idx_mat <- matrix(1:no_w, nrow = no_sp, ncol = no_w, byrow = TRUE)
    mask_below <- w_idx_mat < j_start
    if (any(mask_below)) {
        a[mask_below] <- 0
        b[mask_below] <- 0
        c[mask_below] <- 0
        S[mask_below] <- 0
    }

    list(a = a, b = b, c = c, S = S)
}

# Flux limiter functions chi(r) of the smoothness ratio r. To add another
# limiter, give it a name here and offer that name in the `flux_limiter`
# argument of `project()`.
flux_limiter_funcs <- list(
    van_leer = function(r) (r + abs(r)) / (1 + abs(r))
)

#' Frozen reconstruction weights at the size-bin faces
#'
#' The second-order scheme reconstructs the density at the face below cell `j`
#' from the two flanking cell averages,
#' \deqn{N_j^{face} = N_{j-1} + \tfrac12\,\chi_j\,(N_j - N_{j-1}),}
#' so that \eqn{\chi_j = 1} gives the centred value \eqn{\tfrac12(N_{j-1}+N_j)}
#' and \eqn{\chi_j = 0} gives pure upwind. This helper returns the weights
#' \eqn{\chi}, frozen at a known density field so that the advective flux
#' \eqn{g_j N_j^{face}} is linear in the unknowns and the operator stays
#' tridiagonal. `chi[, j]` is the weight at the face below cell `j`.
#'
#' For `"centred"` the weight is `1` everywhere (pure second-order centred flux,
#' not TVD: it can over/undershoot but is genuinely second order, including at
#' extrema). For `"van_leer"` it is the van Leer limiter \eqn{\chi(r_j)} of the
#' smoothness ratio \eqn{r_j} of successive jumps of the density \eqn{N}, which
#' falls to `0` at extrema and so keeps the update positivity-friendly (TVD). In
#' both cases \eqn{\chi} is forced to `0` at and below the non-smooth recruitment
#' boundary, where the advective flux stays first-order upwind.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param n Density field (species x size) at which to freeze the weights.
#' @param g Growth rate (species x size). Unused; kept for interface symmetry.
#' @param flux_limiter `"none"`, `"centred"` or `"van_leer"`.
#' @return A (species x size) matrix `chi`. Zero everywhere when
#'   `flux_limiter == "none"`.
#' @noRd
flux_limiter_chi <- function(params, n, g, flux_limiter) {
    no_sp <- nrow(n)
    no_w <- ncol(n)
    chi <- matrix(0, nrow = no_sp, ncol = no_w)
    # The high-order term needs an upstream and a downstream bin, so it is only
    # defined once there are at least three size bins.
    if (flux_limiter == "none" || no_w < 3) {
        return(chi)
    }

    if (flux_limiter == "centred") {
        chi[, 3:no_w] <- 1
    } else {
        limiter <- flux_limiter_funcs[[flux_limiter]]
        # Jump of the reconstructed variable N across the face below cell j,
        # stored in column j: delta[, j] = N_j - N_{j-1}.
        delta <- matrix(0, nrow = no_sp, ncol = no_w)
        delta[, 2:no_w] <- n[, 2:no_w] - n[, 1:(no_w - 1)]
        jint <- 3:no_w
        num <- delta[, jint - 1, drop = FALSE]   # delta_{j-1}
        den <- delta[, jint, drop = FALSE]       # delta_j
        r <- ifelse(den == 0, 0, num / den)      # upwind where the jump is 0
        chi[, jint] <- limiter(r)
    }

    # First-order upwind at and below the recruitment boundary, and at the first
    # two faces above it where the smoothness ratio would otherwise reference
    # the inactive region below the boundary or the non-smooth recruitment boundary cell.
    j_start <- params@w_min_idx
    w_idx_mat <- matrix(1:no_w, nrow = no_sp, ncol = no_w, byrow = TRUE)
    chi[w_idx_mat <= j_start + 2] <- 0
    chi
}
