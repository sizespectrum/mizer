#' Helper function to calculate the transport coefficients for the upwind-difference scheme
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param n An array (species x size) with the number density at the current time step.
#' @param g The growth rate.
#' @param mu The mortality rate.
#' @param dt Time step.
#' 
#' This calculates the coefficients A, B, C and S for the linear system
#' A_j * N_{j-1} + B_j * N_j + C_j * N_{j+1} = S_j
#' For details see the [Numerical Details](https://sizespectrum.org/mizer/articles/numerical_details.html#discretised-equation) vignette.
#'
#' @return A list with the coefficients A, B, C and S.
#' @noRd
get_transport_coefs <- function(params, n, g, mu, dt, recruitment_flux, d,
                                flux_limiter = "none") {

    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)
    
    # Pre-calculate some common terms
    # dw_j
    dw <- params@dw
    # Delta t / Delta w_j
    dt_dw <- matrix(dt / dw, nrow = no_sp, ncol = no_w, byrow = TRUE)
    
    # Initialize matrices
    a <- matrix(0, nrow = no_sp, ncol = no_w, 
                dimnames = list(params@species_params$species, NULL))
    b <- matrix(0, nrow = no_sp, ncol = no_w,
                dimnames = list(params@species_params$species, NULL))
    c <- matrix(0, nrow = no_sp, ncol = no_w,
                dimnames = list(params@species_params$species, NULL))
    S <- matrix(0, nrow = no_sp, ncol = no_w,
                dimnames = list(params@species_params$species, NULL))
    
    # We assume d is 0 at the boundaries for simplicity or it is handled by the loop constraints
    # Actually for j=1, A is 0, for j=no_w, C is 0 (boundary condition)

    # Frozen flux limiter. `psi[, j]` is the limiter value at the interface w_j
    # (between cells j-1 and j), computed from the field `n`. With psi == 0 the
    # advective flux is plain first-order upwind; with psi > 0 it blends in the
    # centred high-order flux, J_j = g_{j-1} N_{j-1} + 0.5 psi_j (g_j N_j -
    # g_{j-1} N_{j-1}). Because psi is frozen (lagged), J_j is linear in the
    # unknowns and the operator stays tridiagonal.
    psi <- flux_limiter_psi(params, n, g, flux_limiter)

    # A_j = - dt/dw_j * (g_{j-1} (1 - 0.5 psi_j) + D_{j-1} / (2 * dw_{j-1}))
    # We compute A for j in idx (2:no_w). g_{j-1} corresponds to columns 1:(no_w-1)

    # Indices for j-1 and j
    idx <- 2:no_w
    idx_minus_1 <- idx - 1

    term_diff_minus_1 <- 0.5 * d[, idx_minus_1] /
        matrix(dw[idx_minus_1], nrow = no_sp, ncol = length(idx), byrow = TRUE)

    a[, idx] <- -dt_dw[, idx] *
        (g[, idx_minus_1] * (1 - 0.5 * psi[, idx]) + term_diff_minus_1)

    # C_j gains a high-order advective term +0.5 psi_{j+1} g_{j+1} on top of the
    # diffusive term. At j=no_w we assume N_{j+1}=0, so C_{no_w} is left at 0.
    idx_j <- 1:(no_w - 1)
    term_diff_plus_1 <- 0.5 * d[, idx_j + 1] /
        matrix(dw[idx_j], nrow = no_sp, ncol = length(idx_j), byrow = TRUE)

    c[, idx_j] <- -dt_dw[, idx_j] *
        (term_diff_plus_1 - 0.5 * psi[, idx_j + 1] * g[, idx_j + 1])
    # c[, no_w] is 0 as initialized

    # B_j advective term g_j (1 - 0.5 psi_j - 0.5 psi_{j+1}) plus the diffusive
    # terms. psi_{j+1} for j = no_w is the outflow interface, taken to be 0.
    term_diff_j <- 0.5 * d[, idx] / matrix(dw[idx], nrow = no_sp, ncol = length(idx), byrow = TRUE)
    term_diff_j_minus_1 <- 0.5 * d[, idx] / matrix(dw[idx_minus_1], nrow = no_sp, ncol = length(idx), byrow = TRUE)

    # psi at the upper interface j+1, aligned to idx = 2:no_w (0 at the top).
    psi_upper <- matrix(0, nrow = no_sp, ncol = length(idx))
    if (no_w >= 3) {
        psi_upper[, seq_len(no_w - 2)] <- psi[, 3:no_w]
    }

    b[, idx] <- 1 + dt * mu[, idx] +
        dt_dw[, idx] * (g[, idx] * (1 - 0.5 * psi[, idx] - 0.5 * psi_upper) +
                            term_diff_j + term_diff_j_minus_1)

    # Boundary condition updates
    # We treat the start of the size spectrum (j_start) for each species as a boundary.
    # At j_start, the incoming flux is the recruitment flux R_dd.
    # The boundary condition for B (LHS) reflects that there is no transport from "below" j_start
    # (other than R_dd which is on RHS).

    j_start <- params@w_min_idx
    idxs <- cbind(1:no_sp, j_start)

    # S_j_start = N_old + dt/dw * R_dd
    S[] <- n
    S[idxs] <- S[idxs] + recruitment_flux * dt / params@dw[j_start]

    # b_j_start = 1 + dt*mu + dt/dw * (g (1 - 0.5 psi_{j_start+1}) + D/(2*dw))
    # This formula excludes the upstream diffusion term D/(2*dw_{j-1}) because
    # the flux from below is replaced by the recruitment flux (psi_{j_start} = 0
    # there). It keeps the high-order term at the upper interface j_start+1.
    psi_start_upper <- psi[cbind(1:no_sp, pmin(j_start + 1, no_w))]
    b[idxs] <- 1 + dt * mu[idxs] +
        dt_dw[idxs] * (g[idxs] * (1 - 0.5 * psi_start_upper) +
                           0.5 * d[idxs] / params@dw[j_start])

    # a_j_start = 0
    a[idxs] <- 0

    # Zero out elements for sizes smaller than w_min_idx
    # This ensures that there is no transport or dynamics below the recruitment size
    # Create a logical mask where col index < w_min_idx
    # Using outer is efficient enough for this size
    w_idx_mat <- matrix(1:no_w, nrow = no_sp, ncol = no_w, byrow = TRUE)
    mask_below <- w_idx_mat < j_start
    
    if (any(mask_below)) {
        a[mask_below] <- 0
        b[mask_below] <- 0
        c[mask_below] <- 0
        S[mask_below] <- 0
    }
    
    return(list(a = a, b = b, c = c, S = S))
}

# Flux limiter functions psi(r) of the smoothness ratio r. To add another
# limiter, give it a name here and offer that name in the `flux_limiter`
# argument of `project()`.
flux_limiter_funcs <- list(
    van_leer = function(r) (r + abs(r)) / (1 + abs(r))
)

#' Frozen flux-limiter values at the size-bin interfaces
#'
#' The first-order upwind advective flux \eqn{g_{j-1} N_{j-1}} carries a
#' numerical diffusion of size \eqn{d_{num}\approx g\,w\,\log\beta}, which
#' dominates the error on coarse logarithmic grids. To reduce it,
#' `get_transport_coefs()` uses the flux-limited (TVD) high-order flux
#' \deqn{J_j = g_{j-1}N_{j-1} + \tfrac12\,\psi_j\,(g_j N_j - g_{j-1}N_{j-1}),}
#' whose high-order target is the centred flux
#' \eqn{\tfrac12(g_{j-1}N_{j-1} + g_j N_j)}. This helper returns the limiter
#' values \eqn{\psi_j = \psi(r_j)} at each interface \eqn{w_j}, computed from a
#' known ("frozen") density field, with \eqn{r_j} the ratio of successive jumps
#' of \eqn{gN}. Because \eqn{\psi_j} is frozen, the flux is linear in the
#' unknowns and the implicit operator stays tridiagonal.
#'
#' The van Leer limiter switches the high-order term off at extrema, which keeps
#' the update positivity-friendly, and \eqn{\psi} is forced to zero at and below
#' the non-smooth recruitment boundary so the advective flux there stays
#' first-order upwind.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param n Density field (species x size) at which to freeze the limiter.
#' @param g Growth rate (species x size).
#' @param flux_limiter Name of the flux limiter, or `"none"`.
#' @return A (species x size) matrix `psi`, where `psi[, j]` is the limiter
#'   value at the interface \eqn{w_j}. Zero everywhere when
#'   `flux_limiter == "none"`.
#' @noRd
flux_limiter_psi <- function(params, n, g, flux_limiter) {
    no_sp <- nrow(n)
    no_w <- ncol(n)
    psi <- matrix(0, nrow = no_sp, ncol = no_w)
    # The limiter needs an upstream and a downstream jump, so it is only defined
    # once there are at least three size bins.
    if (flux_limiter == "none" || no_w < 3) {
        return(psi)
    }
    limiter <- flux_limiter_funcs[[flux_limiter]]

    # Advected flux quantity phi_j = g_j N_j and the jump across interface j
    # (the interface at w_j sits between cells j-1 and j), stored at column j:
    # delta[, j] = phi_j - phi_{j-1}.
    phi <- g * n
    delta <- matrix(0, nrow = no_sp, ncol = no_w)
    delta[, 2:no_w] <- phi[, 2:no_w] - phi[, 1:(no_w - 1)]

    # psi_j = psi(r_j) at interior interfaces j = 3..no_w, with
    # r_j = delta_{j-1} / delta_j. psi_2 stays 0 (no upstream jump).
    jint <- 3:no_w
    num <- delta[, jint - 1, drop = FALSE]   # delta_{j-1}
    den <- delta[, jint, drop = FALSE]       # delta_j
    r <- ifelse(den == 0, 0, num / den)      # no high-order term where jump is 0
    psi[, jint] <- limiter(r)

    # First-order upwind at and below the recruitment boundary.
    j_start <- params@w_min_idx
    w_idx_mat <- matrix(1:no_w, nrow = no_sp, ncol = no_w, byrow = TRUE)
    psi[w_idx_mat <= j_start] <- 0
    psi
}
