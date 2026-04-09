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
get_transport_coefs <- function(params, n, g, mu, dt, recruitment_flux) {
    
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)
    
    # Diffusion coefficient D_i(w)
    d <- params@diffusion # species x size
    
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
    
    # A_j = - dt/dw_j * (g_{j-1} + D_{j-1} / (2 * dw_{j-1}))
    # Note: efficient calculation avoiding loop:
    # We compute A for j in idx (2:no_w). g_{j-1} corresponds to columns 1:(no_w-1)
    
    # Indices for j-1 and j
    idx <- 2:no_w
    idx_minus_1 <- idx - 1
    
    term_diff_minus_1 <- 0.5 * d[, idx_minus_1] / 
        matrix(dw[idx_minus_1], nrow = no_sp, ncol = length(idx), byrow = TRUE)
    
    a[, idx] <- -dt_dw[, idx] * (g[, idx_minus_1] + term_diff_minus_1)
    
    # C_j = - dt/dw_j * (D_{j+1} / (2 * dw_j))
    # Note: For j=no_w, we assume N_{j+1}=0, so we don't need C_{no_w} effectively, or it is 0 flux.
    # The equation involves C_j * N_{j+1}. At j=no_w, N_{no_w+1} is 0. So C_{no_w} doesn't matter.
    # We can compute C for j=1:(no_w-1).
    idx_j <- 1:(no_w - 1)
    term_diff_plus_1 <- 0.5 * d[, idx_j + 1] / 
        matrix(dw[idx_j], nrow = no_sp, ncol = length(idx_j), byrow = TRUE)
    
    c[, idx_j] <- -dt_dw[, idx_j] * term_diff_plus_1
    # c[, no_w] is 0 as initialized
    
    # B_j
    # B_j = 1 + dt * mu_j + dt/dw_j * (g_j + D_j / (2 * dw_j) + D_j / (2 * dw_{j-1}))
    # Careful with j=1 term for D_j / (2 * dw_{j-1}). dw_0 is not defined.
    # At j=1 boundary condition comes from recruitment.
    # Standard formula works for j > 1.
    
    term_diff_j <- 0.5 * d[, idx] / matrix(dw[idx], nrow = no_sp, ncol = length(idx), byrow = TRUE)
    term_diff_j_minus_1 <- 0.5 * d[, idx] / matrix(dw[idx_minus_1], nrow = no_sp, ncol = length(idx), byrow = TRUE)
    
    b[, idx] <- 1 + dt * mu[, idx] + dt_dw[, idx] * (g[, idx] + term_diff_j + term_diff_j_minus_1)
    
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
    
    # b_j_start = 1 + dt*mu + dt/dw * (g + D/(2*dw))
    # This formula excludes the upstream diffusion term D/(2*dw_{j-1}) because 
    # flux from below is replaced by recruitment flux.
    b[idxs] <- 1 + dt * mu[idxs] + dt_dw[idxs] * (g[idxs] + 0.5 * d[idxs] / params@dw[j_start])
    
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
