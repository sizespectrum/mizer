#' Helper function to calculate the transport coefficients for the upwind-difference scheme
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param n An array (species x size) with the number density at the current time step.
#' @param g The growth rate.
#' @param mu The mortality rate.
#' @param dt Time step.
#'
#' @return A list with the coefficients A, B, C and S.
#' @noRd
get_transport_coefs <- function(params, n, g, mu, dt) {
    # We solve the system A_j * N_{j-1} + B_j * N_j + C_j * N_{j+1} = S_j
    
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
    c[, no_w] <- 0 # Boundary condition N_{no_w+1} = 0
    
    # B_j
    # B_j = 1 + dt * mu_j + dt/dw_j * (g_j + D_j / (2 * dw_j) + D_j / (2 * dw_{j-1}))
    # Careful with j=1 term for D_j / (2 * dw_{j-1}). dw_0 is not defined.
    # At j=1 boundary condition comes from recruitment.
    # Standard formula works for j > 1.
    
    term_diff_j <- 0.5 * d[, idx] / matrix(dw[idx], nrow = no_sp, ncol = length(idx), byrow = TRUE)
    term_diff_j_minus_1 <- 0.5 * d[, idx] / matrix(dw[idx_minus_1], nrow = no_sp, ncol = length(idx), byrow = TRUE)
    
    b[, idx] <- 1 + dt * mu[, idx] + dt_dw[, idx] * (g[, idx] + term_diff_j + term_diff_j_minus_1)
    
    # Boundary j=1 (approx w_min)
    # Equation: N_1 - N_1^old / dt + (J_{3/2} - J_{1/2}) / dw_1 = -mu_1 N_1
    # J_{1/2} = R_dd (recruitment flux).
    # J_{3/2} follows standard flux definition.
    # result: (1 + dt*mu_1 + dt/dw_1 * (g_1 + D_1/(2*dw_1))) N_1 - dt/dw_1 * (D_2/(2*dw_1)) N_2 = N_1^old + dt/dw_1 * R_dd
    # So for j=1:
    # B_1 = 1 + dt*mu_1 + dt/dw_1 * (g_1 + d_1/(2*dw_1))
    # C_1 = - dt/dw_1 * d_2/(2*dw_1)  (Matches standard formula)
    # A_1 = 0
    
    dw_1 <- dw[1]
    dt_dw_1 <- dt / dw_1
    b[, 1] <- 1 + dt * mu[, 1] + dt_dw_1 * (g[, 1] + 0.5 * d[, 1] / dw_1)
    # c[, 1] already computed correctly above
    a[, 1] <- 0

    # RHS S
    S[] <- n
    
    return(list(a = a, b = b, c = c, S = S))
}
