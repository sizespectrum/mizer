#' Project values for first time step of Euler method
#' 
#' This is an internal function used by the user-facing `project()` function.
#' It is of potential interest only to mizer extension authors.
#' 
#' @details
#' The function calculates the abundance at the next time step using the
#' McKendrick-von Foerster equation:
#' \deqn{\frac{\partial N}{\partial t} + \frac{\partial}{\partial w} \left( g N - \frac{1}{2}\frac{\partial(D N)}{\partial w} \right) = -\mu N}
#' which is solved using a semi-implicit upwind finite volume scheme.
#' 
#' @param params A \linkS4class{MizerParams} object.
#' @param r A list of rates as returned by `mizerRates()`.
#' @param n An array (species x size) with the number density at the current time step.
#' @param dt Time step.
#' @param a A matrix (species x size) used in the solver (transport term).
#' @param b A matrix (species x size) used in the solver (diagonal term).
#' @param c A matrix (species x size) used in the solver (transport term).
#' @param S A matrix (species x size) used in the solver (source term).
#' @param idx Index vector for size bins (excluding the first one).
#' @param w_min_idx_array_ref Index vector for the start of the size spectrum for each species.
#' @param no_sp Number of species.
#' @param no_w Number of size bins.
#' 
#' @return The updated abundance density matrix `n`.
#' @seealso \code{\link{project}}, \code{\link{mizerRates}}
#' @concept helper
#' @export
project_n <- function(params, r, n, dt, a, b, c, S, idx, w_min_idx_array_ref,
                      no_sp, no_w) {
    # We solve the system A_j * N_{j-1} + B_j * N_j + C_j * N_{j+1} = S_j
    
    # Diffusion coefficient D_i(w)
    d <- params@diffusion # species x size
    
    # Pre-calculate some common terms
    # g_i(w_j)
    g <- r$e_growth
    # mu_i(w_j)
    mu <- r$mort
    # dw_j
    dw <- params@dw
    # Delta t / Delta w_j
    dt_dw <- matrix(dt / dw, nrow = no_sp, ncol = no_w, byrow = TRUE)
    
    # We assume d is 0 at the boundaries for simplicity or it is handled by the loop constraints
    # Actually for j=1, A is 0, for j=no_w, C is 0 (boundary condition)
    
    # A_j = - dt/dw_j * (g_{j-1} + D_{j-1} / (2 * dw_{j-1}))
    # Note: efficient calculation avoiding loop:
    # We compute A for j in idx (2:no_w). g_{j-1} corresponds to columns 1:(no_w-1)
    
    # Indices for j-1 and j
    # idx is 2:no_w
    idx_minus_1 <- idx - 1
    
    term_diff_minus_1 <- 0.5 * d[, idx_minus_1] / 
        matrix(dw[idx_minus_1], nrow = no_sp, ncol = length(idx), byrow = TRUE)
    
    a[, idx] <- -dt_dw[, idx] * (g[, idx_minus_1] + term_diff_minus_1)
    
    # C_j = - dt/dw_j * (D_{j+1} / (2 * dw_j))
    # Note: For j=no_w, we assume N_{j+1}=0, so we don't need C_{no_w} effectively, or it is 0 flux.
    # The equation involves C_j * N_{j+1}. At j=no_w, N_{no_w+1} is 0. So C_{no_w} doesn't matter.
    # We can compute C for j=1:(no_w-1).
    # idx_plus_1 is 2:no_w
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
    # S_1 = N_1^old + dt/dw_1 * R_dd
    
    dw_1 <- dw[1]
    dt_dw_1 <- dt / dw_1
    b[, 1] <- 1 + dt * mu[, 1] + dt_dw_1 * (g[, 1] + 0.5 * d[, 1] / dw_1)
    # c[, 1] already computed correctly above
    a[, 1] <- 0
    
    # RHS S
    S[] <- n
    # Add recruitment to S[, 1] for each species
    # We need to distribute R_dd correctly.
    # r$rdd is vector of length no_sp.
    # "r$rdd * dt / params@dw[params@w_min_idx]"
    # Wait, the prompt says "w_min_idx_array_ref" handles the start index.
    # Different species can have different w_min_idx.
    # So the "j=1" above really refers to w_min_idx[i].
    # But currently 'a' and 'b' and 'c' are computed for the whole grid.
    # Species i only exists from w_min_idx[i] to w_max_idx[i] (implicitly).
    # We should iterate the Thomas algorithm for each species from w_min_idx[i] to no_w.
    
    # Let's adjust S for the recruitment term at the start index for each species.
    # S[w_min_idx_array_ref] <- S[w_min_idx_array_ref] + r$rdd * dt / params@dw[params@w_min_idx]
    # But wait, we need to respect the diffusion B_1 term calculation which might be different at the boundary.
    # The B calculation above `b[, idx]` used `idx` which starts at 2.
    # If w_min_idx[i] > 1, then the "standard" formula for B at w_min_idx[i] might be using w_{j-1} which is essentially 0 since N is 0 there?
    # Actually, if we assume N is 0 below w_min_idx, then the flux from below is just the recruitment.
    # So for j = w_min_idx[i], the term A_j should be effectively 0 (or replaced by recruitment boundary condition).
    # The B_j term should not include diffusion from below (or handled as boundary).
    # We can handle this by essentially running the solver from w_min_idx[i].
    
    # Correct B matrix for start indices
    # We need to loop over species to correct B at w_min_idx, because vectorization is hard with variable indices.
    # However, w_min_idx might be same for all species or not.
    
    # Loop over species to apply boundary condition and solve
    
    # Temporary copy of C for modification during Thomas algo
    c_prime <- c 
    # Temporary copy of S (d in Thomas algo context)
    d_prime <- S
    
    for (i in 1:no_sp) {
        # Start index for this species
        j_start <- params@w_min_idx[i]
        
        # Apply boundary condition to S (RHS)
        # S_j_start = N_old + dt/dw * R_dd
        d_prime[i, j_start] <- d_prime[i, j_start] + r$rdd[i] * dt / dw[j_start]
        
        # Apply boundary condition to B (LHS)
        # Remove the influence of "below" diffusion/growth which is replaced by recruitment flux
        # B_j = 1 + ... + dt/dw * (g_j + D_j/(2*dw_j) + D_j/(2*dw_{j-1}))
        # The D_j/(2*dw_{j-1}) term came from flux J_j.
        # At boundary, J_j is explicitly R_dd.
        # So we should remove the D_j/(2*dw_{j-1}) term from B[i, j_start]
        # and A[i, j_start] is 0.
        
        # Original B[i, j_start] contains: ... + dt/dw * ( ... + 0.5 * d[i, j_start] / dw[j_start-1])
        # We need to subtract that term
        if (j_start > 1) {
             # The term added was: dt/dw[j_start] * 0.5 * d[i, j_start] / dw[j_start-1]
             correction <- (dt / dw[j_start]) * 0.5 * d[i, j_start] / dw[j_start - 1]
             b[i, j_start] <- b[i, j_start] - correction
             
             # Also A[i, j_start] should be ignored/0.
             a[i, j_start] <- 0
        } else {
             # j_start == 1. The formula used 'dw[1]' for prev bin width in approximation?
             # My code above: `b[, 1] <- ... + g[, 1] + 0.5 * d[, 1] / dw_1`
             # This formula for j=1 was already "boundary-like" (omitted term from below).
             # So if j_start == 1, B is already correct for the boundary condition derived.
        }
        
        # Thomas Algorithm: Forward Elimination
        # For j = j_start
        # c'[j] = c[j] / b[j]
        # d'[j] = d[j] / b[j]
        
        c_prime[i, j_start] <- c_prime[i, j_start] / b[i, j_start]
        d_prime[i, j_start] <- d_prime[i, j_start] / b[i, j_start]
        
        # Loop j from j_start+1 to no_w
        if (j_start < no_w) {
            for (j in (j_start + 1):no_w) {
                temp <- b[i, j] - a[i, j] * c_prime[i, j - 1]
                if (j < no_w) {
                    c_prime[i, j] <- c_prime[i, j] / temp
                }
                d_prime[i, j] <- (d_prime[i, j] - a[i, j] * d_prime[i, j - 1]) / temp
            }
        }
        
        # Backward Substitution
        n[i, no_w] <- d_prime[i, no_w]
        if (no_w > j_start) {
            for (j in (no_w - 1):j_start) {
                n[i, j] <- d_prime[i, j] - c_prime[i, j] * n[i, j + 1]
            }
        }
        # Species density is 0 below j_start? Mizer usually keeps it 0 or doesn't update.
        # The loop range updates n for j >= j_start.
    }
    
    n
}