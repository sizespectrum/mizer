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
    coefs <- get_transport_coefs(params, n, n_pp, n_other, r, dt)
    a <- coefs$a
    b <- coefs$b
    c <- coefs$c
    S <- coefs$S

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
        d_prime[i, j_start] <- d_prime[i, j_start] + r$rdd[i] * dt / params@dw[j_start]
        
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
             correction <- (dt / params@dw[j_start]) * 0.5 * params@diffusion[i, j_start] / params@dw[j_start - 1]
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