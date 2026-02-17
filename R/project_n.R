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
    coefs <- get_transport_coefs(params, n, r, dt)
    a <- coefs$a
    b <- coefs$b
    c <- coefs$c
    S <- coefs$S
    
    # Boundary condition updates
    # Vectorized implementation of the loop over species
    
    # S_j_start = N_old + dt/dw * R_dd
    # Matrix of indices for (i, j_start)
    # w_min_idx is 1-based index of start size bin for each species
    j_start <- params@w_min_idx
    idxs <- cbind(1:no_sp, j_start)
    
    S[idxs] <- S[idxs] + r$rdd * dt / params@dw[j_start]
    
    # Apply boundary condition to B (LHS)
    # Remove the influence of "below" diffusion/growth which is replaced by recruitment flux
    # if (j_start > 1) { ... }
    
    # Identify species where start index > 1
    mask <- j_start > 1
    if (any(mask)) {
        i_sub <- which(mask)
        j_sub <- j_start[i_sub]
        idxs_sub <- cbind(i_sub, j_sub)
        
        # correction <- (dt / params@dw[j_start]) * 0.5 * params@diffusion[i, j_start] / params@dw[j_start - 1]
        correction <- (dt / params@dw[j_sub]) * 0.5 * 
            params@diffusion[idxs_sub] / params@dw[j_sub - 1]
        
        b[idxs_sub] <- b[idxs_sub] - correction
        a[idxs_sub] <- 0
    }
    
    # Call C++ function to solve tridiagonal system
    n <- project_n_loop(n, a, b, c, S, j_start)
    
    n
}

#' @rdname project_n
project_n_diffusion_R <- function(params, r, n, dt, a, b, c, S, idx, w_min_idx_array_ref,
                                  no_sp, no_w) {
    coefs <- get_transport_coefs(params, n, r, dt)
    a <- coefs$a
    b <- coefs$b
    c <- coefs$c
    S <- coefs$S
    
    # Loop over species to apply boundary condition and solve
    
    for (i in 1:no_sp) {
        # Start index for this species
        j_start <- params@w_min_idx[i]
        
        # Apply boundary condition to S (RHS)
        # S_j_start = N_old + dt/dw * R_dd
        S[i, j_start] <- S[i, j_start] + r$rdd[i] * dt / params@dw[j_start]
        
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
        }
        
        # Thomas Algorithm
        # We need to pass the sub-vectors for the current species i, starting from j_start
        # We are solving for n[i, j_start:no_w]
        
        # Extract the relevant parts of the vectors
        # Note: thomas_solve accepts vectors of length N
        # a, b, c, d are vectors of length N
        
        # Correctly slicing from j_start to no_w
        # a[i, j_start] is effectively 0 or ignored by thomas_solve if it's the first element passed
        # c[i, no_w] is 0 or ignored by thomas_solve
        
        # We solve for the segment of the size spectrum inhabited by the species
        relevant_indices <- j_start:no_w
        
        n[i, relevant_indices] <- thomas_solve(
            a = a[i, relevant_indices],
            b = b[i, relevant_indices],
            c = c[i, relevant_indices],
            d = S[i, relevant_indices]
        )
    }
    
    n
}

project_n_no_diffusion <- function(params, r, n, dt, a, b, S, idx, w_min_idx_array_ref,
                                   no_sp, no_w) {
    # a_{ij} = - g_i(w_{j-1}) / dw_j dt
    a[, idx] <- sweep(
        -r$e_growth[, idx - 1, drop = FALSE] * dt, 2,
        params@dw[idx], "/"
    )
    # b_{ij} = 1 + g_i(w_j) / dw_j dt + \mu_i(w_j) dt
    b[] <- 1 + sweep(r$e_growth * dt, 2, params@dw, "/") + r$mort * dt
    # S_{ij} <- N_i(w_j)
    S[, idx] <- n[, idx, drop = FALSE]
    # Update first size group of n
    n[w_min_idx_array_ref] <-
        (n[w_min_idx_array_ref] + r$rdd * dt /
             params@dw[params@w_min_idx]) /
        b[w_min_idx_array_ref]
    # Update n
    # for (i in 1:no_sp) # number of species assumed small, so no need to
    #                      vectorize this loop over species
    #     for (j in (params@w_min_idx[i]+1):no_w)
    #         n[i,j] <- (S[i,j] - A[i,j]*n[i,j-1]) / B[i,j]
    # This is implemented via Rcpp
    n <- inner_project_loop(
        no_sp = no_sp, no_w = no_w, n = n,
        A = a, B = b, S = S,
        w_min_idx = params@w_min_idx
    )
    n
}