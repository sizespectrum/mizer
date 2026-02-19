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
    coefs <- get_transport_coefs(params, n, r$e_growth, r$mort, dt,
                                 recruitment_flux = r$rdd)
    a <- coefs$a
    b <- coefs$b
    c <- coefs$c
    S <- coefs$S
    
    # Call C++ function to solve tridiagonal system
    # j_start is needed for the C++ loop, we can get it from params
    params@w_min_idx
    
    # Note: project_n_loop takes j_start as argument
    n <- project_n_loop(n, a, b, c, S, params@w_min_idx)
    
    n
}

#' @rdname project_n
project_n_diffusion_R <- function(params, r, n, dt, a, b, c, S, idx, w_min_idx_array_ref,
                                  no_sp, no_w) {
    coefs <- get_transport_coefs(params, n, r$e_growth, r$mort, dt,
                                 recruitment_flux = r$rdd)
    a <- coefs$a
    b <- coefs$b
    c <- coefs$c
    S <- coefs$S
    
    # Loop over species to solve
    
    for (i in 1:no_sp) {
        # Start index for this species
        j_start <- params@w_min_idx[i]
        
        # Boundary conditions are handled in get_transport_coefs
        
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