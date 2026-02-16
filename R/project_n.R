#' Project values for first time step of Euler method
#' 
#' This is an internal function used by the user-facing `project()` function.
#' It is of potential interest only to mizer extension authors.
#' 
#' @details
#' The function calculates the abundance at the next time step using the
#' McKendrick-von Foerster equation:
#' \deqn{\frac{\partial N}{\partial t} + \frac{\partial g N}{\partial w} = -\mu N}
#' which is solved using a semi-implicit upwind finite difference scheme.
#' 
#' @param params A \linkS4class{MizerParams} object.
#' @param r A list of rates as returned by `mizerRates()`.
#' @param n An array (species x size) with the number density at the current time step.
#' @param dt Time step.
#' @param a A matrix (species x size) used in the solver (transport term).
#' @param b A matrix (species x size) used in the solver (diagonal term).
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
project_n <- function(params, r, n, dt, a, b, S, idx, w_min_idx_array_ref,
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