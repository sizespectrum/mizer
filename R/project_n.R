#' Advance species abundance densities by one time step
#'
#' Update the consumer abundance density matrix by one explicit time step of the
#' semi-implicit solver used in [project_simple()] and [project()]. The first
#' occupied size class of each species is updated from the density-dependent
#' reproduction rate `r$rdd`, and all larger size classes are then advanced with
#' the tridiagonal update implemented in `inner_project_loop()`.
#'
#' @param params A [MizerParams-class()] object.
#' @param r A list of rates as returned by [getRates()] or [mizerRates()]. This
#'   function uses the `e_growth`, `mort`, and `rdd` entries.
#' @param n A two-dimensional array (species x size) with the current consumer
#'   abundance densities.
#' @param dt The time step in years.
#' @param a A matrix with the same dimensions as `n`, used as workspace for the
#'   lower diagonal coefficients of the semi-implicit update.
#' @param b A matrix with the same dimensions as `n`, used as workspace for the
#'   diagonal coefficients of the semi-implicit update.
#' @param S A matrix with the same dimensions as `n`, used as workspace for the
#'   right-hand side of the semi-implicit update.
#' @param idx Integer indices of the non-egg size classes to be updated with the
#'   tridiagonal recursion, typically `2:no_w`.
#' @param w_min_idx_array_ref Integer indices for one-dimensional indexing into
#'   `n[, params@w_min_idx]`, one entry per species.
#' @param no_sp Number of species, equal to `nrow(n)`.
#' @param no_w Number of consumer size classes, equal to `ncol(n)`.
#'
#' @return A two-dimensional array (species x size) with the updated consumer
#'   abundance densities after one time step.
#'
#' @keywords internal
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
