#' Advance abundance density by one time step
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