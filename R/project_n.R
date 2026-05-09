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
                                 recruitment_flux = r$rdd, d = r$diffusion)
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

#' Project values with a predictor-corrector method
#'
#' This is an experimental second-order time stepping variant of [project_n()].
#' It first predicts the new consumer densities with [project_n()], optionally
#' recalculates rates from that prediction, and then applies a Crank-Nicolson
#' corrector using midpoint rates.
#'
#' If the rate recalculation arguments are not supplied, the corrector uses the
#' supplied rates as fixed rates. In that case the corrector is second order only
#' for the frozen-rate transport problem, not for the full nonlinear mizer
#' dynamics.
#'
#' @inheritParams project_n
#' @param rates_fns Optional named list of rate functions, as used by
#'   [mizerRates()]. If supplied together with `n_pp`, `n_other` and `effort`,
#'   provisional end-of-step rates are calculated from the predicted densities.
#' @param n_pp Resource abundance used when recalculating provisional rates.
#' @param n_other Other ecosystem components used when recalculating provisional
#'   rates.
#' @param t Current time.
#' @param effort Fishing effort used when recalculating provisional rates.
#' @param r_hat Optional provisional end-of-step rates. If supplied, these are
#'   used instead of recalculating them.
#' @param r_mid Optional midpoint rates. If supplied, these are used directly in
#'   the Crank-Nicolson corrector.
#' @param ... Further arguments passed to the rate functions.
#'
#' @return The updated abundance density matrix `n`.
#' @seealso \code{\link{project_n}}
#' @concept helper
#' @export
project_n_2 <- function(params, r, n, dt, a, b, c, S, idx,
                        w_min_idx_array_ref, no_sp, no_w, rates_fns = NULL,
                        n_pp = NULL, n_other = NULL, t = 0, effort = NULL,
                        r_hat = NULL, r_mid = NULL, ...) {
    if (is.null(r_mid)) {
        if (is.null(r_hat) && !is.null(rates_fns) && !is.null(n_pp) &&
            !is.null(n_other) && !is.null(effort)) {
            n_hat <- project_n(params, r, n, dt, a, b, c, S, idx,
                               w_min_idx_array_ref, no_sp, no_w)
            r_hat <- rates_fns$Rates(
                params,
                n = n_hat, n_pp = n_pp, n_other = n_other,
                t = t + dt, effort = effort, rates_fns = rates_fns, ...
            )
        }

        if (is.null(r_hat)) {
            r_mid <- r
        } else {
            r_mid <- r
            r_mid$e_growth <- 0.5 * (r$e_growth + r_hat$e_growth)
            r_mid$mort <- 0.5 * (r$mort + r_hat$mort)
            r_mid$diffusion <- 0.5 * (r$diffusion + r_hat$diffusion)
            r_mid$rdd <- 0.5 * (r$rdd + r_hat$rdd)
        }
    }

    coefs <- get_transport_coefs(params, n, r_mid$e_growth, r_mid$mort,
                                 dt / 2, recruitment_flux = r_mid$rdd,
                                 d = r_mid$diffusion)
    a <- coefs$a
    b <- coefs$b
    c <- coefs$c
    S <- project_n_2_rhs(params, n, a, b, c, dt, r_mid$rdd)

    project_n_loop(n, a, b, c, S, params@w_min_idx)
}

project_n_2_rhs <- function(params, n, a, b, c, dt, recruitment_flux) {
    no_sp <- nrow(n)
    no_w <- ncol(n)

    matrix_n <- b * n
    idx <- 2:no_w
    matrix_n[, idx] <- matrix_n[, idx] + a[, idx] * n[, idx - 1, drop = FALSE]
    idx <- 1:(no_w - 1)
    matrix_n[, idx] <- matrix_n[, idx] + c[, idx] * n[, idx + 1, drop = FALSE]

    source <- matrix(0, nrow = no_sp, ncol = no_w, dimnames = dimnames(n))
    j_start <- params@w_min_idx
    idxs <- cbind(seq_len(no_sp), j_start)
    source[idxs] <- recruitment_flux / params@dw[j_start]

    2 * n + dt * source - matrix_n
}

#' @rdname project_n
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
