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

# Average the start-of-step rates `r` with the provisional end-of-step rates
# `r_hat` to obtain second-order-accurate midpoint rates. Only the fields that
# are consumed by the density and resource updates are averaged: the consumer
# transport fields and the resource mortality.
average_rates <- function(r, r_hat) {
    for (field in c("e_growth", "mort", "diffusion", "rdd", "resource_mort")) {
        r[[field]] <- 0.5 * (r[[field]] + r_hat[[field]])
    }
    r
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
            r_mid <- average_rates(r, r_hat)
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

#' Project values with the TR-BDF2 method
#'
#' This is an L-stable, second-order time stepping variant of [project_n()]. It
#' takes one TR-BDF2 step, consisting of a trapezoidal (Crank-Nicolson) stage
#' over the first part of the time step followed by a second-order backward
#' differentiation (BDF2) stage over the remainder.
#'
#' The nonlinear rates are handled exactly as in [project_n_2()]: a provisional
#' Euler predictor gives end-of-step rates, which are averaged with the
#' start-of-step rates to obtain second-order-accurate midpoint rates `r_mid`.
#' Both TR-BDF2 stages then use this single frozen operator.
#'
#' With the standard parameter \eqn{\gamma = 2 - \sqrt 2} both stages share the
#' same implicit coefficient \eqn{\alpha\,\Delta t} with
#' \eqn{\alpha = \gamma/2 = 1 - 1/\sqrt 2}, so the operator
#' \eqn{I - \alpha\,\Delta t\,L} is assembled once with `get_transport_coefs()`
#' and each stage is a single tridiagonal solve with `project_n_loop()`, exactly
#' as in [project_n()] and [project_n_2()]. Unlike the Crank-Nicolson corrector
#' in [project_n_2()], TR-BDF2 is L-stable and therefore damps the stiff modes
#' that cause Crank-Nicolson to oscillate at large time steps.
#'
#' If the rate recalculation arguments are not supplied, the step uses the
#' supplied rates as fixed rates. In that case the method is second order only
#' for the frozen-rate transport problem, not for the full nonlinear mizer
#' dynamics, but it remains L-stable.
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
#' @param r_mid Optional midpoint rates. If supplied, these are used directly to
#'   build the TR-BDF2 operator.
#' @param ... Further arguments passed to the rate functions.
#'
#' @return The updated abundance density matrix `n`.
#' @seealso \code{\link{project_n}}, \code{\link{project_n_2}}
#' @concept helper
#' @export
project_n_tr_bdf2 <- function(params, r, n, dt, a, b, c, S, idx,
                              w_min_idx_array_ref, no_sp, no_w,
                              rates_fns = NULL, n_pp = NULL, n_other = NULL,
                              t = 0, effort = NULL, r_hat = NULL, r_mid = NULL,
                              ...) {
    gamma <- 2 - sqrt(2)
    alpha <- 1 - 1 / sqrt(2)   # = gamma / 2; shared implicit coefficient
    c1 <- (sqrt(2) + 1) / 2
    c0 <- (sqrt(2) - 1) / 2    # note c1 - c0 = 1

    # Determine the frozen midpoint rates, exactly as in project_n_2().
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
            r_mid <- average_rates(r, r_hat)
        }
    }

    # Build the shared implicit operator (I - alpha * dt * L) once, then take
    # the two TR-BDF2 stages as tridiagonal solves against it.
    coefs <- get_transport_coefs(params, n, r_mid$e_growth, r_mid$mort,
                                 dt * alpha, recruitment_flux = r_mid$rdd,
                                 d = r_mid$diffusion)

    # Stage 1: trapezoidal step over gamma * dt. Its right-hand side is the
    # Crank-Nicolson form, reusing project_n_2_rhs() with sub-step gamma * dt.
    S1 <- project_n_2_rhs(params, n, coefs$a, coefs$b, coefs$c,
                          dt * gamma, r_mid$rdd)
    n_gamma <- project_n_loop(n, coefs$a, coefs$b, coefs$c, S1,
                              params@w_min_idx)

    # Stage 2: BDF2 step using the same operator.
    S2 <- project_n_tr_bdf2_rhs(params, n, n_gamma, dt * alpha, r_mid$rdd,
                                c1, c0)
    project_n_loop(n_gamma, coefs$a, coefs$b, coefs$c, S2, params@w_min_idx)
}

project_n_tr_bdf2_rhs <- function(params, n, n_gamma, dt_eff, recruitment_flux,
                                  c1, c0) {
    no_sp <- nrow(n)

    # BDF2 right-hand side: c1 * N^{t+gamma} - c0 * N^t, with the recruitment
    # source added at the lower boundary.
    rhs <- c1 * n_gamma - c0 * n

    j_start <- params@w_min_idx
    idxs <- cbind(seq_len(no_sp), j_start)
    rhs[idxs] <- rhs[idxs] + dt_eff * recruitment_flux / params@dw[j_start]

    rhs
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
