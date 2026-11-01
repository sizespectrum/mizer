#' Determine reproduction rate needed for initial egg abundance
#'
#' @param params A MizerParams object
#' @param ... Unused.
#' @return A vector of reproduction rates for all species
#' @seealso [getR0()] for the lifetime reproductive output, which is
#'   `getRDI(params) / getRequiredRDD(params)` at a steady state.
#' @export
getRequiredRDD <- function(params, ...) {
    UseMethod("getRequiredRDD")
}

#' @export
getRequiredRDD.MizerParams <- function(params, ...) {
    # The advective-flux scheme is read from the model's second_order_w slot, so
    # that the egg-cell balance uses the same flux as project() will. The
    # flux-limited correction changes the flux leaving the egg size class and
    # hence the reproduction required to hold the egg density steady.
    flux_limiter <- flux_limiter_scheme(params)
    # Calculate required rdd
    no_sp <- nrow(params@species_params)

    # Calculate transport coefficients
    dt <- 1
    # Compute EGrowth, Mort, and Diffusion together so each upstream rate
    # (Encounter, FeedingLevel, PredMort, ...) is computed only once.
    rates_fns <- projectRateFunctions(params)
    r <- mizer_rates_subset(params,
                             n = params@initial_n,
                             n_pp = params@initial_n_pp,
                             n_other = params@initial_n_other,
                             t = 0, effort = params@initial_effort,
                             rates_fns = rates_fns,
                             targets = c("EGrowth", "Mort", "Diffusion"))
    # We pass a dummy recruitment flux of 0 to trigger the boundary condition
    # corrections for a and b in get_transport_coefs
    coefs <- get_transport_coefs(params, n = params@initial_n,
                                 g = r$e_growth,
                                 mu = r$mort, dt,
                                 recruitment_flux = numeric(no_sp),
                                 d = r$diffusion,
                                 flux_limiter = flux_limiter)
    # The flux limiter is folded into the coefficients a, b, c above, so the
    # egg-cell balance already accounts for the high-order flux leaving the egg
    # size class. This keeps getRequiredRDD() consistent with a project() call
    # on the same params.
    required_recruitment(params, params@initial_n, coefs, dt)
}

#' Recruitment that holds the abundance in the egg size class steady
#'
#' The time-stepper solves, at the egg size class \eqn{j} of each species,
#' \deqn{A_j N_{j-1} + B_j N_j + C_j N_{j+1} = N^{old}_j + R\, dt / \Delta w_j.}
#' In a steady state \eqn{N_j = N^{old}_j}, so the recruitment \eqn{R} that
#' holds the abundance `n` steady there is
#' \deqn{R = (A_j N_{j-1} + (B_j - 1) N_j + C_j N_{j+1})\, \Delta w_j / dt.}
#' This is the one place that this balance is written, shared by
#' [getRequiredRDD()] and [getNPerRecruit()].
#'
#' @param params A MizerParams object.
#' @param n An array (species x size) with the abundance to hold steady.
#' @param coefs The transport coefficients from `get_transport_coefs()`,
#'   calculated at `n` with time step `dt` and zero recruitment flux, so that
#'   the boundary corrections for the egg size class are applied.
#' @param dt The time step that `coefs` were calculated with.
#' @return A named vector with the required recruitment for each species.
#' @noRd
required_recruitment <- function(params, n, coefs, dt = 1) {
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)
    j <- params@w_min_idx
    idx <- cbind(seq_len(no_sp), j)
    # A at the egg size class is zero, but keep the term so that this is the
    # balance exactly as the time-stepper solves it.
    n_prev <- ifelse(j > 1, n[cbind(seq_len(no_sp), pmax(j - 1, 1))], 0)
    n_next <- ifelse(j < no_w, n[cbind(seq_len(no_sp), pmin(j + 1, no_w))], 0)

    total_rate <- coefs$a[idx] * n_prev + (coefs$b[idx] - 1) * n[idx] +
        coefs$c[idx] * n_next
    reproduction <- total_rate * params@dw[j] / dt
    names(reproduction) <- params@species_params$species
    reproduction
}
