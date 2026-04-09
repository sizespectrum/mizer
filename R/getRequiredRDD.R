#' Determine reproduction rate needed for initial egg abundance
#'
#' @param params A MizerParams object
#' @return A vector of reproduction rates for all species
#' @export
getRequiredRDD <- function(params) {
    UseMethod("getRequiredRDD")
}

#' @export
getRequiredRDD.MizerParams <- function(params) {
    # Calculate required rdd
    no_sp <- nrow(params@species_params)

    # Calculate transport coefficients
    dt <- 1
    # We pass a dummy recruitment flux of 0 to trigger the boundary condition
    # corrections for a and b in get_transport_coefs
    coefs <- get_transport_coefs(params, n = params@initial_n,
                                 g = getEGrowth(params),
                                 mu = getMort(params), dt,
                                 recruitment_flux = numeric(no_sp))

    reproduction <- params@species_params$erepro # vector of correct length
    names(reproduction) <- params@species_params$species

    for (i in (1:no_sp)) {
        w_min_idx <- params@w_min_idx[i]

        # Get coefficients for this species at the boundary
        # The equation for the first node is:
        # (N_new - N_old)/dt = -(Flux_matrix * N) + R/dw
        # In steady state N_new = N_old, so:
        # Flux_matrix * N = R/dw
        # The rows of coefs correspond to the linear system A*N_{j-1} + B*N_j + C*N_{j+1} = ...
        # For the first node j=w_min_idx:
        # A*N_{j-1} + (B-1)/dt * N_j + C/dt * N_{j+1} = R/dw / dt ?
        # No, let's look at project_n again.
        # It solves A N_{i-1} + B N_i + C N_{i+1} = N_old + RHS_source
        # In steady state: A N_{i-1} + B N_i + C N_{i+1} = N_i + R * dt / dw
        # So R = ( A N_{i-1} + (B-1) N_i + C N_{i+1} ) * dw / dt

        # Extract coefficients
        a <- coefs$a[i, w_min_idx]
        b <- coefs$b[i, w_min_idx]
        c <- coefs$c[i, w_min_idx]

        # Boundary corrections for a and b are now handled in get_transport_coefs

        n_current <- params@initial_n[i, w_min_idx]
        n_next <- if (w_min_idx < length(params@w)) params@initial_n[i, w_min_idx + 1] else 0
        n_prev <- if (w_min_idx > 1) params@initial_n[i, w_min_idx - 1] else 0 # Should be irrelevant if A=0 or boundary

        # Calculate R
        # R = ( A * n_prev + (B - 1) * n_current + C * n_next ) * dw / dt

        total_rate <- a * n_prev + (b - 1) * n_current + c * n_next
        reproduction[i] <- total_rate * params@dw[w_min_idx] / dt
    }
    reproduction
}
