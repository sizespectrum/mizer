#' Determine reproduction rate needed for initial egg abundance
#'
#' @param params A MizerParams object
#' @param ... Arguments passed to methods.
#' @param rates Optional already evaluated rate list for this model's initial
#'   state and effort, used to share an effective evaluation with a caller.
#'   Must contain growth, total mortality and diffusion. Usually omitted.
#' @return A vector of reproduction rates for all species
#' @export
getRequiredRDD <- function(params, ..., rates = NULL) {
    UseMethod("getRequiredRDD")
}

#' @rdname getRequiredRDD
#' @usage NULL
#' @export
getRequiredRDD.MizerParams <- function(params, ..., rates = NULL) {
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
    r <- rates
    if (is.null(r)) {
        rates_fns <- projectRateFunctions(params)
        method <- NULL
        for (class_name in c(class(params), "default")) {
            method <- getS3method("projectRates", class_name, optional = TRUE)
            if (!is.null(method)) break
        }
        provider <- rate_callback_provider(rates_fns$Rates)
        standard <- identical(provider, mizerRates) ||
            (identical(provider, projectRates) &&
                 identical(method, projectRates.MizerParams))
        args <- list(params, n = params@initial_n,
                     n_pp = params@initial_n_pp,
                     n_other = params@initial_n_other, t = 0,
                     effort = params@initial_effort, rates_fns = rates_fns)
        if (standard) {
            # The default dependency graph can avoid reproduction itself:
            # a custom RDI may legitimately use this egg-balance diagnostic.
            args$targets <- c("EGrowth", "Mort", "Diffusion")
            r <- do.call(mizer_rates_subset, args)
        } else {
            # A full-rate replacement may change these final rates after their
            # individual providers have run. Its effective output is decisive.
            reentered <- sum(vapply(seq_len(sys.nframe()), function(frame) {
                identical(sys.function(frame), getRequiredRDD.MizerParams)
            }, logical(1))) > 1L
            if (reentered) {
                stop("Recursive full-rate evaluation in getRequiredRDD(); ",
                     "supply the already evaluated growth, mortality and ",
                     "diffusion through `rates`.", call. = FALSE)
            }
            r <- do.call(rates_fns$Rates, args)
        }
    }
    for (field in c("e_growth", "mort", "diffusion")) {
        value <- r[[field]]
        if (!is.numeric(value) || is.complex(value) || any(!is.finite(value)) ||
                !identical(dim(value), dim(params@initial_n))) {
            stop("Required recruitment needs a finite species-by-size `",
                 field, "` rate.", call. = FALSE)
        }
    }
    # We pass a dummy recruitment flux of 0 to trigger the boundary condition
    # corrections for a and b in get_transport_coefs
    coefs <- get_transport_coefs(params, n = params@initial_n,
                                 g = r$e_growth,
                                 mu = r$mort, dt,
                                 recruitment_flux = numeric(no_sp),
                                 d = r$diffusion,
                                 flux_limiter = flux_limiter)

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
    # The flux limiter is folded into the coefficients a, b, c above, so the
    # egg-cell balance already accounts for the high-order flux leaving the egg
    # size class. This keeps getRequiredRDD() consistent with a project() call
    # on the same params.
    reproduction
}
