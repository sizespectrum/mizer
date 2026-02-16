#' Set initial abundances to single-species steady state abundances
#'
#' `r lifecycle::badge("experimental")`
#' This first calculates growth and death rates that arise from the current
#' initial abundances. Then it uses these growth and death rates to
#' determine the steady-state abundances of the selected species.
#'
#' The result of applying this function is of course not a multi-species steady
#' state, because after changing the abundances of the selected species the
#' growth and death rates will have changed.
#'
#' @param params A MizerParams object
#' @param species The species to be selected. Optional. By default all target
#'   species are selected. A vector of species names, or a numeric vector with
#'   the species indices, or a logical vector indicating for each species
#'   whether it is to be selected (TRUE) or not.
#' @param keep A string determining which quantity is to be kept constant. The
#'   choices are "egg" which keeps the egg density constant, "biomass" which
#'   keeps the total biomass of the species constant and "number" which keeps
#'   the total number of individuals constant.
#' @return A MizerParams object in which the initial abundances of the selected
#'   species are changed to their single-species steady state abundances.
#' @export
steadySingleSpecies <- function(params, species = NULL,
                                keep = c("egg", "biomass", "number")) {
    UseMethod("steadySingleSpecies")
}
#' @export
steadySingleSpecies.MizerParams <- function(params, species = NULL,
                                keep = c("egg", "biomass", "number")) {
    species <- valid_species_arg(params, species)
    keep <- match.arg(keep)

    biomass <- getBiomass(params, use_cutoff = TRUE)
    number <- getN(params)

    # Use growth and mortality from current abundances
    # Use growth and mortality from current abundances
    growth_all <- getEGrowth(params)
    mort_all <- getMort(params)

    # Loop through all species and calculate their steady state abundances
    # using the current growth and mortality rates
    
    # We can use get_transport_coefs with an arbitrary dt (e.g. 1) and then transform.
    # We need to construct the rates list
    rates <- list(e_growth = growth_all, mort = mort_all)
    
    # We use a dummy dt = 1
    dt <- 1
    
    coefs <- get_transport_coefs(params, params@initial_n, params@initial_n_pp, 
                                 params@initial_n_other, rates, dt)
    
    # Loop over species
    for (sp in species) {
        w_min_idx <- params@w_min_idx[sp]
        w_max_idx <- sum(params@w <= params@species_params[sp, "w_max"])
        idx <- w_min_idx:w_max_idx
        
        # Check that species can grow to maturity at least
        w_mat_idx <- sum(params@w <= params@species_params[sp, "w_mat"])
        
        # Check growth (existing check)
        growth <- growth_all[sp, ]
        zero_growth_idx <- which(growth[w_min_idx:w_max_idx] == 0)
        if (length(zero_growth_idx) > 0) {
            first_zero_idx <- w_min_idx + zero_growth_idx[1] - 1
            if (first_zero_idx < w_mat_idx) {
                stop(sp, " cannot grow to maturity")
            } else {
                warning(sp, " has zero growth rate after maturity size")
            }
        }
        
        # Extract coefficients for this species
        # We need to perform the same boundary adjustments as in project_n?
        # project_n modifies B and S at w_min_idx.
        # "d_prime[i, j_start] <- d_prime[i, j_start] + r$rdd[i] * dt / params@dw[j_start]"
        # Here S_rhs (right hand side) for steady state is just the recruitment flux term?
        # N_old cancels out with "1" in B.
        # Net equation: Trans(N) + mu N = Source
        # Source at j_start is R_dd / dw[j_start].
        
        # Let's extract A, B, C for the specific species range
        # Note: the coefficients in coefs passed cover all w.
        # We only solve for w_min_idx:w_max_idx.
        # At w_max_idx, we assume N_{w_max_idx+1} = 0 (handled by C=0).
        # At w_min_idx, we have the boundary condition.
        
        a_sp <- coefs$a[sp, ]
        b_sp <- coefs$b[sp, ]
        c_sp <- coefs$c[sp, ]
        
        # Adjust B for steady state (subtract 1 because B includes +1 from time stepping N_new)
        b_sp <- b_sp - 1
        
        # Adjust B at boundary w_min_idx same as in project_n?
        # In project_n:
        # if (j_start > 1) {
        #      correction <- (dt / params@dw[j_start]) * 0.5 * params@diffusion[i, j_start] / params@dw[j_start - 1]
        #      b[i, j_start] <- b[i, j_start] - correction
        #      a[i, j_start] <- 0
        # }
        # This correction removes the "influx from below" due to diffusion, 
        # because we implement influx explicitly via R_dd.
        # We must apply the same correction here.
        
        # Keep egg density constant
        N_fixed <- params@initial_n[sp, w_min_idx]

        # Right hand side vector
        rhs <- numeric(length(params@w))
        rhs[w_min_idx] <- N_fixed
        
        # Modify coefficients for fixed boundary condition
        # Row for w_min_idx: 1 * N_{w_min_idx} + 0 * N_{w_min_idx+1} = N_fixed
        b_sp[w_min_idx] <- 1
        c_sp[w_min_idx] <- 0
        a_sp[w_min_idx] <- 0 # Should be 0 already
        
        # Solve system
        idx_solve <- w_min_idx:w_max_idx
        
        a_sub <- a_sp[idx_solve]
        b_sub <- b_sp[idx_solve]
        c_sub <- c_sp[idx_solve]
        rhs_sub <- rhs[idx_solve]
        
        # Solver
        n_opt <- thomas_solve(a_sub, b_sub, c_sub, rhs_sub)
        
        # Update params
        params@initial_n[sp, ] <- 0
        params@initial_n[sp, idx_solve] <- n_opt
    }

    if (any(is.infinite(params@initial_n))) {
        stop("Candidate steady state holds infinities")
    }
    if (any(is.na(params@initial_n) | is.nan(params@initial_n))) {
        stop("Candidate steady state holds non-numeric values")
    }

    if (keep == "biomass") {
        factor <- biomass / getBiomass(params, use_cutoff = TRUE)
        params@initial_n <- params@initial_n * factor
    }
    if (keep == "number") {
        factor <- number / getN(params)
        params@initial_n <- params@initial_n * factor
    }

    params@time_modified <- lubridate::now()
    params
}
