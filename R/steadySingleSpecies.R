#' Set initial abundances to solution of steady-state equation with current rates
#'
#' `r lifecycle::badge("experimental")`
#' This first calculates growth and death rates that arise from the current
#' initial abundances. Then it solves the steady-state equation with these
#' growth and death rates and the current abundance at the smallest size.
#' It sets the initial abundances of the selected species to this solution.
#'
#' The function only changes the initial abundances. It does not adjust the
#' reproduction parameters or any other parameters. Therefore the result of
#' applying this function is of course not a steady state, because after
#' changing the abundances of the selected species the growth, death and
#' reproduction rates will have changed.
#'
#' If the `keep` argument is supplied, the solution for the selected species
#' are rescaled to keep the specified quantity at the value they had before
#' calling this function.
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

    # Loop over species to make checks
    N0_vec <- numeric(nrow(params@species_params))
    names(N0_vec) <- params@species_params$species
    for (sp in species) {
        w_min_idx <- params@w_min_idx[sp]
        w_max_idx <- sum(params@w <= params@species_params[sp, "w_max"])

        # Check that species can grow to maturity at least
        w_mat_idx <- sum(params@w <= params@species_params[sp, "w_mat"])

        # Check growth (existing check)
        growth <- growth_all[sp, ]
        zero_growth_idx <- which(growth[w_min_idx:w_max_idx] == 0)
        if (length(zero_growth_idx) > 0) {
            first_zero_idx <- w_min_idx + zero_growth_idx[1] - 1
            if (first_zero_idx < w_mat_idx) {
                stop(sp, " cannot grow to maturity")
            }
        }

        N0_vec[sp] <- params@initial_n[sp, w_min_idx]
        params@initial_n[sp, ] <- 0
    }

    # Calculate steady state for all species at once
    n_exact_matrix <- get_steady_state_n(params, growth_all, mort_all, N0_vec)

    # Update initial_n for selected species
    for (sp in species) {
        w_min_idx <- params@w_min_idx[sp]

        if (w_min_idx == length(params@w)) {
             params@initial_n[sp, w_min_idx] <- N0_vec[sp]
        } else {
             params@initial_n[sp, w_min_idx:length(params@w)] <-
                 n_exact_matrix[sp, w_min_idx:length(params@w)]
        }
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
