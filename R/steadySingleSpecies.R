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
#' calling this function. The biomass is counted above the `biomass_cutoff`
#' species parameter, if there is one. A species that had none of the
#' quantity is not rescaled, because there is nothing to keep: it is absent,
#' or it lies entirely below its biomass cutoff. If its solution does have
#' some of the quantity, a warning says so. The function stops with an error
#' if the solution for a species has none of a quantity that the species did
#' have.
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
#' @param info_level Controls the amount of information messages and warnings
#'   that are shown. Higher levels lead to more messages, `info_level = 0`
#'   gives silence. The default is taken from the `mizer_info_level` option,
#'   see [default_info_level()].
#' @return A MizerParams object in which the initial abundances of the selected
#'   species are changed to their single-species steady state abundances.
#' @export
#' @examples
#' # Set initial abundance of Cod to its single-species steady state
#' params <- steadySingleSpecies(NS_params, species = "Cod")
steadySingleSpecies <- function(params, species = NULL,
                                keep = c("egg", "biomass", "number"),
                                info_level = default_info_level()) {
    UseMethod("steadySingleSpecies")
}
#' @export
steadySingleSpecies.MizerParams <- function(params, species = NULL,
                                keep = c("egg", "biomass", "number"),
                                info_level = default_info_level()) {
    with_info_level(info_level = info_level, {
    # A species listed twice would on its second pass take its egg abundance
    # from the row that its first pass has already emptied.
    species <- unique(valid_species_arg(params, species))
    keep <- match.arg(keep)

    # The quantity that `keep` holds fixed, measured before `initial_n` changes
    kept_quantity <- function(params) {
        switch(keep,
               biomass = getBiomass(params, use_cutoff = TRUE),
               number = getN(params))[species]
    }
    if (keep != "egg") {
        old <- kept_quantity(params)
    }

    # Use rates from current abundances.
    growth_all <- getEGrowth(params)
    mort_all <- getMort(params)
    diffusion_all <- getDiffusion(params)

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
        # The solution is grown up from the egg size, so an empty egg size
        # class would give an all-zero spectrum. That is the right answer only
        # for a species that is absent altogether.
        if (N0_vec[sp] == 0 && any(params@initial_n[sp, ] > 0)) {
            stop("The abundance of ", sp, " is zero in the size class ",
                 "containing its egg size w_min but non-zero at larger ",
                 "sizes. Because the steady state is determined by growing ",
                 "the egg abundance up through the size classes, this would ",
                 "empty the species. Give ", sp, " a positive abundance in ",
                 "size class ", w_min_idx, ", the class containing w_min = ",
                 params@species_params[sp, "w_min"], ".")
        }
        params@initial_n[sp, ] <- 0
    }

    # Calculate steady state for all species at once
    n_exact_matrix <- get_steady_state_n(params, growth_all, mort_all,
                                         diffusion_all, N0_vec)

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

    if (keep != "egg") {
        new <- kept_quantity(params)
        # `which()` drops the NA of a non-numeric candidate, which the checks
        # below report
        lost <- which(old > 0 & new == 0)
        if (length(lost) > 0) {
            stop("The steady state has zero ", keep, " for ",
                 toString(species[lost]), ", so it cannot be rescaled to ",
                 "keep the ", keep, " it had before.")
        }
        # A species that had none of the quantity is not rescaled, because
        # there is nothing to keep. A factor of 0 would empty a species that
        # lies entirely below its biomass cutoff, and 0 / 0 would fill an
        # absent species with NaN. The steady state can still reach above the
        # cutoff, because it extends a bin beyond `w_max`.
        gained <- which(old == 0 & new > 0)
        if (length(gained) > 0) {
            signal_info("keep",
                        paste0("There was no ", keep, " of ",
                               toString(species[gained]), " to keep, so its ",
                               "steady state has not been rescaled and now ",
                               "has some ", keep, "."),
                        level = 1, severity = "warning", unhandled = "show")
        }
        factor <- ifelse(old == 0, 1, old / new)
        params@initial_n[species, ] <- params@initial_n[species, ] * factor
    }

    # Checked after the rescaling, which can itself overflow
    if (any(is.infinite(params@initial_n))) {
        stop("Candidate steady state holds infinities")
    }
    if (any(is.na(params@initial_n) | is.nan(params@initial_n))) {
        stop("Candidate steady state holds non-numeric values")
    }

    params@time_modified <- lubridate::now()
    params
    })
}
