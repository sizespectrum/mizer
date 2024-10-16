#' Validate species parameter data frame
#' 
#' These functions check the validity of a species parameter frame and, where
#' necessary, make corrections. `validGivenSpeciesParams()` only checks and
#' corrects the given species parameters but does not add default values for
#' species parameters that were not provided. `validSpeciesParams()` first calls
#' `validGivenSpeciesParams()` but then goes further by adding default values
#' for species parameters that were not provided.
#' 
#' `validGivenSpeciesParams()` checks the validity of the given species
#' parameter It throws an error if
#' * the `species` column does not exist or contains duplicates
#' * the maximum size is not specified for all species
#' 
#' If a weight-based parameter is missing but the corresponding length-based
#' parameter is given, as well as the `a` and `b` parameters for length-weight
#' conversion, then the weight-based parameters are added. If both length and
#' weight are given, then weight is used and a warning is issued if the two are
#' inconsistent.
#' 
#' If a `w_inf` column is given but no `w_max` then the value from `w_inf` is
#' used. This is for backwards compatibility. But note that the von Bertalanffy
#' parameter `w_inf` is not the maximum size of the largest individual, but the
#' asymptotic size of an average individual.
#' 
#' Some inconsistencies in the size parameters are resolved as follows:
#' * Any `w_mat` that is not smaller than `w_max` is set to `w_max / 4`.
#' * Any `w_mat25` that is not smaller than `w_mat` is set to NA.
#' * Any `w_min` that is not smaller than `w_mat` is set to `0.001` or 
#'   `w_mat /10`, whichever is smaller.
#' * Any `w_repro_max` that is not larger than `w_mat` is set to `4 * w_mat`.
#' 
#' The row names of the returned data frame will be the species names.
#' If `species_params` was provided as a tibble it is converted back to an
#' ordinary data frame.
#' 
#' The function tests for some typical misspellings of parameter names, like
#' wrong capitalisation or missing underscores and issues a warning if it 
#' detects such a name.
#' 
#' `validSpeciesParams()` first calls `validateGivenSpeciesParams()` but then
#' goes further by adding default values for species parameters that were not
#' provided. The function sets default values if any of the following species
#' parameters are missing or NA:
#' * `w_repro_max` is set to `w_max`
#' * `w_mat` is set to `w_max/4`
#' * `w_min` is set to `0.001`
#' * `alpha` is set to `0.6`
#' * `interaction_resource` is set to `1`
#' * `n` is set to `3/4`
#' 
#' Note that the species parameters returned by these functions are not
#' guaranteed to produce a viable model. More checks of the parameters are
#' performed by the individual rate-setting functions (see [setParams()] for the
#' list of these functions).
#' 
#' @param species_params The user-supplied species parameter data frame
#' @return For `validSpeciesParams()`: A valid species parameter data frame with
#'   additional parameters with default values.
#' 
#' @seealso [species_params()], [validGearParams()], [validParams()], [validSim()]
#' @concept helper
#' @export
validSpeciesParams <- function(species_params) {
    sp <- validGivenSpeciesParams(species_params)
    sp <- set_species_param_default(sp, "w_repro_max", sp$w_max)
    sp <- set_species_param_default(sp, "w_mat", sp$w_max / 4)
    sp <- set_species_param_default(sp, "w_min", 0.001)
    sp <- set_species_param_default(sp, "alpha", 0.6)
    sp <- set_species_param_default(sp, "interaction_resource", 1)
    sp <- set_species_param_default(sp, "n", 3/4)
    return(sp)
}

#' @rdname validSpeciesParams
#' @return For `validGivenSpeciesParams()`: A valid species parameter data frame
#'   without additional parameters.
#' @export
validGivenSpeciesParams <- function(species_params) {
    assert_that(is.data.frame(species_params))
    # Convert a tibble back to an ordinary data frame
    sp <- as.data.frame(species_params,
                        stringsAsFactors = FALSE) # for old versions of R
    sp$species <- as.character(sp$species)
    
    # Check for misspellings ----
    misspellings <- c("wmin", "wmax", "wmat", "wmat25", "w_mat_25", "Rmax",
                      "Species", "Gamma", "Beta", "Sigma", "Alpha",
                      "W_min", "W_max", "W_mat", "e_repro", "Age_mat",
                      "w_max_mat")
    query <- intersect(misspellings, names(sp))
    if (length(query) > 0) {
        warning("Some column names in your species parameter data ",
                "frame are very close to standard parameter names: ",
                paste(query, collapse = ", "),
                ". Did you perhaps mis-spell the names?")
    }
    
    # check species ----
    if (!("species" %in% colnames(sp))) {
        stop("The species params dataframe needs a column 'species' with the species names")
    }
    species_names <- as.character(sp$species)
    sp$species <- species_names
    row.names(sp) <- species_names
    no_sp <- nrow(sp)
    if (length(unique(species_names)) != no_sp) {
        stop("The species parameter data frame has multiple rows for the same species")
    }
    
    ## For backwards compatibility, allow r_max instead of R_max
    if (!("R_max" %in% names(sp)) &&
        "r_max" %in% names(sp)) {
        names(sp)[names(sp) == "r_max"] <- "R_max"
    }
    
    # Convert lengths to weights ----
    if (all(c("a", "b") %in% names(sp))) {
        sp <- sp %>%
            set_species_param_from_length("w_mat", "l_mat") %>%
            set_species_param_from_length("w_mat25", "l_mat25") %>%
            set_species_param_from_length("w_repro_max", "l_repro_max") %>%
            set_species_param_from_length("w_max", "l_max") %>%
            set_species_param_from_length("w_min", "l_min")
    }
    
    # check w_max ----
    if (!("w_max" %in% names(sp))) {
        # If old name `w_inf` is used, then copy over to `w_max`
        if ("w_inf" %in% names(sp)) {
            sp$w_max <- sp$w_inf
            warning("The species parameter data frame is missing a `w_max` column. I am copying over the values from the `w_inf` column. But note that `w_max` should be the maximum size of the largest individual, not the asymptotic size of an average indivdidual.")
        } else {
            sp$w_max <- rep(NA, no_sp)
        }
    }
    missing <- is.na(sp$w_max)
    if (any(missing)) {
        stop("You need to specify maximum sizes for all species.")
    }
    if (!is.numeric(sp$w_max)) {
        stop("`w_max` contains non-numeric values.")
    }
    
    # check w_mat ----
    if ("w_mat" %in% names(sp)) {
        wrong <- !is.na(sp$w_mat) & sp$w_mat >= sp$w_max
        if (any(wrong)) {
            warning("For the species ", 
                    paste(sp$species[wrong], collapse = ", "),
                    " the value for `w_mat` is not smaller than that of `w_max`.",
                    " I have corrected that by setting it to 25% of `w_max.")
            sp$w_mat[wrong] <- sp$w_max[wrong] / 4
        }
        
        # check w_mat25 ----
        if ("w_mat25" %in% names(sp)) {
            wrong <- !is.na(sp$w_mat) & !is.na(sp$w_mat25) & sp$w_mat25 >= sp$w_mat
            if (any(wrong)) {
                warning("For the species ", 
                        paste(sp$species[wrong], collapse = ", "),
                        " the value for `w_mat25` is not smaller than that of `w_mat`.",
                        " I have corrected that by setting it to NA.")
                sp$w_mat25[wrong] <- NA
            }
        }
        
        # check w_min ----
        if ("w_min" %in% names(sp)) {
            wrong <- !is.na(sp$w_min) & !is.na(sp$w_mat) & sp$w_min >= sp$w_mat
            if (any(wrong)) {
                sp$w_min[wrong] <- pmin(0.001, sp$w_mat[wrong] / 10)
                warning("For the species ", 
                        paste(sp$species[wrong], collapse = ", "),
                        " the value for `w_min` is not smaller than that of `w_mat`.",
                        " I have reduced the values.")
            }
        }
    }
    
    # check w_repro_max ----
    if ("w_repro_max" %in% names(sp) &&
        "w_mat" %in% names(sp) &&
        "w_mat" %in% names(sp)) {
        wrong <- !is.na(sp$w_repro_max) &
            !is.na(sp$w_mat) &
            sp$w_repro_max <= sp$w_mat
        if (any(wrong)) {
            warning("For the species ", 
                    paste(sp$species[wrong], collapse = ", "),
                    " the value for `w_repro_max` is smaller than that of `w_mat`.",
                    " I have corrected that by setting it to 4 times `w_mat.")
            sp$w_repro_max[wrong] <- 4 * sp$w_mat[wrong]
        }
    }
    sp
}

# Set weight-based parameter from length-based parameter
set_species_param_from_length <- function(sp, pw, pl) {
    if (pl %in% names(sp)) {
        vw <- l2w(sp[[pl]], sp)
        if (any(vw <= 0, na.rm = TRUE)) {
            stop("All lengths should be positive and non-zero.")
        }
        sp <- set_species_param_default(sp, pw, vw)
        # If both weight and length are given, check that they agree to
        # within 10% at least.
        incons <- !is.na(sp[[pw]]) & !is.na(sp[[pl]]) &
            (abs(sp[[pw]] - vw) / pmax(sp[[pw]], vw) > 0.1)
        if (any(incons)) {
            signal(paste0("For the following species I will ignore your value for ",
                          pl, " because it is not consistent with your value for ",
                          pw, ": ", paste(sp$species[incons], collapse = ", ")),
                   class = "info_about_default", var = pl, level = 3)
        }
    }
    sp
}

#' Alias for `validSpeciesParams()`
#' 
#' @description
#' `r lifecycle::badge("deprecated")`
#' 
#' An alias provided for backward compatibility with mizer version <= 2.5.2
#' @inherit validSpeciesParams
#' @export
#' @concept deprecated
completeSpeciesParams <- validSpeciesParams
