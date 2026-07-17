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
#' parameters. It throws an error if
#' * the `species` column does not exist or contains duplicates
#' * the asymptotic size `w_inf` is not specified for all species (but see the
#'   backwards-compatibility note below)
#' 
#' If a weight-based parameter is missing but the corresponding length-based
#' parameter is given, as well as the `a` and `b` parameters for length-weight
#' conversion, then the weight-based parameters are added. If both length and
#' weight are given, then weight is used and an `info_about_default` condition
#' is signalled if the two are inconsistent.
#' 
#' The required maximum-size parameter is `w_inf`, the von Bertalanffy
#' asymptotic size of an average individual. For backwards compatibility, if no
#' `w_inf` column is given, its values are taken from the `w_repro_max` column
#' if that is present, or otherwise from the `w_max` column, and an
#' informational message is issued. (`w_repro_max` is preferred over `w_max`
#' because in earlier versions of mizer it was the size at which growth stopped
#' and is therefore the closest analogue to the asymptotic size.)
#'
#' Some inconsistencies in the size parameters are resolved as follows:
#' * Any `w_mat` that is not smaller than `w_inf` is set to `w_inf / 4`.
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
#' `validSpeciesParams()` first calls `validGivenSpeciesParams()` but then
#' goes further by adding default values for species parameters that were not
#' provided. It only sets defaults for those species parameters that are not
#' owned by a single rate-setting function, namely those that are read by
#' several of them (`n`), that are used only when projecting (`alpha`), that
#' determine the size grid (`w_min`, `w_max`), that are needed for the
#' length-weight conversion (`a`, `b`) or that are used only for reporting
#' (`is_background`). The function sets default values if any of the following
#' species parameters are missing or NA:
#' * `w_max` is set to `1.5 * w_inf` (it is only a computational boundary)
#' * `w_repro_max` is set to `w_inf`
#' * `w_mat` is set to `w_inf/4`
#' * `w_min` is set to `0.001`
#' * `alpha` is set to `0.6`
#' * `n` is set to `3/4`
#' * `a` is set to `0.01`
#' * `b` is set to `3`
#' * `is_background` is set to `FALSE`
#'
#' All other species parameters are given their default values by the
#' rate-setting function that uses them, so that each default has a single
#' home. For example `p` and `k` are set by [setMetabolicRate()], `z_ext`, `d`
#' and `z0` by [setExtMort()], `E_ext` by [setExtEncounter()], `D_ext` by
#' [setExtDiffusion()], `interaction_resource` by [setInteraction()], `beta`
#' and `sigma` by [setPredKernel()], `q` and `gamma` by [setSearchVolume()],
#' and `erepro`, `m`, `w_mat25` and `R_max` by [setReproduction()]. These
#' columns are therefore absent from the data frame returned by
#' `validSpeciesParams()` but present in the species parameters of a
#' `MizerParams` object, because [setParams()] calls all the rate-setting
#' functions.
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
    species_params(species_params, strict = TRUE)
}

#' @rdname validSpeciesParams
#' @return For `validGivenSpeciesParams()`: A valid species parameter data frame
#'   without additional parameters.
#' @export
validGivenSpeciesParams <- function(species_params) {
    given_species_params(species_params, strict = TRUE)
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
