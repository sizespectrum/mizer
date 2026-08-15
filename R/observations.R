# Shared internals for the functions that compare a model with observations.
#
# The biomass and the number variants of the calibration and matching functions
# differ only in which pair of species parameter columns they read and in
# whether the size integral carries a factor of the weight. Everything they
# have in common lives here, so that a change to how an observation is compared
# with the model is made once rather than once per variant.

#' The species parameter columns holding an observation
#'
#' Internal helper for the calibration and matching functions. Observations of
#' type `to` live in the `<to>_observed` column of the species parameters,
#' alongside an optional `<to>_cutoff` column giving the smallest size that the
#' observation includes.
#'
#' @param to The type of observation, either "biomass" or "number".
#' @return A list with entries `to`, `observed` and `cutoff`, the latter two
#'   giving the names of the corresponding species parameter columns.
#' @concept helper
#' @keywords internal
observation_columns <- function(to = c("biomass", "number")) {
    to <- match.arg(to)
    list(to = to,
         observed = paste0(to, "_observed"),
         cutoff = paste0(to, "_cutoff"))
}

#' The minimum weights given by a cutoff species parameter
#'
#' Internal helper for [getBiomass()] and for the calibration and matching
#' functions. Returns the `<to>_cutoff` column of the species parameters, with
#' any NAs replaced by the smallest weight in the model. If the model has no
#' such column at all, the smallest weight is used for every species, so that
#' the whole size range is counted.
#'
#' @param params A MizerParams object.
#' @param to The type of observation, either "biomass" or "number".
#' @return A numeric vector with one minimum weight for each species.
#' @concept helper
#' @keywords internal
cutoff_min_w <- function(params, to = c("biomass", "number")) {
    cols <- observation_columns(to)
    min_w <- params@species_params[[cols$cutoff]]
    if (is.null(min_w)) {
        min_w <- rep(NA_real_, nrow(params@species_params))
    }
    min_w[is.na(min_w)] <- min(params@w)
    min_w
}

#' The modelled counterpart of an observation
#'
#' Internal helper for the calibration and matching functions. Integrates the
#' initial abundance of each species over the sizes that the observation covers
#' and returns the total biomass in grams for `to = "biomass"` or the total
#' number of individuals for `to = "number"`. The `<to>_cutoff` species
#' parameter sets the smallest size counted for each species; where it is
#' missing the whole size range of the species is counted.
#'
#' The integral is done by [sizeIntegral()], so it follows the model's
#' quadrature scheme and lets the bin straddling the cutoff contribute only the
#' part of it that lies above the cutoff.
#'
#' @param params A MizerParams object.
#' @param to The type of observation, either "biomass" or "number".
#' @return A named vector with one value for each species.
#' @concept helper
#' @keywords internal
model_observation <- function(params, to = c("biomass", "number")) {
    cols <- observation_columns(to)
    weight <- if (cols$to == "biomass") params@w else 1
    sizeIntegral(params, weight = weight,
                 min_w = cutoff_min_w(params, cols$to))
}
