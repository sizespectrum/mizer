#' Match biomasses to observations
#' 
#' `r lifecycle::badge("experimental")`
#' The function adjusts the abundances of the species in the model so that their
#' biomasses match with observations.
#' 
#' The function works by multiplying for each species the abundance density
#' at all sizes by the same factor. This will of course not give a steady
#' state solution, even if the initial abundance densities were at steady state.
#' So after using this function you may want to use `steady()` to run the model 
#' to steady state, after which of course the biomasses will no longer match
#' exactly. You could then iterate this process. This is described in the
#' blog post at \url{https://blog.mizer.sizespectrum.org/posts/2021-08-20-a-5-step-recipe-for-tuning-the-model-steady-state/}.
#' 
#' Before you can use this function you will need to have added a
#' `biomass_observed` column to your model which gives the observed biomass in
#' grams.  For species for which you have no observed biomass, you should set
#' the value in the `biomass_observed` column to 0 or NA.
#'
#' Biomass observations usually only include individuals above a certain size.
#' This size should be specified in a `biomass_cutoff` column of the species
#' parameter data frame. If this is missing, it is assumed that all sizes are
#' included in the observed biomass, i.e., it includes larval biomass.
#' 
#' @param params A MizerParams object
#' @param species The species to be affected. Optional. By default all observed
#'   biomasses will be matched. A vector of species names, or a numeric vector
#'   with the species indices, or a logical vector indicating for each species
#'   whether it is to be affected (TRUE) or not.
#' @param info_level Controls the amount of information messages that are shown.
#'   Higher levels lead to more messages, `info_level = 0` gives silence. The
#'   default is taken from the `mizer_info_level` option, see
#'   [default_info_level()].
#' @param ... Additional arguments passed to the method.
#' @return A MizerParams object
#' @export
#' @examples 
#' params <- NS_params
#' species_params(params)$biomass_observed <- 
#'     c(0.8, 61, 12, 35, 1.6, 20, 10, 7.6, 135, 60, 30, 78)
#' species_params(params)$biomass_cutoff <- 10
#' params <- calibrateBiomass(params)
#' params <- matchBiomasses(params)
#' plotBiomassObservedVsModel(params)
matchBiomasses <- function(params, species = NULL,
                           info_level = default_info_level(), ...) {
    UseMethod("matchBiomasses")
}

#' @export
matchBiomasses.MizerParams <- function(params, species = NULL,
                                       info_level = default_info_level(), ...) {
    with_info_level(info_level = info_level,
                    match_to(params, species = species, to = "biomass",
                             fname = "matchBiomasses"))
}

#' Match numbers to observations
#'
#' `r lifecycle::badge("experimental")`
#' The function adjusts the numbers of the species in the model so that their
#' numbers match with observations.
#'
#' The function works by multiplying for each species the number density
#' at all sizes by the same factor. This will of course not give a steady
#' state solution, even if the initial number densities were at steady state.
#' So after using this function you may want to use `steady()` to run the model
#' to steady state, after which of course the numbers will no longer match
#' exactly. You could then iterate this process. This is described in the
#' blog post at \url{https://blog.mizer.sizespectrum.org/posts/2021-08-20-a-5-step-recipe-for-tuning-the-model-steady-state/}.
#'
#' Before you can use this function you will need to have added a
#' `number_observed` column to your model which gives the observed number of
#' individuals. For species for which you have no observed number, you should set
#' the value in the `number_observed` column to 0 or NA.
#'
#' Number observations usually only include individuals above a certain size.
#' This size should be specified in a `number_cutoff` column of the species
#' parameter data frame. If this is missing, it is assumed that all sizes are
#' included in the observed number, i.e., it includes larval number.
#'
#' @param params A MizerParams object
#' @param species The species to be affected. Optional. By default all observed
#'   numbers will be matched. A vector of species names, or a numeric vector
#'   with the species indices, or a logical vector indicating for each species
#'   whether it is to be affected (TRUE) or not.
#' @param info_level Controls the amount of information messages that are shown.
#'   Higher levels lead to more messages, `info_level = 0` gives silence. The
#'   default is taken from the `mizer_info_level` option, see
#'   [default_info_level()].
#' @param ... Additional arguments passed to the method.
#' @return A MizerParams object
#' @export
#' @examples
#' params <- NS_params
#' species_params(params)$number_observed <-
#'     c(0.8, 61, 12, 35, 1.6, 20, 10, 7.6, 135, 60, 30, 78)
#' species_params(params)$number_cutoff <- 10
#' params <- calibrateNumber(params)
#' params <- matchNumbers(params)
matchNumbers <- function(params, species = NULL,
                         info_level = default_info_level(), ...) {
    UseMethod("matchNumbers")
}

#' @export
matchNumbers.MizerParams <- function(params, species = NULL,
                                     info_level = default_info_level(), ...) {
    with_info_level(info_level = info_level,
                    match_to(params, species = species, to = "number",
                             fname = "matchNumbers"))
}

#' Match a quantity to observations species by species
#'
#' Internal implementation shared by [matchBiomasses()] and [matchNumbers()].
#' Multiplies the abundance density of each selected species at all sizes by
#' the factor that brings the modelled quantity onto the observation. Species
#' that were not selected, or that have no positive observation, are left
#' alone.
#'
#' @param params A MizerParams object.
#' @param species The species to be affected, in any of the forms accepted by
#'   [valid_species_arg()].
#' @param to The type of observation, either "biomass" or "number".
#' @param fname The name of the calling function, used when reporting that the
#'   model has been moved off its steady state.
#' @return A MizerParams object.
#' @concept helper
#' @keywords internal
match_to <- function(params, species = NULL, to = c("biomass", "number"),
                     fname) {
    cols <- observation_columns(to)
    observed <- params@species_params[[cols$observed]]
    if (is.null(observed)) {
        return(params)
    }
    selected <- valid_species_arg(params, species = species,
                                  return.logical = TRUE) &
        !is.na(observed) & observed > 0
    if (!any(selected)) {
        return(params)
    }

    model <- model_observation(params, cols$to)[selected]
    unreachable <- is.na(model) | model <= 0
    if (any(unreachable)) {
        cutoff <- cutoff_min_w(params, cols$to)[selected][unreachable]
        error_species <- params@species_params$species[selected][unreachable]
        stop(paste(
            paste(error_species, "does not grow up to the",
                  cols$cutoff, "size of", cutoff, "grams."),
            collapse = "\n"
        ))
    }
    factors <- observed[selected] / model
    params@initial_n[selected, ] <-
        params@initial_n[selected, , drop = FALSE] * factors

    signal_off_steady(fname)
    setBevertonHolt(params)
}
