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
    with_info_level(info_level = info_level, {
    if (!("biomass_observed" %in% names(params@species_params))) {
        return(params)
    }
    species_sel <- valid_species_arg(params, species = species, 
                                 return.logical = TRUE) &
        !is.na(params@species_params$biomass_observed) &
        params@species_params$biomass_observed > 0
    if (!any(species_sel)) {
        return(params)
    }
    
    model_biomass <- getBiomass(params, use_cutoff = TRUE)
    observed_biomass <- params@species_params$biomass_observed
    
    # Only consider selected species
    selected_idx <- which(species_sel)
    zero_biomass <- model_biomass[selected_idx] <= 0 | is.na(model_biomass[selected_idx])
    if (any(zero_biomass)) {
        cutoff <- params@species_params$biomass_cutoff[selected_idx][zero_biomass]
        error_species <- params@species_params$species[selected_idx][zero_biomass]
        stop(paste(
            paste(error_species, "does not grow up to the biomass_cutoff size of",
                  cutoff, "grams."),
            collapse = "\n"
        ))
    }
    factors <- observed_biomass[selected_idx] / model_biomass[selected_idx]
    params@initial_n[selected_idx, ] <- params@initial_n[selected_idx, , drop = FALSE] * factors

    signal_off_steady("matchBiomasses")
    setBevertonHolt(params)
    })
}

# The following is a copy of the code for `matchBiomasses()` just with
# the text replacements "Biomass" -> "Number" and "biomass" to "number" and
# the removal of the `params@w` factor in the calculations.

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
    with_info_level(info_level = info_level, {
    if (!("number_observed" %in% names(params@species_params))) {
        return(params)
    }
    species <- valid_species_arg(params, species = species,
                                 return.logical = TRUE) &
        !is.na(params@species_params$number_observed) &
        params@species_params$number_observed > 0
    if (length(species) == 0) {
        return(params)
    }
    
    error_message <- ""
    for (sp in seq_len(nrow(params@species_params))[species]) {
        cutoff <- params@species_params$number_cutoff[[sp]]
        if (is.null(cutoff) || is.na(cutoff)) {
            cutoff <- 0
        }
        total <- sum((params@initial_n[sp, ] * params@dw)
                     [params@w >= cutoff])
        if (!(total > 0)) {
            error_message <- paste(
                error_message, params@species_params$species[[sp]],
                "does not grow up to the number_cutoff size of",
                cutoff, "grams.\n")
        }
        factor <- params@species_params$number_observed[[sp]] / total
        params@initial_n[sp, ] <- params@initial_n[sp, ] * factor
    }
    if (error_message != "") {
        stop(error_message)
    }

    signal_off_steady("matchNumbers")
    setBevertonHolt(params)
    })
}
