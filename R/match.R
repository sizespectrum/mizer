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
#' blog post at https://bit.ly/2YqXESV.
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
matchBiomasses <- function(params, species = NULL) {
    if (!("biomass_observed" %in% names(params@species_params))) {
        return(params)
    }
    species <- valid_species_arg(params, species = species, 
                                 return.logical = TRUE) &
        !is.na(params@species_params$biomass_observed) &
        params@species_params$biomass_observed > 0
    if (length(species) == 0) {
        return(params)
    }
    for (sp in seq_len(nrow(params@species_params))[species]) {
        cutoff <- params@species_params$biomass_cutoff[[sp]]
        if (is.null(cutoff) || is.na(cutoff)) {
            cutoff <- 0
        }
        total <- sum((params@initial_n[sp, ] * params@w * params@dw)
                     [params@w >= cutoff])
        factor <- params@species_params$biomass_observed[[sp]] / total
        params@initial_n[sp, ] <- params@initial_n[sp, ] * factor
    }
    
    setBevertonHolt(params)
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
#' blog post at https://bit.ly/2YqXESV.
#'
#' Before you can use this function you will need to have added a
#' `number_observed` column to your model which gives the observed number in
#' grams.  For species for which you have no observed number, you should set
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
#' @return A MizerParams object
#' @export
#' @examples
#' params <- NS_params
#' species_params(params)$number_observed <-
#'     c(0.8, 61, 12, 35, 1.6, 20, 10, 7.6, 135, 60, 30, 78)
#' species_params(params)$number_cutoff <- 10
#' params <- calibrateNumber(params)
#' params <- matchNumbers(params)
matchNumbers <- function(params, species = NULL) {
    if (!("number_observed" %in% names(params@species_params))) {
        return(params)
    }
    species <- valid_species_arg(params, species = species,
                                 return.logical = TRUE) &
        !is.na(params@species_params$number_observed) &
        params@species_params$number_observed > 0
    for (sp in seq_len(nrow(params@species_params))[species]) {
        cutoff <- params@species_params$number_cutoff[[sp]]
        if (is.null(cutoff) || is.na(cutoff)) {
            cutoff <- 0
        }
        total <- sum((params@initial_n[sp, ] * params@dw)
                     [params@w >= cutoff])
        factor <- params@species_params$number_observed[[sp]] / total
        params@initial_n[sp, ] <- params@initial_n[sp, ] * factor
    }
    
    setBevertonHolt(params)
}

#' Match yields to observations
#' 
#' `r lifecycle::badge("deprecated")`
#' This function has been deprecated and will be removed in the future unless
#' you have a usecase for it. If you do have a use case for it, please let the
#' developers know by creating an issue at
#' <https://github.com/sizespectrum/mizer/issues>.
#' 
#' If you want to match the yields to observations, you should use the
#' matchYield() function from the mizerExperimental package instead, which
#' adjusts the catchability to match the yield rather than by adjusting the
#' biomass.
#' 
#' The function adjusts the abundances of the species in the model so that their
#' yearly yields under the given fishing mortalities match with observations.
#' 
#' The function works by multiplying for each species the abundance density
#' at all sizes by the same factor. This will of course not give a steady
#' state solution, even if the initial abundance densities were at steady state.
#' So after using this function you may want to use `steady()` to run the model 
#' to steady state, after which of course the yields will no longer match
#' exactly. You could then iterate this process. This is described in the
#' blog post at https://bit.ly/2YqXESV.
#' 
#' Before you can use this function you will need to have added a
#' `yield_observed` column to your model which gives the observed yields in
#' grams per year.  For species for which you have no observed biomass, you
#' should set the value in the `yield_observed` column to 0 or NA.
#' 
#' @param params A MizerParams object
#' @param species The species to be affected. Optional. By default all observed
#'   yields will be matched. A vector of species names, or a numeric vector
#'   with the species indices, or a logical vector indicating for each species
#'   whether it is to be affected (TRUE) or not.
#' @return A MizerParams object
#' @concept deprecated
#' @export
#' @examples 
#' params <- NS_params
#' species_params(params)$yield_observed <- 
#'     c(0.8, 61, 12, 35, 1.6, 20, 10, 7.6, 135, 60, 30, 78)
#' gear_params(params)$catchability <-
#'     c(1.3, 0.065, 0.31, 0.18, 0.98, 0.24, 0.37, 0.46, 0.18, 0.30, 0.27, 0.39)
#' params <- calibrateYield(params)
#' params <- matchYields(params)
#' plotYieldObservedVsModel(params)
matchYields <- function(params, species = NULL) {
    lifecycle::deprecate_warn(
        "2.6.0", "matchYields()", "mizerExperimental::matchYield()",
        details = "This function has not proven useful. If you do have a use case for it, please let the developers know by creating an issue at https://github.com/sizespectrum/mizer/issues"
    )
    if (!("yield_observed" %in% names(params@species_params))) {
        return(params)
    }
    biomass <- sweep(params@initial_n, 2, params@w * params@dw, "*")
    yield_model <- rowSums(biomass * getFMort(params))
    
    all_species <- params@species_params$species
    include <- valid_species_arg(params, species = species, 
                                 return.logical = TRUE)
    
    # ignore species with no observations
    no_obs <- include &
        (is.na(params@species_params$yield_observed) |
        params@species_params$yield_observed <= 0)
    if (any(no_obs)) {
        message("The following species have no yield observations and ",
                "their abundances will not be changed: ",
                paste0(all_species[no_obs], collapse = ", "), ".")
        include <- include & !no_obs
    }
    
    # ignore species with no model yield
    no_yield <- include & yield_model == 0
    if (any(no_yield)) {
        message("The following species are not being fished in your model ",
                "and their abundances will not be changed: ",
                paste0(all_species[no_yield], collapse = ", "), ".")
        include <- include & !no_yield
    }
    
    factors <- params@species_params$yield_observed[include] /
        yield_model[include]
    params@initial_n[include, ] <- 
        sweep(params@initial_n[include, , drop = FALSE], 1, factors, "*")
    
    setBevertonHolt(params)
}
