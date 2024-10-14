#' Calibrate the model scale to match total observed biomass
#' 
#' `r lifecycle::badge("experimental")`
#' Given a MizerParams object `params` for which biomass observations are
#' available for at least some species via the `biomass_observed` column in the
#' species_params data frame, this function returns an updated MizerParams
#' object which is rescaled with [scaleModel()] so that the total biomass in
#' the model agrees with the total observed biomass.
#' 
#' Biomass observations usually only include individuals above a certain size.
#' This size should be specified in a biomass_cutoff column of the species
#' parameter data frame. If this is missing, it is assumed that all sizes are
#' included in the observed biomass, i.e., it includes larval biomass.
#' 
#' After using this function the total biomass in the model will match the
#' total biomass, summed over all species. However the biomasses of the
#' individual species will not match observations yet, with some species
#' having biomasses that are too high and others too low. So after this
#' function you may want to use [matchBiomasses()]. This is described in the
#' blog post at https://bit.ly/2YqXESV.
#' 
#' If you have observations of the yearly yield instead of biomasses, you can
#' use [calibrateYield()] instead of this function.
#' 
#' @param params A MizerParams object
#' @return A MizerParams object
#' @export
#' @examples 
#' params <- NS_params
#' species_params(params)$biomass_observed <- 
#'     c(0.8, 61, 12, 35, 1.6, 20, 10, 7.6, 135, 60, 30, 78)
#' species_params(params)$biomass_cutoff <- 10
#' params2 <- calibrateBiomass(params)
#' plotBiomassObservedVsModel(params2)
calibrateBiomass <- function(params) {
    if ((!("biomass_observed" %in% names(params@species_params))) ||
        all(is.na(params@species_params$biomass_observed))) {
        return(params)
    }
    no_sp <- nrow(params@species_params)
    cutoff <- params@species_params$biomass_cutoff
    # When no cutoff known, set it to 0
    if (is.null(cutoff)) cutoff <- rep(0, no_sp)
    cutoff[is.na(cutoff)] <- 0
    observed <- params@species_params$biomass_observed
    observed_total <- sum(observed, na.rm = TRUE)
    sp_observed <- which(!is.na(observed))
    model_total <- 0
    for (sp_idx in sp_observed) {
        model_total <- 
            model_total + 
            sum((params@initial_n[sp_idx, ] * params@w * params@dw)
                [params@w >= cutoff[[sp_idx]]])
    }
    scaleModel(params, factor = observed_total / model_total)
}

# The following is a copy of the code for `calibrateBiomass()` just with
# the text replacements "Biomass" -> "Number" and "biomass" to "number" and
# the removal of the `params@w` factor in the calculations.

#' Calibrate the model scale to match total observed number
#'
#' `r lifecycle::badge("experimental")`
#' Given a MizerParams object `params` for which number observations are
#' available for at least some species via the `number_observed` column in the
#' species_params data frame, this function returns an updated MizerParams
#' object which is rescaled with [scaleModel()] so that the total number in
#' the model agrees with the total observed number.
#'
#' Number observations usually only include individuals above a certain size.
#' This size should be specified in a number_cutoff column of the species
#' parameter data frame. If this is missing, it is assumed that all sizes are
#' included in the observed number, i.e., it includes larval number.
#'
#' After using this function the total number in the model will match the
#' total number, summed over all species. However the numbers of the
#' individual species will not match observations yet, with some species
#' having numbers that are too high and others too low. So after this
#' function you may want to use [matchNumbers()]. This is described in the
#' blog post at https://bit.ly/2YqXESV.
#'
#' If you have observations of the yearly yield instead of numbers, you can
#' use [calibrateYield()] instead of this function.
#'
#' @param params A MizerParams object
#' @return A MizerParams object
#' @export
#' @examples
#' params <- NS_params
#' species_params(params)$number_observed <-
#'     c(0.8, 61, 12, 35, 1.6, 20, 10, 7.6, 135, 60, 30, 78)
#' species_params(params)$number_cutoff <- 10
#' params2 <- calibrateNumber(params)
calibrateNumber <- function(params) {
    if ((!("number_observed" %in% names(params@species_params))) ||
        all(is.na(params@species_params$number_observed))) {
        return(params)
    }
    no_sp <- nrow(params@species_params)
    cutoff <- params@species_params$number_cutoff
    # When no cutoff known, set it to 0
    if (is.null(cutoff)) cutoff <- rep(0, no_sp)
    cutoff[is.na(cutoff)] <- 0
    observed <- params@species_params$number_observed
    observed_total <- sum(observed, na.rm = TRUE)
    sp_observed <- which(!is.na(observed))
    model_total <- 0
    for (sp_idx in sp_observed) {
        model_total <-
            model_total +
            sum((params@initial_n[sp_idx, ] * params@dw)
                [params@w >= cutoff[[sp_idx]]])
    }
    scaleModel(params, factor = observed_total / model_total)
}

#' Calibrate the model scale to match total observed yield
#' 
#' `r lifecycle::badge("deprecated")`
#' 
#' This function has been deprecated and will be removed in the future unless
#' you have a use case for it. If you do have a use case for it, please let the
#' developers know by creating an issue at
#' <https://github.com/sizespectrum/mizer/issues>.
#' 
#' Given a MizerParams object `params` for which yield observations are
#' available for at least some species via the `yield_observed` column in the
#' species_params data frame, this function returns an updated MizerParams
#' object which is rescaled with [scaleModel()] so that the total yield in
#' the model agrees with the total observed yield.
#' 
#' After using this function the total yield in the model will match the
#' total observed yield, summed over all species. However the yields of the
#' individual species will not match observations yet, with some species
#' having yields that are too high and others too low. So after this
#' function you may want to use [matchYields()].
#' 
#' If you have observations of species biomasses instead of yields, you can
#' use [calibrateBiomass()] instead of this function.
#' 
#' @param params A MizerParams object
#' @return A MizerParams object
#' @concept deprecated
#' @export
#' @examples 
#' params <- NS_params
#' species_params(params)$yield_observed <-
#'     c(0.8, 61, 12, 35, 1.6, 20, 10, 7.6, 135, 60, 30, 78)
#' gear_params(params)$catchability <-
#'     c(1.3, 0.065, 0.31, 0.18, 0.98, 0.24, 0.37, 0.46, 0.18, 0.30, 0.27, 0.39)
#' params2 <- calibrateYield(params)
#' plotYieldObservedVsModel(params2)
calibrateYield <- function(params) {
    lifecycle::deprecate_warn(
        "2.6.0", "calibrateYield()",
        details = "This function has not proven useful. If you do have a use case for it, please let the developers know by creating an issue at https://github.com/sizespectrum/mizer/issues"
    )
    if ((!("yield_observed" %in% names(params@species_params))) ||
        all(is.na(params@species_params$yield_observed))) {
        return(params)
    }
    observed <- params@species_params$yield_observed
    observed_total <- sum(observed, na.rm = TRUE)
    sp_observed <- which(!is.na(observed))
    biomass <- sweep(params@initial_n, 2, params@w * params@dw, "*")
    yield_model <- rowSums(biomass * getFMort(params))[sp_observed]
    model_total <- sum(yield_model)
    scaleModel(params, factor = observed_total / model_total)
}


#' Change scale of the model
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' The abundances in mizer and some rates depend on the size of the area to
#' which they refer. So they could be given per square meter or per square
#' kilometer or for an entire study area or any other choice of yours. This
#' function allows you to change the scale of the model by automatically
#' changing the abundances and rates accordingly.
#'
#' @details
#' If you rescale the model by a factor \eqn{c} then this function makes the
#' following rescalings in the params object:
#' \itemize{
#' \item The initial abundances are rescaled by \eqn{c}.
#' \item The search volume is rescaled by \eqn{1/c}.
#' \item The resource carrying capacity is rescaled by \eqn{c}
#' \item The maximum reproduction rate \eqn{R_{max}} is rescaled by
#'   \eqn{c}.
#' }
#' The effect of this is that the dynamics of the rescaled model are identical
#' to those of the unscaled model, in the sense that it does not matter whether
#' one first calls [scaleModel()] and then runs a simulation with
#' [project()] or whether one first runs a simulation and then rescales the
#' resulting abundances.
#'
#' Note that if you use non-standard resource dynamics or other components then
#' you may need to rescale additional parameters that appear in those dynamics.
#' 
#' In practice you will need to use some observations to set the scale for your
#' model. If you have biomass observations you can use [calibrateBiomass()],
#' if you have yearly yields you can use [calibrateYield()].
#'
#' @param params A MizerParams object
#' @param factor The factor by which the scale is multiplied
#'
#' @return The rescaled MizerParams object
#' @export
scaleModel <- function(params, factor) {
    params <- validParams(params)
    assert_that(is.number(factor),
                factor > 0)
    
    # Resource replenishment rate
    params@cc_pp <- params@cc_pp * factor
    params@resource_params$kappa <- params@resource_params$kappa * factor
    
    # Rmax
    # r_max is a deprecated spelling of R_max. Get rid of it.
    if ("r_max" %in% names(params@species_params)) {
        params@species_params$R_max <- params@species_params$r_max
        params@species_params$r_max <- NULL
        message("The 'r_max' column has been renamed to 'R_max'.")
    }
    if ("R_max" %in% names(params@species_params)) {
        params@species_params$R_max <- params@species_params$R_max * factor
    }
    
    # Search volume
    params@search_vol <- params@search_vol / factor
    if ("gamma" %in% names(params@species_params)) {
        params@species_params$gamma <- params@species_params$gamma / factor
    }
    
    # Initial values
    initial_n_other <- params@initial_n_other
    for (res in names(initial_n_other)) {
        initial_n_other[[res]] <- initial_n_other[[res]] * factor
    }
    initialN(params) <- params@initial_n * factor
    initialNResource(params) <- params@initial_n_pp * factor
    initialNOther(params) <- initial_n_other
    
    # community
    params@sc <- params@sc * factor
    
    return(params)
}
