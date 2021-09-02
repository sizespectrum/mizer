#' Calibrate model biomass to observed biomass
#' 
#' Given a MizerParams object `params` for which biomass observations are
#' available for at least some species via 
#' `species_params(params)$biomass_observed`, this function returns an
#' updated MizerParams object which is rescaled so that the total biomass
#' of the observed species in the model agrees with the total observed
#' biomass. 
#' 
#' Observed biomasses usually only include individuals above a certain size.
#' This size should either be specified in 
#' `species_params(params)$biomass_cutoff` in grams, or else a cutoff size
#' of `w_mat/20` is assumed.
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


#' Change scale of the model
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' The abundances in mizer and some rates depend on the size of the area to
#' which they refer. So they could be given per square meter or per square
#' kilometer or for an entire study area or any other choice of yours. This
#' function allows you to change the size by automatically changing the
#' abundances and rates accordingly.
#'
#' @details
#' If you rescale the system by a factor \eqn{c} then this function makes the
#' following rescalings in the params object:
#' \itemize{
#' \item The initial abundances `initial_n`, `initial_n_pp` and
#'   `initial_n_other` are rescaled by \eqn{c}.
#' \item The search volume is rescaled by \eqn{1/c}.
#' \item The resource carrying capacity is rescaled by \eqn{c}
#' \item The maximum reproduction rate \eqn{R_{max}}, if used, is rescaled by
#'   \eqn{c}.
#' }
#' The effect of this is that the dynamics of the rescaled system are identical
#' to those of the unscaled system, in the sense that it does not matter whether
#' one first calls [scaleModel()] and then runs a simulation with
#' [project()] or whether one first runs a simulation and then rescales the
#' resulting abundances.
#'
#' Note that if you use non-standard resource dynamics or other components then you
#' may need to rescale additional parameters that appear in those dynamics.
#'
#' @param params A mizer params object
#' @param factor The factor by which the size is rescaled with respect to which
#'   the abundances are given.
#'
#' @return An object of type \linkS4class{MizerParams}
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
    params@search_vol = params@search_vol / factor
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
    initialNOther(params) = initial_n_other
    
    # community
    params@sc <- params@sc * factor
    
    return(params)
}