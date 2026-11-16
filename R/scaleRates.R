#' Rescale all rates in a mizer model
#'
#' `r lifecycle::badge("experimental")`
#' Multiplies all rates in the model by a given factor. Rescaling all rates by
#' a factor \eqn{f} is equivalent to rescaling time by \eqn{f}: it speeds up
#' (or slows down) all dynamics without affecting the steady state of each
#' species and resource. Extension providers must declare how their dynamics
#' transform under this change of time scale.
#'
#' The following rates and their associated species parameters are rescaled:
#' \itemize{
#'   \item Search volume (`search_vol` slot and `gamma` species parameter)
#'   \item Maximum intake rate (`intake_max` slot and `h` species parameter)
#'   \item Metabolic rate (`metab` slot and `ks`, `k` species parameters)
#'   \item External mortality (`mu_b` slot and `z0`, `z_ext`, `z0pre` species
#'     parameters)
#'   \item External encounter rate (`ext_encounter` slot and `E_ext` species
#'     parameter)
#'   \item External diffusion (`ext_diffusion` slot and `D_ext` species
#'     parameter)
#'   \item Catchability (`catchability` slot and `catchability` column in
#'     `gear_params`)
#'   \item Reproduction rate parameters (`R_max` and
#'     `constant_reproduction` multiply by the factor; `ricker_b` and
#'     `sheperd_b` divide by it)
#'   \item Resource growth rate (`rr_pp` slot and `r_pp` resource parameter)
#' }
#'
#' Both the rate arrays stored in the MizerParams slots and the associated
#' species parameters in `species_params` and `given_species_params` are
#' rescaled, so that the parameters remain consistent with the rate arrays.
#' Changed scalar coefficients are recorded as given values, including those
#' previously calculated by mizer. Existing manually set arrays retain their
#' protection from recalculation.
#'
#' Components and custom providers require declared transformation support.
#' Their hooks run after the core fields are scaled and before final validation,
#' with the original model available for reference. Component states are left
#' to their hooks; there is no automatic state rescaling.
#'
#' @param params A MizerParams object
#' @param factor The positive factor by which all rates are multiplied.
#' @param ... Currently unused.
#'
#' @return The MizerParams object with all rates rescaled by `factor`.
#' @export
#' @seealso [scaleModel()]
scaleRates <- function(params, factor, ...) {
    UseMethod("scaleRates")
}

#' @export
scaleRates.MizerParams <- function(params, factor, ...) {
    params <- validParams(params)
    assert_that(is.number(factor), is.finite(factor), factor > 0)
    if (factor == 1) return(params)
    plan <- prepare_transformation(params, "scale_rates",
                                   details = list(factor = factor))

    # Transform the cached arrays without invoking their rebuilding setters.
    # [] assignment retains dimensions and manually supplied-array comments.
    slots <- c("search_vol", "intake_max", "metab", "mu_b", "ext_encounter",
               "ext_diffusion", "catchability", "rr_pp")
    for (slot in slots) params[[slot]][] <- params[[slot]] * factor
    params <- scale_species_coefficients(params, factor, c(
        gamma = 1, h = 1, ks = 1, k = 1, z0 = 1, z_ext = 1, z0pre = 1,
        E_ext = 1, D_ext = 1
    ))
    params <- scale_repro_coefficients(params, factor)

    if ("catchability" %in% names(params@gear_params)) {
        params@gear_params$catchability <-
            params@gear_params$catchability * factor
    }
    if (!is.null(params@resource_params$r_pp)) {
        params@resource_params$r_pp <- params@resource_params$r_pp * factor
    }

    params <- apply_transformation_hooks(params, plan)
    params <- validParams(params)
    params@time_modified <- lubridate::now()
    params
}

# Scalar coefficients must continue to describe the transformed arrays after a
# later setter rebuild. Recording changed entries also protects a coefficient
# that was originally calculated from an unscaled default or observation.
scale_species_coefficients <- function(params, factor, powers) {
    before <- params@species_params
    for (column in intersect(names(powers), names(before))) {
        params@species_params[[column]] <-
            before[[column]] * factor ^ powers[[column]]
    }
    params@given_species_params <- record_given_species_params(
        params@given_species_params, params@species_params, before
    )
    params
}

# These standard columns retain their units even when another built-in law is
# currently active. A custom provider's hook can override their interpretation.
scale_repro_coefficients <- function(params, factor) {
    scale_species_coefficients(params, factor, c(
        R_max = 1, constant_reproduction = 1, ricker_b = -1, sheperd_b = -1
    ))
}
