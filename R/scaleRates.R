#' Rescale all rates in a mizer model
#'
#' Multiplies all rates in the model by a given factor. Rescaling all rates by
#' a factor \eqn{f} is equivalent to rescaling time by \eqn{f}: it speeds up
#' (or slows down) all dynamics without affecting the steady state of each
#' species, provided the resource spectrum is held at its steady-state value.
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
#'   \item Maximum reproduction rate (`R_max` species parameter)
#'   \item Resource growth rate (`rr_pp` slot)
#' }
#'
#' Both the rate arrays stored in the MizerParams slots and the associated
#' species parameters in `species_params` and `given_species_params` are
#' rescaled, so that the parameters remain consistent with the rate arrays.
#'
#' @param params A MizerParams object
#' @param factor The positive factor by which all rates are multiplied.
#' @param ... Currently unused.
#'
#' @return The MizerParams object with all rates rescaled by `factor`.
#' @export
#' @seealso [scaleModel()]
scaleRates <- function(params, factor, ...)
    UseMethod("scaleRates")

#' @export
scaleRates.MizerParams <- function(params, factor, ...) {
    params <- validParams(params)
    assert_that(is.number(factor), factor > 0)

    # Scale rate slots directly
    params@search_vol[] <- params@search_vol * factor
    params@intake_max[] <- params@intake_max * factor
    params@metab[] <- params@metab * factor
    params@mu_b[] <- params@mu_b * factor
    params@ext_encounter[] <- params@ext_encounter * factor
    params@ext_diffusion[] <- params@ext_diffusion * factor
    params@catchability[] <- params@catchability * factor
    params@rr_pp[] <- params@rr_pp * factor

    # Scale species parameters associated with each rate
    sp_cols <- c("gamma", "h", "ks", "k", "z0", "z_ext", "z0pre",
                 "E_ext", "D_ext", "R_max")
    for (col in sp_cols) {
        if (col %in% names(params@species_params)) {
            params@species_params[[col]] <- params@species_params[[col]] * factor
        }
        if (col %in% names(params@given_species_params)) {
            params@given_species_params[[col]] <- params@given_species_params[[col]] * factor
        }
    }

    # Scale catchability in gear_params
    if ("catchability" %in% names(params@gear_params)) {
        params@gear_params$catchability <- params@gear_params$catchability * factor
    }

    params@time_modified <- lubridate::now()
    return(params)
}
