#' Set external encounter rate
#'
#' You will usually not need to call this function directly. Instead change
#' the `E_ext` and `n` species parameters with
#' `given_species_params(params) <-` and let mizer recalculate the external
#' encounter rate for you. Call `setExtEncounter()` directly only if you want
#' to impose a different functional form for the size dependence of the
#' external encounter rate. See `vignette("guide-change-parameters")`
#' for a full explanation of when to reach for which level of the model.
#'
#' @section Setting external encounter rate:
#' The external encounter rate is the rate at which a predator encounters
#' food that is not explicitly modelled. It is a rate with units mass/year.
#'
#' The `ext_encounter` argument allows you to specify an external encounter rate
#' that depends on species and body size. You can see an example of this in
#' the Examples section of the help page for [setExtEncounter()].
#'
#' If the `ext_encounter` argument is not supplied, then the external encounter
#' rate is calculated as a power law:
#' \deqn{E_{ext.i}(w) = E_{ext.i}\, w^{n_i}.}{E_{ext.i}(w) = E_{ext.i} * w^n_i.}
#' The coefficient \eqn{E_{ext.i}} is taken from the `E_ext` column of the
#' species parameter data frame, which defaults to 0. The exponent \eqn{n_i} is
#' taken from the `n` column of the species parameter data frame.
#'
#' With `reset = FALSE`, manually overridden species are preserved while the
#' other species are recalculated. A message is issued if recalculation would
#' change protected values. See [rate_overrides()] to inspect overrides.
#'
#' @param params MizerParams
#' @param ext_encounter Optional. An array (species x size) holding the external
#'   encounter rate.  If not supplied, a default is calculated from the `E_ext`
#'   and `n` species parameters as described in the section "Setting external
#'   encounter rate".
#' @param reset
#'   If set to TRUE, then the external encounter rate will be reset to the value
#'   calculated from the species parameters, even if it was previously
#'   overwritten with a custom value. If set to FALSE (default) then a
#'   recalculation takes place for species without a custom value.
#' @param species Optional character vector of species to override or reset.
#'   Omit to override or reset the whole array. Values may be full arrays or
#'   arrays containing only the selected species. For a single species its
#'   values without the species dimension are also accepted: a vector, or a
#'   gear x size matrix for the selectivity. Other species keep their existing
#'   override status.
#'   See [rate_overrides()] to inspect overrides.
#' @param ... Unused
#'
#' @return `setExtEncounter()`: A MizerParams object with updated external encounter
#'   rate.
#' @export
#' @family functions for setting parameters
#' @examples
#' params <- newMultispeciesParams(NS_species_params)
#'
#' #### Setting allometric encounter rate #######################
#'
#' # Set coefficient for each species. Here we choose 0.1 for each species
#' encounter_pre <- rep(0.1, nrow(species_params(params)))
#'
#' # Multiply by power of size with exponent, here chosen to be 3/4
#' # The outer() function makes it an array species x size
#' allo_encounter <- outer(encounter_pre, w(params)^(3/4))
#'
#' # Change the external encounter rate in the params object
#' ext_encounter(params) <- allo_encounter
setExtEncounter <- function(params, ext_encounter = NULL, reset = FALSE, ..., species = NULL) {
    UseMethod("setExtEncounter")
}
#' @export
setExtEncounter.MizerParams <- function(params, ext_encounter = NULL,
                                        reset = FALSE, ..., species = NULL) {
    assert_that(is.flag(reset))

    if (reset) {
        if (!is.null(ext_encounter)) {
            warning("Because you set `reset = TRUE`, the value you provided ",
                    "for `ext_encounter` will be ignored and a value will be ",
                    "calculated from the species parameters.")
            ext_encounter <- NULL
        }
        params <- reset_rate_override(params, "ext_encounter", species)
    }

    # If ext_encounter array is supplied, check it, store it and return
    ext_encounter <- prepare_rate_value(params, "ext_encounter", ext_encounter, species)

    if (!is.null(ext_encounter)) {
        if (is.null(comment(ext_encounter))) {
            if (is.null(comment(params@ext_encounter))) {
                comment(ext_encounter) <- "set manually"
            } else {
                comment(ext_encounter) <- comment(params@ext_encounter)
            }
        }
        assert_that(is.array(ext_encounter),
                    identical(dim(ext_encounter), dim(params@ext_encounter)))
        params@ext_encounter[] <- ext_encounter
        params@ext_encounter <- copy_override_marker(params@ext_encounter, ext_encounter)

        params@time_modified <- lubridate::now()
        return(params)
    }

    # Else recalculate from species params
    params <- set_species_param_default(params, "E_ext", 0)

    ext_encounter <- sweep(outer(params@species_params[["n"]],
                                 params@w, function(x, y) y^x),
                           1, params@species_params[["E_ext"]], "*")

    ext_encounter <- preserve_rate_overrides(params, "ext_encounter", ext_encounter)
    params@ext_encounter[] <- ext_encounter

    params@time_modified <- lubridate::now()
    return(params)
}

#' @rdname setExtEncounter
#' @return `ext_encounter()`: An `ArraySpeciesBySize` object (species x size)
#'   with the external encounter rate.
#' @export
ext_encounter <- function(params) {
    UseMethod("ext_encounter")
}
#' @export
ext_encounter.MizerParams <- function(params) {
    ArraySpeciesBySize(params@ext_encounter,
                       value_name = "External encounter rate",
                       units = "g/year",
                       params = params)
}

#' @rdname setExtEncounter
#' @param value ext_encounter
#' @export
`ext_encounter<-` <- function(params, species = NULL, value) {
    # Preserve direct two-argument calls to the replacement function.
    if (missing(value)) {
        value <- species
        species <- NULL
    }
    setExtEncounter(params, ext_encounter = value, species = species)
}
