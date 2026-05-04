#' Set external encounter rate
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
#' If the `ext_encounter` slot has a comment and `reset = FALSE`, then a
#' recalculation from the species parameters is suppressed and a message is
#' issued if the recalculated values would differ from the stored ones.
#'
#' @param params MizerParams
#' @param ext_encounter Optional. An array (species x size) holding the external
#'   encounter rate.  If not supplied, a default is calculated from the `E_ext`
#'   and `n` species parameters as described in the section "Setting external
#'   encounter rate".
#' @param reset `r lifecycle::badge("experimental")`
#'   If set to TRUE, then the external encounter rate will be reset to the value
#'   calculated from the species parameters, even if it was previously
#'   overwritten with a custom value. If set to FALSE (default) then a
#'   recalculation from the species parameters will take place only if no
#'   custom value has been set.
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
setExtEncounter <- function(params, ext_encounter = NULL, reset = FALSE, ...) {
    UseMethod("setExtEncounter")
}
#' @export
setExtEncounter.MizerParams <- function(params, ext_encounter = NULL,
                                        reset = FALSE, ...) {
    assert_that(is.flag(reset))

    if (reset) {
        if (!is.null(ext_encounter)) {
            warning("Because you set `reset = TRUE`, the value you provided ",
                    "for `ext_encounter` will be ignored and a value will be ",
                    "calculated from the species parameters.")
            ext_encounter <- NULL
        }
        comment(params@ext_encounter) <- NULL
    }

    # If ext_encounter array is supplied, check it, store it and return
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
        comment(params@ext_encounter) <- comment(ext_encounter)

        params@time_modified <- lubridate::now()
        return(params)
    }

    # Else recalculate from species params
    params <- set_species_param_default(params, "E_ext", 0)

    ext_encounter <- sweep(outer(params@species_params[["n"]],
                                 params@w, function(x, y) y^x),
                           1, params@species_params[["E_ext"]], "*")

    # Prevent overwriting slot if it has been commented
    if (!is.null(comment(params@ext_encounter))) {
        if (different(ext_encounter, params@ext_encounter)) {
            message("The external encounter rate has been commented and therefore will ",
                    "not be recalculated from the species parameters.")
        }
        return(params)
    }
    params@ext_encounter[] <- ext_encounter

    params@time_modified <- lubridate::now()
    return(params)
}

#' @rdname setExtEncounter
#' @return `getExtEncounter()` or equivalently `ext_encounter()`: A
#'   `ArraySpeciesBySize` object (species x size) with the external encounter rate.
#' @export
getExtEncounter <- function(params) {
    UseMethod("getExtEncounter")
}
#' @export
getExtEncounter.MizerParams <- function(params) {
    ArraySpeciesBySize(params@ext_encounter,
                       value_name = "External encounter rate",
                       units = "g/year",
                       params = params)
}

#' @rdname setExtEncounter
#' @export
ext_encounter <- function(params) {
    getExtEncounter(params)
}

#' @rdname setExtEncounter
#' @param value ext_encounter
#' @export
`ext_encounter<-` <- function(params, value) {
    setExtEncounter(params, ext_encounter = value)
}
