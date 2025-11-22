#' Set external mortality rate
#'
#' @section Setting external mortality rate:
#' The external mortality is all the mortality that is not due to fishing or
#' predation by predators included in the model. The external mortality could be
#' due to predation by predators that are not explicitly included in the model
#' (e.g. mammals or seabirds) or due to other causes like illness. It is a rate
#' with units 1/year.
#'
#' The `ext_mort` argument allows you to specify an external mortality rate
#' that depends on species and body size. You can see an example of this in
#' the Examples section of the help page for [setExtMort()].
#'
#' If the `ext_mort` argument is not supplied, then the external mortality is
#' assumed to depend only on the species, not on the size of the individual:
#' \eqn{\mu_{ext.i}(w) = z_{0.i}}. The value of the constant \eqn{z_0} for each
#' species is taken from the `z0` column of the species parameter data frame, if
#' that column exists. Otherwise it is calculated as
#' \deqn{z_{0.i} = {\tt z0pre}_i\, w_{inf}^{\tt z0exp}.}{z_{0.i} = z0pre_i w_{inf}^{z0exp}.}
#'
#' @param params MizerParams
#' @param ext_mort Optional. An array (species x size) holding the external
#'   mortality rate.  If not supplied, a default is set as described in the
#'   section "Setting external mortality rate".
#' @param reset `r lifecycle::badge("experimental")`
#'   If set to TRUE, then the external mortality rate will be reset
#'   to the value calculated from the `z0` parameters, even if it was
#'   previously overwritten with a custom value. If set to FALSE (default) then
#'   a recalculation from the species parameters will take place only if no
#'   custom value has been set.
#' @param z0pre If `z0`, the mortality from other sources, is not a column
#'   in the species data frame, it is calculated as z0pre * w_max ^ z0exp.
#'   Default value is 0.6.
#' @param z0exp If `z0`, the mortality from other sources, is not a column in
#'   the species data frame, it is calculated as \code{z0pre * w_max ^ z0exp}.
#'   Default value is \code{n-1}.
#' @param z0 `r lifecycle::badge("deprecated")` Use `ext_mort` instead. Not to
#'   be confused with the species_parameter `z0`.
#' @param ... Unused
#'
#' @return `setExtMort()`: A MizerParams object with updated external mortality
#'   rate.
#' @export
#' @family functions for setting parameters
#' @examples
#' params <- newMultispeciesParams(NS_species_params)
#'
#' #### Setting allometric death rate #######################
#'
#' # Set coefficient for each species. Here we choose 0.1 for each species
#' z0pre <- rep(0.1, nrow(species_params(params)))
#'
#' # Multiply by power of size with exponent, here chosen to be -1/4
#' # The outer() function makes it an array species x size
#' allo_mort <- outer(z0pre, w(params)^(-1/4))
#'
#' # Change the external mortality rate in the params object
#' ext_mort(params) <- allo_mort
setExtMort <- function(params, ext_mort = NULL, z0pre = 0.6,
                       z0exp = params@resource_params$n - 1, reset = FALSE, ...) {
    UseMethod("setExtMort")
}
#' @export
setExtMort.MizerParams <- function(params, ext_mort = NULL,
                       z0pre = 0.6, z0exp = params@resource_params$n - 1,
                       reset = FALSE, z0 = deprecated(), ...) {
    if (lifecycle::is_present(z0)) {
        lifecycle::deprecate_warn("2.2.3", "setExtMort(z0)", "setExtMort(ext_mort)")
        ext_mort <- z0
    }
    assert_that(is(params, "MizerParams"),
                is.flag(reset),
                is.number(z0pre), is.number(z0exp))

    if (reset) {
        if (!is.null(ext_mort)) {
            warning("Because you set `reset = TRUE`, the value you provided ",
                    "for `ext_mort` will be ignored and a value will be ",
                    "calculated from the species parameters.")
            ext_mort <- NULL
        }
        comment(params@mu_b) <- NULL
    }

    if (!is.null(ext_mort)) {
        if (is.null(comment(ext_mort))) {
            if (is.null(comment(params@mu_b))) {
                comment(ext_mort) <- "set manually"
            } else {
                comment(ext_mort) <- comment(params@mu_b)
            }
        }
        assert_that(is.array(ext_mort),
                    identical(dim(ext_mort), dim(params@mu_b)))
        params@mu_b[] <- ext_mort
        comment(params@mu_b) <- comment(ext_mort)

        params@time_modified <- lubridate::now()
        return(params)
    }

    assert_that(is.number(z0pre), z0pre >= 0,
                is.number(z0exp))
    species_params <- params@species_params
    assert_that(noNA(species_params$w_max))
    # Sort out z0 (external mortality)
    message <- ("Using z0 = z0pre * w_max ^ z0exp for missing z0 values.")
    params <- set_species_param_default(params, "z0",
                                        z0pre * species_params$w_max^z0exp,
                                        message)
    mu_b <- params@mu_b
    mu_b[] <- params@species_params$z0

    # Prevent overwriting slot if it has been commented
    if (!is.null(comment(params@mu_b))) {
        # Issue warning but only if a change was actually requested
        if (different(mu_b, params@mu_b)) {
            message("The external mortality rate has been commented and therefore ",
                    "will not be recalculated from the species parameters.")
        }
        return(params)
    }
    params@mu_b[] <- mu_b

    params@time_modified <- lubridate::now()
    return(params)
}

#' @rdname setExtMort
#' @return `getExtMort()` or equivalently `ext_mort()`: An array (species x
#'   size) with the external mortality.
#' @export
getExtMort <- function(params) {
    UseMethod("getExtMort")
}
#' @export
getExtMort.MizerParams <- function(params) {
    params@mu_b
}

#' @rdname setExtMort
#' @export
ext_mort <- function(params) {
    UseMethod("ext_mort")
}
#' @export
ext_mort.MizerParams <- function(params) {
    params@mu_b
}

#' @rdname setExtMort
#' @param value ext_mort
#' @export
`ext_mort<-` <- function(params, value) {
    UseMethod("ext_mort<-")
}
#' @export
`ext_mort<-.MizerParams` <- function(params, value) {
    setExtMort(params, ext_mort = value)
}
