#' Set external mortality rate
#' 
#' @section Setting external mortality rate:
#' The external mortality is all the mortality that is not due to fishing or
#' predation by predators included in the model. The external mortality could be
#' due to predation by predators that are not explicitly included in the model
#' (e.g. mammals or seabirds) or due to other causes like illness. It is a rate
#' with units 1/year.
#' 
#' The `z0` argument allows you to specify an external mortality rate
#' that depends on species and body size. You can see an example of this in
#' the Examples section of the help page for [setExtMort()].
#' 
#' If the `z0` argument is not supplied, then the external mortality
#' is assumed to depend only on the species, not on the
#' size of the individual: \eqn{\mu_{b.i}(w) = z_{0.i}}. The value of the
#' constant \eqn{z_0} for each species is taken from the `z0` column of the
#' species_params data frame, if that column exists. Otherwise it is calculated
#' as 
#' \deqn{z_{0.i} = {\tt z0pre}_i\, w_{inf}^{\tt z0exp}.}{z_{0.i} = z0pre_i w_{inf}^{z0exp}.}
#' 
#' @param params MizerParams
#' @param z0 Optional. An array (species x size) holding the external
#'   mortality rate.
#' @param z0pre If `z0`, the mortality from other sources, is not a column
#'   in the species data frame, it is calculated as z0pre * w_inf ^ z0exp.
#'   Default value is 0.6.
#' @param z0exp If `z0`, the mortality from other sources, is not a column
#'   in the species data frame, it is calculated as \code{z0pre * w_inf ^ z0exp}.
#'   Default value is \code{n-1}.
#' @param ... Unused
#' 
#' @return MizerParams object with updated external mortality rate. Because of
#'   the way the R language works, `setExtMort()` does not make the changes
#'   to the params object that you pass to it but instead returns a new params
#'   object. So to affect the change you call the function in the form
#'   `params <- setExtMort(params, ...)`.
#' @export
#' @family functions for setting parameters
#' @examples
#' \dontrun{
#' params <- newMultispeciesParams(NS_species_params)
#'
#' #### Setting allometric death rate #######################
#' 
#' # Set coefficient for each species. Here we choose 0.1 for each species
#' z0pre <- rep(0.1, nrow(params@species_params))
#' 
#' # Multiply by power of size with exponent, here chosen to be -1/4
#' # The outer() function makes it an array species x size
#' z0 <- outer(z0pre, params@w^(-1/4))
#' 
#' # Change the external mortality rate in the params object
#' params <- setExtMort(params, z0 = z0)
#' }
setExtMort <- function(params, z0 = NULL, z0pre = 0.6, z0exp = -1/4, ...) {
    assert_that(is(params, "MizerParams"))
    if (!is.null(z0)) {
        assert_that(is.array(z0),
                    identical(dim(z0), dim(params@mu_b)))
        params@mu_b[] <- z0
        comment(params@mu_b) <- comment(z0)
        return(params)
    }
    
    assert_that(is.number(z0pre), z0pre >= 0,
                is.number(z0exp))
    species_params <- params@species_params
    assert_that(noNA(species_params$w_inf))
    # Sort out z0 (external mortality)
    message <- ("Note: Using z0 = z0pre * w_inf ^ z0exp for missing z0 values.")
    params <- set_species_param_default(params, "z0",
                                        z0pre * species_params$w_inf^z0exp,
                                        message)
    mu_b <- params@species_params$z0
    
    # Prevent overwriting slot if it has been commented
    if (!is.null(comment(params@mu_b))) {
        # Issue warning but only if a change was actually requested
        if (any(mu_b != params@mu_b)) {
            message("The external mortality rate has been commented and therefore ",
                    "will not be recalculated from the species parameters.")
        }
        return(params)
    }
    params@mu_b[] <- mu_b
    return(params)
}

#' @rdname setExtMort
#' @export
getExtMort <- function(params) {
    params@mu_b
}
