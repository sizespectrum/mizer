#' Set metabolic rate
#' 
#' Sets the rate at which energy is used for metabolism and activity
#' 
#' @section Setting metabolic rate:
#' The metabolic rate is subtracted from the energy income rate to calculate
#' the rate at which energy is available for growth and reproduction, see
#' [getEReproAndGrowth()]. It is measured in grams/year.
#' 
#' If the `metab` argument is not supplied, then for each species the
#' metabolic rate \eqn{k(w)} for an individual of size \eqn{w} is set to
#' \deqn{k(w) = ks w^p + k w,}
#' where \eqn{ks w^p} represents the rate of standard metabolism and \eqn{k w}
#' is the rate at which energy is expended on activity and movement. The values
#' of \eqn{ks}, \eqn{p} and \eqn{k} are taken from the `ks`, `p` and
#' `k` columns in the species parameter dataframe. If any of these
#' parameters are not supplied, the defaults are \eqn{k = 0}, \eqn{p = n} and
#' \deqn{ks = f_c h \alpha w_{mat}^{n-p},}{ks = f_c * h * alpha * w_mat^(n - p),}
#' where \eqn{f_c} is the critical feeding level taken from the `fc` column
#' in the species parameter data frame. If the critical feeding level is not
#' specified, a default of \eqn{f_c = 0.2} is used.
#' 
#' @param params MizerParams
#' @param metab Optional. An array (species x size) holding the metabolic rate
#'   for each species at size. If not supplied, a default is set as described in
#'   the section "Setting metabolic rate".
#' @param p The allometric metabolic exponent. This is only used if `metab`
#'   is not given explicitly and if the exponent is not specified in a `p`
#'   column in the `species_params`.
#' @param ... Unused
#' 
#' @return MizerParams object with updated metabolic rate. Because of the way
#'   the R language works, `setMetabolicRate()` does not make the changes to the
#'   params object that you pass to it but instead returns a new params object.
#'   So to affect the change you call the function in the form
#'   `params <- setMetabolicRate(params, ...)`.
#' @export
#' @family functions for setting parameters
#' @examples
#' \dontrun{
#' params <- NS_params
#' # Change activity coefficient for species 3
#' params@species_params$k[3] <- 8
#' params <- setMetabolicRate(params)
#' }
setMetabolicRate <- function(params, 
                             metab = NULL, p = NULL, ...) {
    assert_that(is(params, "MizerParams"))
    if (!is.null(p)) {
        assert_that(is.numeric(p))
        params <- set_species_param_default(params, "p", p)
    }
    species_params <- params@species_params
    if (!is.null(metab)) {
        assert_that(is.array(metab),
                    identical(dim(metab), dim(params@metab)))
        if (!is.null(dimnames(metab)) && 
            !all(dimnames(metab)[[1]] == species_params$species)) {
            stop("You need to use the same ordering of species as in the ",
                 "params object: ", toString(species_params$species))
        }
        assert_that(all(metab >= 0))
        params@metab[] <- metab
        comment(params@metab) <- comment(metab)
        return(params)
    }
    
    params <- set_species_param_default(params, "k", 0)
    params@species_params$ks <- get_ks_default(params)
    metab <- 
        sweep(outer(params@species_params$p, params@w,
                    function(x, y) y ^ x),
              1, params@species_params$ks, "*") +
        outer(params@species_params$k, params@w)
    
    # Prevent overwriting slot if it has been commented
    if (!is.null(comment(params@metab))) {
        # Issue warning but only if a change was actually requested
        if (any(metab != params@metab)) {
            message("The metabolic rate has been commented and therefore will ",
                    "not be recalculated from the species parameters.")
        }
        return(params)
    }
    params@metab[] <- metab
    return(params)
}

#' @rdname setMetabolicRate
#' @export
getMetabolicRate <- function(params) {
    params@metab
}
