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
#' @param reset `r lifecycle::badge("experimental")`
#'   If set to TRUE, then the metabolic rate will be reset to the
#'   value calculated from the species parameters, even if it was previously
#'   overwritten with a custom value. If set to FALSE (default) then a
#'   recalculation from the species parameters will take place only if no
#'   custom value has been set.
#' @param ... Unused
#' 
#' @return `setMetabolicRate()`: A MizerParams object with updated metabolic rate.
#' @export
#' @family functions for setting parameters
setMetabolicRate <- function(params, metab = NULL, p = NULL, 
                             reset = FALSE, ...) {
    assert_that(is(params, "MizerParams"),
                is.flag(reset))
    if (!is.null(p)) {
        assert_that(is.numeric(p))
        params <- set_species_param_default(params, "p", p)
    }
    species_params <- params@species_params
    
    if (reset) {
        if (!is.null(metab)) {
            warning("Because you set `reset = TRUE`, the value you provided ", 
                    "for `metab` will be ignored and a value will be ",
                    "calculated from the species parameters.")
            metab <- NULL
        }
        comment(params@metab) <- NULL
    }
    
    if (!is.null(metab)) {
        if (is.null(comment(metab))) {
            if (is.null(comment(params@metab))) {
                comment(metab) <- "set manually"
            } else {
                comment(metab) <- comment(params@metab)
            }
        }
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
        
        params@time_modified <- lubridate::now()
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
        if (different(metab, params@metab)) {
            message("The metabolic rate has been commented and therefore will ",
                    "not be recalculated from the species parameters.")
        }
        return(params)
    }
    params@metab[] <- metab
    
    params@time_modified <- lubridate::now()
    return(params)
}

#' @rdname setMetabolicRate
#' @return `getMetabolicRate()` or equivalently `metab()`: An array
#'   (species x size) with the metabolic rate.
#' @export
getMetabolicRate <- function(params) {
    params@metab
}

#' @rdname setMetabolicRate
#' @export
metab <- function(params) {
    params@metab
}

#' @rdname setMetabolicRate
#' @param value metab
#' @export
`metab<-` <- function(params, value) {
    setMetabolicRate(params, metab = value)
}