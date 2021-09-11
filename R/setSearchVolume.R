#' Set search volume
#' 
#' @section Setting search volume:
#' The search volume \eqn{\gamma_i(w)} of an individual of species \eqn{i}
#' and weight \eqn{w} multiplies the predation kernel when
#' calculating the encounter rate in [getEncounter()] and the 
#' predation rate in [getPredRate()].
#' 
#' The name "search volume" is a bit misleading, because \eqn{\gamma_i(w)} does
#' not have units of volume. It is simply a parameter that determines the rate
#' of predation. Its units depend on your choice, see section "Units in mizer".
#' If you have chosen to work with total abundances, then it is a rate with units
#' 1/year. If you have chosen to work with abundances per m^2 then it has units
#' of m^2/year. If you have chosen to work with abundances per m^3 then it has
#' units of m^3/year.
#' 
#' If the `search_vol` argument is not supplied, then the search volume is 
#' set to
#' \deqn{\gamma_i(w) = \gamma_i w^q_i.} 
#' The values of \eqn{\gamma_i} (the search volume at 1g) and \eqn{q_i} (the
#' allometric exponent of the search volume) are taken from the `gamma` and
#' `q` columns in the species parameter dataframe. If the `gamma`
#' column is not supplied in the species parameter dataframe, a default is
#' calculated by the [get_gamma_default()] function. Note that only
#' for predators of size \eqn{w = 1} gram is the value of the species parameter
#' \eqn{\gamma_i} the same as the value of the search volume \eqn{\gamma_i(w)}.
#' 
#' @param params MizerParams
#' @param search_vol Optional. An array (species x size) holding the search volume
#'   for each species at size. If not supplied, a default is set as described in
#'   the section "Setting search volume". 
#' @param reset `r lifecycle::badge("experimental")`
#'   If set to TRUE, then the search volume will be reset to the
#'   value calculated from the species parameters, even if it was previously
#'   overwritten with a custom value. If set to FALSE (default) then a
#'   recalculation from the species parameters will take place only if no custom
#'   value has been set.
#' @param ... Unused
#' 
#' @return `setSearchVolume()`: A MizerParams object with updated search volume.
#' @export
#' @family functions for setting parameters
setSearchVolume <- function(params, search_vol = NULL, reset = FALSE, ...) {
    assert_that(is(params, "MizerParams"),
                is.flag(reset))
    species_params <- params@species_params
    
    if (reset) {
        if (!is.null(search_vol)) {
            warning("Because you set `reset = TRUE`, the value you provided ", 
                    "for `search_vol` will be ignored and a value will be ",
                    "calculated from the species parameters.")
            search_vol <- NULL
        }
        comment(params@search_vol) <- NULL
    }
    
    # If search_vol array is supplied, check it, store it and return
    if (!is.null(search_vol)) {
        if (is.null(comment(search_vol))) {
            if (is.null(comment(params@search_vol))) {
                comment(search_vol) <- "set manually"
            } else {
                comment(search_vol) <- comment(params@search_vol)
            }
        }
        assert_that(is.array(search_vol))
        assert_that(identical(dim(search_vol), dim(params@search_vol)))
        if (!is.null(dimnames(search_vol)) && 
            !all(dimnames(search_vol)[[1]] == species_params$species)) {
            stop("You need to use the same ordering of species in the ",
                 "search_vol array as in the params object: ", 
                 toString(species_params$species))
        }
        assert_that(all(search_vol >= 0))
        params@search_vol[] <- search_vol
        comment(params@search_vol) <- comment(search_vol)
        
        params@time_modified <- lubridate::now()
        return(params)
    }
    
    # Calculate default for any missing gammas
    params@species_params$gamma <- get_gamma_default(params)
    
    search_vol <- 
        sweep(outer(params@species_params[["q"]], params@w,
                    function(x, y) y ^ x),
              1, params@species_params$gamma, "*")
    
    # Prevent overwriting slot if it has been commented
    if (!is.null(comment(params@search_vol))) {
        # Issue warning but only if a change was actually requested
        if (different(search_vol, params@search_vol)) {
            message("The search volume has been commented and therefore will ",
                    "not be recalculated from the species parameters.")
        }
        return(params)
    }
    params@search_vol[] <- search_vol
    
    params@time_modified <- lubridate::now()
    return(params)
}

#' @rdname setSearchVolume
#' @return `getSearchVolume()` or equivalently `search_vol()`: An array (species
#'   x size) holding the search volume
#' @export
getSearchVolume <- function(params) {
    params@search_vol
}


#' @rdname setSearchVolume
#' @export
search_vol <- function(params) {
    params@search_vol
}

#' @rdname setSearchVolume
#' @param value search_vol
#' @export
`search_vol<-` <- function(params, value) {
    setSearchVolume(params, search_vol = value)
}