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
#' If you have chose to work with total abundances, then it is a rate with units
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
#' @param ... Unused
#' 
#' @return MizerParams with updated search volume. Because of the way the R
#'   language works, `setSearchVolume()` does not make the changes to the params
#'   object that you pass to it but instead returns a new params object. So to
#'   affect the change you call the function in the form
#'   `params <- setSearchVolume(params, ...)`.
#' @export
#' @family functions for setting parameters
#' @examples
#' \dontrun{
#' params <- newTraitParams()
#' params@species_params$gamma[3] <- 1000
#' params <- setSearchVolume(params)
#' }
setSearchVolume <- function(params, 
                            search_vol = NULL, ...) {
    assert_that(is(params, "MizerParams"))
    species_params <- params@species_params
    # If search_vol array is supplied, check it, store it and return
    if (!is.null(search_vol)) {
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
        if (any(search_vol != params@search_vol)) {
            message("The search volume has been commented and therefore will ",
                    "not be recalculated from the species parameters.")
        }
        return(params)
    }
    params@search_vol[] <- search_vol
    return(params)
}

#' @rdname setSearchVolume
#' @export
getSearchVolume <- function(params) {
    params@search_vol
}
