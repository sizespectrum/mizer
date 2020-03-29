#' Set maximum intake rate
#'
#' @section Setting maximum intake rate:
#' The maximum intake rate \eqn{h_i(w)} of an individual of species \eqn{i} and
#' weight \eqn{w} determines the feeding level, calculated with
#' [getFeedingLevel()]. It is measured in grams/year.
#'
#' If the `intake_max` argument is not supplied, then the maximum intake
#' rate is set to \deqn{h_i(w) = h_i w^n_i.} 
#' The values of \eqn{h_i} (the maximum intake rate of an individual of size 1
#' gram) and \eqn{n_i} (the allometric exponent for the intake rate) are taken
#' from the `h` and `n` columns in the species parameter dataframe. If
#' the `h` column is not supplied in the species parameter dataframe, it is
#' calculated by the [get_h_default()] function, using `f0` and
#' the `k_vb` column, if they are supplied.
#' 
#' If \eqn{h_i} is set to `Inf`, fish will consume all encountered food.
#'
#' @param params MizerParams
#' @param intake_max Optional. An array (species x size) holding the maximum
#'   intake rate for each species at size. If not supplied, a default is set as
#'   described in the section "Setting maximum intake rate".
#' @param ... Unused
#' 
#' @return A `MizerParams` object with updated maximum intake rate. Because
#'   of the way the R language works, `setMaxIntakeRate()` does not make the
#'   changes to the params object that you pass to it but instead returns a new
#'   params object. So to affect the change you call the function in the form
#'   `params <- setMaxIntakeRate(params, ...)`.
#' @export
#' @family functions for setting parameters
#' @examples
#' \dontrun{
#' params <- NS_params
#' params@species_params$h[3] <- 35
#' params <- setMaxIntakeRate(params)
#' }
setMaxIntakeRate <- function(params, 
                             intake_max = NULL, ...) {
    assert_that(is(params, "MizerParams"))
    species_params <- params@species_params
    
    # If intake_max array is supplied, check it, store it and return
    if (!is.null(intake_max)) {
        assert_that(is.array(intake_max),
                    identical(dim(intake_max), dim(params@intake_max)))
        if (!is.null(dimnames(intake_max)) && 
            !all(dimnames(intake_max)[[1]] == species_params$species)) {
            stop("You need to use the same ordering of species in the ",
                 "intake_max array as in the params object: ", 
                 toString(species_params$species))
        }
        assert_that(all(intake_max >= 0))
        params@intake_max[] <- intake_max
        comment(params@intake_max) <- comment(intake_max)
        return(params)
    }
    
    params@species_params$h <- get_h_default(params)
    
    intake_max <- sweep(outer(params@species_params[["n"]], 
                              params@w, function(x, y) y^x),
                        1, params@species_params[["h"]], "*") 
    
    # Prevent overwriting slot if it has been commented
    if (!is.null(comment(params@intake_max))) {
        # Issue warning but only if a change was actually requested
        if (any(intake_max != params@intake_max)) {
            message("The max intake rate has been commented and therefore will ",
                    "not be recalculated from the species parameters.")
        }
        return(params)
    }
    params@intake_max[] <- intake_max
    return(params)
}

#' @rdname setMaxIntakeRate
#' @export
getMaxIntakeRate <- function(params) {
    params@intake_max
}
