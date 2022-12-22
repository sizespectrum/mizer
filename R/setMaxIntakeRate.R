#' Set maximum intake rate
#'
#' @section Setting maximum intake rate:
#' The maximum intake rate \eqn{h_i(w)} of an individual of species \eqn{i} and
#' weight \eqn{w} determines the feeding level, calculated with
#' [getFeedingLevel()]. It is measured in grams/year.
#'
#' If the `intake_max` argument is not supplied, then the maximum intake
#' rate is set to \deqn{h_i(w) = h_i w^{n_i}.} 
#' The values of \eqn{h_i} (the maximum intake rate of an individual of size 1
#' gram) and \eqn{n_i} (the allometric exponent for the intake rate) are taken
#' from the `h` and `n` columns in the species parameter dataframe. If
#' the `h` column is not supplied in the species parameter dataframe, it is
#' calculated by the [get_h_default()] function.
#' 
#' If \eqn{h_i} is set to `Inf`, fish of species i will consume all encountered
#' food.
#'
#' @param params MizerParams
#' @param intake_max Optional. An array (species x size) holding the maximum
#'   intake rate for each species at size. If not supplied, a default is set as
#'   described in the section "Setting maximum intake rate".
#' @param reset `r lifecycle::badge("experimental")`
#'   If set to TRUE, then the intake rate will be reset to the value
#'   calculated from the species parameters, even if it was previously 
#'   overwritten with a custom value. If set to FALSE (default) then a
#'   recalculation from the species parameters will take place only if no
#'   custom value has been set.
#' @param ... Unused
#' 
#' @return `setReproduction()`: A MizerParams object with updated maximum
#'   intake rate.
#' @export
#' @family functions for setting parameters
setMaxIntakeRate <- function(params, intake_max = NULL, reset = FALSE, ...) {
    assert_that(is(params, "MizerParams"),
                is.flag(reset))
    species_params <- params@species_params
    
    if (reset) {
        if (!is.null(intake_max)) {
            warning("Because you set `reset = TRUE`, the value you provided ", 
                    "for `intake_max` will be ignored and a value will be ",
                    "calculated from the species parameters.")
            intake_max <- NULL
        }
        comment(params@intake_max) <- NULL
    }
    
    # If intake_max array is supplied, check it, store it and return
    if (!is.null(intake_max)) {
        if (is.null(comment(intake_max))) {
            if (is.null(comment(params@intake_max))) {
                comment(intake_max) <- "set manually"
            } else {
                comment(intake_max) <- comment(params@intake_max)
            }
        }
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
        
        params@time_modified <- lubridate::now()
        return(params)
    }
    
    # Else recalculate from species params
    params@species_params[["h"]] <- get_h_default(params)
    
    intake_max <- sweep(outer(params@species_params[["n"]], 
                              params@w, function(x, y) y^x),
                        1, params@species_params[["h"]], "*") 
    
    # Prevent overwriting slot if it has been commented
    if (!is.null(comment(params@intake_max))) {
        # Issue warning but only if a change was actually requested
        if (different(intake_max, params@intake_max)) {
            message("The max intake rate has been commented and therefore will ",
                    "not be recalculated from the species parameters.")
        }
        return(params)
    }
    params@intake_max[] <- intake_max
    
    params@time_modified <- lubridate::now()
    return(params)
}

#' @rdname setMaxIntakeRate
#' @return `getMaxIntakeRate()` or equivalently `intake_max()`: An array
#'   (species x size) with the maximum intake rate.
#' @export
getMaxIntakeRate <- function(params) {
    params@intake_max
}


#' @rdname setMaxIntakeRate
#' @export
intake_max <- function(params) {
    params@intake_max
}

#' @rdname setMaxIntakeRate
#' @param value intake_max
#' @export
`intake_max<-` <- function(params, value) {
    setMaxIntakeRate(params, intake_max = value)
}
