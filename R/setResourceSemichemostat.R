#' Tune the density dependence in reproduction
#'
#' `r lifecycle::badge("experimental")`
#' Takes a MizerParams object `params` with arbitrary density dependence in
#' the resource dynamics and
#' returns a MizerParams object with the resource parameters chosen in such a
#' way that the 
#' Hence if you have tuned your `params` object to describe a
#' particular steady state, then setting the resource parameters
#' with this function will leave you with the exact same steady state.
#'
#' If you do not provide a value for any of the resource parameter
#' arguments, then the resource rate will be kept at its current value. 
#' If you do provide one of the reproduction
#' parameters, this can be either a vector with the same length as 
#' `w_full(params)`, giving a
#' value for every size class, or it can be a single value that is taken as the
#' value at 1g with the values at other sizes calculated according to a 
#' power-law.
#'
#' The values for `resource_capacity` must be larger than the current resource number
#' density. If a smaller value is requested a warning is issued and the
#' value is increased to twice the current resource number density.
#' 
#' The values for the `resource_level` must be positive and less than 1. 
#' Otherwise an error is raised.
#' 
#' The values for `resource_rate` must be large enough to allow the required
#' resource production rate. If a smaller value is requested a warning is issued
#' and the value is increased to the smallest possible value.
#'
#' As can be seen in the graph above, choosing a lower value for
#' `resource_capacity` or a higher value for `resource_rate` means that near the
#' steady state the resource abundance will be less sensitive to a change in the
#' predation mortality and hence less sensitive to changes in the fish
#' abundances.
#'
#' @inheritParams setResource
#' @param resource_level Sets `resource_capacity` to the current resource
#'   number density divided by `resource_level`.
#'
#' @return A MizerParams object with updated resource parameters.
#' @export
setResourceSemichemostat <- function(params, resource_rate, resource_capacity,
                                     resource_level) {
    assert_that(is(params, "MizerParams"))
    
    # check number of arguments
    num_args <- hasArg("resource_rate") + hasArg("resource_capacity") + 
        hasArg("resource_level")
    if (num_args > 1) {
        stop("You should only provide `params` and one other argument.")
    }
    if (num_args == 0) {
        # no values given, so use previous resource_rate
        resource_rate <- resource_rate(params)
    }

    mu <- getResourceMort(params)
    NR <- initialNResource(params)
    n <- params@resource_params[["n"]]
    lambda <- params@resource_params[["lambda"]]
    w_full <- w_full(params)
    no_w_full <- length(w_full)
    
    params@resource_dynamics <- "resource_semichemostat"
    
    if (!missing(resource_level)) {
        assert_that(is.numeric(resource_level))
        if (length(resource_level) != 1 || length(resource_level) != no_w_full) {
            stop("The 'resource_level' should have length 1 or length ",
                 no_w_full, ".")
        }
        # The resource level is allowed to be NaN only where the resource is 0
        if (any(NR > 0 & is.nan(resource_level))) {
            stop("The resource level must be defined everywhere where the current resource is non-vanishing.")
        }
        if (any(NR > 0 &
                (resource_level <= 0 | resource_level >= 1))) {
            stop("The 'resource_level' must always be strictly between 0 and 1.")
        }
        capacity <- NR / resource_level
        capacity[is.nan(resource_level)] <- 0
        rate <- mu / (1 / resource_level - 1)
        rate[is.nan(resource_level)] <- 0
    }
    
    if (!missing(resource_rate)) {
        assert_that(is.numeric(resource_rate))
        if (length(resource_rate) == 1) {
            resource_rate <- resource_rate * w_full ^ (n - 1)
        } else if (length(resource_rate) != no_w_full) {
            stop("The 'resource_rate' should have length 1 or length ",
                 no_w_full, ".")
        }
        if (any(resource_rate <= 0)) {
            stop("The 'resource_rate' must always be strictly positive.")
        }
        capacity <- NR * (resource_rate + mu) / resource_rate
        rate <- resource_rate
    }
    
    if (!missing(resource_capacity)) {
        assert_that(is.numeric(resource_capacity))
        if (length(resource_capacity) == 1) {
            resource_capacity <- resource_capacity * w_full ^ (-lambda)
            resource_capacity[w_full > params@resource_params$w_pp_cutoff] <- 0
        } else if (length(resource_capacity) != no_w_full) {
            stop("The 'resource_rate' should have length 1 or length ",
                 no_w_full, ".")
        }
        if (any(resource_capacity < 0)) {
            stop("The 'resource_capacity' must never be negative.")
        }
        # resource capacity must be above current abundance except where
        # both are zero
        if (any(resource_capacity <= NR & NR > 0)) {
            stop("The 'resource_capacity' must always be greater than current resource number density.")
        }
        capacity <- resource_capacity
        rate <- mu * NR / (resource_capacity - NR)
    }
    
    params@rr_pp <- rate
    params@cc_pp <- capacity
    params@time_modified <- lubridate::now()
    return(params)
}


#' Get resource level
#'
#' `r lifecycle::badge("experimental")` The resource level is the ratio between
#' the current resource number density and the resource carrying capacity.
#' Where both are zero, the resource level is NaN.
#'
#' @param params A MizerParams object
#'
#' @return A vector of the same length as `w_full(params)` with the resource
#'   level in each size class.
#' @export
#' @examples
#' getResourceLevel(NS_params)
#'
#' # The resource level can be changed without changing the steady state:
#' params <- setResourceSemichemostat(NS_params, resource_level = 0.9)
#' getResourceLevel(params)
getResourceLevel <- function(params) {
    assert_that(is(params, "MizerParams"))
    initialNResource(params) / getResourceCapacity(params)
}

