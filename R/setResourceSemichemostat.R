#' Set the parameters for the resource dynamics without changing the steady state
#'
#' @description
#' `r lifecycle::badge("experimental")`
#' 
#' Sets the resource capacity and the resource rate in such a way that at the
#' current abundances the resource replenishes as quickly as it is consumed.
#' Hence if you have tuned your `params` object to describe a particular steady
#' state, then setting the resource parameters with this function will leave you
#' with the exact same steady state.
#' 
#' @details
#' This function sets the resource dynamics to use `resource_semichemostat()`.
#' This means that the resource rate \eqn{r_R} and the resource capacity 
#' \eqn{c_R} enter the equation for the resource abundance density \eqn{N_R} as
#' \deqn{\frac{\partial N_R(w,t)}{\partial t} = r_R(w) \Big[ c_R (w) - N_R(w,t) \Big] - \mu_R(w, t) N_R(w,t)}{dN_R(w,t)/d t  = r_R(w) ( c_R (w) - N_R(w,t) ) - \mu_R(w,t ) N_R(w,t)}
#' where the mortality \eqn{\mu_R(w, t)} is
#' due to predation by consumers and is calculate with [getResourceMort()].
#' 
#' Because of the requirement that the resource is replenished at the same rate
#' at which it is consumed, there is only one free parameter left to be chosen.
#' Therefore you should provide at most one of the three resource arguments
#' `resource_rate`, `resource_capcacity` or `resource_level`. 
#' 
#' The argument be either a vector with the same length as `w_full(params)`,
#' giving a value for every size class, or it can be a single value that is
#' taken as the value at 1g. This is then completed to a vector of values for
#' all sizes as follows:
#' 
#' * The `resource_rate` argument is completed to a power law with exponent 
#'   `n-1` where `n` is taken from the [resource_params()] list.
#' * The `resource_capacity` argument is completed to a power low with exponent
#'   `-lambda`, with `lamba` taken from the [resource_params()] list. The 
#'   resulting values must be larger than the current resource number density, 
#'   otherwise an error is raised.
#' * The `resource_level` argument is completed to the same value at all sizes. 
#'   It must be a number strictly between 0 and 1.
#' 
#' You can change the exponents `n` and `lambda` using [resource_params()].
#' 
#' You can however set up different resource dynamics with
#' [resource_dynamics<-()].
#'
#' @inheritParams setResource
#' @param resource_rate Optional. Vector of resource intrinsic birth rates.
#'   Must be strictly positive.
#' @param resource_capacity Optional. Vector of resource intrinsic carrying
#'   capacities.
#' @param resource_level Optional. Vector with the ratio between the current
#'   resource number density and the resource capacity. Must be strictly
#'   between 0 and 1, except at sizes where the resource is zero, where it can 
#'   be `NaN`.
#'
#' @return A MizerParams object with updated resource rate and resource
#'   capacity.
#' @concept resource dynamics
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
        resource_rate <- getResourceRate(params)
    }

    mu <- getResourceMort(params)
    NR <- initialNResource(params)
    n <- params@resource_params[["n"]]
    lambda <- params@resource_params[["lambda"]]
    w_full <- w_full(params)
    no_w_full <- length(w_full)
    
    if (!missing(resource_level)) {
        assert_that(is.numeric(resource_level))
        if (length(resource_level) != 1 && length(resource_level) != no_w_full) {
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
    
    params@resource_dynamics <- "resource_semichemostat"
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
#' @concept resource dynamics
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
