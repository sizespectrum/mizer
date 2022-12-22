#' Project resource using semichemostat model
#' 
#' If you set your resource dynamics to use this function then the time
#' evolution of the resource spectrum is described by a semi-chemostat equation
#' \deqn{\frac{\partial N_R(w,t)}{\partial t} = r_R(w) \Big[ c_R (w) - N_R(w,t) \Big] - \mu_R(w, t) N_R(w,t)}{dN_R(w,t)/d t  = r_R(w) ( c_R (w) - N_R(w,t) ) - \mu_R(w,t ) N_R(w,t)}
#' 
#' Here \eqn{r_R(w)} is the resource regeneration rate and \eqn{c_R(w)} is the
#' carrying capacity in the absence of predation. These parameters are changed
#' with [setResource()]. The mortality \eqn{\mu_R(w, t)} is
#' due to predation by consumers and is calculate with [getResourceMort()].
#' 
#' To set your model to use semichemostat dynamics for the resource you do
#' ```
#' resource_dynamics(params) <- "resource_semichemostat"
#' ```
#' where you should replace `params` with the name of the variable holding your
#' MizerParams object.
#' 
#' This function uses the analytic solution of the above equation to calculate
#' the resource abundance at time `t + dt` from all abundances and rates at time
#' `t`, keeping the mortality fixed during the timestep.
#' 
#' @param params A [MizerParams] object
#' @param n A matrix of species abundances (species x size)
#' @param n_pp A vector of the resource abundance by size
#' @param n_other A list with the abundances of other components
#' @param rates A list of rates as returned by [mizerRates()]
#' @param t The current time
#' @param dt Time step
#' @param resource_rate Resource replenishment rate
#' @param resource_capacity Resource carrying capacity
#' @param ... Unused
#'   
#' @return Vector containing resource spectrum at next timestep
#' @export
#' @family resource dynamics
#' @examples
#' params <- NS_params
#' resource_dynamics(params) <- "resource_semichemostat"
resource_semichemostat <- function(params, n, n_pp, n_other, rates, t, dt,
                                   resource_rate, resource_capacity, ...) {
    # We use the exact solution under the assumption of constant mortality 
    # during timestep
    mur <- resource_rate + rates$resource_mort
    n_steady <- resource_rate * resource_capacity / mur
    n_pp_new <- n_steady + (n_pp - n_steady) * exp(-mur * dt)
    
    # Here is an alternative expression that looks as if it might be more
    # precise when the sum of the rates is small due to the use of expm1.
    # However the above has the advantage of preserving the steady state
    # n_steady exactly.
    # n_pp_new <- n_pp * exp(-mur * dt) + n_steady * expm1(-mur * dt)
    
    # if growth rate and death rate are zero then the above would give NaN
    # whereas the value should simply not change
    sel <- !is.finite(n_pp_new)
    n_pp_new[sel] <- n_pp[sel]
    
    n_pp_new
}


#' @rdname resource_semichemostat
#' @details 
#' The [balance_resource_semichemostat()] function is called by [setResource()]
#' to determine the values of the resource parameters that are needed to make
#' the replenishment rate at each size equal the consumption rate at that size,
#' as calculated by [getResourceMort()]. It should be called with only one of
#' `resource_rate` or `resource_capacity` should and will return a named list
#' with the values for both.
#' @export
balance_resource_semichemostat <- function(params,
                                           resource_rate, resource_capacity) {
    assert_that(is(params, "MizerParams"))
    mu <- getResourceMort(params)
    NR <- initialNResource(params)
    
    if (!is.null(resource_rate)) {
        assert_that(is.numeric(resource_rate),
                    length(resource_rate) == length(params@rr_pp))
        if (any(resource_rate[mu != 0] == 0)) {
            stop("I can't balance the resource if the resource rate is zero while the resource mortality is not.")
        }
        cc <- NR * (resource_rate + mu) / resource_rate
        # If there is no death then the capacity must match the abundance
        cc[mu == 0] <- NR[mu == 0]
        # unless the rate is also zero in which case the capacity is not
        # determined, so we'll keep it at the current value.
        cc[mu == 0 & resource_rate == 0] <- 
            params@cc_pp[mu == 0 & resource_rate == 0]
        balance <- list(
            resource_rate = resource_rate,
            resource_capacity = cc
        )
    } else if (!is.null(resource_capacity)) {
        assert_that(is.numeric(resource_capacity),
                    length(resource_capacity) == length(params@cc_pp))
        if (any(resource_capacity < NR)) {
            "I can't balance the resource if the capacity is less than the current abundance."
        }
        # If the capacity equals the abundance then there is no 
        # replenishment, so this is only allowed when there is no death either.
        has_death <- (mu * NR) != 0
        if (any(resource_capacity[has_death] <= NR[has_death])) {
            stop("I can't balance the resource unless the capacity is greater than the current abundance wherever there is consumption.")
        }
        rr <- mu * NR / (resource_capacity - NR)
        # If the capacity equals the abundance then the rate is not determined,
        # so we'll keep it at the current level
        rr[resource_capacity == NR] <- params@rr_pp[resource_capacity == NR]
        balance <- list(
            resource_rate = rr,
            resource_capacity = resource_capacity
        )
    } else {
        stop("Both `resource_capacity` and `resource_rate` were NULL.")
    }
    
    # The following should never happen if the tests above worked as desired.
    if (!all(is.finite(balance$resource_rate)) ||
        !all(is.finite(balance$resource_capacity))) {
        stop("Non-finite results obtained.")
    }
    
    balance
}
