#' Set resource dynamics
#' 
#' Sets the intrinsic resource birth rate and the intrinsic resource carrying
#' capacity as well as the name of the function used to simulate the resource
#' dynamics. By default, the birth rate and the carrying capacity are changed
#' together in such a way that the resource replenishes at the same rate at
#' which it is consumed. So you should only provide either the 
#' `resource_rate` or the `resource_capacity` (or `resource_level`) because
#' the other is determined by the requirement that the resource replenishes
#' at the same rate at which it is consumed.
#' 
#' @section Setting resource dynamics:
#' You would usually set the resource dynamics only after having finished the 
#' calibration of the steady state. Then setting the resource dynamics with
#' this function will preserve that steady state, unless you explicitly 
#' choose to set `balance = FALSE`. Your choice of the resource dynamics only
#' affects the dynamics around the steady state. The higher the resource rate
#' or the lower the resource capacity the less sensitive the model will be to
#' changes in the competition for resource.
#' 
#' The `resource_dynamics` argument allows you to choose the resource dynamics
#' function. By default, mizer uses a semichemostat model to describe the
#' resource dynamics in each size class independently. This semichemostat
#' dynamics is implemented by the function [resource_semichemostat()]. You can
#' change that to use a logistic model implemented by [resource_logistic()] or
#' you can use [resource_constant()] which keeps the resource constant or you
#' can write your own function.
#' 
#' Both the [resource_semichemostat()] and the [resource_logistic()] dynamics
#' are parametrised in terms of a size-dependent rate \eqn{r_R(w)} and a 
#' size-dependent capacity \eqn{c_R}. The help pages of these functions give
#' the details.
#' 
#' The `resource_rate` argument can be a vector (with the same length as
#' `w_full(params)`) specifying the intrinsic resource birth rate for each size
#' class. Alternatively it can be a single number, which is then used as the
#' coefficient in a power law: then the intrinsic birth rate \eqn{r_R(w)} at
#' size \eqn{w} is set to
#' \deqn{r_R(w) = r_R w^{n-1}.}
#' The power-law exponent \eqn{n} is taken from the `n` argument.
#' 
#' The `resource_capacity` argument can be a vector specifying the intrinsic
#' resource carrying capacity for each size class. Alternatively it can be a
#' single number, which is then used as the coefficient in a truncated power
#' law: then the intrinsic carrying capacity \eqn{c_R(w)} at size \eqn{w}
#' is set to
#' \deqn{c(w) = \kappa\, w^{-\lambda}}{c(w) = \kappa w^{-\lambda}}
#' for all \eqn{w} less than `w_pp_cutoff` and zero for larger sizes.
#' The power-law exponent \eqn{\lambda} is taken from the `lambda` argument.
#'
#' The values for `kappa`, `lambda`, `n` and `w_pp_cutoff` are stored in a list
#' in the `resource_params` slot of the MizerParams object so that they can be
#' re-used automatically in the future. That list can be accessed with
#' [resource_params()].
#' 
#' @param params A MizerParams object
#' @param resource_rate Optional. A vector of per-capita resource birth
#'   rate for each size class or a single number giving the coefficient in the
#'   power-law for this rate, see Details. Must be strictly positive.
#' @param resource_capacity Optional. Vector of resource intrinsic carrying
#'   capacities or coefficient in the power-law for the capacity, see Details.
#'   The resource capacity must be larger than the resource abundance.
#' @param resource_level Optional. The ratio between the current resource number
#'   density and the resource capacity. Either a number used at all sizes or a
#'   vector specifying a value for each size. Must be strictly between 0 and 1,
#'   except at sizes where the resource is zero, where it can be `NaN`. This
#'   determines the resource capacity, so do not specify both this and
#'   `resource_capacity`.
#' @param resource_dynamics Optional. Name of the function that determines the
#'   resource dynamics by calculating the resource spectrum at the next time
#'   step from the current state.
#' @param balance By default, if possible, the resource parameters are 
#'   set so that the resource replenishes at the same rate at which it is 
#'   consumed. In this case you should only specify either the resource rate
#'   or the resource capacity (or resource level) because the other is then
#'   determined automatically. Set to FALSE if you do not want the balancing.
#' @param n Used to set power-law exponent for resource rate if the
#'   `resource_rate` argument is given as a single number.
#' @param lambda Used to set power-law exponent for resource capacity if the
#'   `resource_capacity` argument is given as a single number.
#' @param w_pp_cutoff The upper cut off size of the resource spectrum power law
#'   used only if `resource_capacity` is given as a single number.
#' @param r_pp `r lifecycle::badge("deprecated")`. Use `resource_rate` argument
#'   instead.
#' @param kappa `r lifecycle::badge("deprecated")`. Use `resource_capacity`
#'   argument instead.
#' @param ... Unused
#' 
#' @return `setResource`: A MizerParams object with updated resource parameters
#' @family resource parameters
#' @export
setResource <- function(params,
                        resource_rate = NULL,
                        resource_capacity = NULL,
                        resource_level = NULL,
                        resource_dynamics = NULL,
                        balance = NULL,
                        lambda = resource_params(params)[["lambda"]],
                        n = resource_params(params)[["n"]],
                        w_pp_cutoff = resource_params(params)[["w_pp_cutoff"]],
                        r_pp = deprecated(),
                        kappa = deprecated(),
                        ...) {
    
    if (lifecycle::is_present(r_pp)) {
        lifecycle::deprecate_warn("1.0.0", "setParams(r_pp)", 
                                  "setParams(resource_rate)")
        resource_rate <- r_pp
    }
    if (lifecycle::is_present(kappa)) {
        lifecycle::deprecate_warn("1.0.0", "setParams(kappa)", 
                                  "setParams(resource_capacity)")
        resource_capacity <- kappa
    } 
    assert_that(is(params, "MizerParams"),
                is.number(lambda),
                is.number(w_pp_cutoff), w_pp_cutoff > 0,
                is.number(n))
    params@resource_params[["lambda"]] <- lambda
    params@resource_params[["n"]] <- n
    params@resource_params[["w_pp_cutoff"]] <- w_pp_cutoff
    
    if (!is.null(resource_capacity) && !is.null(resource_level)) {
        stop("You should specify only either 'resource_level' or 'resource_capacity'.")
    }
    
    # Check and set dynamics function ----
    if (!is.null(resource_dynamics)) {
        assert_that(is.character(resource_dynamics))
        if (!is.function(get0(resource_dynamics))) {
            stop('The resource dynamics function "', resource_dynamics, '" is not defined.')
        }
        params@resource_dynamics <- resource_dynamics
    }
    
    w_full <- w_full(params)
    no_w_full <- length(w_full)
    mu <- getResourceMort(params)
    NR <- initialNResource(params)
    
    # Check resource level ----
    if (!is.null(resource_level)) {
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
        resource_capacity <- NR / resource_level
        resource_capacity[is.nan(resource_level)] <- 0
        comment(resource_capacity) <- comment(resource_level)
    }
    
    # Check growth rate ----
    if (!is.null(resource_rate)) {
        assert_that(is.numeric(resource_rate))
        if (length(resource_rate) == 1) {
            co <- comment(resource_rate)
            resource_rate <- resource_rate * w_full ^ (n - 1)
            comment(resource_rate) <- co
        } else if (length(resource_rate) != no_w_full) {
            stop("The 'resource_rate' should have length 1 or length ",
                 no_w_full, ".")
        }
        if (any(resource_rate < 0)) {
            stop("The 'resource_rate' must always be non-negative.")
        }
    }
    
    # Check capacity ----
    if (!is.null(resource_capacity)) {
        assert_that(is.numeric(resource_capacity))
        if (length(resource_capacity) == 1) {
            co <- comment(resource_capacity)
            resource_capacity <- resource_capacity * w_full ^ (-lambda)
            resource_capacity[w_full >= params@resource_params$w_pp_cutoff] <- 0
            comment(resource_capacity) <- co
        } else if (length(resource_capacity) != no_w_full) {
            stop("The 'resource_rate' should have length 1 or length ",
                 no_w_full, ".")
        }
        if (any(resource_capacity < 0)) {
            stop("The 'resource_capacity' must never be negative.")
        }
    }
    
    # Balance ----
    balance_fn <- get0(paste0("balance_", params@resource_dynamics))
    if (is.null(balance)) {
        balance <- is.function(balance_fn)
    }
    if (balance) {
        # check number of arguments
        num_args <- (!is.null(resource_rate)) +
            (!is.null(resource_capacity))
        if (num_args > 1) {
            stop("You should only provide either the `resource_rate` or `resource_capacity` (or `resource_level`) because the other is determined by the requirement that the resource replenishes at the same rate at which it is consumed.")
        }
        if (num_args == 0) {
            # no values given, so use previous resource_rate
            resource_rate <- params@rr_pp
        }
        
        # For balancing the resource capacity must be above current abundance 
        # except where both are zero
        if (!is.null(resource_capacity) &&
            any(resource_capacity <= NR & NR > 0)) {
            stop("The 'resource_capacity' must always be greater than current resource number density.")
        }
        
        balance_fn <- get0(paste0("balance_", params@resource_dynamics))
        if (!is.function(balance_fn)) {
            stop("There is no balancing function available for ",
                 params@resource_dynamics, 
                 ". You should not set `balance = TRUE`.")
        }
        balance <- balance_fn(params,
                              resource_rate = resource_rate,
                              resource_capacity = resource_capacity)
        resource_rate <- balance$resource_rate
        resource_capacity <- balance$resource_capacity
    }
    
    # Set rates
    if (!is.null(resource_rate)) {
        params@rr_pp[] <- resource_rate
        comment(params@rr_pp) <- comment(resource_rate)
    }
    if (!is.null(resource_capacity)) {
        params@cc_pp[] <- resource_capacity
        comment(params@cc_pp) <- comment(resource_capacity)
    }
    
    params@time_modified <- lubridate::now()
    return(params)
}

#' @rdname setResource
#' @export
resource_rate <- function(params) {
    params@rr_pp
}

#' @rdname setResource
#' @param value The desired new value for the respective parameter.
#' @export
`resource_rate<-` <- function(params, value) {
    setResource(params, resource_rate = value)
}

#' @rdname setResource
#' @export
resource_capacity <- function(params) {
    params@cc_pp
}

#' @rdname setResource
#' @export
`resource_capacity<-` <- function(params, value) {
    setResource(params, resource_capacity = value)
}


#' @rdname setResource
#' @export
resource_level <- function(params) {
    params@initial_n_pp / params@cc_pp
}

#' @rdname setResource
#' @export
`resource_level<-` <- function(params, value) {
    setResource(params, resource_level = value)
}


#' @rdname setResource
#' @export
resource_dynamics <- function(params) {
    params@resource_dynamics
}


#' @rdname setResource
#' @export
#' @examples
#' params <- NS_params
#' resource_dynamics(params)
#' resource_dynamics(params) <- "resource_constant"
`resource_dynamics<-` <- function(params, value) {
    setResource(params, resource_dynamics = value)
}
