#' Set up resource
#' 
#' Sets the intrinsic resource growth rate and the intrinsic resource carrying
#' capacity as well as the name of the function used to simulate the resource
#' dynamics
#' 
#' @section Setting resource dynamics:
#' By default, mizer uses a semichemostat model to describe the resource
#' dynamics in each size class independently. This semichemostat dynamics is implemented
#' by the function [resource_semichemostat()]. You can change the
#' resource dynamics by writing your own function, modelled on
#' [resource_semichemostat()], and then passing the name of your
#' function in the `resource_dynamics` argument.
#' 
#' The `resource_rate` argument is a vector specifying the intrinsic resource
#' growth rate for each size class. If it is not supplied, then the intrinsic growth
#' rate \eqn{r(w)} at size \eqn{w}
#' is set to \deqn{r(w) = r_{pp}\, w^{n-1}.}{r(w) = r_pp w^{n-1}}
#' The values of \eqn{r_{pp}} and \eqn{n} are taken from the `r_pp`
#' and `n` arguments.
#' 
#' The `resource_capacity` argument is a vector specifying the intrinsic resource
#' carrying capacity for each size class. If it is not supplied, then the 
#' intrinsic carrying capacity \eqn{c(w)} at size \eqn{w}
#' is set to \deqn{c(w) = \kappa\, w^{-\lambda}}{c(w) = \kappa w^{-\lambda}}
#' for all \eqn{w} less than `w_pp_cutoff` and zero for larger sizes.
#' The values of \eqn{\kappa} and \eqn{\lambda} are taken from the `kappa`
#' and `lambda` arguments.
#' 
#' @param params A MizerParams object
#' @param resource_rate Optional. Vector of resource intrinsic birth rates 
#' @param resource_capacity Optional. Vector of resource intrinsic carrying 
#'   capacity 
#' @param reset If set to TRUE, then both `resource_rate` and
#'   `resource_capacity` will be reset to the value calculated from the resource
#'   parameters, even if they were previously overwritten with custom values. If
#'   set to FALSE (default) then a recalculation from the resource parameters
#'   will take place only if no custom values have been set.
#' @param r_pp Coefficient of the intrinsic resource birth rate
#' @param n Allometric growth exponent for resource
#' @param kappa Coefficient of the intrinsic resource carrying capacity
#' @param lambda Scaling exponent of the intrinsic resource carrying capacity
#' @param w_pp_cutoff The upper cut off size of the resource spectrum.  The
#'   carrying capacity will be set to 0 above this size.
#'   Default is 10 g.
#' @param resource_dynamics Optional. Name of the function that determines the
#'   resource dynamics by calculating the resource spectrum at the next time
#'   step from the current state. You only need to specify this if you do not
#'   want to use the default [resource_semichemostat()].
#' @param ... Unused
#' 
#' @return `setResource`: A MizerParams object with updated resource parameters
#' @export
#' @seealso [resource_params()]
#' @family functions for setting parameters
setResource <- function(params,
                        resource_rate = NULL,
                        resource_capacity = NULL,
                        reset = FALSE,
                        r_pp = resource_params(params)[["r_pp"]],
                        kappa = resource_params(params)[["kappa"]],
                        lambda = resource_params(params)[["lambda"]],
                        n = resource_params(params)[["n"]],
                        w_pp_cutoff = resource_params(params)[["w_pp_cutoff"]],
                        resource_dynamics = NULL,
                        ...) {
    assert_that(is(params, "MizerParams"),
                is.flag(reset),
                is.number(kappa), kappa > 0,
                is.number(lambda),
                is.number(r_pp), r_pp > 0,
                is.number(w_pp_cutoff), w_pp_cutoff > 0,
                is.number(n))
    params@resource_params[["kappa"]] <- kappa
    params@resource_params[["lambda"]] <- lambda
    params@resource_params[["r_pp"]] <- r_pp
    params@resource_params[["n"]] <- n
    params@resource_params[["w_pp_cutoff"]] <- w_pp_cutoff
    
    if (reset) {
        if (!is.null(resource_rate)) {
            warning("Because you set `reset = TRUE`, the value you provided ", 
                    "for `resource_rate` will be ignored and a value will be ",
                    "calculated from the resource parameters.")
            resource_rate <- NULL
        }
        comment(params@rr_pp) <- NULL
        if (!is.null(resource_capacity)) {
            warning("Because you set `reset = TRUE`, the value you provided ", 
                    "for `resource_capacity` will be ignored and a value will be ",
                    "calculated from the resource parameters.")
            resource_capacity <- NULL
        }
        comment(params@cc_pp) <- NULL
    }
    
    # weight specific resource growth rate
    if (!is.null(resource_rate)) {
        if (is.null(comment(resource_rate))) {
            if (is.null(comment(params@rr_pp))) {
                comment(resource_rate) <- "set manually"
            } else {
                comment(resource_rate) <- comment(params@rr_pp)
            }
        }
        assert_that(is.numeric(resource_rate),
                    identical(length(resource_rate), length(params@rr_pp)))
        params@rr_pp[] <- resource_rate
        comment(params@rr_pp) <- comment(resource_rate)
    } else {
        rr_pp <- r_pp * params@w_full^(n - 1)
        if (!is.null(comment(params@rr_pp)) &&
            different(params@rr_pp, rr_pp)) {
            message("The resource intrinsic growth rate has been commented and therefore will ",
                    "not be recalculated from the resource parameters.")
        } else {
            params@rr_pp[] <- rr_pp
        }
    }
    # the resource carrying capacity
    if (!is.null(resource_capacity)) {
        if (is.null(comment(resource_capacity))) {
            if (is.null(comment(params@cc_pp))) {
                comment(resource_capacity) <- "set manually"
            } else {
                comment(resource_capacity) <- comment(params@cc_pp)
            }
        }
        assert_that(is.numeric(resource_capacity),
                    identical(length(resource_capacity), length(params@cc_pp)))
        params@cc_pp[] <- resource_capacity
        comment(params@cc_pp) <- comment(resource_capacity)
    } else {
        cc_pp <- kappa*params@w_full^(-lambda)
        cc_pp[params@w_full > w_pp_cutoff] <- 0
        if (!is.null(comment(params@cc_pp)) &&
            different(params@cc_pp, cc_pp)) {
            message("The resource carrying capacity has been commented and therefore will ",
                    "not be recalculated from the resource parameters.")
        } else {
            params@cc_pp[] <- cc_pp
        }
    }
    if (!is.null(resource_dynamics)) {
        assert_that(is.character(resource_dynamics))
        if (!is.function(get0(resource_dynamics))) {
            stop('The resource dynamics function "', resource_dynamics, '" is not defined.')
        }
        params@resource_dynamics <- resource_dynamics
    }
    
    params@time_modified <- lubridate::now()
    return(params)
}

#' @rdname setResource
#' @export
getResourceRate <- function(params) {
    params@rr_pp
}

#' @rdname setResource
#' @export
resource_rate <- function(params) {
    params@rr_pp
}

#' @rdname setResource
#' @param value .
#' @export
`resource_rate<-` <- function(params, value) {
    setResource(params, resource_rate = value)
}

#' @rdname setResource
#' @export
getResourceCapacity <- function(params) {
    params@cc_pp
}

#' @rdname setResource
#' @export
resource_capacity <- function(params) {
    params@cc_pp
}

#' @rdname setResource
#' @param value .
#' @export
`resource_capacity<-` <- function(params, value) {
    setResource(params, resource_capacity = value)
}

#' @rdname setResource
#' @export
getResourceDynamics <- function(params) {
    params@resource_dynamics
}

#' @rdname setResource
#' @export
resource_dynamics <- function(params) {
    params@resource_dynamics
}

#' @rdname setResource
#' @param value .
#' @export
`resource_dynamics<-` <- function(params, value) {
    setResource(params, resource_dynamics = value)
}

#' Resource parameters
#' 
#' These functions allow you to get or set the resource parameters stored in a
#' MizerParams object. The resource parameters are stored as a named list with
#' the slot names `r_pp`, `kappa`, `lambda`, `n`, `w_pp_cutoff`. For their
#' meaning see Details below. If you change these parameters then this will
#' recalculate the resource rate and the resource capacity, unless you have set
#' custom values for these. If you have specified a different resource dynamics
#' function that requires additional parameters, then these should also be added
#' to the `resource_params` list.
#' 
#' The resource parameters `r_pp` and `n` are used to set the intrinsic
#' replenishment rate \eqn{r_R(w)} for the resource at size \eqn{w} to 
#' \deqn{r_R(w) = r_{pp}\, w^{n-1}.}{r_R(w) = r_pp w^{n-1}}
#' 
#' The resource paramters `kappa`, `lambda` and `w_pp_cutoff` are used to set
#' the intrinsic resource carrying capacity capacity \eqn{c_R(w)} at size \eqn{w}
#' is set to
#' \deqn{c_R(w) = \kappa\, w^{-\lambda}}{c_R(w) = \kappa w^{-\lambda}}
#' for all \eqn{w} less than `w_pp_cutoff` and zero for larger sizes.
#' 
#' If you use the default semichemostat dynamics for the resource then these
#' rates enter the equation for the resource abundance density as
#' \deqn{\frac{\partial N_R(w,t)}{\partial t} = r_R(w) \Big[ c_R (w) - N_R(w,t) \Big] - \mu_R(w, t) N_R(w,t)}{dN_R(w,t)/d t  = r_R(w) ( c_R (w) - N_R(w,t) ) - \mu_R(w,t ) N_R(w,t)}
#' where the mortality \eqn{\mu_R(w, t)} is
#' due to predation by consumers and is calculate with [getResourceMort()].
#' 
#' You can however set up different resource dynamics with
#' [resource_dynamics<-()].
#' 
#' @param params A MizerParams object
#' @export
#' @family functions for setting parameters
#' @examples 
#' resource_params(NS_params)
resource_params <- function(params) {
    params@resource_params
}

#' @rdname resource_params
#' @param value A named list of resource parameters.
#' @export
#' @examples
#' # Doubling the replenishment rate
#' params <- NS_params
#' resource_params(params)$r_pp <- 2 * resource_params(params)$r_pp
`resource_params<-` <- function(params, value) {
    assert_that(
        is(params, "MizerParams"),
        setequal(names(value), names(params@resource_params)),
        is.number(value$lambda),
        value$lambda >= 0,
        is.number(value$kappa),
        value$kappa >= 0,
        is.number(value$r_pp),
        value$r_pp >= 0,
        is.number(value$n),
        value$n >= 0,
        is.number(value$w_pp_cutoff),
        value$w_pp_cutoff > min(params@w_full),
        value$w_pp_cutoff < max(params@w_full)
    )
    params@resource_params <- value
    setResource(params)
}
