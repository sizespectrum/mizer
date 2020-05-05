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
#' rate \eqn{r_p(w)} at size \eqn{w}
#' is set to \deqn{r_p(w) = r_p\, w^{n-1}.}{r_p(w) = r_p w^{n-1}}
#' The values of \eqn{r_p} and \eqn{n} are taken from the `resource_rate`
#' and `n` arguments.
#' 
#' The `resource_capacity` argument is a vector specifying the intrinsic resource
#' carrying capacity for each size class. If it is not supplied, then the 
#' intrinsic carrying capacity \eqn{c_p(w)} at size \eqn{w}
#' is set to \deqn{c_p(w) = \kappa\, w^{-\lambda}}{c_p(w) = \kappa w^{-\lambda}}
#' for all \eqn{w} less than `w_pp_cutoff` and zero for larger sizes.
#' The values of \eqn{\kappa} and \eqn{\lambda} are taken from the `kappa`
#' and `lambda` arguments.
#' 
#' @param params A MizerParams object
#' @param resource_rate Optional. Vector of resource intrinsic birth rates
#' @param resource_capacity Optional. Vector of resource intrinsic carrying capacity
#' @param r_pp Coefficient of the intrinsic resource birth rate
#' @param n Allometric growth exponent for resource
#' @param kappa Coefficient of the intrinsic resource carrying capacity
#' @param lambda Scaling exponent of the intrinsic resource carrying capacity
#' @param w_pp_cutoff The upper cut off size of the resource spectrum. 
#'   Default is 10 g.
#' @param resource_dynamics Function that determines resource dynamics by
#'   calculating the resource spectrum at the next time step from the current
#'   state.
#' @param ... Unused
#' 
#' @return A MizerParams object with updated resource parameters. Because of the
#'   way the R language works, `setResource()` does not make the changes to the
#'   params object that you pass to it but instead returns a new params object.
#'   So to affect the change you call the function in the form
#'   `params <- setResource(params, ...)`.
#' @export
#' @family functions for setting parameters
setResource <- function(params,
                        resource_rate = NULL,
                        resource_capacity = NULL,
                        r_pp = resource_params(params)[["r_pp"]],
                        kappa = resource_params(params)[["kappa"]],
                        lambda = resource_params(params)[["lambda"]],
                        n = resource_params(params)[["n"]],
                        w_pp_cutoff = resource_params(params)[["w_pp_cutoff"]],
                        resource_dynamics = NULL,
                        ...) {
    assert_that(is(params, "MizerParams"),
                is.number(kappa), kappa > 0,
                is.number(lambda),
                is.number(r_pp), r_pp > 0,
                is.number(w_pp_cutoff),
                is.number(n))
    params@resource_params[["kappa"]] <- kappa
    params@resource_params[["lambda"]] <- lambda
    params@resource_params[["r_pp"]] <- r_pp
    params@resource_params[["n"]] <- n
    params@resource_params[["w_pp_cutoff"]] <- w_pp_cutoff
    # weight specific resource growth rate
    if (!is.null(resource_rate)) {
        assert_that(is.numeric(resource_rate),
                    identical(length(resource_rate), length(params@rr_pp)))
        params@rr_pp[] <- resource_rate
        comment(params@rr_pp) <- comment(resource_rate)
    } else {
        rr_pp <- r_pp * params@w_full^(n - 1)
        if (!is.null(comment(params@rr_pp)) &&
            any(params@rr_pp != rr_pp)) {
            message("The resource intrinsic growth rate has been commented and therefore will ",
                    "not be recalculated from the resource parameters.")
        } else {
            params@rr_pp[] <- rr_pp
        }
    }
    # the resource carrying capacity
    if (!is.null(resource_capacity)) {
        assert_that(is.numeric(resource_capacity),
                    identical(length(resource_capacity), length(params@cc_pp)))
        params@cc_pp[] <- resource_capacity
        comment(params@cc_pp) <- comment(resource_capacity)
    } else {
        cc_pp <- kappa*params@w_full^(-lambda)
        cc_pp[params@w_full > w_pp_cutoff] <- 0
        if (!is.null(comment(params@cc_pp)) &&
            any(params@cc_pp != cc_pp)) {
            message("The resource carrying capacity has been commented and therefore will ",
                    "not be recalculated from the resource parameters.")
        } else {
            params@cc_pp[] <- cc_pp
        }
    }
    if (!is.null(resource_dynamics)) {
        assert_that(is.character(resource_dynamics))
        if (!is.function(get0(resource_dynamics))) {
            stop("The function ", resource_dynamics, "is not defined.")
        }
        params@resource_dynamics <- resource_dynamics
    }
    
    return(params)
}

#' @rdname setResource
#' @export
getResourceRate <- function(params) {
    params@rr_pp
}

#' @rdname setResource
#' @export
getResourceCapacity <- function(params) {
    params@cc_pp
}

#' @rdname setResource
#' @export
getResourceDynamics <- function(params) {
    params@resource_dynamics
}

#' @rdname setResource
#' @export
resource_params <- function(params) {
    params@resource_params
}

#' @rdname setResource
#' @param value List of resource parameters
#' @export
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
