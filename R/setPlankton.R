#' Set up plankton
#' 
#' Sets the intrinsic plankton growth rate and the intrinsic plankton carrying
#' capacity as well as the name of the function used to simulate the plankton
#' dynamics
#' 
#' @section Setting plankton dynamics:
#' By default, mizer uses a semichemostat model to describe the plankton
#' dynamics in each size class independently. This semichemostat dynamics is implemented
#' by the function [plankton_semichemostat()]. You can change the
#' plankton dynamics by writing your own function, modelled on
#' [plankton_semichemostat()], and then passing the name of your
#' function in the `plankton_dynamics` argument.
#' 
#' The `r_plankton` argument is a vector specifying the intrinsic plankton
#' growth rate for each size class. If it is not supplied, then the intrinsic growth
#' rate \eqn{r_p(w)} at size \eqn{w}
#' is set to \deqn{r_p(w) = r_p\, w^{n-1}.}{r_p(w) = r_p w^{n-1}}
#' The values of \eqn{r_p} and \eqn{n} are taken from the `r_plankton`
#' and `n` arguments.
#' 
#' The `K_plankton` argument is a vector specifying the intrinsic plankton
#' carrying capacity for each size class. If it is not supplied, then the 
#' intrinsic carrying capacity \eqn{c_p(w)} at size \eqn{w}
#' is set to \deqn{c_p(w) = \kappa\, w^{-\lambda}}{c_p(w) = \kappa w^{-\lambda}}
#' for all \eqn{w} less than `w_pp_cutoff` and zero for larger sizes.
#' The values of \eqn{\kappa} and \eqn{\lambda} are taken from the `kappa`
#' and `lambda` arguments.
#' 
#' @param params A MizerParams object
#' @param r_plankton Optional. Vector of plankton intrinsic birth rates
#' @param K_plankton Optional. Vector of plankton intrinsic carrying capacity
#' @param r_pp Coefficient of the intrinsic plankton birth rate
#' @param n Allometric growth exponent for plankton
#' @param kappa Coefficient of the intrinsic plankton carrying capacity
#' @param lambda Scaling exponent of the intrinsic plankton carrying capacity
#' @param w_pp_cutoff The upper cut off size of the plankton spectrum. 
#'   Default is 10 g.
#' @param plankton_dynamics Function that determines plankton dynamics by
#'   calculating the plankton spectrum at the next time step from the current
#'   state.
#' @param ... Unused
#' 
#' @return A MizerParams object with updated plankton parameters. Because of the
#'   way the R language works, `setPlankton()` does not make the changes to the
#'   params object that you pass to it but instead returns a new params object.
#'   So to affect the change you call the function in the form
#'   `params <- setPlankton(params, ...)`.
#' @export
#' @family functions for setting parameters
setPlankton <- function(params,
                        r_plankton = NULL,
                        K_plankton = NULL,
                        r_pp = plankton_params(params)[["r_pp"]],
                        kappa = plankton_params(params)[["kappa"]],
                        lambda = plankton_params(params)[["lambda"]],
                        n = plankton_params(params)[["n"]],
                        w_pp_cutoff = plankton_params(params)[["w_pp_cutoff"]],
                        plankton_dynamics = NULL,
                        ...) {
    assert_that(is(params, "MizerParams"),
                is.number(kappa), kappa > 0,
                is.number(lambda),
                is.number(r_pp), r_pp > 0,
                is.number(w_pp_cutoff),
                is.number(n))
    params@plankton_params[["kappa"]] <- kappa
    params@plankton_params[["lambda"]] <- lambda
    params@plankton_params[["r_pp"]] <- r_pp
    params@plankton_params[["n"]] <- n
    params@plankton_params[["w_pp_cutoff"]] <- w_pp_cutoff
    # weight specific plankton growth rate
    if (!is.null(r_plankton)) {
        assert_that(is.numeric(r_plankton),
                    identical(length(r_plankton), length(params@rr_pp)))
        params@rr_pp[] <- r_plankton
        comment(params@rr_pp) <- comment(r_plankton)
    } else {
        rr_pp <- r_pp * params@w_full^(n - 1)
        if (!is.null(comment(params@rr_pp)) &&
            any(params@rr_pp != rr_pp)) {
            message("The plankton intrinsic growth rate has been commented and therefore will ",
                    "not be recalculated from the plankton parameters.")
        } else {
            params@rr_pp[] <- rr_pp
        }
    }
    # the plankton carrying capacity
    if (!is.null(K_plankton)) {
        assert_that(is.numeric(K_plankton),
                    identical(length(K_plankton), length(params@cc_pp)))
        params@cc_pp[] <- K_plankton
        comment(params@cc_pp) <- comment(K_plankton)
    } else {
        cc_pp <- kappa*params@w_full^(-lambda)
        cc_pp[params@w_full > w_pp_cutoff] <- 0
        if (!is.null(comment(params@cc_pp)) &&
            any(params@cc_pp != cc_pp)) {
            message("The plankton carrying capacity has been commented and therefore will ",
                    "not be recalculated from the plankton parameters.")
        } else {
            params@cc_pp[] <- cc_pp
        }
    }
    if (!is.null(plankton_dynamics)) {
        assert_that(is.character(plankton_dynamics))
        if (!is.function(get0(plankton_dynamics))) {
            stop("The function ", plankton_dynamics, "is not defined.")
        }
        params@plankton_dynamics <- plankton_dynamics
    }
    
    return(params)
}

#' @rdname setPlankton
#' @export
getPlanktonRate <- function(params) {
    params@rr_pp
}

#' @rdname setPlankton
#' @export
getPlanktonCapacity <- function(params) {
    params@cc_pp
}

#' @rdname setPlankton
#' @export
getPlanktonDynamics <- function(params) {
    params@plankton_dynamics
}

#' @rdname setPlankton
#' @export
plankton_params <- function(params) {
    params@plankton_params
}

#' @rdname setPlankton
#' @param value List of plankton parameters
#' @export
`plankton_params<-` <- function(params, value) {
    assert_that(
        is(params, "MizerParams"),
        setequal(names(value) == names(params@plankton_params)),
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
    params@plankton_params <- value
    setPlankton(params)
}
