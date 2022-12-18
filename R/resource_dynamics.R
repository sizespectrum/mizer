#' Keep resource abundance constant
#' 
#' This function can be used instead of the standard 
#' [resource_semichemostat()] in order to keep the resource
#' spectrum constant over time.
#' 
#' @inheritParams resource_semichemostat
#' @param ... Unused
#'   
#' @return Vector containing resource spectrum at next timestep
#' @export
#' @family resource dynamics
#' @examples
#' \dontrun{
#' params <- NS_params
#' resource_dynamics(params) <- "resource_constant"
#' }
resource_constant <- function(params, n_pp, ...) {
    return(n_pp)
}



#' Resource parameters
#' 
#' Functions for working with the resource parameters.
#' 
#' Both the [resource_semichemostat()] and the [resource_logistic()] dynamics
#' are parametrised in terms of a size-dependent rate \eqn{r_R(w)} and a 
#' size-dependent capacity \eqn{c_R}. You can see their current values with
#' [getResourceRate()] and [getResourceCapacity()].
#' 
#' Due to the predation mortality, the actual resource abundance \eqn{N_R} is
#' always lower than the capacity \eqn{c_R}. The ratio \eqn{N_R / c_R} is
#' denoted as the resource level. This can be NaN when both \eqn{N_R} and
#' \eqn{c_R} are zero. You can see the current values with [getResourceLevel()].
#' 
#' The recommended way to change these paramters is with the functions
#' [setResource()]. The `resource_params` list only contains values that are
#' helpful in setting up the actual size-dependent parameters. Also If you have
#' specified a different resource dynamics function that requires additional
#' parameters, then these should also be added to the `resource_params` list.
#' 
#' The `resource_params` list will at least contain the slots
#' `kappa`, `lambda`, `w_pp_cutoff` and `n`.
#' 
#' The resource parameter `n` is the 
#' exponent for the power-law form for the replenishment rate \eqn{r_R(w)}:
#' \deqn{r_R(w) = r_R\, w^{n-1}.}{r_R(w) = r_R w^{n-1}.}
#' 
#' The resource parameter `lambda` (\eqn{\lambda}) is the exponent for the
#' power-law form for the carrying capacity \eqn{c_R(w)} and `w_pp_cutoff` is
#' its cutoff value:
#' \deqn{c_R(w) = c_R w^{-\lambda}}
#' for all \eqn{w} less than `w_pp_cutoff` and zero for larger sizes.
#' 
#' The resource parameter `kappa` (\eqn{\kappa}) determines the initial 
#' resource abundance:
#' \deqn{N_R(w) = \kappa\, w^{-\lambda}}{c_R(w) = \kappa w^{-\lambda}}
#' for all \eqn{w} less than `w_pp_cutoff` and zero for larger sizes. Of course
#' you can overrule this with [initialNResource()].
#' 
#' @param params A MizerParams object
#' @export
#' @family resource parameters
#' @examples 
#' resource_params(NS_params)
resource_params <- function(params) {
    params@resource_params
}

#' @rdname resource_params
#' @param value A named list of resource parameters.
#' @export
`resource_params<-` <- function(params, value) {
    assert_that(
        is(params, "MizerParams"),
        is.number(value$lambda),
        value$lambda >= 0,
        is.number(value$kappa),
        value$kappa >= 0,
        is.number(value$n),
        value$n >= 0,
        is.number(value$w_pp_cutoff),
        value$w_pp_cutoff > min(params@w_full),
        value$w_pp_cutoff < max(params@w_full)
    )
    params@resource_params <- value
    params
}



#' Deprecated functions for getting resource parameters
#' 
#' `r lifecycle::badge("deprecated")` Use [resource_dynamics()],
#' [resource_level()], [resource_rate()] and [resource_capacity()] instead.
#' 
#' @param params A MizerParams object
#' @keywords internal
#' @export
getResourceDynamics <- function(params) {
    lifecycle::deprecate_warn("2.4.0", "getResourceDynamics()", 
                              "resource_dynamics()")
    resource_dynamics(params)
}

#' @rdname getResourceDynamics
#' @keywords internal
#' @export
getResourceLevel <- function(params) {
    lifecycle::deprecate_warn("2.4.0", "getResourceLevel()", 
                              "resource_level()")
    resource_level(params)
}

#' @rdname getResourceDynamics
#' @keywords internal
#' @export
getResourceRate <- function(params) {
    lifecycle::deprecate_warn("2.4.0", "getResourceRate()", 
                              "resource_rate()")
    resource_rate(params)
}

#' @rdname getResourceDynamics
#' @keywords internal
#' @export
getResourceCapacity <- function(params) {
    lifecycle::deprecate_warn("2.4.0", "getResourceCapacity()", 
                              "resource_capacity()")
    resource_capacity(params)
}
