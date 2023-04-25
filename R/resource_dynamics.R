#' Keep resource abundance constant
#' 
#' If you set your resource dynamics to use this function then the resource
#' abundances are kept constant over time.
#' 
#' 
#' To set your model to keep the resource constant over time you do
#' ```
#' resource_dynamics(params) <- "resource_constant"
#' ```
#' where you should replace `params` with the name of the variable holding your
#' MizerParams object.
#' 
#' @inheritParams resource_semichemostat
#' @param ... Unused
#'   
#' @return Vector containing resource spectrum at next timestep
#' @export
#' @family resource dynamics
#' @examples
#' params <- NS_params
#' resource_dynamics(params) <- "resource_constant"
resource_constant <- function(params, n_pp, ...) {
    return(n_pp)
}



#' Resource parameters
#' 
#' The recommended way to change the resource dynamics parameters is to use
#' [setResource()]. The `resource_params` list contains values that are helpful
#' in setting up the actual size-dependent parameters with [setResource()]. If
#' you have specified a custom resource dynamics function that requires
#' additional parameters, then these should also be added to the
#' `resource_params` list.
#' 
#' The `resource_params` list will at least contain the slots `kappa`, `lambda`,
#' `w_pp_cutoff` and `n`.
#' 
#' The resource parameter `n` is the exponent for the power-law form for the
#' replenishment rate \eqn{r_R(w)}: \deqn{r_R(w) = r_R\, w^{n-1}.}{r_R(w) = r_R
#' w^{n-1}.}
#'
#' The resource parameter `lambda` (\eqn{\lambda}) is the exponent for the
#' power-law form for the carrying capacity \eqn{c_R(w)} and `w_pp_cutoff` is
#' its cutoff value: \deqn{c_R(w) = c_R w^{-\lambda}} for all \eqn{w} less than
#' `w_pp_cutoff` and zero for larger sizes.
#'
#' The resource parameter `kappa` (\eqn{\kappa}) determines the initial resource
#' abundance: \deqn{N_R(w) = \kappa\, w^{-\lambda}}{c_R(w) = \kappa
#' w^{-\lambda}} for all \eqn{w} less than `w_pp_cutoff` and zero for larger
#' sizes.
#' 
#' @param params A MizerParams object
#' @return A named list of resource parameters.
#' @export
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
        is.number(value[["n"]]),
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
#' @return Name of the function that determines the resource dynamics
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
