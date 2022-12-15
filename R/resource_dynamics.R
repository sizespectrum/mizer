
#' Set or get resource dynamics function
#' 
#' Mizer can model the dynamics of the resource spectrum in different ways,
#' depending on the choice of the resource dynamics function. This function is
#' called at each time step to calculate the value of the resource spectrum at
#' the next time step. Mizer has three built-in resource dynamics functions:
#' [resource_semichemostat()], [resource_logistic()] and [resource_constant()].
#' But you can write your own function. 
#' 
#' @param params A MizerParams object
#' @export
#' @family resource parameters
resource_dynamics <- function(params) {
    params@resource_dynamics
}

#' @rdname resource_dynamics
#' @param value A string with the name of the resource dynamics function.
#' @export
#' @examples
#' params <- NS_params
#' resource_dynamics(params)
#' resource_dynamics(params) <- "resource_constant"
`resource_dynamics<-` <- function(params, value) {
    assert_that(is(params, "MizerParams"),
                is.character(value))
    if (!is.function(get0(value))) {
        stop('The resource dynamics function "', value, '" is not defined.')
    }
    params@resource_dynamics <- value
    
    params@time_modified <- lubridate::now()
    return(params)
    setResource(params, resource_dynamics = value)
}

#' @rdname resource_dynamics
#' @export
getResourceDynamics <- function(params) {
    params@resource_dynamics
}


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
#' [setResourceSemichemostat()] or [setResourceLogistic()] (dependent on your
#' choice of dynamics). The `resource_params` list only contains values that are
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



#' @rdname resource_params
#' @param params A MizerParams object
#'
#' @return A vector of the same length as `w_full(params)` with the value of the
#' resource dynamics parameter in each size class.
#' @export
getResourceLevel <- function(params) {
    assert_that(is(params, "MizerParams"))
    initialNResource(params) / getResourceCapacity(params)
}

#' @rdname resource_params
#' @export
getResourceRate <- function(params) {
    params@rr_pp
}

#' @rdname resource_params
#' @export
getResourceCapacity <- function(params) {
    params@cc_pp
}
