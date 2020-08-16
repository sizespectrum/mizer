#' Project resource using semichemostat model
#' 
#' This function calculates the resource abundance at time `t + dt` from all
#' abundances and rates at time `t`. 
#' 
#' The time evolution of the resource spectrum is described by a 
#' semi-chemostat equation
#' \deqn{\frac{\partial N_R(w,t)}{\partial t} = r_R(w) \Big[ c_R (w) - N_R(w,t) \Big] - \mu_R(w, t) N_R(w,t)}{dN_R(w,t)/d t  = r_R(w) ( c_R (w) - N_R(w,t) ) - \mu_R(w,t ) N_R(w,t)}
#' 
#' Here \eqn{r_R(w)} is the resource regeneration rate and \eqn{c_R(w)} is the
#' carrying capacity in the absence of predation. These parameters are changed
#' with [setResource()]. The mortality \eqn{\mu_R(w, t)} is
#' due to predation by consumers and is calculate with [getResourceMort()].
#' 
#' This function uses the analytic solution of the above equation, keeping the
#' mortality fixed during the timestep.
#' 
#' It is also possible to implement other resource dynamics, as
#' described in the help page for [setResource()].
#' 
#' @param params A [MizerParams] object
#' @param n A matrix of species abundances (species x size)
#' @param n_pp A vector of the resource abundance by size
#' @param n_other A list with the abundances of other components
#' @param rates A list of rates as returned by [mizerRates()]
#' @param t The current time
#' @param dt Time step
#' @param ... Unused
#'   
#' @return Vector containing resource spectrum at next timestep
#' @export
#' @examples
#' \dontrun{
#' params <- newMultispeciesParams(NS_species_params_gears, inter,
#'                                 resource_dynamics = "resource_semichemostat")
#' }
resource_semichemostat <- function(params, n, n_pp, n_other, rates, t, dt, ...) {
    # We use the exact solution under the assumption of constant mortality 
    # during timestep
    tmp <- params@rr_pp * params@cc_pp / (params@rr_pp + rates$resource_mort)
    return(tmp - (tmp - n_pp) * exp(-(params@rr_pp + rates$resource_mort) * dt))
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
#' @examples
#' \dontrun{
#' params <- newMultispeciesParams(NS_species_params_gears, inter,
#'                                 resource_dynamics = "resource_constant")
#' }
resource_constant <- function(params, n, n_pp, n_other, rates, t, dt, ...) {
    return(n_pp)
}
