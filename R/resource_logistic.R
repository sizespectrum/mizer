#' Project resource using logistic model
#' 
#' This function calculates the resource abundance at time `t + dt` from all
#' abundances and rates at time `t`. 
#' 
#' The time evolution of the resource spectrum is described by a 
#' logistic equation
#' \deqn{\frac{\partial N_R(w,t)}{\partial t} = r_R(w) N_R(w)\Big[ 1 - \frac{N_R(w,t)}{c_R (w)} \Big] - \mu_R(w, t) N_R(w,t)}{dN_R(w,t)/d t  = r_R(w) N_r(w)( 1 - N_R(w,t) / c_R (w)) - \mu_R(w,t ) N_R(w,t)}
#' 
#' Here \eqn{r_R(w)} is the resource regeneration rate and \eqn{c_R(w)} is the
#' carrying capacity in the absence of predation. These parameters are changed
#' with [setResourceLogistic()]. The mortality \eqn{\mu_R(w, t)} is
#' due to predation by consumers and is calculate with [getResourceMort()].
#' 
#' This function uses the analytic solution of the above equation, keeping the
#' mortality fixed during the timestep.
#' 
#' It is also possible to implement other resource dynamics, as
#' described in the help page for [resource_dynamics()<-].
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
#' \dontrun{
#' params <- newMultispeciesParams(NS_species_params_gears, NS_interaction,
#'                                 resource_dynamics = "resource_logistic")
#' }
resource_logistic <- function(params, n, n_pp, n_other, rates, t, dt,
                              resource_rate, resource_capacity, ...) {
    # We use the exact solution under the assumption of constant mortality 
    # during timestep
    mur <- resource_rate - rates$resource_mort
    n_steady <- mur * resource_capacity / resource_rate
    n_pp_new <- n_steady / ((n_steady - n_pp) / n_pp * exp(-mur * dt) + 1)
    
    # Here is an alternative expression that looks as if it might be more
    # precise when the sum of the rates is small due to the use of expm1.
    # However the above has the advantage of preserving the steady state
    # n_steady exactly.
    # n_pp_new <- n_steady / (n_steady / n_pp * exp(-mur * dt) +  expm1(-mur * dt))
    
    # if growth rate and death rate are zero then the above would give NaN
    # whereas the value should simply not change
    sel <- is.nan(n_pp_new)
    n_pp_new[sel] <- n_pp[sel]
    
    n_pp_new
}