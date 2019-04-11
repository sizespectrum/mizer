#' Detritus dynamics
#' 
#' Calculates the detritus biomass at the next timestep from the current
#' detritus biomass.
#' 
#' The equation for the time evolution of the detritus biomass \eqn{B} is
#' assumed to be of the form
#' \deqn{dB/dt = inflow - consumption * B + external} 
#' where 
#' * `inflow` comes from feces, calculated as a proportion 
#'   `detritus_proportion` of the biomass consumed by all consumers.
#' * `consumption` is by detritivorous species, where the encounter rate is 
#'   specified by `params@rho[, "detritus", ]`.
#' * `external` is an influx from external sources. It can be negative in which
#'   case it represents a loss to external sources.
#' 
#' This equation is solved analytically to
#' \deqn{B(t+dt) = B(t)\exp(-consumption \cdot dt)
#'   +\frac{inflow + external}{consumption}
#'   (1-\exp(-consumption \cdot dt)).}{B(t+dt) 
#'   = B(t) exp(-consumption * dt)
#'   +(inflow + external)/(consumption) * (1 - exp(-consumption * dt)).}
#' This avoids the stability problems that would arise if we used the Euler
#' method to solve the equation numerically.
#' 
#' @param params A [MizerParams] object
#' @param n A matrix of current species abundances (species x size)
#' @param n_pp A vector of current plankton abundance by size
#' @param B A vector of current resource biomasses
#' @param rates A list of rates as returned by [getRates()]
#' @param dt Time step
#' @param detritus_external Rate of change from external sources
#' @param detritus_proportion Proportion of consumption by fish that flows into
#'   the detritus component.
#' @param ... Unused
#'   
#' @return A single number giving the biomass of detritus at next time step
#' @export
#' @family resource dynamics functions
#' @md
detritus_dynamics <- 
    function(params, n, n_pp, B, rates, dt,
             detritus_external = params@resource_params$detritus_external,
             detritus_proportion = params@resource_params$detritus_proportion,
             ...) {
    
    # Total consumption is obtained by multiplying the encounter rate per
    # consumer by the number density of consumers and their feeding rate
    # and then integrating over all consumer weights and summing over all
    # consumer species. This should not be changed because it needs to be in
    # agreement with the calculation of consumption in getEGrowth().
    consumption <- sum((params@rho[, "detritus",] * n *
                            (1 - rates$feeding_level)) %*% params@dw)
    inflow <-
        detritus_proportion *
          sum((rates$feeding_level * params@intake_max * n) %*% params@dw)
    
    et <- exp(consumption * dt)
    return(B["detritus"] * et + 
               (inflow  + detritus_external) / consumption  * (1 - et))
    }


#' Carrion dynamics
#' 
#' Calculates the biomass of carrion (dead animals) at the next timestep from
#' the current biomass.
#' 
#' The equation for the time evolution of the carrion biomass \eqn{B} is
#' assumed to be of the form
#' \deqn{dB/dt = inflow - consumption * B + external} 
#' where 
#' * `inflow` comes from
#'     + Discards from fishing.
#'     + Animals killed by fishing gear.
#'     + Animals that have died by causes other than predation. 
#'   `detritus_proportion` of the biomass consumed by all consumers.
#' * `consumption` is by scavenger species, where the encounter rate is 
#'   specified by `params@rho[, "carrion", ]`.
#' * `external` is an influx from external sources. It can be negative in which
#'   case it represents a loss to external sources.
#' 
#' This equation is solved analytically to
#' \deqn{B(t+dt) = B(t)\exp(-consumption \cdot dt)
#'   +\frac{inflow + external}{consumption}
#'   (1-\exp(-consumption \cdot dt)).}{B(t+dt) 
#'   = B(t) exp(-consumption * dt)
#'   +(inflow + external)/(consumption) * (1 - exp(-consumption * dt)).}
#' This avoids the stability problems that would arise if we used the Euler
#' method to solve the equation numerically.
#' 
#' @param params A [MizerParams] object
#' @param n A matrix of current species abundances (species x size)
#' @param n_pp A vector of current plankton abundance by size
#' @param B A vector of current resource biomasses
#' @param rates A list of rates as returned by [getRates()]
#' @param dt Time step
#' @param carrion_external External inflow rate of carrion biomass
#' @param ... Unused
#'   
#' @return A single number giving the biomass of carrion at next time step
#' @export
#' @family resource dynamics functions
#' @md
carrion_dynamics <- 
    function(params, n, n_pp, B, rates, dt,
             carrion_external = params@resource_params$carrion_external,
             ...) {
        
        consumption <- sum((params@rho[, "carrion", ] * n *
                                (1 - rates$feeding_level)) %*% params@dw)
        inflow <- 
            # still need to be written
            0
        
        et <- exp(consumption * dt)
        return(B["carrion"] * et + 
                   (inflow  + carrion_external) / consumption  * (1 - et))
}

