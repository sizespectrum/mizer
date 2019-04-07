#' Description of the dynamics of unstructured resources
#'
#' Besides the size-structured planktonic resource, mizer can also model
#' unstructured resource components. Such unstructured components are
#' appropriate whenever the predation on these componente is not size based.
#' Examples include detritus as a resource for detritovores, carrion as a
#' resource for scavengers, or macroflora on which fish can graze. 
#' 
#' Mizer allows an arbitrary number of unstructured components. The biomasses of
#' all components are collected into a named vector and stored in the `B` slot
#' of the `MizerSim` object. The initial value for `B` is stored in the
#' `initial_B` slot of the `MizerParams` object but can also be specified
#' explicitly as an argument to [project()].
#' 
#' @section Setting up a Model with Resources:
#' You can set up a Mizer model with unstructured resource components using
#' the [set_multispecies_model()] function. The arguments to that function that
#' relate to the resources are:
#' * `rho` The encounter rate, see "Feeding on Resources" section below.
#' * `resource_dynamics` Functions to update resource biomass during 
#'     simulations, see "Resource Dynamics" section below.
#' * `resource_params` Model parameters for the resource dynamics, see 
#'     "Resource Dynamics" section below.
#' 
#' You can also set up the resource dynamics by hand by overwriting the
#' following slots in an existing `MizerParams` object:
#' * `rho` A 3-dim array described in "Feeding on Resources" dynamics below
#' * `resource_dynamics`
#' * `resource_params`
#'
#' @section Feeding on Resources:
#' We denote by \eqn{\rho_{id}(w)} the rate at which a predator of species
#' \eqn{i} and size \eqn{w} encounters biomass from the d-th unstructured
#' resource component. The contribution of the unstructured resources to the
#' total rate at which biomass is encountered by an individual of species \eqn{i}
#' and size \eqn{w} is thus
#' \deqn{\sum_d \rho_{id}(w) B_d.}
#' Resource consumption is subject to satiation in the same way as other food,
#' so that a consumer only consumes a fraction \eqn{1-f_i(w)} of the encountered
#' resource biomass, where \eqn{f_i(w)} is the feeding level.
#' 
#' A sensible assumption for the size dependence of the encounter rate is that
#' it scales allometrically with the same exponent as the maximum intake rate,
#' so that
#' \deqn{\rho_{id}(w) = \rho_{id} w^n.}
#' This is the default set by the [set_multispecies_model()]
#' function. However you can always overwrite it with your own choice using the
#' [changeResourceEncounter()] function.
#' 
#' @section Resource Dynamics:
#' During a simulation using [project()], the biomasses of the resources is
#' updated at each time step by calling the functions specified in the
#' `resource_dynamics` slot of the `MizerParams` object. This slot holds a
#' vector of functions, one for each unstructured resource component. Mizer
#' provides two functions that you can use: [detritus_dynamics()] and
#' [carrion_dynamics()], but you can easily implement others by following those
#' templates.
#' 
#' As you can see in the documentation of these functions, their arguments are:
#' the `MizerParams` object `params`, the current fish size spectra `n`, the
#' current plankton spectrum `n_pp`, the current resource biomasses `B` and the current
#' rates calculated by the [getRates()] function. This last argument is passed
#' for efficiency reasons, so that the rates do not have to be recomputed. See
#' the documentation of [getRates()] for a list of what it contains.
#' 
#' The other arguments to the resource dynamics functions are model parameters, 
#' like for example growth rates.
#' These need to be stored in the `resource_params` slot of `params`. One model
#' parameter that should always be present is the rate of change due to external
#' causes. This should be given a name of the form `resource_external` where
#' `resource` should be replaced by the name of the resource, see for example
#' `detritus_external` in [detritus_dynamics()]
#' 
#' When writing your own resource dynamics functions, you can choose any names
#' for your other model parameters, but you must make sure not to use the same
#' name in the function for another resource component. One way to ensure this
#' is to prefix all parameter namse with your resource name.
#' 
#' The dynamics for a resource should always have a loss term accounting for
#' the consumption of the resource. This should always have the form used in the
#' example function [detritus_dynamics()], in order to be in agreement with the
#' feeding by consumers described above in the "Feeding on Resources" section.
#' 
#' @name resource_dynamics
#' @md
NULL

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

