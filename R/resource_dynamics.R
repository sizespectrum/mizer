#' Unstructured resources
#'
#' Besides the size-structured planktonic resource, mizer can also model
#' unstructured resource components. Such unstructured components are
#' appropriate whenever the predation on this component is not size based.
#' Examples include detritus as a resource for detritovores, carrion as a
#' resource for scavengers, or macroflora on which fish can graze. 
#' 
#' Mizer allows an arbitrary number of unstructured components. The biomasses of
#' all components are collected into a vector and stored in the `B` slot of the
#' `MizerSim` object. The initial value `B` is stored in the `initial_B` slot of
#' the `MizerParams` object but can also be specified explicitly as an argument
#' to [project()].
#' 
#' @section Setting up a Model with Resources:
#' You can set up a Mizer model with unstructured resource components using
#' the [multispeciesParams()] function. The arguments to that function that
#' relate to the resources are:
#' * `rho`
#' * `resource_names`
#' * `resource_dyn`
#' * `resource_params`
#' 
#' See the documentation of [multispeciesParams()] for details.
#'
#' @section Feeding on Resources:
#' We denote by \eqn{\rho_{id}(w)} the rate at which a predator of species
#' \eqn{i} and weight \eqn{w} encounters biomass from the d-th unstructured
#' resource component. The contribution of the unstructured resources to the
#' total rate at which biomass encountered is thus 
#' \deqn{\sum_d \rho_{id}(w) B_d.} 
#' The values of \eqn{\rho{id}(w)} are stored in the `rho` slot of the
#' `MizerParams` object, which is a 3-dim array (predator species x resource x
#' predator weight).
#' 
#' @section Resource Dynamics:
#' During a simulation using [project()], the biomasses of the resources is
#' updated at each time step by calling the function specified in the
#' `resource_dyn` slot of the `MizerParams` object. Mizer provides some
#' functions that can be used: [detritus_dyn()], [carrion_dyn()], and 
#' [dead_matter_dyn()]. 
#' 
#' As you can see in the documentation of these functions,
#' their arguments are: the `MizerParams` object `params`, the current fish size
#' spectra `n`, the plankton spectrum `n_pp`, the current resource biomasses
#' `B` and the current rates calculated by the [get_rates()] function. This
#' last argument is passed for efficiency reasons, so that the rates do not
#' have to be recomputed. See the documentation of [get_rates()] for a list of 
#' what it contains. 
#' 
#' The resource dynamics will depend on some model parameters, like for
#' example growth rates. These 
#' in the `species_params` slot of `params`.
#' 
#' @name resources
#' @md
NULL

#' Detritus dynamics
#' 
#' Calculates the detritus biomass at the next timestep from the current
#' detritus biomass.
#' 
#' The inflow of biomass consists of
#' * an external inflow at a rate given by 
#'   `params@resource_params$detritus_external` and
#' * an inflow of feces, calculated as a proportion 
#'   `params@resource_params$detritus_proportion` of the biomass consumed.
#' 
#' The outflow of biomass is due to 
#' * consumption by detritivorous species, where the encounter rate is 
#'   specified by `params@rho[, "detritus", ]`.
#' 
#' This function can be used as the project_resources function in the case
#' where there is only a single unstructured resource that represents 
#' detritus, or as part of a more comprehensive function when there are also
#' other unstructured resource components. It is used for example in the
#' [dead_matter_dyn()] function.
#' 
#' @param params A [MizerParams] object
#' @param n A matrix of species abundances (species x size)
#' @param n_pp A vector of the plankton abundance by size
#' @param B The biomass of detritus
#' @param rates A list of rates as returned by [getRates()]
#' @param dt Time step
#' @param ... Unused
#'   
#' @return Biomass of detritus at next time step
#' @export
#' @md
detritus_dyn <- function(params, n, n_pp, B, rates, dt, ...) {
    consumption <- 
        B * sum((params@rho[, "detritus", ] * n * (1 - rates$feeding_level)) %*%
               params@dw)
    creation <- params@resource_params$detritus_external +
        params@resource_params$detritus_proportion *
          sum((rates$feeding_level * params@intake_max * n) %*% params@dw)
    return(B + (creation - consumption) * dt)
}


#' Carrion dynamics
#' 
#' Calculates the biomass of carrion (dead animals) at the next timestep from
#' the current biomass.
#' 
#' The inflow of carrion biomass consists of 
#' * An external inflow at a rate given by 
#'   `params@resource_params$carrion_external`.
#' * Discards from fishing.
#' * Animals killed by fishing gear.
#' * Animals that have died by causes other than predation.
#' 
#' The outflow of biomass is due to 
#' * consumption by scavenger species, where the encounter rate is 
#'   specified by `params@rho[, "carrion", ]`.
#' 
#' This function can be used as the project_resources function in the case
#' where there is only a single unstructured resource that represents 
#' carrion, or as part of a more comprehensive function when there are also
#' other unstructured resource components. It is used for example in the
#' [dead_matter_dyn()] function.
#' 
#' @param params A [MizerParams] object
#' @param n A matrix of species abundances (species x size)
#' @param n_pp A vector of the plankton abundance by size
#' @param B The biomass of carrion
#' @param rates A list of rates as returned by [getRates()]
#' @param dt Time step
#' @param ... Unused
#'   
#' @return Biomass of carrion at next time step
#' @export
#' @md
carrion_dyn <- function(params, n, n_pp, B, rates, dt, ...) {
    consumption <- 
        B * sum((params@rho[, "carrion", ] * n * (1 - rates$feeding_level)) %*%
                   params@dw)
    creation <- params@resource_params$carrion_external +
        # the other parts still need to be written
        0
    return(B + (creation - consumption) * dt)
}

#' Project two-component dead-matter biomass
#' 
#' This function can be used as the `resource_dyn` function in the case
#' where there are two unstructured resource components named "detritus"
#' and "dead_animals". It gives the biomass of these components at the next
#' time step.
#' 
#' This function in turn calls [detritus_dyn()] and
#' [carrion_dyn()]. Hence the parameter list in `params@resource_params`
#' needs to contain the entries required by these functions.
#' 
#' @param params A [MizerParams] object
#' @param n A matrix of species abundances (species x size)
#' @param n_pp A vector of the plankton abundance by size
#' @param B Vector of biomasses with named components "detritus" and
#'   "carrion".
#' @param rates A list of rates as returned by [getRates()]
#' @param dt Time step
#' @param ... Unused
#'   
#' @return Named vector of biomasses at next time step
#' @export
#' @md
dead_matter_dyn <- function(params, n, n_pp, B, rates, dt, ...) {
    Bn <- B
    Bn["detritus"] <- detritus_dyn(params, n, n_pp, B["detritus"],
                                   rates, dt)
    Bn["carrion"] <- carrion_dyn(params, n, n_pp, B["carrion"],
                                 rates, dt)
    return(Bn)
}