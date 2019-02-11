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
        B* sum((params@rho[, "detritus", ] * n * (1 - rates$feeding_level)) %*%
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
#' This function can be used as the `project_resources` function in the case
#' where there are two unstructured resource components named "detritus"
#' and "dead_animals". It gives the biomass of these components at the next
#' time step.
#' 
#' This function in turn calls [project_detritus()] and
#' [project_carrion()]. Hence the parameter list in `params@resource_params`
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