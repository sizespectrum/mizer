#' Project resource using semichemostat model
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
#' [resource_semichemostat()] in order to keep the Resource
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
