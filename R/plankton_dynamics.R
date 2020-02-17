#' Project plankton using semichemostat model
#' 
#' @param params A [MizerParams] object
#' @param n A matrix of species abundances (species x size)
#' @param n_pp A vector of the plankton abundance by size
#' @param n_other A list with the abundances of other components
#' @param rates A list of rates as returned by [getRates()]
#' @param t The current time
#' @param dt Time step
#' @param ... Unused
#'   
#' @return Vector containing plankton spectrum at next timestep
#' @export
#' @md
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- newMultispeciesParams(NS_species_params_gears, inter,
#'                                 plankton_dynamics = "plankton_semichemostat")
#' }
plankton_semichemostat <- function(params, n, n_pp, n_other, rates, t, dt, ...) {
    # We use the exact solution under the assumption of constant mortality 
    # during timestep
    tmp <- params@rr_pp * params@cc_pp / (params@rr_pp + rates$plankton_mort)
    return(tmp - (tmp - n_pp) * exp(-(params@rr_pp + rates$plankton_mort) * dt))
}


#' Keep plankton abundance constant
#' 
#' This function can be used instead of the standard 
#' \code{\link{plankton_semichemostat}} in order to keep the Plankton
#' spectrum constant over time.
#' 
#' @inheritParams plankton_semichemostat
#' @param ... Unused
#'   
#' @return Vector containing plankton spectrum at next timestep
#' @export
#' @md
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- newMultispeciesParams(NS_species_params_gears, inter,
#'                                 plankton_dynamics = "plankton_constant")
#' }
plankton_constant <- function(params, n, n_pp, n_other, rates, t, dt, ...) {
    return(n_pp)
}
