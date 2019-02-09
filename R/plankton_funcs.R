#' Project plankton using semichemostat model
#' 
#' @param params A [MizerParams] object
#' @param n A matrix of species abundances (species x size)
#' @param n_pp A vector of the plankton abundance by size
#' @param B The biomass of detritus
#' @param rates A list of rates as returned by [getRates()]
#' @param dt Time step
#' @param ... Unused
#'   
#' @return Vector containing plankton spectrum at next timestep
#' @export
#' @md
plankton_semichemostat <- function(params, n, n_pp, B, rates, dt, ...) {
    # Dynamics of plankton spectrum uses a semi-chemostat model (de Roos - ask Ken)
    # We use the exact solution under the assumption of constant mortality during timestep
    tmp <- params@rr_pp * params@cc_pp / (params@rr_pp + rates$plankton_mort)
    return(tmp - (tmp - n_pp) * exp(-(params@rr_pp + rates$plankton_mort) * dt))
}