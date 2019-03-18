#' Lognormal predation kernel
#' 
#' This is the most commonly-used predation kernel. The log of the predator/prey
#' mass ratio is normally distributed.
#' 
#' Writing the predator mass as \eqn{w} and the prey mass as \eqn{w_p},
#' the feeding kernel is given as
#' \deqn{\phi_i(w, w_p) = 
#' \exp \left[ \dfrac{-(\ln(w / w_p / \beta_i))^2}{2\sigma_i^2} \right]
#' }{\phi_i(w/w_p) = exp(-(ln(w/w_p/\beta_i))^2/(2\sigma_i^2))}
#' if \eqn{w/w_p} is between 1 and \eqn{\beta_i\exp(3\sigma_i)}{\beta_i exp(3\sigma_i)}
#' and zero otherwise.
#' Here \eqn{\beta_i} is the preferred predator-prey mass ratio and \eqn{\sigma_i}
#' determines the width of the kernel.
#' These two parameters need to be given in the species parameter dataframe in
#' the colunns \code{beta} and \code{sigma}.
#' 
#' This function is called from \code{\link{defaultPredKernel}} to set up the
#' predation kernel slots in a \code{MizerParams} object. Associated to this
#' function there is the function \code{\link{lognormal_max_ppmr}} that 
#' returns the largest allowed predator/prey mass ratio, i.e.,
#' \eqn{\beta_i\exp(3\sigma_i)}{\beta_i exp(3\sigma_i)}.
#' 
#' @param ppmr A vector of predator/prey size ratios
#' @param sp The index of the predator species
#' @param params A MizerParams object
#' 
#' @return A vector giving the value of the predation kernel at each of the
#'   predator/prey mass ratios in the \code{ppmr} argument.
#' @md
#' @export
lognormal_pred_kernel <- function(ppmr, sp, params) {
    Beta <- log(params@species_params$beta[sp])
    sigma <- params@species_params$sigma[sp]
    phi <- exp(-(log(ppmr) - Beta)^2 / (2 * sigma^2))
    # rr is the maximal log predator/prey mass ratio
    rr <- exp(Beta + 3 * sigma)
    phi[ppmr > rr] <- 0
    # Do not allow feeding at own size
    phi[1] <- 0
    return(phi)
}

#' Maximum predator/prey mass ratio for each species
#' 
#' Gives the maximum predator/prey mass ratio at which the predation kernel
#' computed by \code{\link{lognormal_pred_kernel}} is non-zero.
#' This is used in \code{\link{set_multispecies_model}} to determine the 
#' smallest relevant plankton size if the \code{min_w_pp} argument is not
#' specified there.
#' 
#' Because of how we cut off the Gaussian feeding kernel in
#' \code{\link{lognormal_pred_kernel}}, the maximum predator/prey mass ratio is
#' \eqn{\beta  \exp(3 \sigma)}{\beta exp(3 \sigma)}.
#' 
#' @param species_params A species parameter dataframe from which this function
#'   uses the \code{beta} and \code{sigma} slots.
#' 
#' @return Vector with the maximum predator/prey ratio for each species
#' @export
lognormal_max_ppmr <- function(species_params) {
    #TODO: check validity of beta and sigma
    species_params$beta * exp(3 * species_params$sigma)
}