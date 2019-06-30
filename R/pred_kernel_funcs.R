#' Lognormal predation kernel
#' 
#' This is the most commonly-used predation kernel. The log of the predator/prey
#' mass ratio is normally distributed.
#' 
#' Writing the predator mass as \eqn{w} and the prey mass as \eqn{w_p},
#' the feeding kernel is given as
#' \deqn{\phi_i(w, w_p) = 
#' \exp \left[ \frac{-(\ln(w / w_p / \beta_i))^2}{2\sigma_i^2} \right]
#' }{\phi_i(w/w_p) = exp(-(ln(w/w_p/\beta_i))^2/(2\sigma_i^2))}
#' if \eqn{w/w_p} is between 1 and \eqn{\beta_i\exp(3\sigma_i)}{\beta_i exp(3\sigma_i)}
#' and zero otherwise.
#' Here \eqn{\beta_i} is the preferred predator-prey mass ratio and \eqn{\sigma_i}
#' determines the width of the kernel.
#' These two parameters need to be given in the species parameter dataframe in
#' the columns \code{beta} and \code{sigma}.
#' 
#' This function is called from \code{\link{setPredKernel}} to set up the
#' predation kernel slots in a \code{MizerParams} object. 
#' 
#' @param ppmr A vector of predator/prey size ratios
#' @param beta The preferred predator/prey size ratio
#' @param sigma The width parameter of the log-normal kernel
#' 
#' @return A vector giving the value of the predation kernel at each of the
#'   predator/prey mass ratios in the \code{ppmr} argument.
#' @md
#' @export
lognormal_pred_kernel <- function(ppmr, beta, sigma) {
    Beta <- log(beta)
    phi <- exp(-(log(ppmr) - Beta)^2 / (2 * sigma^2))
    # rr is the maximal log predator/prey mass ratio
    rr <- exp(Beta + 3 * sigma)
    phi[ppmr > rr] <- 0
    # Do not allow feeding at own size
    phi[1] <- 0
    return(phi)
}

#' Box predation kernel
#' 
#' A predation kernel where the predator/prey mass ratio is uniformly
#' distributed on an interval.
#' 
#' Writing the predator mass as \eqn{w} and the prey mass as \eqn{w_p}, the
#' feeding kernel is 1 if \eqn{w/w_p} is between \code{ppmr_min} and
#' \code{ppmr_max} and zero otherwise. The parameters need to be given in the
#' species parameter dataframe in the columns \code{ppmr_min} and
#' \code{ppmr_max}.
#' 
#' @param ppmr A vector of predator/prey size ratios
#' @param ppmr_min Minimum predator/prey mass ratio
#' @param ppmr_max Maximum predator/prey mass ratio
#' 
#' @return A vector giving the value of the predation kernel at each of the
#'   predator/prey mass ratios in the \code{ppmr} argument.
#' @md
#' @export
box_pred_kernel <- function(ppmr, ppmr_min, ppmr_max) {
    assert_that(ppmr_min < ppmr_max)
    phi <- rep(1, length(ppmr))
    phi[ppmr > ppmr_max] <- 0
    phi[ppmr < ppmr_min] <- 0
    # Do not allow feeding at own size
    phi[1] <- 0
    return(phi)
}
