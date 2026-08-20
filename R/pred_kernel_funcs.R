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
#' if \eqn{w/w_p} is larger than 1 and zero otherwise. Here \eqn{\beta_i} is the
#' preferred predator-prey mass ratio and \eqn{\sigma_i} determines the width of
#' the kernel. These two parameters need to be given in the species parameter
#' dataframe in the columns \code{beta} and \code{sigma}.
#'
#' This function is called from [setPredKernel()] to set up the
#' predation kernel slots in a MizerParams object.
#'
#' @param ppmr A vector of predator/prey size ratios
#' @param beta The preferred predator/prey size ratio
#' @param sigma The width parameter of the log-normal kernel
#'
#' @return A vector giving the value of the predation kernel at each of the
#'   predator/prey mass ratios in the \code{ppmr} argument.
#' @export
#' @family predation kernel
#' @seealso [setPredKernel()]
#' @examples
#' params <- NS_params
#' plot(w_full(params), pred_kernel(params)["Cod", 10, ], type="l", log="x")
#' # The restriction that the kernel is zero for w/w_p < 1 is more
#' # noticeable for larger sigma
#' species_params(params)$sigma <- 4
#' plot(w_full(params), pred_kernel(params)["Cod", 10, ], type="l", log="x")
lognormal_pred_kernel <- function(ppmr, beta, sigma) {
    Beta <- log(beta)
    phi <- exp(-(log(ppmr) - Beta)^2 / (2 * sigma^2))
    phi[ppmr < 1] <- 0
    return(phi)
}

#' Truncated lognormal predation kernel
#'
#' This is like the [lognormal_pred_kernel()] but with an imposed maximum
#' predator/prey mass ratio
#'
#' Writing the predator mass as \eqn{w} and the prey mass as \eqn{w_p},
#' the feeding kernel is given as
#' \deqn{\phi_i(w, w_p) =
#' \exp \left[ \frac{-(\ln(w / w_p / \beta_i))^2}{2\sigma_i^2} \right]
#' }{\phi_i(w/w_p) = exp(-(ln(w/w_p/\beta_i))^2/(2\sigma_i^2))}
#' if \eqn{w/w_p} is between 1 and
#' \eqn{\beta_i\exp(3\sigma_i)}{\beta_i exp(3\sigma_i)}
#' and zero otherwise. Here \eqn{\beta_i} is the preferred predator-prey mass
#' ratio and \eqn{\sigma_i} determines the width of the kernel. These two
#' parameters need to be given in the species parameter dataframe in the columns
#' `beta` and `sigma`.
#'
#' This function is called from [setPredKernel()] to set up the
#' predation kernel slots in a MizerParams object.
#'
#' @param ppmr A vector of predator/prey size ratios
#' @param beta The preferred predator/prey size ratio
#' @param sigma The width parameter of the log-normal kernel
#'
#' @return A vector giving the value of the predation kernel at each of the
#'   predator/prey mass ratios in the `ppmr` argument.
#' @export
#' @family predation kernel
#' @seealso [setPredKernel()]
#' @examples
#' params <- NS_params
#' species_params(params)$pred_kernel_type <- "truncated_lognormal"
#' plot(w_full(params), pred_kernel(params)["Cod", 10, ], type="l", log="x")
truncated_lognormal_pred_kernel <- function(ppmr, beta, sigma) {
    phi <- lognormal_pred_kernel(ppmr, beta, sigma)
    rr <- exp(log(beta) + 3 * sigma)
    phi[ppmr > rr] <- 0
    return(phi)
}

#' Gaussian-mixture predation kernel
#'
#' `r lifecycle::badge("experimental")`
#' A predation kernel for which the log predator/prey mass ratio follows a
#' mixture of Gaussian distributions.
#'
#' Writing the predator mass as \eqn{w}, the prey mass as \eqn{w_p}, and
#' \eqn{x = \ln(w / w_p)}, the feeding kernel is
#' \deqn{
#' \phi_i(w, w_p) = \sum_j a_{ij}
#' \exp\left[-\frac{(x - \mu_{ij})^2}{2\sigma_{ij}^2}\right],
#' \qquad
#' a_{ij} = \frac{p_{ij}/\sigma_{ij}}
#' {\sum_k p_{ik}/\sigma_{ik}}.
#' }
#' for predator/prey mass ratios greater than or equal to one, and zero for
#' smaller ratios.
#'
#' This is proportional to the Gaussian-mixture probability density with
#' mixing proportions \eqn{p_{ij}}, means \eqn{\mu_{ij}}, and standard
#' deviations \eqn{\sigma_{ij}}. The scaling makes the sum of the component
#' peak heights equal to one. Consequently the kernel is at most one, and a
#' one-component mixture is identical to [lognormal_pred_kernel()] with
#' `beta = exp(kernel_mean)` and `sigma = kernel_sd`.
#'
#' The three component parameters are vectors of equal length. When this
#' function is selected in a species parameter data frame, they should be held
#' in the list-columns `kernel_p`, `kernel_mean`, and `kernel_sd`. The values in
#' `kernel_p` must be non-negative with at least one positive value, but they do
#' not need to sum to one because they are normalised by the function.
#'
#' @param ppmr A vector of predator/prey mass ratios.
#' @param kernel_p A numeric vector of relative component proportions.
#' @param kernel_mean A numeric vector of component means on the log
#'   predator/prey mass-ratio scale.
#' @param kernel_sd A numeric vector of positive component standard deviations.
#'
#' @return A vector giving the value of the predation kernel at each of the
#'   predator/prey mass ratios in the `ppmr` argument.
#' @export
#' @family predation kernel
#' @seealso [setPredKernel()]
#' @examples
#' ppmr <- exp(seq(0, 12, length.out = 200))
#' phi <- gaussian_mixture_pred_kernel(
#'     ppmr,
#'     kernel_p = c(0.3, 0.7),
#'     kernel_mean = c(4, 8),
#'     kernel_sd = c(0.8, 1.5)
#' )
#' plot(ppmr, phi, type = "l", log = "x")
gaussian_mixture_pred_kernel <- function(ppmr, kernel_p,
                                         kernel_mean, kernel_sd) {
    if (!is.numeric(ppmr) || any(!is.finite(ppmr)) || any(ppmr <= 0)) {
        stop("`ppmr` must contain only positive finite values.")
    }
    if (!is.numeric(kernel_p) || !is.numeric(kernel_mean) ||
            !is.numeric(kernel_sd)) {
        stop("Gaussian-mixture parameters must be numeric vectors.")
    }
    no_components <- length(kernel_p)
    if (no_components == 0 || length(kernel_mean) != no_components ||
            length(kernel_sd) != no_components) {
        stop("`kernel_p`, `kernel_mean`, and `kernel_sd` must have the same ",
             "positive length.")
    }
    if (any(!is.finite(kernel_p)) || any(!is.finite(kernel_mean)) ||
            any(!is.finite(kernel_sd))) {
        stop("Gaussian-mixture parameters must contain only finite values.")
    }
    if (any(kernel_p < 0) || sum(kernel_p) <= 0) {
        stop("`kernel_p` must be non-negative with at least one positive ",
             "value.")
    }
    if (any(kernel_sd <= 0)) {
        stop("`kernel_sd` must contain only positive values.")
    }

    # A Gaussian density has peak height proportional to p / sd. Normalising
    # these component heights preserves the mixture-density shape, keeps the
    # kernel at or below one, and agrees with lognormal_pred_kernel() when
    # there is only one component.
    component_height <- kernel_p / sum(kernel_p) / kernel_sd
    component_height <- component_height / sum(component_height)
    phi <- numeric(length(ppmr))
    feeding <- ppmr >= 1
    log_ppmr <- log(ppmr[feeding])
    for (j in seq_len(no_components)) {
        phi[feeding] <- phi[feeding] + component_height[j] *
            exp(-(log_ppmr - kernel_mean[j])^2 / (2 * kernel_sd[j]^2))
    }
    phi
}

#' Box predation kernel
#'
#' A predation kernel where the predator/prey mass ratio is uniformly
#' distributed on an interval.
#'
#' Writing the predator mass as \eqn{w} and the prey mass as \eqn{w_p}, the
#' feeding kernel is 1 if \eqn{w/w_p} is between `ppmr_min` and
#' `ppmr_max` inclusive and zero otherwise. `ppmr_min` must be strictly smaller
#' than `ppmr_max`. The parameters need to be given in the species parameter
#' dataframe in the columns `ppmr_min` and `ppmr_max`.
#'
#' @param ppmr A vector of predator/prey size ratios
#' @param ppmr_min Minimum predator/prey mass ratio
#' @param ppmr_max Maximum predator/prey mass ratio
#'
#' @return A vector giving the value of the predation kernel at each of the
#'   predator/prey mass ratios in the `ppmr` argument.
#' @export
#' @family predation kernel
#' @seealso [setPredKernel()]
#' @examples
#' params <- NS_params
#' # Set all required paramters before changing kernel type
#' species_params(params)$ppmr_max <- 4000
#' species_params(params)$ppmr_min <- 200
#' species_params(params)$pred_kernel_type <- "box"
#' plot(w_full(params), pred_kernel(params)["Cod", 10, ], type="l", log="x")
box_pred_kernel <- function(ppmr, ppmr_min, ppmr_max) {
    assert_that(ppmr_min < ppmr_max)
    phi <- rep(1, length(ppmr))
    phi[ppmr > ppmr_max] <- 0
    phi[ppmr < ppmr_min] <- 0
    return(phi)
}

#' Power-law predation kernel
#'
#' This predation kernel is a power-law, with sigmoidal cut-offs at large and
#' small predator/prey mass ratios.
#'
#' The return value is calculated as
#'
#' \code{
#' ppmr^kernel_exp /
#'   (1 + (exp(kernel_l_l) / ppmr)^kernel_u_l) /
#'   (1 + (ppmr / exp(kernel_l_r))^kernel_u_r)
#' }
#'
#' The parameters need to be given as columns in the species parameter
#' dataframe.
#'
#' @param ppmr A vector of predator/prey size ratios at which to evaluate the
#'   predation kernel.
#' @param kernel_exp The exponent of the power law
#' @param kernel_l_l The location of the left, rising sigmoid
#' @param kernel_u_l The shape of the left, rising sigmoid
#' @param kernel_l_r The location of the right, falling sigmoid
#' @param kernel_u_r The shape of the right, falling sigmoid
#'
#' @return A vector giving the value of the predation kernel at each of the
#'   predator/prey mass ratios in the `ppmr` argument.
#' @export
#' @family predation kernel
#' @seealso [setPredKernel()]
#' @examples
#' params <- NS_params
#' # Set all required paramters before changing kernel type
#' species_params(params)["Cod", "kernel_exp"] <- -0.8
#' species_params(params)["Cod", "kernel_l_l"] <- 4.6
#' species_params(params)["Cod", "kernel_u_l"] <- 3
#' species_params(params)["Cod", "kernel_l_r"] <- 12.5
#' species_params(params)["Cod", "kernel_u_r"] <- 4.3
#' species_params(params)["Cod", "pred_kernel_type"] <- "power_law"
#' plot(w_full(params), pred_kernel(params)["Cod", 10, ], type="l", log="x")
power_law_pred_kernel <- function(ppmr, kernel_exp,
                                  kernel_l_l, kernel_u_l,
                                  kernel_l_r, kernel_u_r) {
    ppmr^kernel_exp /
        (1 + (exp(kernel_l_l) / ppmr)^kernel_u_l) /
        (1 + (ppmr / exp(kernel_l_r))^kernel_u_r)
}
