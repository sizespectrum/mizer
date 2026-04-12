#' Get diffusion rate from predation
#'
#' @description
#' Calculates the diffusion rate \eqn{D_i(w)} (grams^2/year) for each species.
#' This is the rate at which the abundance density is diffused along the
#' size axis due to the variability in prey sizes. This is the diffusion
#' term from the jump-growth equation.
#'
#' @param params A MizerParams object
#' @param n A matrix of species abundances (species x size). Defaults to the
#'   initial abundances in `params`.
#' @param n_pp A vector of the resource abundance by size. Defaults to the
#'   initial resource abundances in `params`.
#' @param n_other A list of abundances for other dynamical components.
#' @param t The time for which to do the calculation.
#' @param ... Unused
#'
#' @return An array of dimensions species x size holding the diffusion rates.
#' @export
#' @family rate functions
#' @references
#' Datta, S., Delius, G. W. and Law, R. (2010). A jump-growth model for
#' predator-prey dynamics: derivation and application to marine ecosystems.
#' Bulletin of Mathematical Biology, 72(6):1361–1382
getDiffusion <- function(params, n = initialN(params),
                         n_pp = initialNResource(params),
                         n_other = initialNOther(params),
                         t = 0,
                         ...) {
    UseMethod("getDiffusion")
}
#' @export
getDiffusion.MizerParams <- function(params, n = initialN(params),
                                     n_pp = initialNResource(params),
                                     n_other = initialNOther(params),
                                     t = 0,
                                     ...) {
    params <- validParams(params)
    f <- get(params@rates_funcs$Diffusion)
    d <- f(params, n = n, n_pp = n_pp, n_other = n_other, t = t,
           feeding_level = getFeedingLevel(params, n = n, n_pp = n_pp,
                                           n_other = n_other, t = t), ...)
    ArraySpeciesBySize(d, value_name = "Diffusion rate",
                       units = "g^2/year", params = params)
}


#' Calculate diffusion rate
#'
#' @description
#' Calculates the diffusion rate \eqn{D_i(w)} (grams^2/year) for each species.
#' This diffusion rate has two components: 
#' 1. Thethe diffusion due due to the variability in prey sizes. This is the
#'    diffusion term from the jump-growth equation.
#' 2. Any externally specified diffusion, which is added via [setExtDiffusion()]
#'
#' You would not usually call this function directly but instead use
#' [getDiffusion()], which then calls this function unless an alternative
#' diffusion rate function has been registered, see [setRateFunction()].
#'
#' @param params A MizerParams object
#' @param n A matrix of species abundances (species x size).
#' @param n_pp A vector of the resource abundance by size
#' @param n_other A list of abundances for other dynamical components
#' @param t The time for which to do the calculation (Not used by standard
#'   mizer rate functions but useful for extensions.)
#' @param feeding_level An array (species x size) with the feeding level.
#'   If not provided, it is calculated from the given abundances.
#' @param ... Unused
#'
#' @return A two dimensional array (species x size) holding the diffusion rate.
mizerDiffusion <- function(params, n, n_pp, n_other, t, feeding_level, ...) {

    if (missing(feeding_level)) {
        feeding_level <- getFeedingLevel(params, n = n, n_pp = n_pp,
                                         n_other = n_other, t = t)
    }

    # idx_sp are the indices into w_full that correspond to consumer sizes w
    idx_sp <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)

    if (isTRUE(params@use_predation_diffusion)) {
        # Calculate the diffusion integral
        # I_d(w) = sum_prey theta_i * N_prey(w_p) * w_p^2 * dw_p
        # This is the same convolution as in mizerEncounter but with w_p^2 * dw_p
        # weighting instead of w_p * dw_p.
        prey_sq <- outer(params@species_params$interaction_resource, n_pp)
        prey_sq[, idx_sp] <- prey_sq[, idx_sp] + params@interaction %*% n
        prey_sq <- sweep(prey_sq, 2, params@w_full^2 * params@dw_full, "*")

        # Convolve with the predation kernel via FFT.
        # mvfft() transforms each column, so we transpose to get row-wise FFTs,
        # following the same pattern as mizerEncounter.
        integral_d <- Re(base::t(mvfft(base::t(params@ft_pred_kernel_e) *
                                           mvfft(base::t(prey_sq)),
                                       inverse = TRUE))) / length(params@w_full)
        # Keep only the consumer sizes
        integral_d <- integral_d[, idx_sp, drop = FALSE]
        # Remove numerical noise
        integral_d[integral_d < 0] <- 0

        # D(w) = (1 - f(w)) * gamma(w) * alpha^2 * I_d(w)
        alpha <- params@species_params$alpha
        D <- (1 - feeding_level) * params@search_vol * alpha^2 * integral_d
        dimnames(D) <- dimnames(params@metab)
    } else {
        D <- matrix(0, nrow = nrow(feeding_level), ncol = ncol(feeding_level),
                    dimnames = dimnames(params@metab))
    }

    # Add any externally specified diffusion
    D <- D + params@ext_diffusion

    return(D)
}
