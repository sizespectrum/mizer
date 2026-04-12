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
    
    # The integral of prey biomass squared weighted by predation kernel
    # This is a convolution in log-space.
    # We calculate it in Fourier space.
    
    # Extend n to full grid
    n_full <- matrix(0, nrow = nrow(n), ncol = length(params@w_full))
    n_full[, (length(params@w_full) - length(params@w) + 1):length(params@w_full)] <- n
    
    # Prey abundance is a combination of resource and consumers
    prey <- rbind(n_pp, n_full)
    
    # We want to calculate a convolution with prey biomass squared density
    # prey_density * w_p^2
    prey_sq_biomass <- sweep(prey, 2, params@w_full^2, "*")
    prey_ft <- mvfft(prey_sq_biomass)
    
    # The fft of the predation kernel
    kernel_ft <- params@ft_pred_kernel_e
    
    # The interaction matrix is predator x prey. The first prey is the resource.
    interaction <- cbind(params@species_params$interaction_resource,
                         params@interaction)
    
    # The convolution theorem says that the Fourier transform of a
    # convolution is the product of the Fourier transforms.
    integral_d_ft <- kernel_ft * (interaction %*% prey_ft)
    
    # Inverse FFT
    integral_d <- Re(mvfft(integral_d_ft, inverse = TRUE)) / length(params@w_full)
    
    # We are only interested in the values for the consumer sizes
    # and we want the result as predator x size
    integral_d <- integral_d[, (length(params@w_full) - length(params@w) + 1):length(params@w_full),
                             drop = FALSE]
    
    # Get assimilation efficiency
    alpha <- params@species_params$alpha

    # Calculate diffusion rate
    # D(w) = (1-f(w)) * gamma(w) * alpha^2 * I_d(w)
    D <- (1 - feeding_level) * params@search_vol *
        alpha^2 * integral_d
    
    dimnames(D) <- dimnames(params@metab)

    # Add any externally specified diffusion
    D <- D + params@ext_diffusion

    return(D)
}
