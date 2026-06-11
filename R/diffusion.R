#' Get diffusion rate from predation
#'
#' @description
#' Calculates the diffusion rate \eqn{D_i(w)} (grams^2/year) for each species.
#' This diffusion rate has two components:
#' 1. The diffusion due due to the variability in prey sizes. This is the
#'    diffusion term from the jump-growth equation.
#' 2. Any externally specified diffusion, which is added via [setExtDiffusion()]
#'
#' @details
#' The diffusion due due to the variability in prey sizes
#' is determined by summing over all prey
#' species and the resource spectrum and then integrating over all prey sizes
#' \eqn{w_p}, weighted by predation kernel \eqn{\phi(w,w_p)}:
#' \deqn{
#' d_i(w) = (1-f_i(w))(\alpha_i(1-\psi_i(w)))^2\gamma_i(w) \int
#' \left( \theta_{ip} N_R(w_p) + \sum_{j} \theta_{ij} N_j(w_p) \right)
#' \phi_i(w,w_p) w_p^2 \, dw_p.
#' }{(1-f_i(w))(\alpha_i(1-\psi_i(w)))^2\gamma_i(w) \int
#' ( \theta_{ip} N_R(w_p) + \sum_{j} \theta_{ij} N_j(w_p) )
#' \phi_i(w,w_p) w_p^2 dw_p.}
#' Here \eqn{N_j(w)} is the abundance density of species \eqn{j} and
#' \eqn{N_R(w)} is the abundance density of resource.
#' The overall prefactor \eqn{\gamma_i(w)} determines the predation power of the
#' predator. It could be interpreted as a search volume and is set with the
#' [setSearchVolume()] function. The predation kernel
#' \eqn{\phi(w,w_p)} is set with the [setPredKernel()] function. The
#' species interaction matrix \eqn{\theta_{ij}} is set with [setInteraction()]
#' and the resource interaction vector \eqn{\theta_{ip}} is taken from the
#' `interaction_resource` column in [species_params()].
#' \eqn{f(w)} is the feeding level calculated with
#' [getFeedingLevel()]. \eqn{\psi(w)} is the proportion of the available energy
#' that is invested in reproduction instead of growth, obtained with [psi()].
#'
#' @template param_object_dots
#'
#' @return
#' * `MizerParams`: An `ArraySpeciesBySize` object (predator species x predator
#'   size) with the diffusion rates.
#' * `MizerSim`: An `ArrayTimeBySpeciesBySize` object (time step x predator
#'   species x predator size) with the diffusion rates at every time step.
#'   If `drop = TRUE` then dimensions of length 1 will be removed.
#' @export
#' @family rate functions
#' @references
#' Datta, S., Delius, G. W. and Law, R. (2010). A jump-growth model for
#' predator-prey dynamics: derivation and application to marine ecosystems.
#' Bulletin of Mathematical Biology, 72(6):1361–1382
getDiffusion <- function(object, ...) {
    UseMethod("getDiffusion")
}
#' @export
getDiffusion.MizerParams <- function(object, n = initialN(object),
                                     n_pp = initialNResource(object),
                                     n_other = initialNOther(object),
                                     t = 0,
                                     ...) {
    params <- object
    params <- validParams(params)
    feeding_level <- getFeedingLevel(params, n = n, n_pp = n_pp,
                                     n_other = n_other, time_range = t)
    if (usesExtensionDispatch(params)) {
        d <- projectDiffusion(params, n = n, n_pp = n_pp, n_other = n_other,
                              t = t, feeding_level = feeding_level, ...)
    } else {
        f <- get(params@rates_funcs$Diffusion)
        d <- f(params, n = n, n_pp = n_pp, n_other = n_other, t = t,
               feeding_level = feeding_level, ...)
    }
    ArraySpeciesBySize(d, value_name = "Diffusion rate",
                       units = "g^2/year", params = params)
}

#' @export
getDiffusion.MizerSim <- function(object, n, n_pp, n_other, t = 0,
                                  time_range, drop = FALSE, ...) {
    sim <- object
    sim_size_rate(sim, time_range, drop, target = "Diffusion",
                  slot = "diffusion", value_name = "Diffusion rate",
                  units = "g^2/year", ...)
}

#' @name mizerDiffusion
#' @rdname mizerDiffusion
#' @export
projectDiffusion <- function(params, n, n_pp, n_other, t = 0,
                             feeding_level, ...) {
    UseMethod("projectDiffusion")
}

#' Calculate diffusion rate
#'
#' @description
#' Calculates the diffusion rate \eqn{D_i(w)} (grams^2/year) for each species.
#' This diffusion rate has two components:
#' 1. The diffusion due due to the variability in prey sizes. This is the
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
#' @rdname mizerDiffusion
#' @export
projectDiffusion.MizerParams <- function(params, n, n_pp, n_other, t = 0,
                                         feeding_level, ...) {

    if (missing(feeding_level)) {
        feeding_level <- getFeedingLevel(params, n = n, n_pp = n_pp,
                                         n_other = n_other, time_range = t)
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
        # following the same pattern as mizerEncounter. We use the dedicated
        # diffusion kernel `ft_pred_kernel_d`: the diffusion integrand carries
        # w_p^2 dw_p (one more power of prey size than the encounter's
        # w_p dw_p), so under `second_order_w` its bin-integral needs the e^{3t}
        # Jacobian rather than the encounter's e^{2t}. In the default first-order
        # scheme `ft_pred_kernel_d` equals `ft_pred_kernel_e`, so the result is
        # byte-identical to before.
        integral_d <- Re(base::t(mvfft(base::t(params@ft_pred_kernel_d) *
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

#' @rdname mizerDiffusion
#' @export
mizerDiffusion <- projectDiffusion.MizerParams


#' Get or set the use_predation_diffusion flag
#'
#' Controls whether predation-induced diffusion is included when calculating
#' rates with [mizerDiffusion()]. When `FALSE` (the default), the
#' predation-driven diffusion term is omitted, preserving the behaviour of
#' previous mizer versions. Set to `TRUE` to enable the diffusion term from
#' the jump-growth equation.
#'
#' @param params A MizerParams object.
#' @return `use_predation_diffusion()`: A single logical value.
#' @export
#' @family functions for setting parameters
use_predation_diffusion <- function(params) {
    params@use_predation_diffusion
}

#' @rdname use_predation_diffusion
#' @param value A single logical value (`TRUE` or `FALSE`).
#' @return `use_predation_diffusion<-`: A MizerParams object with the
#'   `use_predation_diffusion` flag updated.
#' @export
`use_predation_diffusion<-` <- function(params, value) {
    assert_that(is.flag(value))
    params@use_predation_diffusion <- value
    params@time_modified <- lubridate::now()
    params
}
