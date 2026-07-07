#' Set predation kernel
#'
#' You will usually not need to call this function directly. Instead change
#' the relevant species parameters (`pred_kernel_type`, and, depending on its
#' value, `beta`/`sigma` or `ppmr_min`/`ppmr_max`) with
#' `given_species_params(params) <-` and let mizer recalculate the predation
#' kernel for you. Call `setPredKernel()` directly only if you want to supply
#' the full kernel array yourself. See
#' `vignette("cheatsheet-changing-parameters")` for a full explanation of when
#' to reach for which level of the model.
#'
#' The predation kernel determines the distribution of prey sizes that a
#' predator feeds on. It is used in [getEncounter()] when calculating
#' the rate at which food is encountered and in [getPredRate()] when
#' calculating the rate at which a prey is predated upon. The predation kernel
#' can be a function of the predator/prey size ratio or it can be a function of
#' the predator size and the prey size separately. Both types can be set up with
#' this function.
#'
#' @section Setting predation kernel:
#' \strong{Kernel dependent on predator to prey size ratio}
#'
#' If the `pred_kernel` argument is not supplied, then this function sets a
#' predation kernel that depends only on the ratio of predator mass to prey
#' mass, not on the two masses independently. The shape of that kernel is then
#' determined by the `pred_kernel_type` column in species_params.
#'
#' The default for `pred_kernel_type` is "lognormal". This will call the function
#' [lognormal_pred_kernel()] to calculate the predation kernel.
#' An alternative pred_kernel type is "box", implemented by the function
#' [box_pred_kernel()], and "power_law", implemented by the function
#' [power_law_pred_kernel()]. These functions require certain species
#' parameters in the species_params data frame. For the lognormal kernel these
#' are `beta` and `sigma`, for the box kernel they are `ppmr_min`
#' and `ppmr_max`. They are explained in the help pages for the kernel
#' functions. Except for `beta` and `sigma`, no defaults are set for
#' these parameters. If they are missing from the species_params data frame then
#' mizer will issue an error message.
#'
#' You can use any other string for `pred_kernel_type`. If for example you
#' choose "my" then you need to define a function `my_pred_kernel` that you can
#' model on the existing functions like [lognormal_pred_kernel()].
#'
#' When using a kernel that depends on the predator/prey size ratio only, mizer
#' does not need to store the entire three dimensional array in the MizerParams
#' object. Such an array can be very big when there is a large number of size
#' bins. Instead, mizer only needs to store two two-dimensional arrays that hold
#' Fourier transforms of the feeding kernel function that allow the encounter
#' rate and the predation rate to be calculated very efficiently. However, if
#' you need the full three-dimensional array you can calculate it with the
#' [getPredKernel()] function.
#'
#' \strong{Kernel dependent on both predator and prey size}
#'
#' If you want to work with a feeding kernel that depends on predator mass and
#' prey mass independently, you can specify the full feeding kernel as a
#' three-dimensional array (predator species x predator size x prey size).
#'
#' You should use this option only if a kernel dependent only on the
#' predator/prey mass ratio is not appropriate. Using a kernel dependent on
#' predator/prey mass ratio only allows mizer to use fast Fourier transform
#' methods to significantly reduce the running time of simulations.
#'
#' The order of the predator species in `pred_kernel` should be the same
#' as the order in the species params dataframe in the `params` object. If you
#' supply a named array then the function will check the order and warn if it is
#' different.
#'
#' @section Higher-order quadrature:
#' When the predation kernel depends only on the predator/prey mass ratio, the
#' encounter and predation rates are evaluated as convolutions using the fast
#' Fourier transform. By default mizer point-samples the kernel at the grid
#' nodes, which is a first-order (rectangle-rule) quadrature. When the
#' `bin_average` entry of the `second_order_w` slot is `TRUE`, the kernel is
#' instead integrated over each logarithmic size bin when building the
#' Fourier-transformed kernels. This finite-volume consistent quadrature lifts
#' the encounter and predation rates towards second order at no extra runtime
#' cost, because the integration is performed once here and the rate functions
#' themselves are unchanged. The predation kernel is additionally averaged over
#' the prey bin (a trapezoid fold), so that the predation rate returned by
#' [getPredRate()] is the prey-bin average that the predation-mortality sink
#' needs to be second order. The default remains the first-order scheme so that
#' existing models are unaffected. Enable it with `second_order_w(params) <-
#' TRUE` (see [second_order_w()]), which also turns on the other bin-averaged
#' rate quadratures so the whole model stays consistent. See the
#' `vignette("fft")` for the mathematical details.
#'
#' @param params A MizerParams object
#' @param pred_kernel Optional. An array (species x predator size x prey size)
#'   that holds the predation coefficient of each predator at size on each prey
#'   size. If not supplied, a default is set as described in section "Setting
#'   predation kernel".
#' @param reset
#'   If set to TRUE, then the predation kernel will be reset to the
#'   value calculated from the species parameters, even if it was previously
#'   overwritten with a custom value. If set to FALSE (default) then a
#'   recalculation from the species parameters will take place only if no custom
#'   value has been set.
#' @param ... Unused
#'
#' @return `setPredKernel()`: A MizerParams object with updated predation kernel.
#' @export
#' @family functions for setting parameters
#' @examples
#' ## Set up a MizerParams object
#' params <-  NS_params
#'
#' ## If you change predation kernel parameters after setting up a model,
#' #  this will be used to recalculate the kernel
#' species_params(params)["Cod", "beta"] <- 200
#'
#' ## You can change to a different predation kernel type
#' species_params(params)$ppmr_max <- 4000
#' species_params(params)$ppmr_min <- 200
#' species_params(params)$pred_kernel_type <- "box"
#' plot(w_full(params), getPredKernel(params)["Cod", 100, ], type="l", log="x")
#'
#' ## If you need a kernel that depends also on prey size you need to define
#' # it yourself.
#' pred_kernel <- getPredKernel(params)
#' pred_kernel["Herring", , ] <- sweep(pred_kernel["Herring", , ], 2,
#'                                     params@w_full, "*")
#' params<- setPredKernel(params, pred_kernel = pred_kernel)
setPredKernel <- function(params, pred_kernel = NULL, reset = FALSE, ...) {
    UseMethod("setPredKernel")
}
#' @export
setPredKernel.MizerParams <- function(params,
                          pred_kernel = NULL,
                          reset = FALSE, ...) {
    assert_that(is.flag(reset))

    # Which quadrature to use for the Fourier-transformed kernels is controlled
    # by the `bin_average` entry of the `second_order_w` slot, the single switch
    # for all the bin-averaged rate quadratures. Enable it with
    # `second_order_w(params) <- TRUE`, which re-runs setParams() so that the
    # kernels and the other bin-averaged quantities stay mutually consistent.
    high_order <- isTRUE(params@second_order_w[["bin_average"]])

    if (reset) {
        if (!is.null(pred_kernel)) {
            warning("Because you set `reset = TRUE`, the value you provided ",
                    "for `pred_kernel` will be ignored and a value will be ",
                    "calculated from the species parameters.")
            pred_kernel <- NULL
        }
        comment(params@pred_kernel) <- NULL
    }

    if (!is.null(pred_kernel)) {
        if (is.null(comment(pred_kernel))) {
            if (is.null(comment(params@pred_kernel))) {
                comment(pred_kernel) <- "set manually"
            } else {
                comment(pred_kernel) <- comment(params@pred_kernel)
            }
        }
        # A pred kernel was supplied, so check it and store it
        assert_that(is.array(pred_kernel))
        # psi is used in the next line just because it has the right dimension
        assert_that(identical(dim(pred_kernel),
                              c(dim(params@psi), length(params@w_full))))
        if (!is.null(dimnames(pred_kernel)) &&
            !all(dimnames(pred_kernel)[[1]] == params@species_params$species)) {
            stop(paste0("You need to use the same ordering of species as in the ",
                        "params object: ", toString(params@species_params$species)))
        }
        assert_that(all(pred_kernel >= 0))
        dimnames(pred_kernel) <-
            list(sp = params@species_params$species,
                 w_pred = signif(params@w, 3),
                 w_prey = signif(params@w_full, 3))
        params@pred_kernel <- pred_kernel
        params@time_modified <- lubridate::now()
        return(params)
    }

    ## Set a pred kernel dependent on predator/prey size ratio only

    # If pred_kernel_type is not supplied use "lognormal"
    params <- default_pred_kernel_params(params)

    species_params <- params@species_params
    pred_kernel_type <- species_params$pred_kernel_type
    no_sp <- nrow(species_params)
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)
    ft_pred_kernel_e <-
        array(NA, dim = c(no_sp, no_w_full),
              dimnames = list(sp = species_params$species, k = 1:no_w_full))
    ft_pred_kernel_p <- ft_pred_kernel_e
    ft_pred_kernel_d <- ft_pred_kernel_e
    # The kernel weights phi_e[i, ], phi_p[i, ] and phi_d[i, ] give, for each
    # offset m = 0, ..., no_w_full - 1 between the predator and prey bins, the
    # value that enters the encounter, predation and predation-diffusion
    # convolutions respectively. With the default first-order scheme all three
    # are the kernel point-sampled at the grid node, tilde_phi(beta^m). With the
    # higher-order scheme they are the kernel integrated over the prey bin
    # (encounter and diffusion) and predator bin (predation); see below and the
    # vignette("fft"). The diffusion integrand carries one more power of prey
    # size than the encounter (w_p^2 dw_p vs w_p dw_p), so it uses the same
    # prey-bin ratios but the e^{3t} Jacobian instead of e^{2t}.
    if (!high_order) {
        # First-order (rectangle-rule) quadrature: point-sample the kernel.
        # The smallest predator/prey mass ratio is 1.
        ppmr <- params@w_full / params@w_full[1]
        phi_e <- get_phi(species_params, ppmr)
        # Do not allow feeding at own size
        phi_e[, 1] <- 0
        phi_p <- phi_e
        # The diffusion convolution reuses the point-sampled kernel; the extra
        # power of prey size lives entirely in the prey vector in this scheme.
        phi_d <- phi_e
    } else {
        # Higher-order (bin-integrated) quadrature. The grid is geometric, so
        # the bin ratio beta = w_full[k+1] / w_full[k] is constant; Delta is the
        # natural-log bin width. We integrate the kernel over each log-bin by a
        # composite-midpoint rule with Q sub-cells, parametrised by
        # t = s * Delta with s in [0, 1].
        beta_grid <- params@w_full[2] / params@w_full[1]
        Delta <- log(beta_grid)
        m <- 0:(no_w_full - 1)
        Q <- 100L
        tt <- (seq_len(Q) - 0.5) / Q * Delta
        # Ratios at which to evaluate the kernel. Encounter integrates over the
        # prey bin, so the ratio w/w_p = beta^m e^{-t} decreases across the bin;
        # predation integrates over the predator bin, so w/w_p = beta^m e^{t}
        # increases across the bin.
        ppmr_e <- exp(outer(m * Delta, tt, `-`))
        ppmr_p <- exp(outer(m * Delta, tt, `+`))
        kernel_e <- get_phi(species_params, as.vector(ppmr_e))
        kernel_p <- get_phi(species_params, as.vector(ppmr_p))
        # Quadrature weights including the bin-integral Jacobian. The prey
        # integrand carries w_p dw_p (two powers of w_p -> e^{2t}); the predator
        # integrand carries dw (one power of w -> e^{t}). Dividing by (beta - 1)
        # cancels the factor w * dw = (beta - 1) w^2 that mizerEncounter() and
        # mizerPredRate() already fold into the prey and predator vectors, so
        # those rate functions need no change.
        # The diffusion integrand is the encounter integrand with one more power
        # of prey size (w_p^2 dw_p instead of w_p dw_p), so it integrates over
        # the same prey bin (ppmr_e ratios) but with the e^{3t} Jacobian.
        weight_e <- exp(2 * tt) * Delta / Q / (beta_grid - 1)
        weight_p <- exp(tt)     * Delta / Q / (beta_grid - 1)
        weight_d <- exp(3 * tt) * Delta / Q / (beta_grid - 1)
        phi_e <- matrix(0, nrow = no_sp, ncol = no_w_full)
        phi_p <- phi_e
        phi_d <- phi_e
        for (i in 1:no_sp) {
            phi_e[i, ] <- matrix(kernel_e[i, ], nrow = no_w_full) %*% weight_e
            phi_p[i, ] <- matrix(kernel_p[i, ], nrow = no_w_full) %*% weight_p
            phi_d[i, ] <- matrix(kernel_e[i, ], nrow = no_w_full) %*% weight_d
        }
        # Prey-bin average of the predation kernel (issue #381). The predation
        # convolution outputs at the prey node, but the predation-mortality sink
        # integrates that rate against the prey density over the prey bin, so it
        # wants the prey-bin average (1/Delta w_p) int pred_rate(w_p) dw_p, not a
        # point value. Bin-averaging the output over the prey bin is a trapezoid
        # fold of the kernel over adjacent offsets: with pred_rate at prey node j
        # using kernel offset m = (predator index) - j, the cell average
        # 1/2 (pred_rate[j] + pred_rate[j+1]) replaces the offset-m weight by
        # 1/2 (phi_p(beta^m) + phi_p(beta^{m-1})). This costs nothing at runtime
        # and completes #374's predator-bin integral with the prey-bin integral,
        # so the encounter integrates the prey bins and the predation integrates
        # both predator and prey bins. The lowest offset is one-sided (the
        # own-size term is dropped from the convolution anyway). Encounter
        # (phi_e) is left untouched: it feeds the growth flux, a point/face
        # quantity, and is bin-averaged only at the reproduction integral.
        phi_p[, -1] <- 0.5 * (phi_p[, -1, drop = FALSE] +
                                  phi_p[, -no_w_full, drop = FALSE])
        # Predator-bin average of the diffusion kernel. The diffusion
        # convolution outputs at the predator node, and the transport solver
        # needs the cell-averaged diffusion coefficient over the predator bin,
        # so it wants the predator-bin average. Bin-averaging the output over
        # the predator bin is a trapezoid fold of the kernel over adjacent offsets:
        # with diffusion at predator node j using kernel offset m = j - (prey index),
        # the cell average 1/2 (integral_d[j] + integral_d[j+1]) replaces the
        # offset-m weight by 1/2 (phi_d(beta^m) + phi_d(beta^{m+1})). The largest
        # offset is one-sided.
        phi_d[, -no_w_full] <- 0.5 * (phi_d[, -no_w_full, drop = FALSE] +
                                          phi_d[, -1, drop = FALSE])
    }

    for (i in 1:no_sp) {
        # Fourier transform of feeding kernel for evaluating available energy
        ft_pred_kernel_e[i, ] <- fft(phi_e[i, ])
        # Fourier transform of feeding kernel for evaluating the predation
        # diffusion integral (same convolution as encounter, extra power of w_p)
        ft_pred_kernel_d[i, ] <- fft(phi_d[i, ])
        # Fourier transform of feeding kernel for evaluating predation rate
        ri <- min(max(which(phi_p[i, ] > 0)), no_w_full - 1)  # index of largest ppmr
        phi_p_rev <- rep(0, no_w_full)
        phi_p_rev[(no_w_full - ri + 1):no_w_full] <- phi_p[i, (ri + 1):2]
        ft_pred_kernel_p[i, ] <- fft(phi_p_rev)
    }

    # Prevent resetting if full slot has been commented
    if (!is.null(comment(params@pred_kernel))) {
        # Issue warning but only if a change was actually requested
        if (different(ft_pred_kernel_e, params@ft_pred_kernel_e) ||
            different(ft_pred_kernel_p, params@ft_pred_kernel_p) ||
            different(ft_pred_kernel_d, params@ft_pred_kernel_d)) {
            message("You have set a custom predation kernel and so it ",
                    "will not be recalculated from the species parameters ",
                    "unless you set `reset = TRUE`.")
        }
        return(params)
    }
    params@ft_pred_kernel_e[] <- ft_pred_kernel_e
    params@ft_pred_kernel_p[] <- ft_pred_kernel_p
    params@ft_pred_kernel_d[] <- ft_pred_kernel_d

    params@time_modified <- lubridate::now()
    return(params)
}

#' @rdname setPredKernel
#' @return `getPredKernel()` or equivalently `pred_kernel()`: An array (predator
#'   species x predator_size x prey_size)
#' @export
getPredKernel <- function(params) {
    UseMethod("getPredKernel")
}
#' @export
getPredKernel.MizerParams <- function(params) {
    # This function is more complicated than you might have thought because
    # usually the predation kernel is not stored in the MizerParams object,
    # but rather only the Fourier coefficients needed for fast calculation of
    # the convolution integrals.
    if (length(dim(params@pred_kernel)) > 1) {
        return(params@pred_kernel)
    }
    species_params <- default_pred_kernel_params(params@species_params)
    pred_kernel_type <- species_params$pred_kernel_type
    no_sp <- nrow(species_params)
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)
    # Vector of predator/prey mass ratios
    # The smallest predator/prey mass ratio is 1
    ppmr <- params@w_full / params@w_full[1]
    phis <- get_phi(species_params, ppmr)
    # Do not allow feeding at own size
    phis[, 1] <- 0
    pred_kernel <-
        array(0,
              dim = c(no_sp, no_w, no_w_full),
              dimnames = list(sp = species_params$species,
                              w_pred = signif(params@w, 3),
                              w_prey = signif(params@w_full, 3)))
    for (i in 1:no_sp) {
        min_w_idx <- no_w_full - no_w + 1
        for (k in seq_len(no_w)) {
            pred_kernel[i, k, (min_w_idx - 1 + k):1] <-
                phis[i, 1:(min_w_idx - 1 + k)]
        }
    }
    return(pred_kernel)
}

#' @rdname setPredKernel
#' @export
pred_kernel <- function(params) {
    getPredKernel(params)
}

#' @rdname setPredKernel
#' @param value pred_kernel
#' @export
`pred_kernel<-` <- function(params, value) {
    setPredKernel(params, pred_kernel = value)
}

#' Set defaults for predation kernel parameters
#'
#' If the predation kernel type has not been specified for a species, then it
#' is set to "lognormal" and the default values are set for the parameters
#' `beta` and `sigma`.
#' @param object Either a MizerParams object or a species parameter data frame
#' @return The `object` with updated columns in the species params data frame.
#' @export
#' @concept helper
default_pred_kernel_params <- function(object) {
    if (is(object, "MizerParams")) {
        # Nothing to do if full pred kernel has been specified
        if (length(dim(object@pred_kernel)) > 1) {
            return(object)
        }
        species_params <- object@species_params
    } else {
        species_params <- object
    }

    species_params <- set_species_param_default(species_params,
                                                "pred_kernel_type",
                                                "lognormal")
    # For species where the pred_kernel_type is lognormal, set defaults for
    # sigma and beta if none are supplied
    if (any(species_params$pred_kernel_type == "lognormal")) {
        species_params <- set_species_param_default(species_params,
                                                    "beta", 30)
        species_params <- set_species_param_default(species_params,
                                                    "sigma", 2)
    }
    if  (is(object, "MizerParams")) {
        object@species_params <- species_params
        return(object)
    } else {
        return(species_params)
    }
}

#' Get values from feeding kernel function
#'
#' This involves finding the feeding kernel function for each species, using the
#' pred_kernel_type parameter in the species_params data frame, checking that it
#' is valid and all its arguments are contained in the species_params data
#' frame, and then calling this function with the ppmr vector.
#'
#' @param species_params A species parameter data frame
#' @param ppmr Values of the predator/prey mass ratio at which to evaluate the
#'   predation kernel function
#' @return An array (species x ppmr) with the values of the predation kernel
#'   function
#' @export
#' @concept helper
get_phi <- function(species_params, ppmr) {
    assert_that(is.data.frame(species_params))
    no_sp <- nrow(species_params)
    species_params <- default_pred_kernel_params(species_params)
    phis <- array(dim = c(no_sp, length(ppmr)))
    for (i in 1:no_sp) {
        pred_kernel_func_name <- paste0(species_params$pred_kernel_type[i],
                                        "_pred_kernel")
        pred_kernel_func <- get0(pred_kernel_func_name)
        assert_that(is.function(pred_kernel_func))
        args <- names(formals(pred_kernel_func))
        if (!("ppmr" %in% args)) {
            stop("The predation kernel function ",
                 pred_kernel_func_name,
                 "needs the argument 'ppmr'.")
        }
        # lop off the compulsory arg
        args <- args[!(args %in% "ppmr")]
        missing <- !(args %in% colnames(species_params))
        if (any(missing)) {
            stop("The following arguments for the predation kernel function ",
                 pred_kernel_func_name,
                 " are missing from the parameter dataframe: ",
                 toString(args[missing]))
        }
        if (any(is.na(species_params[i, args]))) {
            stop("For species ",
                 species_params$species[i],
                 " the following arguments for the predation kernel function ",
                 pred_kernel_func_name,
                 " are NA in the parameter dataframe: ",
                 toString(args[is.na(species_params[i, args])]))
        }
        pars <- c(ppmr = list(ppmr), as.list(species_params[i, args]))
        phi <- do.call(pred_kernel_func_name, args = pars)

        if (any(is.na(phi))) {
            stop("The function ", pred_kernel_func_name,
                 " returned NA. Did you correctly specify all required",
                 " parameters in the species parameter dataframe?")
        }
        if (any(phi < 0)) {
            stop("The function ", pred_kernel_func_name,
                 " returned negative values.")
        }
        if (all(phi == 0)) {
            stop("The function ", pred_kernel_func_name,
                 " returned a zero predation kernel.")
        }
        phis[i, ] <- phi
    }
    return(phis)
}
