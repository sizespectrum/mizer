#' Set external diffusion rate
#'
#' You will usually not need to call this function directly. Instead change
#' the `D_ext` and `n` species parameters with
#' `given_species_params(params) <-` and let mizer recalculate the external
#' diffusion rate for you. Call `setExtDiffusion()` directly only if you want
#' to impose a different functional form for the size dependence of the
#' external diffusion rate. See `vignette("cheatsheet-changing-parameters")`
#' for a full explanation of when to reach for which level of the model.
#'
#' @section Setting external diffusion rate:
#' The external diffusion rate allows you to impose additional diffusion
#' beyond the predation-driven diffusion that can be internally modelled by mizer.
#'
#' The `ext_diffusion` argument allows you to specify a diffusion rate that
#' depends on species and body size.
#'
#' If the `ext_diffusion` argument is not supplied, then the external diffusion
#' rate is calculated as a power law:
#' \deqn{D_{ext.i}(w) = D_{ext.i}\, w^{n_i+1}.}{
#'   D_{ext.i}(w) = D_{ext.i} * w^(n_i+1).}
#' The coefficient \eqn{D_{ext.i}} is taken from the `D_ext` column of the
#' species parameter data frame, which defaults to 0. The exponent
#' \eqn{n_i + 1} uses the `n` column of the species parameter data frame.
#'
#' If the `ext_diffusion` slot has a comment and `reset = FALSE`, then a
#' recalculation from the species parameters is suppressed and a message is
#' issued if the recalculated values would differ from the stored ones.
#'
#' @param params MizerParams
#' @param ext_diffusion Optional. An array (species x size) holding the
#'   external diffusion rate. If not supplied, a default is calculated from the
#'   `D_ext` and `n` species parameters as described in the section "Setting
#'   external diffusion rate".
#' @param reset
#'   If set to TRUE, then the external diffusion rate will be reset to the
#'   value calculated from the species parameters, even if it was previously
#'   overwritten with a custom value. If set to FALSE (default) then a
#'   recalculation from the species parameters will take place only if no
#'   custom value has been set.
#' @param ... Unused
#'
#' @return `setExtDiffusion()`: A MizerParams object with updated external
#'   diffusion rate.
#' @export
#' @family functions for setting parameters
setExtDiffusion <- function(params, ext_diffusion = NULL, reset = FALSE, ...) {
    UseMethod("setExtDiffusion")
}
#' @rdname setExtDiffusion
#' @usage NULL
#' @export
setExtDiffusion.MizerParams <- function(params, ext_diffusion = NULL,
                                        reset = FALSE, ...) {
    assert_that(is.flag(reset))

    if (reset) {
        if (!is.null(ext_diffusion)) {
            warning("Because you set `reset = TRUE`, the value you provided ",
                    "for `ext_diffusion` will be ignored and a value will be ",
                    "calculated from the species parameters.")
            ext_diffusion <- NULL
        }
        comment(params@ext_diffusion) <- NULL
    }

    # If ext_diffusion array is supplied, check it, store it and return
    if (!is.null(ext_diffusion)) {
        if (is.null(comment(ext_diffusion))) {
            if (is.null(comment(params@ext_diffusion))) {
                comment(ext_diffusion) <- "set manually"
            } else {
                comment(ext_diffusion) <- comment(params@ext_diffusion)
            }
        }
        assert_that(is.array(ext_diffusion),
                    identical(dim(ext_diffusion), dim(params@ext_diffusion)))
        params@ext_diffusion[] <- ext_diffusion
        comment(params@ext_diffusion) <- comment(ext_diffusion)

        params@time_modified <- lubridate::now()
        return(params)
    }

    # Else recalculate from species params
    params <- set_species_param_default(params, "D_ext", 0)

    n_exp <- params@species_params[["n"]] + 1
    if (isTRUE(params@second_order_w[["bin_average"]])) {
        # Use the exact bin average of the power law D_ext * w^(n+1) over each
        # bin instead of point-sampling at the left bin boundary, exactly as
        # setExtMort() does for the external mortality. This makes the diffusion
        # coefficient the bin average d_j that the second-order transport flux
        # needs. The point-vs-average choice is gated on bin_average, like the
        # other rates; flux only governs the advective reconstruction.
        ext_diffusion <- matrix(0, nrow = nrow(params@species_params),
                                ncol = length(params@w))
        for (i in seq_len(nrow(params@species_params))) {
            ext_diffusion[i, ] <- params@species_params[["D_ext"]][i] *
                power_law_bin_average(params@w, params@dw, n_exp[i])
        }
    } else {
        ext_diffusion <- sweep(outer(n_exp, params@w, function(x, y) y^x),
                               1, params@species_params[["D_ext"]], "*")
    }

    # Prevent overwriting slot if it has been commented
    if (!is.null(comment(params@ext_diffusion))) {
        if (different(ext_diffusion, params@ext_diffusion)) {
            message("The external diffusion rate has been commented and ",
                    "therefore will not be recalculated from the species ",
                    "parameters.")
        }
        return(params)
    }
    params@ext_diffusion[] <- ext_diffusion

    params@time_modified <- lubridate::now()
    return(params)
}

#' @rdname setExtDiffusion
#' @return `ext_diffusion()`: An `ArraySpeciesBySize` object (species x size)
#'   with the external diffusion rate.
#' @export
ext_diffusion <- function(params) {
    ArraySpeciesBySize(params@ext_diffusion,
                       value_name = "Diffusion rate",
                       params = params)
}

#' @rdname setExtDiffusion
#' @param value ext_diffusion
#' @export
`ext_diffusion<-` <- function(params, value) {
    setExtDiffusion(params, ext_diffusion = value)
}
