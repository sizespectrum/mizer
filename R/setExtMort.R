#' Set external mortality rate
#'
#' You will usually not need to call this function directly. Instead change
#' the `z0`, `z_ext` and `d` species parameters with
#' `given_species_params(params) <-` and let mizer recalculate the external
#' mortality rate for you. Call `setExtMort()` directly only if you want to
#' impose a different functional form for the size dependence of the external
#' mortality. See `vignette("cheatsheet-changing-parameters")` for a full
#' explanation of when to reach for which level of the model.
#'
#' @section Setting external mortality rate:
#' The external mortality is all the mortality that is not due to fishing or
#' predation by predators included in the model. The external mortality could be
#' due to predation by predators that are not explicitly included in the model
#' (e.g. mammals or seabirds) or due to other causes like illness. It is a rate
#' with units 1/year.
#'
#' The `ext_mort` argument allows you to specify an external mortality rate
#' that depends on species and body size. You can see an example of this in
#' the Examples section of the help page for [setExtMort()].
#'
#' If the `ext_mort` argument is not supplied, then the external mortality is
#' taken from the species parameters as
#' \deqn{\mu_{ext.i}(w) = z_{0.i} + z_{ext.i} w^{d_i}.}{mu_ext,i(w) = z0_i + z_ext,i w ^ d_i.}
#' The value of the constant \eqn{z_0} for each species is taken from the `z0`
#' column of [given_species_params()] if it is present there. Otherwise it is
#' recalculated, even if a value from an earlier calculation is still present
#' in `species_params`, as
#' \deqn{z_{0.i} = {\tt z0pre}_i\, w_{inf}^{\tt z0exp}.}{z_{0.i} = z0pre_i w_{inf}^{z0exp}.}
#' When `z0pre` or `z0exp` is supplied explicitly and used to calculate
#' non-given `z0`, the resulting values are recorded in
#' [given_species_params()]. Values calculated from the defaults
#' `z0pre = 0.6` and `z0exp = n - 1` are not recorded there. If either argument
#' is supplied but cannot be used because `z0` is given for every species or
#' because `ext_mort` was supplied, a warning is issued.
#' Missing values of `z_ext` are set to 0 and missing values of `d` are set to
#' `n - 1`.
#'
#' By default the power law is evaluated at the left bin edges \eqn{w_j}
#' (point sampling). If the `bin_average` entry of the `second_order_w` slot is
#' `TRUE` (see [second_order_w()]), then the \eqn{z_{ext} w^d} term is instead
#' replaced by its exact average over each bin \eqn{[w_j, w_{j+1}]},
#' \deqn{\frac{z_{ext}}{\Delta w_j}\int_{w_j}^{w_{j+1}} w^d\, dw
#'   = z_{ext}\,\frac{w_{j+1}^{d+1} - w_j^{d+1}}{(d+1)\,\Delta w_j},}
#' (with the limiting form \eqn{z_{ext}\ln(w_{j+1}/w_j)/\Delta w_j} when
#' \eqn{d = -1}). This is the consistent choice in the finite-volume scheme,
#' where the external mortality multiplies the bin-averaged abundance. The
#' bin-averaging is applied only to the auto-calculated power-law default; a
#' user-supplied `ext_mort` array is left untouched.
#'
#' @param params MizerParams
#' @param ext_mort Optional. An array (species x size) holding the external
#'   mortality rate.  If not supplied, a default is set as described in the
#'   section "Setting external mortality rate".
#' @param reset
#'   If set to TRUE, then the external mortality rate will be reset
#'   to the value calculated from the species parameters, even if it was
#'   previously overwritten with a custom value. If set to FALSE (default) then
#'   a recalculation from the species parameters will take place only if no
#'   custom value has been set.
#' @param z0pre If `z0`, the mortality from other sources, is not present in
#'   [given_species_params()], it is calculated as
#'   `z0pre * w_inf ^ z0exp`. Default value is 0.6.
#' @param z0exp The exponent used with `z0pre` to calculate non-given `z0`.
#'   Default value is `n - 1`.
#' @param z0 `r lifecycle::badge("deprecated")` Use `ext_mort` instead. Not to
#'   be confused with the species_parameter `z0`.
#' @param ... Unused
#'
#' @return `setExtMort()`: A MizerParams object with updated external mortality
#'   rate.
#' @export
#' @family functions for setting parameters
#' @examples
#' params <- newMultispeciesParams(NS_species_params)
#'
#' #### Setting allometric death rate #######################
#'
#' # Set coefficient for each species. Here we choose 0.1 for each species
#' z0pre <- rep(0.1, nrow(species_params(params)))
#'
#' # Multiply by power of size with exponent, here chosen to be -1/4
#' # The outer() function makes it an array species x size
#' allo_mort <- outer(z0pre, w(params)^(-1/4))
#'
#' # Change the external mortality rate in the params object
#' ext_mort(params) <- allo_mort
setExtMort <- function(params, ext_mort = NULL, z0pre = 0.6,
                       z0exp = params@resource_params$n - 1, reset = FALSE,
                       z0 = deprecated(), ...) {
    UseMethod("setExtMort")
}
#' @export
setExtMort.MizerParams <- function(params, ext_mort = NULL,
                                   z0pre = 0.6,
                                   z0exp = params@resource_params$n - 1,
                                   reset = FALSE,
                                   z0 = deprecated(), ...) {
    z0pre_given <- !missing(z0pre)
    z0exp_given <- !missing(z0exp)
    z0_args_given <- z0pre_given || z0exp_given
    if (lifecycle::is_present(z0)) {
        lifecycle::deprecate_warn("2.2.3", "setExtMort(z0)", "setExtMort(ext_mort)")
        ext_mort <- z0
    }
    assert_that(is.flag(reset),
                is.number(z0pre), z0pre >= 0,
                is.number(z0exp))

    if (reset) {
        if (!is.null(ext_mort)) {
            warning("Because you set `reset = TRUE`, the value you provided ",
                    "for `ext_mort` will be ignored and a value will be ",
                    "calculated from the species parameters.")
            ext_mort <- NULL
        }
        comment(params@mu_b) <- NULL
    }

    if (!is.null(ext_mort)) {
        if (z0_args_given) {
            signal_ignored_z0_args(z0pre_given, z0exp_given,
                                   "`ext_mort` was supplied")
        }
        if (is.null(comment(ext_mort))) {
            if (is.null(comment(params@mu_b))) {
                comment(ext_mort) <- "set manually"
            } else {
                comment(ext_mort) <- comment(params@mu_b)
            }
        }
        assert_that(is.array(ext_mort),
                    identical(dim(ext_mort), dim(params@mu_b)))
        params@mu_b[] <- ext_mort
        comment(params@mu_b) <- comment(ext_mort)

        params@time_modified <- lubridate::now()
        return(params)
    }

    assert_that(noNA(params@species_params$w_inf))
    # A z0 value is fixed only if it is present in given_species_params. A
    # value found only in species_params was calculated previously and must be
    # recalculated so that changes to its inputs propagate.
    z0_is_given <- if ("z0" %in% names(params@given_species_params)) {
        !is.na(params@given_species_params$z0)
    } else {
        rep(FALSE, nrow(params@species_params))
    }
    calculate_z0 <- !z0_is_given
    params <- set_z0_default(params, z0pre = z0pre, z0exp = z0exp)
    if (any(calculate_z0)) {
        if (z0_args_given) {
            if (!"z0" %in% names(params@given_species_params)) {
                params@given_species_params$z0 <-
                    rep(NA_real_, nrow(params@given_species_params))
            }
            params@given_species_params$z0[calculate_z0] <-
                params@species_params$z0[calculate_z0]
        }
    } else if (z0_args_given) {
        signal_ignored_z0_args(z0pre_given, z0exp_given,
                               paste0("`z0` is already present in ",
                                      "`given_species_params` for every ",
                                      "species"))
    }
    params <- set_species_param_default(params, "z_ext", 0)
    params <- set_species_param_default(params, "d",
                                        params@species_params$n - 1)
    mu_b <- params@mu_b
    mu_b[] <- params@species_params$z0
    has_power_law <- params@species_params$z_ext != 0
    if (any(has_power_law)) {
        if (isTRUE(params@second_order_w[["bin_average"]])) {
            # Use the exact bin average of the power law z_ext * w^d over each
            # bin instead of point-sampling at the left bin edge. This makes
            # the external mortality sink second-order (in fact exact) in the
            # finite-volume scheme where mu_b multiplies the bin-averaged
            # abundance N_j. See the "Point values and bin averages" section
            # of the numerical-details vignette.
            sp <- which(has_power_law)
            for (i in sp) {
                w_pow <- power_law_bin_average(params@w, params@dw,
                                               params@species_params$d[i])
                mu_b[i, ] <- mu_b[i, ] +
                    params@species_params$z_ext[i] * w_pow
            }
        } else {
            mu_b[has_power_law, ] <- mu_b[has_power_law, ] +
                sweep(
                    outer(params@species_params$d[has_power_law], params@w,
                          function(d, w) w^d),
                    1, params@species_params$z_ext[has_power_law], "*"
                )
        }
    }

    # Prevent overwriting slot if it has been commented
    if (!is.null(comment(params@mu_b))) {
        # Issue warning but only if a change was actually requested
        if (different(mu_b, params@mu_b)) {
            signal_not_recalculated("mu_b", "external mortality rate",
                               "setExtMort(params, reset = TRUE)")
        }
        return(params)
    }
    params@mu_b[] <- mu_b

    params@time_modified <- lubridate::now()
    return(params)
}

signal_ignored_z0_args <- function(z0pre_given, z0exp_given, reason) {
    args <- c(if (z0pre_given) "`z0pre`", if (z0exp_given) "`z0exp`")
    argument_names <- paste(args, collapse = " and ")
    verb <- if (length(args) == 1) " was" else " were"
    message <- paste0(argument_names, verb, " ignored because ", reason, ".")
    signal_info("z0", message, level = 1, severity = "warning",
                unhandled = "show")
}

# Recalculate non-given constant external mortality parameters. This remains in
# the external-mortality setter's file because `z0` is owned by `setExtMort()`.
set_z0_default <- function(params, z0pre = 0.6,
                           z0exp = params@resource_params$n - 1) {
    assert_that(is(params, "MizerParams"),
                is.number(z0pre), z0pre >= 0,
                is.number(z0exp))
    no_sp <- nrow(params@species_params)
    given_z0 <- if ("z0" %in% names(params@given_species_params)) {
        params@given_species_params$z0
    } else {
        rep(NA_real_, no_sp)
    }
    calculate_z0 <- is.na(given_z0)
    default <- z0pre * params@species_params$w_inf^z0exp
    z0_was_missing <- !"z0" %in% names(params@species_params)
    if (z0_was_missing) {
        params@species_params$z0 <- rep(NA_real_, no_sp)
    }
    z0_is_given <- !calculate_z0
    params@species_params$z0[z0_is_given] <- given_z0[z0_is_given]
    params@species_params$z0[calculate_z0] <- default[calculate_z0]
    if (z0_was_missing && any(calculate_z0)) {
        message <- paste0("Using z0 = z0pre * w_inf ^ z0exp for calculated ",
                          "z0 values.")
        signal_info("z0", message)
    }
    params
}

#' @rdname setExtMort
#' @return `ext_mort()`: An `ArraySpeciesBySize` object (species x size) with
#'   the external mortality.
#' @export
ext_mort <- function(params) {
    UseMethod("ext_mort")
}
#' @export
ext_mort.MizerParams <- function(params) {
    # External mortality is a sink integrated against the abundance over the
    # bin; under second-order bin-averaging mu_b is the exact bin average
    # (see setExtMort()), so it is plotted at the geometric bin centre.
    ArraySpeciesBySize(params@mu_b,
                       value_name = "External mortality",
                       units = "1/year",
                       params = params,
                       representation = "average")
}

#' @rdname setExtMort
#' @param value ext_mort
#' @export
`ext_mort<-` <- function(params, value) {
    setExtMort(params, ext_mort = value)
}
