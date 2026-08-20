#' Keep resource abundance constant
#' 
#' If you set your resource dynamics to use this function then the resource
#' abundances are kept constant over time.
#' 
#' 
#' To set your model to keep the resource constant over time you do
#' ```
#' resource_dynamics(params) <- "resource_constant"
#' ```
#' where you should replace `params` with the name of the variable holding your
#' MizerParams object.
#' 
#' @inheritParams resource_semichemostat
#' @param ... Unused
#'   
#' @return Vector containing the resource number density in each size class at
#'   the next timestep
#' @export
#' @family resource dynamics functions
#' @seealso [setResource()]
#' @examples
#' params <- NS_params
#' resource_dynamics(params) <- "resource_constant"
resource_constant <- function(params, n_pp, ...) {
    return(n_pp)
}



#' Resource parameters
#' 
#' The recommended way to change the resource dynamics parameters is to use
#' [setResource()]. The `resource_params` list contains values that are helpful
#' in setting up the actual size-dependent parameters with [setResource()]. If
#' you have specified a custom resource dynamics function that requires
#' additional parameters, then these should also be added to the
#' `resource_params` list.
#' 
#' The `resource_params` list will at least contain the slots `kappa`, `lambda`,
#' `w_pp_cutoff` and `n`.
#' 
#' The resource parameter `n` is the exponent for the power-law form for the
#' replenishment rate \eqn{r_R(w)}: \deqn{r_R(w) = r_R\, w^{n-1}.}{r_R(w) = r_R
#' w^{n-1}.}
#'
#' The resource parameter `lambda` (\eqn{\lambda}) is the exponent for the
#' power-law form for the carrying capacity \eqn{c_R(w)} and `w_pp_cutoff` is
#' its cutoff value: \deqn{c_R(w) = c_R w^{-\lambda}} for all \eqn{w} less than
#' `w_pp_cutoff` and zero for larger sizes.
#'
#' The resource parameter `kappa` (\eqn{\kappa}) is the coefficient \eqn{c_R} of
#' the carrying capacity in the power law above, so
#' \deqn{c_R(w) = \kappa\, w^{-\lambda}}{c_R(w) = \kappa w^{-\lambda}}
#' for all \eqn{w} less than `w_pp_cutoff` and zero for larger sizes. Changing
#' `kappa` therefore rescales the carrying capacity. It has a second role in that
#' the same expression also set the initial resource abundance when the model was
#' created: \deqn{N_R(w) = \kappa\, w^{-\lambda}.}{N_R(w) = \kappa w^{-\lambda}.}
#' Unlike the carrying capacity, however, the initial resource abundance is
#' **not** updated when you subsequently change `kappa` (or call [setResource()]).
#'
#' The resource parameters `a` and `b` give the allometric weight-length
#' relationship \eqn{w = a l^b} of the resource, with \eqn{w} in grams and
#' \eqn{l} in centimetres. They feed none of the rates; they exist so that the
#' resource can be shown on the length-based plots (`size_axis = "l"`) alongside
#' the species. They default to the equivalent spherical diameter of an organism
#' with the density of water, \eqn{a = \pi/6} and \eqn{b = 3}, which is the
#' convention plankton ecology uses for a composite of many taxa. This is a
#' different convention from the one the species use, so the resource and the
#' species each sit on the length axis at their own; see
#' [resource_length_defaults].
#'
#' Assigning to `resource_params` only rebuilds the size-dependent resource rate
#' and capacity arrays from these scalars (leaving any arrays you have set
#' manually untouched). Changing `lambda` also recalculates any `q` and `gamma`
#' species parameters that mizer calculated, and changing `kappa` recalculates
#' any calculated `gamma`; values you supplied explicitly are preserved. It
#' does **not** balance the resource, i.e. it does not adjust one of the rate or
#' capacity to keep the resource at the steady state where it replenishes at
#' the rate at which it is consumed. This mirrors the way the species
#' parameters feed the species rates. If you want to preserve the steady state
#' after changing a resource scalar, call [setResource()] with the appropriate
#' argument (which balances by default).
#'
#' @param params A MizerParams object
#' @return A named list of resource parameters.
#' @seealso [setResource()]
#' @export
resource_params <- function(params) {
    params@resource_params
}

#' @rdname resource_params
#' @param value A named list of resource parameters.
#' @export
`resource_params<-` <- function(params, value) {
    assert_that(
        is.number(value$lambda),
        value$lambda >= 0,
        is.number(value$kappa),
        value$kappa >= 0,
        is.number(value[["n"]]),
        value$n >= 0,
        is.number(value$w_pp_cutoff),
        value$w_pp_cutoff > min(params@w_full),
        value$w_pp_cutoff < max(params@w_full)
    )
    if (!is.null(value$r_pp)) {
        assert_that(is.number(value$r_pp), value$r_pp >= 0)
    }
    if (!is.null(value[["a"]])) {
        assert_that(is.number(value[["a"]]), value[["a"]] > 0)
    }
    if (!is.null(value[["b"]])) {
        assert_that(is.number(value[["b"]]), value[["b"]] > 0)
    }
    
    scalars <- c("kappa", "lambda", "n", "w_pp_cutoff", "r_pp")
    changed <- scalars[vapply(scalars, function(scalar) {
        !identical(params@resource_params[[scalar]], value[[scalar]])
    }, logical(1))]

    # Warn about the changes that cannot take effect because the array they
    # feed has been set by hand, before the change is applied.
    with_info_level(signal_frozen_changes(params, changed))

    params@resource_params <- value
    params@time_modified <- lubridate::now()
    # Setting the scalar resource parameters only rebuilds the size-dependent
    # rate and capacity arrays from those scalars (respecting any manual
    # overrides). It does not balance the resource: that is a deliberate action
    # reserved for `setResource()`. This mirrors how `species_params<-` feeds
    # the species-parameter rates.
    setResource(params, resource_params_changed = changed, balance = FALSE)
}
