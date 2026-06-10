#' Get or set the high_order flag
#'
#' Controls whether second-order bin-averaged rate quadratures are used when
#' calculating size-dependent rates. When `FALSE` (the default), mizer uses
#' point-sampled (left-edge) rates, preserving the behaviour of previous mizer
#' versions. Set to `TRUE` to enable bin-averaged / bin-integrated quantities
#' consistent with the finite-volume representation.
#'
#' Setting `high_order` to `TRUE` switches on all second-order components
#' together for a consistently second-order model. The setter re-runs
#' [setParams()] to rebuild precomputed arrays.
#'
#' @param params A MizerParams object.
#' @return `high_order()`: A single logical value.
#' @export
#' @family functions for setting parameters
high_order <- function(params) {
    params@high_order
}

#' @rdname high_order
#' @param value A single logical value (`TRUE` or `FALSE`).
#' @return `high_order<-`: A MizerParams object with the `high_order` flag
#'   updated and all model parameters recalculated via [setParams()].
#' @export
`high_order<-` <- function(params, value) {
    assert_that(is.flag(value))
    params@high_order <- value
    setParams(params)
}
