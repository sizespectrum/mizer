#' Get or set the second_order_w flag
#'
#' Controls whether second-order bin-averaged rate quadratures are used when
#' calculating size-dependent rates. When `FALSE` (the default), mizer uses
#' point-sampled (left-edge) rates, preserving the behaviour of previous mizer
#' versions. Set to `TRUE` to enable bin-averaged / bin-integrated quantities
#' consistent with the finite-volume representation.
#'
#' Setting `second_order_w` to `TRUE` switches on all second-order components
#' together for a consistently second-order model. The setter re-runs
#' [setParams()] to rebuild precomputed arrays.
#'
#' @param params A MizerParams object.
#' @return `second_order_w()`: A single logical value.
#' @export
#' @family functions for setting parameters
second_order_w <- function(params) {
    params@second_order_w
}

#' @rdname second_order_w
#' @param value A single logical value (`TRUE` or `FALSE`).
#' @return `second_order_w<-`: A MizerParams object with the `second_order_w` flag
#'   updated and all model parameters recalculated via [setParams()].
#' @export
`second_order_w<-` <- function(params, value) {
    assert_that(is.flag(value))
    params@second_order_w <- value
    setParams(params)
}
