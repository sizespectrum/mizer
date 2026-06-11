#' Get or set the second_order_w flags
#'
#' Controls whether second-order numerical methods are used when calculating
#' size-dependent rates. The slot is a named logical vector with entries:
#'
#' \describe{
#'   \item{`flux_limiter`}{Controls whether a second-order advective flux
#'     (with flux limiter) is used for the growth term. When `FALSE`, a
#'     first-order upwind scheme is used.}
#'   \item{`bin_average`}{Controls whether bin-averaging is used for rates
#'     (predation kernel, selectivity, etc.). When `FALSE`, rates are
#'     point-sampled at the left bin edge.}
#' }
#'
#' When both are `FALSE` (the default), mizer preserves the behaviour of
#' previous mizer versions. Setting both to `TRUE` gives a consistently
#' second-order model.
#'
#' The setter accepts either a single logical value (which sets both entries)
#' or a named logical vector to set individual entries. The setter re-runs
#' [setParams()] to rebuild precomputed arrays.
#'
#' @param params A MizerParams object.
#' @return `second_order_w()`: A named logical vector with entries
#'   `flux_limiter` and `bin_average`.
#' @export
#' @family functions for setting parameters
second_order_w <- function(params) {
    params@second_order_w
}

#' @rdname second_order_w
#' @param value A single logical value (`TRUE` or `FALSE`) which sets both
#'   entries, or a named logical vector with entries `flux_limiter` and/or
#'   `bin_average`.
#' @return `second_order_w<-`: A MizerParams object with the `second_order_w`
#'   flags updated and all model parameters recalculated via [setParams()].
#' @export
`second_order_w<-` <- function(params, value) {
    if (is.logical(value) && length(value) == 1 && !is.na(value)) {
        params@second_order_w[] <- value
    } else if (is.logical(value) && !is.null(names(value))) {
        valid_names <- c("flux_limiter", "bin_average")
        unknown <- setdiff(names(value), valid_names)
        if (length(unknown) > 0) {
            stop("Unknown second_order_w entries: ",
                 paste(unknown, collapse = ", "),
                 ". Valid entries are: flux_limiter, bin_average")
        }
        if (any(is.na(value))) {
            stop("second_order_w entries must not be NA")
        }
        params@second_order_w[names(value)] <- value
    } else {
        stop("second_order_w must be a single logical value or a named ",
             "logical vector with entries 'flux_limiter' and/or 'bin_average'")
    }
    setParams(params)
}
