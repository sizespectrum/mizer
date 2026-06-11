#' Get or set the second_order_w flags
#'
#' Controls whether second-order bin-averaged rate quadratures are used when
#' calculating size-dependent rates. The slot is a named logical vector with
#' entries `flux_limiter` and `bin_average`. When both are `FALSE` (the
#' default), mizer uses point-sampled (left-edge) rates, preserving the
#' behaviour of previous mizer versions.
#'
#' The setter accepts either a single logical value (which sets both entries)
#' or a named logical vector to set individual entries. Setting both to `TRUE`
#' gives a consistently second-order model. The setter re-runs [setParams()]
#' to rebuild precomputed arrays.
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
