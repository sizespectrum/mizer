#' Get or set the second_order_w flags
#'
#' Controls whether second-order numerical methods are use.
#'
#' The slot is a named logical vector with entries:
#'
#' \describe{
#'   \item{`flux_limiter`}{Controls whether a second-order advective flux
#'     (with flux limiter) is used in the numerical solver. When `FALSE`, a
#'     first-order upwind scheme is used.}
#'   \item{`bin_average`}{Controls whether bin-averaging is used for quantities
#'     that need it in order to be second-order precise in bin size.
#'     When `FALSE`, point-sampling at the left bin edge is used.}
#' }
#'
#' When both are `FALSE` (the default), mizer preserves the behaviour of
#' previous mizer versions. Setting both to `TRUE` gives a consistently
#' second-order model.
#'
#' The setter accepts either a single logical value (which sets both entries)
#' or a named logical vector to set individual entries. The setter re-runs
#' [setParams()] to rebuild precomputed arrays when `bin_average` is changed.
#'
#' @param params A MizerParams object.
#' @return `second_order_w()`: A named logical vector with entries
#'   `flux_limiter` and `bin_average`.
#' @export
second_order_w <- function(params) {
    params@second_order_w
}

#' @rdname second_order_w
#' @param value A single logical value (`TRUE` or `FALSE`) which sets both
#'   entries, or a named logical vector with entries `flux_limiter` and/or
#'   `bin_average`.
#' @return `second_order_w<-`: A MizerParams object with the `second_order_w`
#'   flags updated and, when `bin_average` is changed, all model parameters
#'   recalculated via [setParams()].
#' @export
`second_order_w<-` <- function(params, value) {
    old_bin_average <- params@second_order_w[["bin_average"]]
    if (is.logical(value) && length(value) == 1 && !is.na(value) && is.null(names(value))) {
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
    new_bin_average <- params@second_order_w[["bin_average"]]
    if (old_bin_average != new_bin_average) {
        params <- setParams(params)
    }
    params
}

#' Trapezoidal bin-average of a per-bin weight
#'
#' Internal helper for the second-order summary integrals. A summary
#' diagnostic \eqn{\int N(w) K(w)\, dw} is discretised on the finite-volume
#' grid as \eqn{\sum_j N_j \bar K_j \Delta w_j}, where \eqn{N_j} is the cell
#' average of the density over bin \eqn{[w_j, w_{j+1}]}. To be second order in
#' the bin width the point weight \eqn{K(w_j)} must be replaced by the bin
#' average
#' \deqn{\bar K_j = \frac{1}{\Delta w_j}\int_{w_j}^{w_{j+1}} K(w)\,dw
#'   \approx \tfrac12\big(K(w_j) + K(w_{j+1})\big).}
#' The trapezoidal average \eqn{\tfrac12(K_j + K_{j+1})} is uniformly second
#' order and exact whenever \eqn{K} is linear in \eqn{w} (e.g. the first
#' moment \eqn{K = w}, for which it equals \eqn{(w_{j+1}^2 - w_j^2)/(2\Delta
#' w_j)}).
#'
#' The weight `K` is supplied already evaluated on the size grid (a vector
#' indexed over the bins, or a matrix with the size dimension running along the
#' columns). The top bin has no right-hand neighbour on the grid, so its weight
#' is left unaveraged (one-sided); the density there is negligible, so this
#' does not affect the second-order accuracy of the totals.
#'
#' This helper is shared with the reproduction integrals (issue #376), which
#' also need the trapezoidal bin-average of a composite weight.
#'
#' @param K A numeric vector of weights indexed over the size grid, or a
#'   numeric array whose last dimension runs over the size grid (e.g. a
#'   species-by-size matrix or a gear-by-species-by-size array).
#' @return An object of the same shape as `K` containing the trapezoidal
#'   bin-averaged weights.
#' @concept helper
#' @keywords internal
bin_average_weight <- function(K) {
    d <- dim(K)
    if (is.null(d)) {
        n <- length(K)
        if (n < 2) return(K)
        Kbar <- K
        Kbar[-n] <- 0.5 * (K[-n] + K[-1])
        return(Kbar)
    }
    # Average along the last (size) dimension.
    n <- d[length(d)]
    if (n < 2) return(K)
    idx_lo <- slice.index(K, length(d)) <= n - 1L
    idx_hi <- slice.index(K, length(d)) >= 2L
    Kbar <- K
    Kbar[idx_lo] <- 0.5 * (K[idx_lo] + K[idx_hi])
    Kbar
}

#' Bin-average a summary-integral weight when second-order is enabled
#'
#' Convenience wrapper around [bin_average_weight()] that is gated on the
#' `bin_average` entry of the model's `second_order_w` slot. When second-order
#' bin-averaging is switched off (the default), the weight `K` is returned
#' unchanged so that the summary functions reproduce the previous left-edge
#' Riemann sums byte-for-byte. When it is switched on, the trapezoidal
#' bin-average of the weight is returned.
#'
#' @param K A numeric vector of weights indexed over the size grid, or a
#'   numeric matrix with the size dimension running along the columns.
#' @param params A MizerParams object whose `second_order_w` slot controls the
#'   gating.
#' @return The weight `K`, bin-averaged when `params@second_order_w[["bin_average"]]`
#'   is `TRUE`, otherwise returned unchanged.
#' @concept helper
#' @keywords internal
bin_average_summary_weight <- function(K, params) {
    if (isTRUE(params@second_order_w[["bin_average"]])) {
        return(bin_average_weight(K))
    }
    K
}
