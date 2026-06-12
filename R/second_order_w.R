#' Get or set the second_order_w flags
#'
#' Controls whether second-order numerical methods are use.
#'
#' The slot is a named list with entries:
#'
#' \describe{
#'   \item{`flux`}{The advective-flux reconstruction scheme used in the
#'     numerical solver. `"upwind"` is the first-order upwind scheme.
#'     `"van_leer"` is the second-order scheme with the total-variation-
#'     diminishing van Leer limiter, which keeps abundances non-negative.
#'     `"centred"` is the second-order scheme with the unlimited centred flux,
#'     which is genuinely second order even at extrema but is not
#'     monotonicity-preserving (it can produce small over/undershoots and is
#'     best used with some physical diffusion).}
#'   \item{`bin_average`}{Logical. Controls whether bin-averaging is used for
#'     quantities that need it in order to be second-order precise in bin size.
#'     When `FALSE`, point-sampling at the left bin edge is used.}
#' }
#'
#' When `flux` is `"upwind"` and `bin_average` is `FALSE` (the defaults),
#' mizer preserves the behaviour of previous mizer versions. Setting both to
#' their second-order values gives a consistently second-order model.
#'
#' The setter accepts a single logical value (which sets both entries), a single
#' scheme name (which sets only `flux`), or a named vector to set individual
#' entries. The setter re-runs [setParams()] to rebuild precomputed arrays when
#' `bin_average` is changed.
#'
#' @param params A MizerParams object.
#' @return `second_order_w()`: A named list with entries `flux` (character) and
#'   `bin_average` (logical).
#' @export
second_order_w <- function(params) {
    params@second_order_w
}

#' @rdname second_order_w
#' @param value A single logical value (`TRUE` or `FALSE`) which sets both
#'   entries, a single flux scheme name (`"upwind"`, `"van_leer"` or
#'   `"centred"`) which sets only `flux`, or a named vector with entries
#'   `flux` (logical or scheme name) and/or `bin_average` (logical).
#' @return `second_order_w<-`: A MizerParams object with the `second_order_w`
#'   flags updated and, when `bin_average` is changed, all model parameters
#'   recalculated via [setParams()].
#' @export
`second_order_w<-` <- function(params, value) {
    old_bin_average <- params@second_order_w[["bin_average"]]

    # Translate a flux value (logical or scheme name) into a scheme string.
    flux_scheme <- function(v) {
        if (is.logical(v)) {
            if (length(v) != 1 || is.na(v)) {
                stop("second_order_w flux entry must not be NA")
            }
            if (v) "van_leer" else "upwind"
        } else {
            match.arg(as.character(v), c("upwind", "van_leer", "centred"))
        }
    }
    as_flag <- function(v, what) {
        v <- as.logical(v)
        if (length(v) != 1 || is.na(v)) {
            stop("second_order_w ", what, " entry must be TRUE or FALSE")
        }
        v
    }

    if (is.null(names(value)) && length(value) == 1) {
        scheme <- flux_scheme(value)
        params@second_order_w[["flux"]] <- scheme
        # A single logical also sets bin_average; a single scheme name does not.
        if (is.logical(value)) {
            params@second_order_w[["bin_average"]] <- as.logical(value)
        }
    } else if (!is.null(names(value))) {
        unknown <- setdiff(names(value), c("flux", "bin_average"))
        if (length(unknown) > 0) {
            stop("Unknown second_order_w entries: ",
                 paste(unknown, collapse = ", "),
                 ". Valid entries are: flux, bin_average")
        }
        if ("flux" %in% names(value)) {
            params@second_order_w[["flux"]] <- flux_scheme(value[["flux"]])
        }
        if ("bin_average" %in% names(value)) {
            params@second_order_w[["bin_average"]] <-
                as_flag(value[["bin_average"]], "bin_average")
        }
    } else {
        stop("second_order_w must be a single logical value, a single flux ",
             "scheme name, or a named vector with entries 'flux' ",
             "and/or 'bin_average'")
    }

    new_bin_average <- params@second_order_w[["bin_average"]]
    if (!identical(old_bin_average, new_bin_average)) {
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
