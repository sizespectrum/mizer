# Integrals over the size spectrum
#
# Copyright 2026 Gustav Delius.
# Distributed under the GPL 3 or later.

#' Integrate a quantity over the size spectrum
#'
#' `r lifecycle::badge("experimental")`
#' Calculates
#' \deqn{\int_{w_{min}}^{w_{max}} N_i(w)\, K_i(w)\, dw}
#' for each species \eqn{i}, using the quadrature scheme that the model is
#' actually using. This is the recommended way to write your own summary or
#' indicator function: it selects the size range, applies the bin-averaging
#' appropriate to the model's [second_order_w()] setting and wraps the result in
#' the appropriate mizer array class, so that none of those rules need to be
#' remembered. The built-in summary functions like [getBiomass()], [getN()],
#' [getSSB()] and [getYield()] are all implemented with it.
#'
#' @section The weighting factor:
#' The weighting factor \eqn{K(w)} is supplied already evaluated on the size
#' grid. It can be
#' \itemize{
#'   \item a single number (the default `weighting = 1` integrates the abundance
#'     density itself, giving numbers),
#'   \item a vector with one value for each size bin, which is then used for all
#'     species,
#'   \item a matrix (species x size), for example `params@maturity`,
#'   \item an array with further dimensions in front, for example the gear x
#'     species x size array returned by [getFMortGear()] or the time x species x
#'     size array returned by `getFMort(sim)`. Those extra dimensions are
#'     carried through to the result.
#' }
#' If the weighting factor is a product of several size-dependent factors, pass
#' the whole product: bin-averaging is applied to the product as a single
#' weighting factor, which is not the same as averaging the factors separately.
#'
#' Do **not** include the bin widths `params@dw` in the weighting factor and do
#' not bin-average it yourself; `sizeIntegral()` does both.
#'
#' @section Shape of the result:
#' The size dimension is integrated out. The remaining dimensions are those of
#' `n` together with any extra dimensions of `weighting`, so
#' \itemize{
#'   \item with a `MizerParams` object and a species x size weighting the result is
#'     a named vector with one value per species,
#'   \item with a `MizerSim` object it is an [ArrayTimeBySpecies] object (time x
#'     species),
#'   \item with a gear x species x size weighting the extra `gear` dimension is
#'     kept, giving a gear x species array (or time x gear x species for a
#'     `MizerSim`).
#' }
#' Dimensions of `weighting` other than the last two are matched to the dimensions
#' of `n` by the names of their dimnames, so a weighting whose first dimension is
#' named `"time"` is lined up with the times of the simulation rather than
#' producing an outer product.
#'
#' @section The size range:
#' The integral runs over the size range given by `min_w`, `max_w`, `min_l` and
#' `max_l`, with lengths converted to weights as described in
#' [get_size_range_array()]. Bins that lie entirely inside the size range
#' contribute in full and bins that lie entirely outside it not at all. How the
#' bins at the two ends of the range are treated is decided by `exact_bounds`:
#' \itemize{
#'   \item With `exact_bounds = TRUE`, the default for a model of
#'     [defaults_edition()] 3 or later, each bin contributes exactly the part of
#'     it that lies inside the size range. When the model bin-averages its
#'     integrals (see [second_order_w()]), the weighting factor is interpolated
#'     linearly across the bin and integrated over the part inside the range by
#'     the trapezoidal rule. Otherwise the weighting factor at the left edge of
#'     the bin is multiplied by the width of that part. The integral therefore
#'     changes continuously as a limit moves through a bin. An upper limit at or
#'     above the largest grid point `max(w(params))` includes the whole of the
#'     last bin, so the default limits cover the whole size spectrum.
#'   \item With `exact_bounds = FALSE`, the default for models of editions 1
#'     and 2, the size range selects whole bins by their left edge, as
#'     [get_size_range_array()] does. When the model bin-averages its
#'     integrals, that selection is averaged together with the weighting
#'     factor, so that a bin at an end of the range counts about half,
#'     wherever in the bin the limit lies.
#' }
#'
#' @param object A `MizerParams` or a `MizerSim` object.
#' @param weighting The weighting factor \eqn{K(w)} of the integral, evaluated
#'   on the size grid. See the section "The weighting factor" below. Defaults to 1.
#' @param n The abundance density. Either a species x size matrix or a time x
#'   species x size array. Defaults to the initial abundance `initialN(object)`
#'   for a `MizerParams` object and to the saved abundances `object@n` for a
#'   `MizerSim` object.
#' @param ... Arguments selecting the size range to integrate over, i.e.
#'   `min_w`, `max_w`, `min_l` and `max_l`, as described in
#'   [get_size_range_array()]. See the section "The size range" below.
#' @param exact_bounds Whether the bins at the ends of the size range contribute
#'   exactly the part of them that lies inside the range (`TRUE`) or are
#'   included or excluded whole according to their left edge (`FALSE`). The
#'   default `NULL` uses `TRUE` for a model of [defaults_edition()] 3 or later
#'   and `FALSE` for earlier editions. See the section "The size range" in
#'   [sizeIntegral()].
#' @param value_name A string giving a human-readable name for the value, used
#'   when the result is wrapped in a mizer array class.
#' @param units A string giving the units of the result, used when the result is
#'   wrapped in a mizer array class.
#'
#' @return The value of the integral, see the section "Shape of the result"
#'   above.
#' @export
#' @seealso [get_size_range_array()], [bin_average_weight()],
#'   [second_order_w()]
#' @examples
#' # The biomass of each species, i.e. what getBiomass() does
#' sizeIntegral(NS_params, weighting = NS_params@w)
#'
#' # ... restricted to a size range
#' sizeIntegral(NS_params, weighting = NS_params@w, min_w = 10, max_w = 1000)
#'
#' # The numbers of individuals larger than 10g
#' sizeIntegral(NS_params, min_w = 10)
#'
#' # ... counting only the part above 10g of the bin that contains 10g
#' sizeIntegral(NS_params, min_w = 10, exact_bounds = TRUE)
#'
#' # Spawning stock biomass: the weighting is the product maturity * w
#' K <- sweep(NS_params@maturity, 2, NS_params@w, "*")
#' sizeIntegral(NS_params, weighting = K)
#'
#' # An indicator through time, ready to plot
#' biomass <- sizeIntegral(NS_sim, weighting = NS_params@w,
#'                         value_name = "Biomass", units = "g")
#' biomass[c("1972", "2010"), c("Herring", "Cod")]
sizeIntegral <- function(object, weighting = 1, n = NULL, ...,
                         exact_bounds = NULL,
                         value_name = NULL, units = NULL) {
    if (is(object, "MizerSim")) {
        params <- object@params
        if (is.null(n)) n <- object@n
    } else if (is(object, "MizerParams")) {
        params <- object
        if (is.null(n)) n <- object@initial_n
    } else {
        stop("`object` must be a MizerParams or a MizerSim object.")
    }
    no_w <- length(params@w)
    no_sp <- nrow(params@species_params)

    n_labels <- size_dim_labels(n, "n", no_sp, no_w)
    if (!("sp" %in% n_labels)) {
        stop("`n` must have a dimension running over the species.")
    }
    weighting_labels <- size_dim_labels(weighting, "weighting", no_sp, no_w)

    # The quadrature weights hold everything except the abundance: the
    # weighting factor, bin-averaged as the model requires, restricted to the
    # size range and multiplied by the bin widths. The abundance is already a
    # bin average and the bin widths are exact, so neither of them is ever
    # bin-averaged.
    qw <- size_integral_weights(params, weighting = weighting, ...,
                                weighting_labels = weighting_labels,
                                exact_bounds = exact_bounds)
    K <- qw$weights
    K_labels <- qw$labels
    extent <- dim_extents(list(K, n), list(K_labels, n_labels))

    # Contract against the abundance over the size dimension. The sum is done
    # as a matrix multiplication so that it agrees to the last bit with the
    # matrix multiplications that mizer has always used for these integrals.
    labels <- merge_dim_labels(n_labels, K_labels)
    integrand <- broadcast_dims(n, n_labels, labels, extent) *
        broadcast_dims(K, K_labels, labels, extent)
    result_labels <- labels[-length(labels)]
    result <- matrix(integrand, ncol = no_w) %*% rep(1, no_w)
    dn <- collect_dimnames(result_labels, list(n, K), list(n_labels, K_labels))
    if (length(result_labels) == 1) {
        result <- drop(result)
        if (!is.null(dn)) names(result) <- dn[[1]]
    } else {
        dim(result) <- as.integer(extent[result_labels])
        dimnames(result) <- dn
    }

    if (identical(result_labels, c("time", "sp"))) {
        return(ArrayTimeBySpecies(result, value_name = value_name,
                                  units = units, params = params))
    }
    result
}

#' Quadrature weights of an integral over a size range
#'
#' Internal helper that builds the weights \eqn{W_{ij}} with which
#' [sizeIntegral()] discretises the integral
#' \deqn{\int_{w_{min}}^{w_{max}} N_i(w)\, K_i(w)\, dw \approx
#'   \sum_j N_{ij}\, W_{ij}.}
#' The weights hold everything except the abundance: the weighting factor
#' \eqn{K}, bin-averaged according to the model's [second_order_w()] setting,
#' the restriction to the size range, including the treatment of the bins at
#' its ends that is described in the section "The size range" of
#' [sizeIntegral()], and the bin widths.
#'
#' The weights are built either on the size grid of the species, `w(params)`,
#' with one row for each species, or with `resource = TRUE` on the size grid of
#' the resource, `w_full(params)`, with a single row. Length limits are
#' converted to weight limits by [size_range_limits()], with the species
#' parameters `a` and `b` for the species and [resource_length_params()] for
#' the resource.
#'
#' A bin that lies entirely inside the size range gets the same weight, to the
#' last bit, whether or not the bounds are exact: \eqn{K_j\,\Delta w_j}, or
#' with bin averaging \eqn{\tfrac12(K_j + K_{j+1})\,\Delta w_j}, with the
#' factors in the same order as in the calculations that mizer has always
#' done. With exact bounds, a bin \eqn{[w_j, w_{j+1}]} of which only the part
#' between the fractions \eqn{s_0 < s_1} of its width lies inside the size
#' range gets the weight
#' \eqn{(s_1 - s_0)\,\Delta w_j\,\big(K_j + (K_{j+1} - K_j)\,\tfrac{s_0 +
#' s_1}{2}\big)} with bin averaging, which is the trapezoidal rule on that part
#' with \eqn{K} interpolated linearly across the bin, and \eqn{(s_1 -
#' s_0)\,\Delta w_j\,K_j} without. The last bin has no right-hand neighbour, so
#' \eqn{K} is taken to be constant across it, as [trapezoidal_bin_average()]
#' does, and an upper limit at or above the last grid point covers it whole.
#'
#' @param params A MizerParams object.
#' @param weighting The weighting factor \eqn{K(w)}, evaluated on the size
#'   grid, as described for [sizeIntegral()]. On the resource grid it runs over
#'   `w_full(params)`.
#' @param min_w,max_w,min_l,max_l The limits of the size range, as for
#'   [get_size_range_array()]. `NULL` leaves the range open on that side.
#' @param ... Ignored, so that the summary functions can pass on all of their
#'   extra arguments.
#' @param weighting_labels The dimension labels of `weighting`, as returned by
#'   [size_dim_labels()]. If `NULL` they are worked out here.
#' @param exact_bounds Whether bins at the ends of the size range contribute
#'   exactly the part that lies inside it. `NULL` means `TRUE` for a model of
#'   [defaults_edition()] 3 or later and `FALSE` otherwise.
#' @param resource If `TRUE`, build the weights on the resource grid
#'   `w_full(params)` instead of the species grid `w(params)`.
#' @return A list with two entries: `weights`, an array whose last dimension
#'   runs over the size grid and whose other dimensions are those of
#'   `weighting` together with a dimension `sp` over the rows of the grid, and
#'   `labels`, the dimension labels of `weights`.
#' @concept helper
#' @keywords internal
size_integral_weights <- function(params, weighting = 1,
                                  min_w = NULL, max_w = NULL,
                                  min_l = NULL, max_l = NULL, ...,
                                  weighting_labels = NULL,
                                  exact_bounds = NULL, resource = FALSE) {
    if (resource) {
        w <- params@w_full
        dw <- params@dw_full
        rows <- "Resource"
    } else {
        w <- params@w
        dw <- params@dw
        rows <- as.character(params@species_params$species)
    }
    no_w <- length(w)
    no_row <- length(rows)
    if (is.null(exact_bounds)) {
        exact_bounds <- defaults_edition(params) >= 3
    }
    if (!isTRUE(exact_bounds) && !isFALSE(exact_bounds)) {
        stop("`exact_bounds` must be TRUE, FALSE or NULL.")
    }
    if (is.null(weighting_labels)) {
        weighting_labels <- size_dim_labels(weighting, "weighting",
                                            no_row, no_w)
    }

    # The weight limits for each row. A missing limit is the end of the grid,
    # which is also what get_size_range_array() uses.
    limits <- size_range_limits(params,
                                min_w = if (is.null(min_w)) min(w) else min_w,
                                max_w = if (is.null(max_w)) max(w) else max_w,
                                min_l = min_l, max_l = max_l,
                                resource = resource)
    mask <- outer(limits$min_w, w, "<=") & outer(limits$max_w, w, ">=")
    storage.mode(mask) <- "double"
    row_labels <- c("sp", "w")

    K_labels <- merge_dim_labels(row_labels, weighting_labels)
    extent <- dim_extents(list(mask, weighting),
                          list(row_labels, weighting_labels))
    K <- broadcast_dims(weighting, weighting_labels, K_labels, extent)
    sweep_dw <- function(x) sweep(x, length(K_labels), dw, "*")

    if (!exact_bounds) {
        # The size-range mask is part of the weighting: bin-averaging the mask
        # together with the rest of the weighting is what makes the bin
        # straddling the boundary of the size range contribute partially.
        K <- broadcast_dims(mask, row_labels, K_labels, extent) * K
        K <- sweep_dw(bin_average_weight(K, params))
    } else {
        K <- exact_range_weights(K, K_labels, extent, w, dw, limits, params)
    }

    dn <- collect_dimnames(K_labels, list(weighting), list(weighting_labels))
    if (is.null(dn)) {
        dn <- stats::setNames(vector("list", length(K_labels)), K_labels)
    }
    if (is.null(dn[["sp"]])) dn["sp"] <- list(rows)
    dimnames(K) <- dn
    list(weights = K, labels = K_labels)
}

#' Quadrature weights with exactly integrated boundary bins
#'
#' Internal helper for [size_integral_weights()] that does its work when
#' `exact_bounds = TRUE`. See there for the formulas.
#'
#' @param K The weighting factor, already broadcast to the labels `K_labels`.
#' @param K_labels The dimension labels of `K`, containing `"sp"` and ending
#'   in `"w"`.
#' @param extent The extent of each label.
#' @param w,dw The size grid and its bin widths.
#' @param limits The weight limits of each row, as returned by
#'   [size_range_limits()].
#' @param params The MizerParams object, whose `second_order_w` slot decides
#'   whether the weighting factor is bin-averaged.
#' @return The quadrature weights, an array of the same shape as `K`.
#' @concept helper
#' @keywords internal
exact_range_weights <- function(K, K_labels, extent, w, dw, limits, params) {
    no_w <- length(w)
    row_labels <- c("sp", "w")
    w_right <- c(w[-1], w[no_w] + dw[no_w])
    bin_width <- w_right - w
    lower <- limits$min_w
    upper <- limits$max_w
    # An upper limit at or above the last grid point covers the whole last
    # bin. This is what makes the default limit, max(w), and also the explicit
    # max_w = max(w) that the plotting functions pass, cover the whole spectrum.
    upper[which(upper >= w[no_w])] <- Inf

    # The part of each bin inside the range, as fractions s0 <= s1 of the bin
    # width (rows x size).
    s0 <- outer(lower, w, "-")
    s0 <- pmin(pmax(sweep(s0, 2, bin_width, "/"), 0), 1)
    s1 <- outer(upper, w, "-")
    s1 <- pmin(pmax(sweep(s1, 2, bin_width, "/"), 0), 1)
    s1 <- pmax(s1, s0)
    # A missing limit is not a whole bin either, so that it gives NA.
    partial <- is.na(s0) | is.na(s1) | s0 > 0 | s1 < 1

    bin_average <- isTRUE(params@second_order_w[["bin_average"]])
    # Whole bins get exactly the weights of the rule with a mask of ones.
    Kbar <- bin_average_weight(K, params)
    if (any(partial)) {
        partial_K <- broadcast_dims(partial, row_labels, K_labels, extent)
        covered <- broadcast_dims(s1 - s0, row_labels, K_labels, extent)
        if (bin_average) {
            # Interpolate K linearly across the bin and take the trapezoidal
            # rule on the covered part, which is the value at its midpoint
            # times its width. The last bin has no right-hand neighbour, so K
            # is constant across it, as in trapezoidal_bin_average().
            size_index <- slice.index(K, length(dim(K)))
            K_right <- K
            K_right[size_index < no_w] <- K[size_index > 1]
            mid <- broadcast_dims((s0 + s1) / 2, row_labels, K_labels, extent)
            K_part <- K + (K_right - K) * mid
        } else {
            K_part <- K
        }
        Kbar[partial_K] <- K_part[partial_K] * covered[partial_K]
    }
    sweep(Kbar, length(K_labels), dw, "*")
}

#' Identify the dimensions of an array over the size grid
#'
#' Internal helper for [sizeIntegral()]. Returns a label for each dimension of
#' `x`. The last dimension must run over the size grid and is labelled `"w"`.
#' Other dimensions are labelled from the names of their dimnames, if they have
#' any, with `"species"` and `"size"` normalised to mizer's `"sp"` and `"w"`.
#' An unnamed second-to-last dimension is labelled `"sp"` if its extent is the
#' number of species. Any remaining unnamed dimension gets a unique label of its
#' own, so that it is carried through to the result rather than matched against
#' a dimension of the abundance.
#'
#' A scalar has no dimensions and gets no labels.
#'
#' @param x The array to label.
#' @param arg The name of the argument holding `x`, used in error messages.
#' @param no_sp The number of species in the model.
#' @param no_w The number of size bins in the model.
#' @return A character vector with one label for each dimension of `x`.
#' @concept helper
#' @keywords internal
size_dim_labels <- function(x, arg, no_sp, no_w) {
    if (!is.numeric(x) && !is.logical(x)) {
        stop("`", arg, "` must be numeric.")
    }
    d <- dim(x)
    if (is.null(d)) {
        if (length(x) == 1) return(character(0))
        d <- length(x)
    }
    nd <- length(d)
    if (d[[nd]] != no_w) {
        stop("The last dimension of `", arg, "` must run over the ", no_w,
             " size bins of the model, but it has length ", d[[nd]], ".")
    }
    labels <- rep(NA_character_, nd)
    nm <- names(dimnames(x))
    if (!is.null(nm)) {
        nm[is.na(nm)] <- ""
        known <- c(sp = "sp", species = "sp", Species = "sp",
                   time = "time", Time = "time", t = "time",
                   w = "w", size = "w", weight = "w")
        labels <- ifelse(nm %in% names(known), known[nm], nm)
        labels[!nzchar(nm)] <- NA_character_
    }
    labels[[nd]] <- "w"
    if (nd >= 2 && is.na(labels[[nd - 1]])) {
        if (d[[nd - 1]] != no_sp) {
            stop("The second-to-last dimension of `", arg, "` must run over ",
                 "the ", no_sp, " species of the model, but it has length ",
                 d[[nd - 1]], ". If it represents something else, give `", arg,
                 "` dimnames whose names say what its dimensions are.")
        }
        labels[[nd - 1]] <- "sp"
    }
    unnamed <- is.na(labels)
    if (any(unnamed)) {
        labels[unnamed] <- paste0(arg, ".dim", which(unnamed))
    }
    if (anyDuplicated(labels)) {
        stop("The dimensions of `", arg, "` do not have distinct names.")
    }
    labels
}

#' Merge two ordered sets of dimension labels
#'
#' Internal helper for [sizeIntegral()]. Interleaves the labels of two arrays
#' into the labels of the array holding their product, keeping the relative
#' order of the labels within each of the two inputs. Labels that occur in both
#' inputs occur once in the result, which is how a weighting array with a
#' `"time"` dimension is lined up with the times of the abundance rather than
#' multiplied out against them.
#'
#' @param a,b Character vectors of dimension labels.
#' @return A character vector of labels containing each label of `a` and `b`
#'   once.
#' @concept helper
#' @keywords internal
merge_dim_labels <- function(a, b) {
    out <- character(0)
    i <- 1L
    j <- 1L
    while (i <= length(a) && j <= length(b)) {
        if (a[[i]] == b[[j]]) {
            out <- c(out, a[[i]])
            i <- i + 1L
            j <- j + 1L
        } else if (!(a[[i]] %in% b[j:length(b)])) {
            # a[[i]] is not waiting for anything in b, so it can go first
            out <- c(out, a[[i]])
            i <- i + 1L
        } else {
            out <- c(out, b[[j]])
            j <- j + 1L
        }
    }
    if (i <= length(a)) out <- c(out, a[i:length(a)])
    if (j <= length(b)) out <- c(out, b[j:length(b)])
    out
}

#' Collect the extent of each labelled dimension
#'
#' Internal helper for [sizeIntegral()] that checks that arrays sharing a
#' dimension label agree on its extent.
#'
#' @param arrays A list of arrays.
#' @param labels A list of the corresponding label vectors.
#' @return A named integer vector giving the extent of each label.
#' @concept helper
#' @keywords internal
dim_extents <- function(arrays, labels) {
    extent <- integer(0)
    for (k in seq_along(arrays)) {
        d <- dim_from_labels(arrays[[k]], labels[[k]])
        for (m in seq_along(d)) {
            label <- labels[[k]][[m]]
            if (label %in% names(extent)) {
                if (extent[[label]] != d[[m]]) {
                    stop("The '", label, "' dimensions do not have the same ",
                         "length in all arguments: ", extent[[label]], " and ",
                         d[[m]], ".")
                }
            } else {
                extent[[label]] <- d[[m]]
            }
        }
    }
    extent
}

#' The dimensions of an array, given its labels
#'
#' Internal helper for [sizeIntegral()]. A scalar has no labels and no
#' dimensions, a vector has one.
#'
#' @param x An array, vector or scalar.
#' @param labels Its dimension labels.
#' @return An integer vector of dimensions, of the same length as `labels`.
#' @concept helper
#' @keywords internal
dim_from_labels <- function(x, labels) {
    if (length(labels) == 0) return(integer(0))
    d <- dim(x)
    if (is.null(d)) length(x) else d
}

#' Expand an array to a larger set of labelled dimensions
#'
#' Internal helper for [sizeIntegral()]. Replicates `x` along the dimensions it
#' does not have and permutes its dimensions into the order given by
#' `target_labels`, so that arrays with different dimensions can be multiplied
#' together elementwise.
#'
#' @param x An array, vector or scalar.
#' @param labels The dimension labels of `x`.
#' @param target_labels The dimension labels of the result. Must contain all of
#'   `labels`.
#' @param extent A named vector giving the extent of each label.
#' @return An array with dimensions `extent[target_labels]`.
#' @concept helper
#' @keywords internal
broadcast_dims <- function(x, labels, target_labels, extent) {
    target_dim <- as.integer(extent[target_labels])
    if (identical(labels, target_labels)) {
        # Nothing to replicate or permute; only make sure it is a bare array
        # of the right shape.
        x <- unclass(x)
        attributes(x) <- NULL
        dim(x) <- target_dim
        return(x)
    }
    idx <- match(labels, target_labels)
    rest <- setdiff(seq_along(target_labels), idx)
    # Replicating along dimensions appended at the end is what `array()` does
    # for free, so put the new dimensions last and permute afterwards.
    tmp <- array(as.vector(x),
                 dim = c(dim_from_labels(x, labels), target_dim[rest]))
    tmp <- aperm(tmp, order(c(idx, rest)))
    dim(tmp) <- target_dim
    tmp
}

#' Assemble the dimnames of a broadcast array
#'
#' Internal helper for [sizeIntegral()] that takes the dimnames of each labelled
#' dimension from the first of the given arrays that has them.
#'
#' @param target_labels The dimension labels of the result.
#' @param arrays A list of arrays.
#' @param labels A list of the corresponding label vectors.
#' @return A named list of dimnames, or `NULL` if none of the arrays has any.
#' @concept helper
#' @keywords internal
collect_dimnames <- function(target_labels, arrays, labels) {
    dn <- stats::setNames(vector("list", length(target_labels)), target_labels)
    for (k in seq_along(arrays)) {
        dnk <- dimnames(arrays[[k]])
        if (is.null(dnk)) next
        for (m in seq_along(labels[[k]])) {
            label <- labels[[k]][[m]]
            if (label %in% target_labels && is.null(dn[[label]])) {
                dn[[label]] <- dnk[[m]]
            }
        }
    }
    if (all(vapply(dn, is.null, logical(1)))) return(NULL)
    dn
}
