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
#' @section The weight:
#' The weight \eqn{K(w)} is supplied already evaluated on the size grid. It can
#' be
#' \itemize{
#'   \item a single number (the default `weight = 1` integrates the abundance
#'     density itself, giving numbers),
#'   \item a vector with one value for each size bin, which is then used for all
#'     species,
#'   \item a matrix (species x size), for example `params@maturity`,
#'   \item an array with further dimensions in front, for example the gear x
#'     species x size array returned by [getFMortGear()] or the time x species x
#'     size array returned by `getFMort(sim)`. Those extra dimensions are
#'     carried through to the result.
#' }
#' If the weight is a product of several size-dependent factors, pass the whole
#' product: bin-averaging is applied to the product as a single weight, which is
#' not the same as averaging the factors separately.
#'
#' Do **not** include the bin widths `params@dw` in the weight and do not
#' bin-average the weight yourself; `sizeIntegral()` does both.
#'
#' @section Shape of the result:
#' The size dimension is integrated out. The remaining dimensions are those of
#' `n` together with any extra dimensions of `weight`, so
#' \itemize{
#'   \item with a `MizerParams` object and a species x size weight the result is
#'     a named vector with one value per species,
#'   \item with a `MizerSim` object it is an [ArrayTimeBySpecies] object (time x
#'     species),
#'   \item with a gear x species x size weight the extra `gear` dimension is
#'     kept, giving a gear x species array (or time x gear x species for a
#'     `MizerSim`).
#' }
#' Dimensions of `weight` other than the last two are matched to the dimensions
#' of `n` by the names of their dimnames, so a weight whose first dimension is
#' named `"time"` is lined up with the times of the simulation rather than
#' producing an outer product.
#'
#' @param object A `MizerParams` or a `MizerSim` object.
#' @param weight The weight \eqn{K(w)} of the integral, evaluated on the size
#'   grid. See the section "The weight" below. Defaults to 1.
#' @param n The abundance density. Either a species x size matrix or a time x
#'   species x size array. Defaults to the initial abundance `initialN(object)`
#'   for a `MizerParams` object and to the saved abundances `object@n` for a
#'   `MizerSim` object.
#' @param ... Arguments passed to [get_size_range_array()] to select the size
#'   range to integrate over, i.e. `min_w`, `max_w`, `min_l` and `max_l`.
#' @param value_name A string giving a human-readable name for the value, used
#'   when the result is wrapped in a mizer array class.
#' @param units A string giving the units of the result, used when the result is
#'   wrapped in a mizer array class.
#'
#' @return The value of the integral, see the section "Shape of the result"
#'   above.
#' @export
#' @family summary functions
#' @concept summary_function
#' @seealso [get_size_range_array()], [bin_average_weight()],
#'   [second_order_w()]
#' @examples
#' # The biomass of each species, i.e. what getBiomass() does
#' sizeIntegral(NS_params, weight = NS_params@w)
#'
#' # ... restricted to a size range
#' sizeIntegral(NS_params, weight = NS_params@w, min_w = 10, max_w = 1000)
#'
#' # The numbers of individuals larger than 10g
#' sizeIntegral(NS_params, min_w = 10)
#'
#' # Spawning stock biomass: the weight is the product maturity * w
#' K <- sweep(NS_params@maturity, 2, NS_params@w, "*")
#' sizeIntegral(NS_params, weight = K)
#'
#' # An indicator through time, ready to plot
#' biomass <- sizeIntegral(NS_sim, weight = NS_params@w,
#'                         value_name = "Biomass", units = "g")
#' biomass[c("1972", "2010"), c("Herring", "Cod")]
sizeIntegral <- function(object, weight = 1, n = NULL, ...,
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
    weight_labels <- size_dim_labels(weight, "weight", no_sp, no_w)

    # The size-range mask is part of the weight: bin-averaging the mask
    # together with the rest of the weight is what makes the bin straddling
    # the boundary of the size range contribute only partially.
    mask <- get_size_range_array(params, ...)
    storage.mode(mask) <- "double"
    mask_labels <- c("sp", "w")

    # Build the full weight, then do the bin integral on it. The abundance is
    # already a bin average and the bin widths are exact, so neither of them is
    # ever bin-averaged.
    K_labels <- merge_dim_labels(mask_labels, weight_labels)
    extent <- dim_extents(list(mask, weight, n),
                          list(mask_labels, weight_labels, n_labels))
    K <- broadcast_dims(mask, mask_labels, K_labels, extent) *
        broadcast_dims(weight, weight_labels, K_labels, extent)
    K <- bin_average_weight(K, params)
    K <- sweep(K, length(K_labels), params@dw, "*")

    # Contract against the abundance over the size dimension. The sum is done
    # as a matrix multiplication so that it agrees to the last bit with the
    # matrix multiplications that mizer has always used for these integrals.
    labels <- merge_dim_labels(n_labels, K_labels)
    integrand <- broadcast_dims(n, n_labels, labels, extent) *
        broadcast_dims(K, K_labels, labels, extent)
    result_labels <- labels[-length(labels)]
    result <- matrix(integrand, ncol = no_w) %*% rep(1, no_w)
    dn <- collect_dimnames(result_labels, list(n, weight, mask),
                           list(n_labels, weight_labels, mask_labels))
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
#' inputs occur once in the result, which is how a weight with a `"time"`
#' dimension is lined up with the times of the abundance rather than multiplied
#' out against them.
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
