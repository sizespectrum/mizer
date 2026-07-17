# ArrayResourceBySize and ArrayTimeByResourceBySize S3 classes for resource
# size spectra
#
# Copyright 2026 Gustav Delius.
# Distributed under the GPL 3 or later.

#' S3 class for resource size spectra
#'
#' Several functions in mizer return a vector over the full size grid holding
#' a resource-related quantity such as the resource number density, the
#' resource mortality, the intrinsic resource birth rate or carrying capacity.
#' The `ArrayResourceBySize` class wraps these vectors to provide convenient
#' `print()`, `summary()`, `plot()`, and `as.data.frame()` methods.
#'
#' An `ArrayResourceBySize` object behaves just like a regular numeric vector
#' for arithmetic operations and subsetting. It carries three lightweight
#' attributes:
#' \itemize{
#'   \item `value_name` – a human-readable name for the value
#'       (e.g. "Resource mortality").
#'   \item `units` – the units of the value (e.g. "1/year").
#'   \item `params` – the `MizerParams` object that the value was computed from.
#' }
#'
#' @param x A numeric vector over the full size grid. For
#'   `is.ArrayResourceBySize()`, any object to test.
#' @param value_name A string giving the human-readable name for the value.
#' @param units A string giving the units (e.g. "1/year").
#' @param params A `MizerParams` object. Used for the resource colour and the
#'   size grid in the `plot()` method.
#'
#' @return An `ArrayResourceBySize` object (inherits from `numeric`).
#' @seealso [print()], [summary()], [as.data.frame()], [plot()]
#' @export
#' @examples
#' \donttest{
#' mort <- getResourceMort(NS_params)
#' is.ArrayResourceBySize(mort)
#' summary(mort)
#' plot(mort)
#' }
ArrayResourceBySize <- function(x, value_name = NULL, units = NULL,
                                params = NULL) {
    if (!is.numeric(x) || !is.null(dim(x))) {
        stop("`x` must be a numeric vector.")
    }
    if (!is.null(params) && length(x) == length(params@initial_n_pp) &&
            is.null(names(x))) {
        names(x) <- names(params@initial_n_pp)
    }
    structure(x,
        class = c("ArrayResourceBySize", "numeric"),
        value_name = value_name,
        units = units,
        params = params
    )
}

#' @rdname ArrayResourceBySize
#' @return `is.ArrayResourceBySize()` returns `TRUE` if `x` is an
#'   `ArrayResourceBySize` object, `FALSE` otherwise.
#' @export
is.ArrayResourceBySize <- function(x) {
    inherits(x, "ArrayResourceBySize")
}

#' @export
print.ArrayResourceBySize <- function(x, ...) {
    value_name <- attr(x, "value_name") %||% "ArrayResourceBySize"
    units_str <- attr(x, "units")
    header <- paste0(value_name, " (", length(x), " sizes)")
    if (!is.null(units_str) && nzchar(units_str)) {
        header <- paste0(header, " [", units_str, "]")
    }
    cat(header, "\n")

    w <- get_ArrayResourceBySize_w(x)
    vec <- unclass_resource(x)
    n <- length(vec)

    size_k <- fit_log_spaced_k(
        n, mizer_print_defaults$size_max, mizer_print_defaults$size_min,
        width_fn = function(k) {
            idx <- pick_log_spaced_indices(n, k)
            vector_display_width(vec[idx])
        })
    sz_idx <- pick_log_spaced_indices(n, size_k)
    print(vec[sz_idx])

    if (length(sz_idx) < n) {
        cat(format_truncation_note(length(sz_idx), n, "sizes",
                                   format_size_range_detail(w)), "\n")
    }
    invisible(x)
}

#' @export
summary.ArrayResourceBySize <- function(object, ...) {
    value_name <- attr(object, "value_name") %||% "ArrayResourceBySize"
    units_str <- attr(object, "units")
    vals <- unclass(object)

    result <- list(
        value_name = value_name,
        units = units_str,
        length = length(object),
        stats = data.frame(
            Min = min(vals, na.rm = TRUE),
            Mean = mean(vals, na.rm = TRUE),
            Max = max(vals, na.rm = TRUE),
            row.names = NULL
        )
    )
    class(result) <- "summary.ArrayResourceBySize"
    result
}

#' @export
print.summary.ArrayResourceBySize <- function(x, ...) {
    header <- x$value_name
    if (!is.null(x$units) && nzchar(x$units)) {
        header <- paste0(header, " [", x$units, "]")
    }
    cat(header, "\n")
    cat(x$length, "sizes\n\n")
    print(x$stats, row.names = FALSE)
    invisible(x)
}

#' @rdname plot
#' @usage NULL
#' @export
#' @examples
#' \donttest{
#' plot(getResourceMort(NS_params))
#' plot(initialNResource(NS_params))
#' }
plot.ArrayResourceBySize <- function(x, return_data = FALSE,
                                     log_x = TRUE, log_y = TRUE, log = NULL,
                                     wlim = c(NA, NA), ylim = c(NA, NA),
                                     y_ticks = 6, ...) {
    log_axes <- parsePlotLog(log, log_x = log_x, log_y = log_y)
    log_x <- log_axes$log_x
    log_y <- log_axes$log_y

    assert_that(length(wlim) == 2,
                length(ylim) == 2)
    value_name <- attr(x, "value_name") %||% "value"
    units_str <- attr(x, "units")
    params <- attr(x, "params")

    plot_dat <- prepare_ArrayResourceBySize_plot_data(x, wlim = wlim)

    if (return_data) return(plot_dat)

    y_label <- value_name
    if (!is.null(units_str) && nzchar(units_str)) {
        y_label <- paste0(value_name, " [", units_str, "]")
    }

    plotDataFrame(plot_dat, params, xlab = "Weight (g)",
                  ylab = y_label,
                  xtrans = if (log_x) "log10" else "identity",
                  ytrans = if (log_y) "log10" else "identity",
                  xlim = wlim, ylim = ylim,
                  y_ticks = y_ticks, legend_var = "Legend")
}

prepare_ArrayResourceBySize_plot_data <- function(x, wlim = c(NA, NA)) {
    w <- get_ArrayResourceBySize_w(x)
    value_name <- attr(x, "value_name") %||% "value"

    plot_dat <- data.frame(
        w = w,
        value = as.numeric(unclass(x)),
        Species = "Resource",
        Legend = "Resource",
        stringsAsFactors = FALSE
    )
    plot_dat <- apply_wlim(plot_dat, wlim)
    names(plot_dat)[2] <- value_name
    plot_dat
}

#' @rdname plotHover
#' @usage NULL
#' @examples
#' \donttest{
#' plotHover(getResourceMort(NS_params))
#' }
#' @export
plotHover.ArrayResourceBySize <- function(x, ...) {
    plotHover(plot(x, ...), ...)
}

#' @export
as.data.frame.ArrayResourceBySize <- function(x, row.names = NULL,
                                              optional = FALSE, ...) {
    w <- get_ArrayResourceBySize_w(x)
    data.frame(
        w = w,
        value = as.numeric(unclass(x)),
        row.names = row.names,
        check.names = !optional,
        stringsAsFactors = FALSE
    )
}

#' Get the size grid for an ArrayResourceBySize object
#'
#' Internal helper that returns the full prey/resource size grid
#' `params@w_full`, or the numeric vector parsed from the names of `x` if no
#' `params` is attached.
#'
#' @param x An `ArrayResourceBySize` object.
#'
#' @return A numeric vector giving the size represented by each element.
#' @keywords internal
get_ArrayResourceBySize_w <- function(x) {
    params <- attr(x, "params")
    if (!is.null(params) && length(x) == length(params@w_full)) {
        return(params@w_full)
    }
    w <- as.numeric(names(x))
    if (any(is.na(w))) {
        w <- seq_along(x)
    }
    w
}

#' @export
`[.ArrayResourceBySize` <- function(x, ...) {
    result <- NextMethod()
    attr(result, "value_name") <- attr(x, "value_name")
    attr(result, "units") <- attr(x, "units")
    attr(result, "params") <- attr(x, "params")
    class(result) <- c("ArrayResourceBySize", "numeric")
    result
}

#' @export
Ops.ArrayResourceBySize <- function(e1, e2) {
    # Strip ArrayResourceBySize class so that arithmetic returns a plain vector.
    if (is.ArrayResourceBySize(e1)) e1 <- unclass_resource(e1)
    if (!missing(e2) && is.ArrayResourceBySize(e2)) e2 <- unclass_resource(e2)
    op <- match.fun(.Generic)
    if (missing(e2)) op(e1) else op(e1, e2)
}

# Helper to strip all ArrayResourceBySize attributes
unclass_resource <- function(x) {
    x <- unclass(x)
    attr(x, "value_name") <- NULL
    attr(x, "units") <- NULL
    attr(x, "params") <- NULL
    x
}

#' @export
str.ArrayResourceBySize <- function(object, ...) {
    params <- attr(object, "params")
    attr(object, "params") <- NULL
    class(object) <- "numeric"
    out <- utils::capture.output(utils::str(object, ...))
    out[1] <- paste0(" 'ArrayResourceBySize' ", sub("^ ", "", out[1]))
    cat(paste0(out, collapse = "\n"), "\n", sep = "")
    if (!is.null(params)) {
        cat(" - attr(*, \"params\")=Formal class 'MizerParams' [package \"mizer\"] with ",
            length(slotNames(params)), " slots\n", sep = "")
    }
    invisible(NULL)
}


# ArrayTimeByResourceBySize ----------------------------------------------------

#' S3 class for time x resource-size arrays
#'
#' The [NResource()] function returns a two-dimensional array (time x size)
#' holding the resource number density through time. The
#' `ArrayTimeByResourceBySize` class wraps this array to provide convenient
#' `print()`, `summary()`, `plot()`, and `as.data.frame()` methods.
#'
#' An `ArrayTimeByResourceBySize` object behaves just like a regular matrix for
#' arithmetic operations and subsetting. It carries these lightweight
#' attributes:
#' \itemize{
#'   \item `value_name` – a human-readable name for the value
#'       (e.g. "Number density").
#'   \item `units` – the units of the value (e.g. "1/g").
#'   \item `params` – the `MizerParams` object that the value was computed from.
#' }
#'
#' @param x A matrix (time x size). For `is.ArrayTimeByResourceBySize()`, any
#'   object to test.
#' @param value_name A string giving the human-readable name for the value.
#' @param units A string giving the units (e.g. "1/g").
#' @param params A `MizerParams` object. Used for the resource colour and the
#'   size grid in the `plot()` method.
#'
#' @return An `ArrayTimeByResourceBySize` object (inherits from `matrix` and
#'   `array`).
#' @seealso [print()], [summary()], [as.data.frame()], [plot()]
#' @export
#' @examples
#' \donttest{
#' nr <- NResource(NS_sim)
#' is.ArrayTimeByResourceBySize(nr)
#' summary(nr)
#' plot(nr)
#' }
ArrayTimeByResourceBySize <- function(x, value_name = NULL, units = NULL,
                                      params = NULL) {
    if (!is.matrix(x)) {
        stop("`x` must be a matrix.")
    }
    structure(x,
        class = c("ArrayTimeByResourceBySize", "matrix", "array"),
        value_name = value_name,
        units = units,
        params = params
    )
}

#' @rdname ArrayTimeByResourceBySize
#' @return `is.ArrayTimeByResourceBySize()` returns `TRUE` if `x` is an
#'   `ArrayTimeByResourceBySize` object, `FALSE` otherwise.
#' @export
is.ArrayTimeByResourceBySize <- function(x) {
    inherits(x, "ArrayTimeByResourceBySize")
}

#' @export
print.ArrayTimeByResourceBySize <- function(x, ...) {
    value_name <- attr(x, "value_name") %||% "ArrayTimeByResourceBySize"
    units_str <- attr(x, "units")
    dims <- dim(x)
    header <- paste0(value_name, " (", dims[1], " times x ", dims[2],
                     " sizes)")
    if (!is.null(units_str) && nzchar(units_str)) {
        header <- paste0(header, " [", units_str, "]")
    }
    cat(header, "\n")

    mat <- unclass_resource(x)
    n_time <- nrow(mat)
    n_sizes <- ncol(mat)
    times <- parse_numeric_labels(rownames(mat), n_time)
    w <- parse_numeric_labels(colnames(mat), n_sizes)

    time_idx <- pick_log_spaced_indices(n_time, mizer_print_defaults$time_max,
                                        mizer_print_defaults$time_threshold)
    size_k <- fit_log_spaced_k(
        n_sizes, mizer_print_defaults$size_max, mizer_print_defaults$size_min,
        width_fn = function(k) {
            sz_idx <- pick_log_spaced_indices(n_sizes, k)
            matrix_display_width(mat[time_idx, sz_idx, drop = FALSE])
        })
    sz_idx <- pick_log_spaced_indices(n_sizes, size_k)

    print(mat[time_idx, sz_idx, drop = FALSE])

    if (length(time_idx) < n_time) {
        cat(format_truncation_note(length(time_idx), n_time, "times",
                                   format_time_range_detail(times)), "\n")
    }
    if (length(sz_idx) < n_sizes) {
        cat(format_truncation_note(length(sz_idx), n_sizes, "sizes",
                                   format_size_range_detail(w)), "\n")
    }
    invisible(x)
}

#' @export
summary.ArrayTimeByResourceBySize <- function(object, ...) {
    value_name <- attr(object, "value_name") %||% "ArrayTimeByResourceBySize"
    units_str <- attr(object, "units")
    vals <- unclass(object)

    result <- list(
        value_name = value_name,
        units = units_str,
        dims = dim(object),
        stats = data.frame(
            Min = min(vals, na.rm = TRUE),
            Mean = mean(vals, na.rm = TRUE),
            Max = max(vals, na.rm = TRUE),
            row.names = NULL
        )
    )
    class(result) <- "summary.ArrayTimeByResourceBySize"
    result
}

#' @export
print.summary.ArrayTimeByResourceBySize <- function(x, ...) {
    header <- x$value_name
    if (!is.null(x$units) && nzchar(x$units)) {
        header <- paste0(header, " [", x$units, "]")
    }
    cat(header, "\n")
    cat(x$dims[1], "times x", x$dims[2], "sizes\n\n")
    print(x$stats, row.names = FALSE)
    invisible(x)
}

#' @rdname plot
#' @usage NULL
#' @export
#' @examples
#' \donttest{
#' plot(NResource(NS_sim))
#' }
plot.ArrayTimeByResourceBySize <- function(x, time = NULL, ...) {
    slice <- ArrayTimeByResourceBySize_slice(x, time = time)
    plot.ArrayResourceBySize(slice, ...)
}

ArrayTimeByResourceBySize_slice <- function(x, time = NULL) {
    params <- attr(x, "params")
    value_name <- attr(x, "value_name")
    units <- attr(x, "units")

    times <- as.numeric(dimnames(x)[[1]])
    if (is.null(time)) {
        tidx <- dim(x)[1]
    } else {
        tidx <- which.min(abs(times - time))
    }

    vec <- unclass(x)[tidx, ]
    ArrayResourceBySize(vec, value_name = value_name,
                        units = units, params = params)
}

#' @rdname plotHover
#' @usage NULL
#' @examples
#' \donttest{
#' plotHover(NResource(NS_sim))
#' }
#' @export
plotHover.ArrayTimeByResourceBySize <- function(x, ...) {
    plotHover(plot(x, ...), ...)
}

#' @export
as.data.frame.ArrayTimeByResourceBySize <- function(x, row.names = NULL,
                                                    optional = FALSE, ...) {
    times <- as.numeric(dimnames(x)[[1]])
    if (any(is.na(times))) times <- seq_len(dim(x)[1])
    w <- as.numeric(dimnames(x)[[2]])
    if (any(is.na(w))) w <- seq_len(dim(x)[2])

    data.frame(
        expand.grid(time = times, w = w, stringsAsFactors = FALSE),
        value = c(unclass(x)),
        row.names = row.names,
        check.names = !optional
    )
}

#' @export
`[.ArrayTimeByResourceBySize` <- function(x, i, j, ..., drop = TRUE) {
    result <- NextMethod()
    if (is.matrix(result) && length(dim(result)) == 2) {
        attr(result, "value_name") <- attr(x, "value_name")
        attr(result, "units") <- attr(x, "units")
        attr(result, "params") <- attr(x, "params")
        class(result) <- c("ArrayTimeByResourceBySize", "matrix", "array")
    } else if (is.null(dim(result)) && !is.null(names(result)) &&
                   identical(names(result), colnames(x))) {
        # A single time step was selected, leaving a resource size spectrum.
        result <- ArrayResourceBySize(result,
                                      value_name = attr(x, "value_name"),
                                      units = attr(x, "units"),
                                      params = attr(x, "params"))
    }
    result
}

#' @export
Ops.ArrayTimeByResourceBySize <- function(e1, e2) {
    if (is.ArrayTimeByResourceBySize(e1)) e1 <- unclass_resource(e1)
    if (!missing(e2) && is.ArrayTimeByResourceBySize(e2)) {
        e2 <- unclass_resource(e2)
    }
    op <- match.fun(.Generic)
    if (missing(e2)) op(e1) else op(e1, e2)
}

#' @export
str.ArrayTimeByResourceBySize <- function(object, ...) {
    params <- attr(object, "params")
    attr(object, "params") <- NULL
    class(object) <- c("matrix", "array")
    out <- utils::capture.output(utils::str(object, ...))
    out[1] <- paste0(" 'ArrayTimeByResourceBySize' ", sub("^ ", "", out[1]))
    cat(paste0(out, collapse = "\n"), "\n", sep = "")
    if (!is.null(params)) {
        cat(" - attr(*, \"params\")=Formal class 'MizerParams' [package \"mizer\"] with ",
            length(slotNames(params)), " slots\n", sep = "")
    }
    invisible(NULL)
}
