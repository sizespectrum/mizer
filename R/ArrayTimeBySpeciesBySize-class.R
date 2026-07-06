# ArrayTimeBySpeciesBySize S3 class for time x species x size arrays
#
# Copyright 2026 Gustav Delius.
# Distributed under the GPL 3 or later.

#' S3 class for time x species x size arrays
#'
#' Some functions in mizer return three-dimensional arrays (time x species x
#' size) holding quantities like fishing mortality, feeding level, or predation
#' mortality through time. The `ArrayTimeBySpeciesBySize` class wraps these
#' arrays to provide convenient `print()`, `summary()`, `plot()`,
#' `animate()`, and `as.data.frame()` methods.
#'
#' An `ArrayTimeBySpeciesBySize` object behaves just like a regular array for
#' arithmetic operations and subsetting. It carries these lightweight attributes:
#' \itemize{
#'   \item `value_name` – a human-readable name for the value
#'       (e.g. "Fishing mortality").
#'   \item `units` – the units of the value (e.g. "1/year").
#'   \item `params` – the `MizerParams` object that the value was computed from.
#' }
#'
#' @param x A 3D array (time x species x size).
#' @param value_name A string giving the human-readable name for the value.
#' @param units A string giving the units (e.g. "1/year").
#' @param params A `MizerParams` object. Used for species colours, linetypes,
#'   and size ranges in the `plot()` and `animateSpectra()` methods.
#' @param representation Either `"point"` (the default) for a quantity sampled
#'   at the grid nodes, or `"average"` for a finite-volume bin average, which is
#'   then drawn at the geometric bin centre when the model uses second-order
#'   bin-averaging (`second_order_w[["bin_average"]]`). See
#'   [ArraySpeciesBySize()].
#'
#' @return An `ArrayTimeBySpeciesBySize` object (inherits from `array`).
#' @seealso [print()], [summary()], [as.data.frame()], [plot()],
#'   [animateSpectra()]
#' @export
#' @examples
#' \donttest{
#' fmort <- getFMort(NS_sim)
#' is.ArrayTimeBySpeciesBySize(fmort)
#' summary(fmort)
#' plot(fmort, time = 2007)
#' }
ArrayTimeBySpeciesBySize <- function(x, value_name = NULL, units = NULL,
                                     params = NULL,
                                     representation = c("point", "average")) {
    if (!is.array(x) || length(dim(x)) != 3) {
        stop("`x` must be a 3D array.")
    }
    representation <- match.arg(representation)
    structure(x,
        class = c("ArrayTimeBySpeciesBySize", "array"),
        value_name = value_name,
        units = units,
        params = params,
        representation = representation
    )
}

#' Get the size grid for an ArrayTimeBySpeciesBySize object
#'
#' Internal helper, the three-dimensional analogue of
#' [get_ArraySpeciesBySize_w()]. Returns the geometric bin centres (see
#' [bin_midpoints()]) when the array is tagged as a bin average and the model
#' uses second-order bin-averaging, otherwise the grid nodes read from the size
#' dimension names. Falls back to the dimension names when no `params` is
#' attached.
#'
#' @param x An `ArrayTimeBySpeciesBySize` object.
#' @return A numeric vector giving the size represented by each size slice.
#' @keywords internal
get_ArrayTimeBySpeciesBySize_w <- function(x) {
    w <- as.numeric(dimnames(x)[[3]])
    if (any(is.na(w))) w <- seq_len(dim(x)[3])
    params <- attr(x, "params")
    average <- !is.null(params) &&
        identical(attr(x, "representation"), "average") &&
        isTRUE(params@second_order_w[["bin_average"]])
    if (!average) return(w)
    if (dim(x)[3] == length(params@w)) return(bin_midpoints(params))
    if (dim(x)[3] == length(params@w_full)) return(bin_midpoints(params, w_full = TRUE))
    w
}

#' Test if an object is an ArrayTimeBySpeciesBySize
#'
#' @param x An object to test.
#' @return `TRUE` if `x` is an `ArrayTimeBySpeciesBySize` object, `FALSE`
#'   otherwise.
#' @export
#' @examples
#' is.ArrayTimeBySpeciesBySize(getFMort(NS_sim))
#' is.ArrayTimeBySpeciesBySize(array(1:8, dim = c(2, 2, 2)))
is.ArrayTimeBySpeciesBySize <- function(x) {
    inherits(x, "ArrayTimeBySpeciesBySize")
}

#' @export
print.ArrayTimeBySpeciesBySize <- function(x, ...) {
    value_name <- attr(x, "value_name") %||% "ArrayTimeBySpeciesBySize"
    units_str <- attr(x, "units")
    dims <- dim(x)
    header <- paste0(value_name, " (", dims[1], " times x ", dims[2],
                     " species x ", dims[3], " sizes)")
    if (!is.null(units_str) && nzchar(units_str)) {
        header <- paste0(header, " [", units_str, "]")
    }
    cat(header, "\n")

    time_labels <- dimnames(x)[[1]]
    shown_label <- if (!is.null(time_labels)) time_labels[dims[1]] else dims[1]
    cat("Showing final time step (t = ", shown_label, ") of ", dims[1],
        " time steps; use as.data.frame() or animate() for the full data.\n",
        sep = "")

    slice <- ArrayTimeBySpeciesBySize_slice(x)
    print_ArraySpeciesBySize_body(slice)
    invisible(x)
}

#' @export
summary.ArrayTimeBySpeciesBySize <- function(object, ...) {
    value_name <- attr(object, "value_name") %||% "ArrayTimeBySpeciesBySize"
    units_str <- attr(object, "units")
    sp_names <- dimnames(object)[[2]]
    arr <- unclass(object)

    df <- data.frame(
        Species = sp_names,
        Min = apply(arr, 2, min, na.rm = TRUE),
        Mean = apply(arr, 2, mean, na.rm = TRUE),
        Max = apply(arr, 2, max, na.rm = TRUE),
        row.names = NULL,
        stringsAsFactors = FALSE
    )

    result <- list(
        value_name = value_name,
        units = units_str,
        dims = dim(object),
        per_species = df
    )
    class(result) <- "summary.ArrayTimeBySpeciesBySize"
    result
}

#' @export
print.summary.ArrayTimeBySpeciesBySize <- function(x, ...) {
    header <- x$value_name
    if (!is.null(x$units) && nzchar(x$units)) {
        header <- paste0(header, " [", x$units, "]")
    }
    cat(header, "\n")
    cat(x$dims[1], "times x", x$dims[2], "species x", x$dims[3], "sizes\n\n")
    print(x$per_species, row.names = FALSE)
    invisible(x)
}

#' @rdname plot
#' @usage NULL
#' @export
#' @examples
#' \donttest{
#' plot(getFMort(NS_sim), time = 2010)
#' }
plot.ArrayTimeBySpeciesBySize <- function(x, species = NULL, time = NULL,
                                          all.sizes = FALSE, highlight = NULL,
                                          return_data = FALSE, log_x = TRUE,
                                          log_y = FALSE, log = NULL,
                                          wlim = c(NA, NA), llim = c(NA, NA),
                                          ylim = c(NA, NA),
                                          size_axis = c("w", "l"),
                                          total = FALSE, background = TRUE,
                                          y_ticks = 6, ...) {
    params <- attr(x, "params")
    value_name <- attr(x, "value_name")
    units <- attr(x, "units")

    times <- as.numeric(dimnames(x)[[1]])
    if (is.null(time)) {
        tidx <- dim(x)[1]
    } else {
        tidx <- which.min(abs(times - time))
    }

    arr <- unclass(x)
    slice <- matrix(arr[tidx, , , drop = FALSE],
                    nrow = dim(arr)[2],
                    dimnames = dimnames(arr)[2:3])
    slice <- ArraySpeciesBySize(slice, value_name = value_name,
                                units = units, params = params)

    plot.ArraySpeciesBySize(slice, species = species, all.sizes = all.sizes,
                            highlight = highlight, return_data = return_data,
                            log_x = log_x, log_y = log_y, log = log,
                            wlim = wlim, ylim = ylim, llim = llim,
                            size_axis = size_axis,
                            total = total, background = background,
                            y_ticks = y_ticks, ...)
}

#' @rdname plot2
#' @usage NULL
#' @export
plot2.ArrayTimeBySpeciesBySize <- function(x, y, name1 = "First",
                                           name2 = "Second",
                                           species = NULL,
                                           log_x = TRUE, log_y = FALSE,
                                           log = NULL,
                                           ylim = c(NA, NA),
                                           total = FALSE,
                                           background = TRUE,
                                           y_ticks = 6,
                                           time = NULL,
                                           all.sizes = FALSE,
                                           wlim = c(NA, NA),
                                           llim = c(NA, NA),
                                           size_axis = c("w", "l"), ...) {
    check_plot2_compatible(x, y, "ArrayTimeBySpeciesBySize")
    slice1 <- ArrayTimeBySpeciesBySize_slice(x, time = time)
    slice2 <- ArrayTimeBySpeciesBySize_slice(y, time = time)

    plot2.ArraySpeciesBySize(slice1, slice2, name1 = name1, name2 = name2,
                             species = species, all.sizes = all.sizes,
                             log_x = log_x, log_y = log_y, log = log,
                             wlim = wlim, ylim = ylim, llim = llim,
                             size_axis = size_axis, total = total,
                             background = background, y_ticks = y_ticks, ...)
}

#' @rdname plotRelative
#' @usage NULL
#' @export
plotRelative.ArrayTimeBySpeciesBySize <- function(x, y, species = NULL,
                                                  log_x = TRUE,
                                                  ylim = c(NA, NA),
                                                  total = FALSE,
                                                  background = TRUE,
                                                  time = NULL,
                                                  all.sizes = FALSE,
                                                  wlim = c(NA, NA),
                                                  llim = c(NA, NA),
                                                  size_axis = c("w", "l"), ...) {
    check_plot2_compatible(x, y, "ArrayTimeBySpeciesBySize")
    slice1 <- ArrayTimeBySpeciesBySize_slice(x, time = time)
    slice2 <- ArrayTimeBySpeciesBySize_slice(y, time = time)

    plotRelative.ArraySpeciesBySize(slice1, slice2, species = species,
                                    all.sizes = all.sizes, log_x = log_x,
                                    wlim = wlim, ylim = ylim, llim = llim,
                                    size_axis = size_axis,
                                    total = total, background = background,
                                    ...)
}

ArrayTimeBySpeciesBySize_slice <- function(x, time = NULL) {
    params <- attr(x, "params")
    value_name <- attr(x, "value_name")
    units <- attr(x, "units")
    representation <- attr(x, "representation") %||% "point"

    times <- as.numeric(dimnames(x)[[1]])
    if (is.null(time)) {
        tidx <- dim(x)[1]
    } else {
        tidx <- which.min(abs(times - time))
    }

    arr <- unclass(x)
    slice <- matrix(arr[tidx, , , drop = FALSE],
                    nrow = dim(arr)[2],
                    dimnames = dimnames(arr)[2:3])
    ArraySpeciesBySize(slice, value_name = value_name,
                       units = units, params = params,
                       representation = representation)
}

#' @rdname plotHover
#' @usage NULL
#' @examples
#' \donttest{
#' plotHover(getFMort(NS_sim))
#' }
#' @export
plotHover.ArrayTimeBySpeciesBySize <- function(x, ...) {
    plotHover(plot(x, ...), ...)
}

#' @rdname animate
#' @usage NULL
#' @export
animate.ArrayTimeBySpeciesBySize <- function(x, species = NULL,
                                             log_x = TRUE,
                                             log_y = TRUE,
                                             log = NULL,
                                             wlim = c(NA, NA),
                                             llim = c(NA, NA),
                                             ylim = c(NA, NA),
                                             tlim = c(NA, NA),
                                             size_axis = c("w", "l"),
                                             total = FALSE,
                                             background = TRUE,
                                             frame_duration = 500,
                                             transition_duration = frame_duration,
                                             easing = "linear",
                                             ...) {
    assert_that(is.flag(total), is.flag(background),
                is.number(frame_duration), frame_duration >= 0,
                is.number(transition_duration), transition_duration >= 0,
                is.string(easing),
                length(wlim) == 2, length(llim) == 2, length(ylim) == 2)
    size_axis <- plot_size_axis(size_axis)
    log_axes <- parsePlotLog(log, log_x = log_x, log_y = log_y)
    log_x <- log_axes$log_x
    log_y <- log_axes$log_y

    params <- attr(x, "params")
    value_name <- attr(x, "value_name") %||% "Value"
    units_str <- attr(x, "units")

    all_species <- dimnames(x)[[2]]

    bkgrd_sp <- character(0)
    if (!is.null(params) && isTRUE(any(params@species_params$is_background))) {
        bkgrd_sp <- params@species_params$species[params@species_params$is_background]
    }

    if (is.null(species)) {
        species <- all_species
    } else {
        species <- intersect(species, all_species)
        if (length(species) == 0) {
            stop("None of the selected species are in the array.")
        }
    }

    times <- as.numeric(dimnames(x)[[1]])
    arr <- unclass(x)
    if (!is.na(tlim[1])) {
        arr <- arr[times >= tlim[1], , , drop = FALSE]
        times <- times[times >= tlim[1]]
    }
    if (!is.na(tlim[2])) {
        arr <- arr[times <= tlim[2], , , drop = FALSE]
        times <- times[times <= tlim[2]]
    }

    w <- get_ArrayTimeBySpeciesBySize_w(x)
    sub <- arr[, species, , drop = FALSE]

    # Build long-format data frame; time varies fastest to match c(sub)
    df <- expand.grid(time = times, Species = species, w = w,
                      stringsAsFactors = FALSE)
    df$value <- c(sub)

    # Compute total across ALL selected species (including background) before
    # any background filtering, matching the behaviour of plot.ArraySpeciesBySize
    if (total) {
        total_sums <- stats::aggregate(value ~ time + w, data = df, FUN = sum,
                                       na.rm = TRUE)
        total_sums$Species <- "Total"
    }

    # Now handle background: exclude rows or group under "Background" legend
    df$legend_name <- df$Species
    if (length(bkgrd_sp) > 0) {
        if (background) {
            df$legend_name[df$Species %in% bkgrd_sp] <- "Background"
        } else {
            df <- df[!df$Species %in% bkgrd_sp, ]
        }
    }

    if (total) {
        total_sums$legend_name <- "Total"
        df <- rbind(df, total_sums[, names(df)])
    }

    y_label <- value_name
    if (!is.null(units_str) && nzchar(units_str)) {
        y_label <- paste0(value_name, " [", units_str, "]")
    }

    animate_plotly(df, params, log_x, log_y, y_label, wlim, llim,
                   ylim,
                   size_axis = size_axis,
                   frame_duration = frame_duration,
                   transition_duration = transition_duration,
                   easing = easing)
}

#' @export
as.data.frame.ArrayTimeBySpeciesBySize <- function(x, row.names = NULL,
                                                   optional = FALSE, ...) {
    times <- as.numeric(dimnames(x)[[1]])
    if (any(is.na(times))) times <- seq_len(dim(x)[1])
    sp_names <- dimnames(x)[[2]]
    w <- get_ArrayTimeBySpeciesBySize_w(x)

    data.frame(
        expand.grid(time = times, Species = sp_names, w = w,
                    stringsAsFactors = FALSE),
        value = c(unclass(x)),
        row.names = row.names,
        check.names = !optional
    )
}

#' @export
`[.ArrayTimeBySpeciesBySize` <- function(x, i, j, k, ..., drop = TRUE) {
    result <- NextMethod()
    if (is.array(result) && length(dim(result)) == 3) {
        attr(result, "value_name") <- attr(x, "value_name")
        attr(result, "units") <- attr(x, "units")
        attr(result, "params") <- attr(x, "params")
        attr(result, "representation") <- attr(x, "representation")
        class(result) <- c("ArrayTimeBySpeciesBySize", "array")
    } else if (is.matrix(result)) {
        dim_names <- names(dimnames(result))
        attrs <- list(value_name = attr(x, "value_name"),
                      units = attr(x, "units"),
                      params = attr(x, "params"),
                      representation = attr(x, "representation") %||% "point")
        if (identical(dim_names, c("sp", "w"))) {
            result <- ArraySpeciesBySize(result,
                                         value_name = attrs$value_name,
                                         units = attrs$units,
                                         params = attrs$params,
                                         representation = attrs$representation)
        } else if (identical(dim_names, c("time", "sp"))) {
            result <- ArrayTimeBySpecies(result,
                                         value_name = attrs$value_name,
                                         units = attrs$units,
                                         params = attrs$params)
        }
    }
    result
}

#' @export
Ops.ArrayTimeBySpeciesBySize <- function(e1, e2) {
    if (is.ArrayTimeBySpeciesBySize(e1)) e1 <- unclass_tss(e1)
    if (!missing(e2) && is.ArrayTimeBySpeciesBySize(e2)) e2 <- unclass_tss(e2)
    op <- match.fun(.Generic)
    if (missing(e2)) op(e1) else op(e1, e2)
}

# Helper to strip all ArrayTimeBySpeciesBySize attributes
unclass_tss <- function(x) {
    x <- unclass(x)
    attr(x, "value_name") <- NULL
    attr(x, "units") <- NULL
    attr(x, "params") <- NULL
    attr(x, "representation") <- NULL
    x
}

#' @export
str.ArrayTimeBySpeciesBySize <- function(object, ...) {
    params <- attr(object, "params")
    attr(object, "params") <- NULL
    class(object) <- "array"
    out <- utils::capture.output(utils::str(object, ...))
    out[1] <- paste0(" 'ArrayTimeBySpeciesBySize' ", sub("^ ", "", out[1]))
    cat(paste0(out, collapse = "\n"), "\n", sep = "")
    if (!is.null(params)) {
        cat(" - attr(*, \"params\")=Formal class 'MizerParams' [package \"mizer\"] with ",
            length(slotNames(params)), " slots\n", sep = "")
    }
    invisible(NULL)
}
