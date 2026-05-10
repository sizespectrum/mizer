# ArrayTimeBySpeciesBySize S3 class for time x species x size arrays
#
# Copyright 2026 Gustav Delius.
# Distributed under the GPL 3 or later.

#' S3 class for time x species x size rate arrays
#'
#' Some functions in mizer return three-dimensional arrays (time x species x
#' size) holding rates like fishing mortality, feeding level, or predation
#' mortality through time. The `ArrayTimeBySpeciesBySize` class wraps these
#' arrays to provide convenient `print()`, `summary()`, `plot()`,
#' `animateSpectra()`, and `as.data.frame()` methods.
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
#'
#' @return An `ArrayTimeBySpeciesBySize` object (inherits from `array`).
#' @export
#' @examples
#' \donttest{
#' fmort <- getFMort(NS_sim)
#' is.ArrayTimeBySpeciesBySize(fmort)
#' summary(fmort)
#' plot(fmort, time = 2007)
#' }
ArrayTimeBySpeciesBySize <- function(x, value_name = NULL, units = NULL,
                                     params = NULL) {
    if (!is.array(x) || length(dim(x)) != 3) {
        stop("`x` must be a 3D array.")
    }
    structure(x,
        class = c("ArrayTimeBySpeciesBySize", "array"),
        value_name = value_name,
        units = units,
        params = params
    )
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
    sp_names <- dimnames(x)[[2]]
    if (!is.null(sp_names)) {
        arr <- unclass(x)
        vals <- apply(arr, 2, function(col) {
            col <- col[is.finite(col)]
            if (length(col) == 0) return("all NA/Inf")
            paste0("min=", signif(min(col), 3),
                   " mean=", signif(mean(col), 3),
                   " max=", signif(max(col), 3))
        })
        for (i in seq_along(sp_names)) {
            cat("  ", sp_names[i], ": ", vals[i], "\n", sep = "")
        }
    }
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
#'
#' @param time The time to display. Default (`NULL`) is the final time step.
#'   Only applies to `ArrayTimeBySpeciesBySize`.
#' @export
#' @examples
#' \donttest{
#' plot(getFMort(NS_sim))
#' plot(getFMort(NS_sim), time = 2010)
#' }
plot.ArrayTimeBySpeciesBySize <- function(x, species = NULL, time = NULL,
                                          all.sizes = FALSE, highlight = NULL,
                                          return_data = FALSE, log_x = TRUE,
                                          log_y = FALSE,
                                          wlim = c(NA, NA), ylim = c(NA, NA),
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
    slice <- arr[tidx, , ]
    slice <- ArraySpeciesBySize(slice, value_name = value_name,
                                units = units, params = params)

    plot.ArraySpeciesBySize(slice, species = species, all.sizes = all.sizes,
                            highlight = highlight, return_data = return_data,
                            log_x = log_x, log_y = log_y, wlim = wlim,
                            ylim = ylim, total = total, background = background,
                            y_ticks = y_ticks, ...)
}

#' @rdname plot
#' @exportS3Method plotly::ggplotly
#' @examples
#' \donttest{
#' ggplotly(getFMort(NS_sim))
#' }
ggplotly.ArrayTimeBySpeciesBySize <- function(x, ...) {
    ggplotly(plot(x, ...))
}

#' @rdname animate
#' @export
animate.ArrayTimeBySpeciesBySize <- function(x, species = NULL,
                                             time_range = NULL,
                                             log_x = TRUE,
                                             wlim = c(NA, NA),
                                             ylim = c(NA, NA),
                                             log_y = TRUE,
                                             total = FALSE,
                                             background = TRUE,
                                             interpolate = TRUE,
                                             ...) {
    assert_that(is.flag(total), is.flag(background), is.flag(interpolate),
                length(wlim) == 2, length(ylim) == 2)
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
    if (!is.null(time_range)) {
        if (length(time_range) == 2 && !all(time_range %in% times)) {
            keep <- times >= time_range[1] & times <= time_range[2]
        } else {
            keep <- times %in% time_range
        }
        arr <- arr[keep, , , drop = FALSE]
        times <- times[keep]
    }

    w <- as.numeric(dimnames(arr)[[3]])
    sub <- arr[, species, , drop = FALSE]

    # Build long-format data frame; time varies fastest to match c(sub)
    df <- expand.grid(time = times, Species = species, w = w,
                      stringsAsFactors = FALSE)
    df$value <- c(sub)

    # Compute total across ALL selected species (including background) before
    # any background filtering, matching the behaviour of plot.ArraySpeciesBySize
    if (total) {
        total_sums <- aggregate(value ~ time + w, data = df, FUN = sum,
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

    animate_plotly(df, params, log_x, log_y, y_label, wlim, ylim, interpolate)
}

#' @export
as.data.frame.ArrayTimeBySpeciesBySize <- function(x, row.names = NULL,
                                                   optional = FALSE, ...) {
    times <- as.numeric(dimnames(x)[[1]])
    if (any(is.na(times))) times <- seq_len(dim(x)[1])
    sp_names <- dimnames(x)[[2]]
    w <- as.numeric(dimnames(x)[[3]])
    if (any(is.na(w))) w <- seq_len(dim(x)[3])

    df <- expand.grid(time = times, Species = sp_names, w = w,
                      stringsAsFactors = FALSE)
    df$value <- c(unclass(x))
    df
}

#' @export
`[.ArrayTimeBySpeciesBySize` <- function(x, i, j, k, ..., drop = TRUE) {
    result <- NextMethod()
    if (is.array(result) && length(dim(result)) == 3) {
        attr(result, "value_name") <- attr(x, "value_name")
        attr(result, "units") <- attr(x, "units")
        attr(result, "params") <- attr(x, "params")
        class(result) <- c("ArrayTimeBySpeciesBySize", "array")
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
    x
}
