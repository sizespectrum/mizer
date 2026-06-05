# ArraySpeciesBySize S3 class for species x size arrays
#
# Copyright 2026 Gustav Delius.
# Distributed under the GPL 3 or later.

#' S3 class for species x size rate arrays
#'
#' Many functions in mizer return two-dimensional arrays (species x size)
#' holding rates like encounter rate, feeding level, growth rate, mortality etc.
#' The `ArraySpeciesBySize` class wraps these arrays to provide convenient
#' `print()`, `summary()`, `plot()`, and `as.data.frame()` methods.
#'
#' An `ArraySpeciesBySize` object behaves just like a regular matrix for
#' arithmetic operations and subsetting. It carries two lightweight attributes:
#' \itemize{
#'   \item `value_name` – a human-readable name for the value
#'       (e.g. "Encounter rate").
#'   \item `units` – the units of the rate (e.g. "g/year").
#' }
#'
#' @param x A matrix (species x size).
#' @param value_name A string giving the human-readable name for the value.
#' @param units A string giving the units (e.g. "g/year", "1/year").
#' @param params A `MizerParams` object. Used for species colours, linetypes,
#'   and size ranges in the `plot()` method.
#'
#' @return An `ArraySpeciesBySize` object (inherits from `matrix` and `array`).
#' @seealso [print()], [summary()], [as.data.frame()], [plot()]
#' @export
#' @examples
#' \donttest{
#' enc <- getEncounter(NS_params)
#' is.ArraySpeciesBySize(enc)
#' summary(enc)
#' }
ArraySpeciesBySize <- function(x, value_name = NULL, units = NULL,
                               params = NULL) {
    if (!is.matrix(x)) {
        stop("`x` must be a matrix.")
    }
    if (!is.null(params) && identical(dim(x), dim(params@metab))) {
        dimnames(x) <- dimnames(params@metab)
    }
    structure(x,
        class = c("ArraySpeciesBySize", "matrix", "array"),
        value_name = value_name,
        units = units,
        params = params
    )
}

#' Test if an object is a ArraySpeciesBySize
#'
#' @param x An object to test.
#' @return `TRUE` if `x` is an `ArraySpeciesBySize` object, `FALSE` otherwise.
#' @export
#' @examples
#' is.ArraySpeciesBySize(getEncounter(NS_params))
#' is.ArraySpeciesBySize(matrix(1:4, nrow = 2))
is.ArraySpeciesBySize <- function(x) {
    inherits(x, "ArraySpeciesBySize")
}

#' @export
print.ArraySpeciesBySize <- function(x, ...) {
    value_name <- attr(x, "value_name") %||% "ArraySpeciesBySize"
    units_str <- attr(x, "units")
    dims <- dim(x)
    header <- paste0(value_name, " (", dims[1], " species x ", dims[2],
                     " sizes)")
    if (!is.null(units_str)) {
        header <- paste0(header, " [", units_str, "]")
    }
    cat(header, "\n")
    # Print a compact summary per species
    sp_names <- rownames(x)
    if (!is.null(sp_names)) {
        vals <- apply(unclass(x), 1, function(row) {
            row <- row[is.finite(row)]
            if (length(row) == 0) return("all NA/Inf")
            paste0("min=", signif(min(row), 3),
                   " mean=", signif(mean(row), 3),
                   " max=", signif(max(row), 3))
        })
        for (i in seq_along(sp_names)) {
            cat("  ", sp_names[i], ": ", vals[i], "\n", sep = "")
        }
    }
    invisible(x)
}

#' @export
summary.ArraySpeciesBySize <- function(object, ...) {
    value_name <- attr(object, "value_name") %||% "ArraySpeciesBySize"
    units_str <- attr(object, "units")
    sp_names <- rownames(object)
    mat <- unclass(object)

    df <- data.frame(
        Species = sp_names,
        Min = apply(mat, 1, min, na.rm = TRUE),
        Mean = apply(mat, 1, mean, na.rm = TRUE),
        Max = apply(mat, 1, max, na.rm = TRUE),
        row.names = NULL,
        stringsAsFactors = FALSE
    )

    result <- list(
        value_name = value_name,
        units = units_str,
        dims = dim(object),
        per_species = df
    )
    class(result) <- "summary.ArraySpeciesBySize"
    result
}

#' @export
print.summary.ArraySpeciesBySize <- function(x, ...) {
    header <- x$value_name
    if (!is.null(x$units)) {
        header <- paste0(header, " [", x$units, "]")
    }
    cat(header, "\n")
    cat(x$dims[1], "species x", x$dims[2], "sizes\n\n")
    print(x$per_species, row.names = FALSE)
    invisible(x)
}

#' Plot mizer arrays
#'
#' Many mizer functions return values that depend on species and either size or
#' time. `plot()` creates a ggplot2 figure with one line for each species
#' showing the values against size or against time (depending on the type of
#' output). [plotHover()] creates an interactive version of the same figure.
#'
#' This works because the mizer functions that give values that depend on
#' species and size return an `ArraySpeciesBySize` object and those that
#' give values that depend on species and time return an `ArrayTimeBySpecies`
#' object. These objects have attributes that store the name of the value,
#' its units, and a reference to the `MizerParams` object that the value was
#' computed from. This allows the plots to be automatically labelled and
#' coloured appropriately.
#'
#' To compare two mizer arrays in a single plot, use [plot2()]. To show the
#' relative difference between two arrays, use [plotRelative()].
#'
#' @param x An `ArraySpeciesBySize`, `ArrayTimeBySpecies`, or
#'   `ArrayTimeBySpeciesBySize` object.
#' @param ...
#'   **Arguments used by all methods:**
#'   \describe{
#'     \item{`species`}{Character vector of species to include. `NULL`
#'       (default) means all species.}
#'     \item{`highlight`}{Name or vector of names of the species to be
#'       highlighted.}
#'     \item{`total`}{A boolean value that determines whether the total over
#'       all selected species is plotted as well. Default is `FALSE`.}
#'     \item{`background`}{A boolean value that determines whether background
#'       species are included. Ignored if the model does not contain background
#'       species. Default is `TRUE`.}
#'     \item{`return_data`}{If `TRUE`, return the data frame instead of the
#'       plot.}
#'     \item{`log_x`}{If `TRUE`, use a log10 x-axis. Default is `TRUE` for size
#'       spectra and `FALSE` for time series.}
#'     \item{`log_y`}{If `TRUE`, use a log10 y-axis. Default is `FALSE` for
#'       `ArraySpeciesBySize` and `TRUE` for `ArrayTimeBySpecies`.}
#'     \item{`log`}{Character string specifying which axes should use log10
#'       scales, in the same form as the base [plot()] argument. For example,
#'       `"x"`, `"y"`, `"xy"` or `""`. If supplied, this overrides `log_x` and
#'       `log_y`.}
#'     \item{`ylim`}{A numeric vector of length two providing lower and upper
#'       limits for the value (y) axis. Use `NA` to refer to the existing
#'       minimum or maximum.}
#'     \item{`y_ticks`}{The approximate number of ticks desired on the y axis.}
#'   }
#'
#'   **For `ArraySpeciesBySize` and `ArrayTimeBySpeciesBySize` methods:**
#'   \describe{
#'     \item{`all.sizes`}{If `FALSE` (default), values outside a species' size
#'       range (`w_min` to `w_max`) are removed.}
#'     \item{`wlim`}{A numeric vector of length two providing lower and upper
#'       limits for the weight (x) axis. Use `NA` to refer to the existing
#'       minimum or maximum.}
#'     \item{`llim`}{A numeric vector of length two providing lower and upper
#'       limits for the length (x) axis when `size_axis = "l"`. Use `NA` to
#'       refer to the existing minimum or maximum.}
#'     \item{`size_axis`}{Whether to plot size as weight (`"w"`, default) or
#'       length (`"l"`), using the allometric weight-length relationship.}
#'   }
#'
#'   **For `ArrayTimeBySpecies` methods:**
#'   \describe{
#'     \item{`tlim`}{A numeric vector of length two providing lower and upper
#'       limits for the time axis, e.g. `c(1980, 2000)`. Use `NA` to apply no
#'       limit at that end. Default is `c(NA, NA)`.}
#'   }
#'
#'   **For `ArrayTimeBySpeciesBySize` methods:**
#'   \describe{
#'     \item{`time`}{The time to display. Default (`NULL`) is the final time
#'       step.}
#'   }
#'
#' @return A ggplot2 object, unless `return_data = TRUE`, in which case a data
#'   frame is returned. [plotHover()] returns a plotly object.
#'
#' @name plot
#' @family plotting functions
#' @usage NULL
#' @export
#' @examples
#' \donttest{
#' plot(getEncounter(NS_params))
#' plot(getFeedingLevel(NS_params), species = c("Cod", "Herring"))
#' plot(getPredMort(NS_params), species = c("Cod", "Herring"),
#'      size_axis = "l")
#' }
plot.ArraySpeciesBySize <- function(x, species = NULL,
                            all.sizes = FALSE, highlight = NULL,
                            return_data = FALSE, log_x = TRUE, log_y = FALSE,
                            log = NULL,
                            wlim = c(NA, NA), llim = c(NA, NA),
                            ylim = c(NA, NA),
                            size_axis = c("w", "l"),
                            total = FALSE, background = TRUE,
                            y_ticks = 6, ...) {
    size_axis <- plot_size_axis(size_axis)
    log_axes <- parsePlotLog(log, log_x = log_x, log_y = log_y)
    log_x <- log_axes$log_x
    log_y <- log_axes$log_y

    assert_that(length(wlim) == 2,
                length(llim) == 2,
                length(ylim) == 2)
    value_name <- attr(x, "value_name") %||% "Rate"
    units_str <- attr(x, "units")
    params <- attr(x, "params")
    plot_dat <- prepare_ArraySpeciesBySize_plot_data(
        x, species = species, all.sizes = all.sizes, wlim = wlim,
        total = total, background = background)
    plot_dat <- convert_plot_size_axis(plot_dat, params, size_axis)
    if (identical(size_axis, "l")) {
        plot_dat <- filter_plot_length_limits(plot_dat, llim)
    }

    if (return_data) return(plot_dat)

    y_label <- value_name
    if (!is.null(units_str)) {
        y_label <- paste0(value_name, " [", units_str, "]")
    }

    plotDataFrame(plot_dat, params, xlab = plot_size_xlab(size_axis),
                  ylab = y_label,
                  xtrans = if (log_x) "log10" else "identity",
                  ytrans = if (log_y) "log10" else "identity",
                  xlim = plot_size_xlim(wlim, size_axis, llim), ylim = ylim,
                  highlight = highlight, y_ticks = y_ticks,
                  legend_var = "Legend")
}

parsePlotLog <- function(log, log_x = FALSE, log_y = FALSE) {
    if (is.null(log)) {
        return(list(log_x = log_x, log_y = log_y))
    }
    # Backward compatibility: legacy logical `log` toggles only the y-axis.
    if (is.logical(log)) {
        if (length(log) != 1 || is.na(log)) {
            stop("`log` must be a single logical value or a character string ",
                 "containing only \"x\" and/or \"y\".")
        }
        return(list(log_x = FALSE, log_y = isTRUE(log)))
    }
    if (!is.character(log) || length(log) != 1 || is.na(log) ||
        grepl("[^xy]", log)) {
        stop("`log` must be a single logical value or a character string ",
             "containing only \"x\" and/or \"y\".")
    }
    list(
        log_x = grepl("x", log, fixed = TRUE),
        log_y = grepl("y", log, fixed = TRUE)
    )
}

#' Compare two mizer arrays in a single plot
#'
#' `plot2()` compares two compatible mizer array objects in a single ggplot.
#' Colours identify species or groups, and linetype identifies which object
#' the values came from.
#'
#' @param x The first of two compatible mizer array objects to compare.
#'   Can be an `ArraySpeciesBySize`, `ArrayTimeBySpecies`, or
#'   `ArrayTimeBySpeciesBySize` object.
#' @param y The second mizer array object, compatible with `x`.
#' @param name1,name2 Labels for the two objects, used in the linetype legend.
#' @param species Character vector of species to include. `NULL` (default) means
#'   all species.
#' @param log_x If `TRUE`, use a log10 x-axis. Default is `TRUE` for size
#'   spectra and `FALSE` for time series.
#' @param log_y If `TRUE`, use a log10 y-axis. Default is `FALSE` for
#'   `ArraySpeciesBySize` and `TRUE` for `ArrayTimeBySpecies`.
#' @param log Character string specifying which axes should use log10 scales,
#'   in the same form as the base [plot()] argument. For example, `"x"`,
#'   `"y"`, `"xy"` or `""`. If supplied, this overrides `log_x` and `log_y`.
#' @param ylim A numeric vector of length two providing lower and upper limits
#'   for the value (y) axis. Use `NA` to refer to the existing minimum or
#'   maximum.
#' @param total A boolean value that determines whether the total over all
#'   selected species is plotted as well. Default is `FALSE`.
#' @param background A boolean value that determines whether background species
#'   are included. Ignored if the model does not contain background species.
#'   Default is `TRUE`.
#' @param y_ticks The approximate number of ticks desired on the y axis.
#' @param ... Further arguments used by only some of the methods:
#'
#'   **For `ArraySpeciesBySize` and `ArrayTimeBySpeciesBySize` methods:**
#'   \describe{
#'     \item{`all.sizes`}{If `FALSE` (default), values outside a species' size
#'       range (`w_min` to `w_max`) are removed.}
#'     \item{`wlim`}{A numeric vector of length two providing lower and upper
#'       limits for the weight (x) axis. Use `NA` to refer to the existing
#'       minimum or maximum.}
#'     \item{`llim`}{A numeric vector of length two providing lower and upper
#'       limits for the length (x) axis when `size_axis = "l"`. Use `NA` to
#'       refer to the existing minimum or maximum.}
#'     \item{`size_axis`}{Whether to plot size as weight (`"w"`, default) or
#'       length (`"l"`), using the allometric weight-length relationship.}
#'   }
#'
#'   **For `ArrayTimeBySpecies` methods:**
#'   \describe{
#'     \item{`tlim`}{A numeric vector of length two providing lower and upper
#'       limits for the time axis, e.g. `c(1980, 2000)`. Use `NA` to apply no
#'       limit at that end. Default is `c(NA, NA)`.}
#'   }
#'
#'   **For `ArrayTimeBySpeciesBySize` methods:**
#'   \describe{
#'     \item{`time`}{The time to display. Default (`NULL`) is the final time
#'       step.}
#'   }
#'
#' @return A ggplot2 object.
#'
#' @family plotting functions
#' @export
#' @examples
#' \donttest{
#' plot2(getEncounter(NS_params), getEncounter(NS_params))
#' }
plot2 <- function(x, y, name1 = "First", name2 = "Second",
                  species = NULL, log_x, log_y, log = NULL,
                  ylim = c(NA, NA), total = FALSE, background = TRUE,
                  y_ticks = 6, ...) {
    UseMethod("plot2", x)
}

#' @rdname plot2
#' @usage NULL
#' @export
plot2.ArraySpeciesBySize <- function(x, y, name1 = "First", name2 = "Second",
                                     species = NULL,
                                     log_x = TRUE, log_y = FALSE, log = NULL,
                                     ylim = c(NA, NA),
                                     total = FALSE, background = TRUE,
                                     y_ticks = 6,
                                     all.sizes = FALSE,
                                     wlim = c(NA, NA), llim = c(NA, NA),
                                     size_axis = c("w", "l"), ...) {
    check_plot2_compatible(x, y, "ArraySpeciesBySize")
    compare_array_metadata(x, y)
    size_axis <- plot_size_axis(size_axis)
    log_axes <- parsePlotLog(log, log_x = log_x, log_y = log_y)
    log_x <- log_axes$log_x
    log_y <- log_axes$log_y
    assert_that(length(wlim) == 2,
                length(llim) == 2,
                length(ylim) == 2)

    params <- attr(x, "params")
    y_label <- array_y_label(x, default = "Rate")
    plot_dat1 <- prepare_ArraySpeciesBySize_plot_data(
        x, species = species, all.sizes = all.sizes, wlim = wlim,
        total = total, background = background)
    plot_dat2 <- prepare_ArraySpeciesBySize_plot_data(
        y, species = species, all.sizes = all.sizes, wlim = wlim,
        total = total, background = background)

    plotComparisonDataFrame(plot_dat1, plot_dat2, params,
                            name1 = name1, name2 = name2,
                            xlab = plot_size_xlab(size_axis), ylab = y_label,
                            xtrans = if (log_x) "log10" else "identity",
                            ytrans = if (log_y) "log10" else "identity",
                            xlim = plot_size_xlim(wlim, size_axis, llim),
                            ylim = ylim,
                            y_ticks = y_ticks, legend_var = "Legend",
                            size_axis = size_axis)
}

#' Plot relative difference between two mizer arrays
#'
#' `plotRelative()` plots the difference between two compatible mizer array
#' objects relative to their average. If the values in the first object are
#' \eqn{N_1} and the values in the second are \eqn{N_2}, it plots
#' \deqn{2 (N_2 - N_1) / (N_1 + N_2).}
#'
#' @param x The first of two compatible mizer array objects to compare.
#'   Can be an `ArraySpeciesBySize`, `ArrayTimeBySpecies`, or
#'   `ArrayTimeBySpeciesBySize` object.
#' @param y The second mizer array object, compatible with `x`.
#' @param species Character vector of species to include. `NULL` (default) means
#'   all species.
#' @param log_x If `TRUE`, use a log10 x-axis. Default is `TRUE` for size
#'   spectra and `FALSE` for time series.
#' @param ylim A numeric vector of length two providing lower and upper limits
#'   for the value (y) axis.
#' @param total A boolean value that determines whether the total over all
#'   selected species is plotted as well. Default is `FALSE`.
#' @param background A boolean value that determines whether background species
#'   are included. Ignored if the model does not contain background species.
#'   Default is `TRUE`.
#' @param ... Further arguments used by only some of the methods:
#'
#'   **For `ArraySpeciesBySize` and `ArrayTimeBySpeciesBySize` methods:**
#'   \describe{
#'     \item{`all.sizes`}{If `FALSE` (default), values outside a species' size
#'       range (`w_min` to `w_max`) are removed.}
#'     \item{`wlim`}{A numeric vector of length two providing lower and upper
#'       limits for the weight (x) axis. Use `NA` to refer to the existing
#'       minimum or maximum.}
#'     \item{`llim`}{A numeric vector of length two providing lower and upper
#'       limits for the length (x) axis when `size_axis = "l"`. Use `NA` to
#'       refer to the existing minimum or maximum.}
#'     \item{`size_axis`}{Whether to plot size as weight (`"w"`, default) or
#'       length (`"l"`), using the allometric weight-length relationship.}
#'   }
#'
#'   **For `ArrayTimeBySpecies` methods:**
#'   \describe{
#'     \item{`tlim`}{A numeric vector of length two providing lower and upper
#'       limits for the time axis, e.g. `c(1980, 2000)`. Use `NA` to apply no
#'       limit at that end. Default is `c(NA, NA)`.}
#'   }
#'
#'   **For `ArrayTimeBySpeciesBySize` methods:**
#'   \describe{
#'     \item{`time`}{The time to display. Default (`NULL`) is the final time
#'       step.}
#'   }
#'
#' @return A ggplot2 object.
#'
#' @family plotting functions
#' @export
#' @examples
#' \donttest{
#' params <- NS_params
#' given_species_params(params)["Cod", "w_mat"] <- 1200
#' plotRelative(getEGrowth(NS_params), getEGrowth(params),
#'              wlim = c(500, 2000), log_x = FALSE, species = "Cod")
#' }
plotRelative <- function(x, y, species = NULL, log_x,
                         ylim = c(NA, NA), total = FALSE,
                         background = TRUE, ...) {
    UseMethod("plotRelative", x)
}

#' @rdname plotRelative
#' @usage NULL
#' @export
plotRelative.ArraySpeciesBySize <- function(x, y, species = NULL,
                                            log_x = TRUE,
                                            ylim = c(NA, NA),
                                            total = FALSE,
                                            background = TRUE,
                                            all.sizes = FALSE,
                                            wlim = c(NA, NA),
                                            llim = c(NA, NA),
                                            size_axis = c("w", "l"), ...) {
    check_plot2_compatible(x, y, "ArraySpeciesBySize")
    compare_array_metadata(x, y)
    size_axis <- plot_size_axis(size_axis)
    assert_that(length(wlim) == 2,
                length(llim) == 2,
                length(ylim) == 2)
    params <- attr(x, "params")
    plot_dat1 <- prepare_ArraySpeciesBySize_plot_data(
        x, species = species, all.sizes = all.sizes, wlim = wlim,
        total = total, background = background)
    plot_dat2 <- prepare_ArraySpeciesBySize_plot_data(
        y, species = species, all.sizes = all.sizes, wlim = wlim,
        total = total, background = background)

    plotRelativeDataFrame(plot_dat1, plot_dat2, params,
                          xlab = plot_size_xlab(size_axis),
                          xtrans = if (log_x) "log10" else "identity",
                          xlim = plot_size_xlim(wlim, size_axis, llim),
                          ylim = ylim,
                          legend_var = "Legend", size_axis = size_axis)
}

check_plot2_compatible <- function(x, y, class) {
    if (!inherits(y, class)) {
        stop("Both objects must be of class `", class, "`.")
    }
}

compare_array_metadata <- function(x, y) {
    value_name1 <- attr(x, "value_name")
    value_name2 <- attr(y, "value_name")
    if (!is.null(value_name1) && !is.null(value_name2) &&
            !identical(value_name1, value_name2)) {
        warning("The first array has value name `", value_name1,
                "`, but the second array has value name `", value_name2, "`.")
    }
    units1 <- attr(x, "units")
    units2 <- attr(y, "units")
    if (!is.null(units1) && !is.null(units2) &&
            nzchar(units1) && nzchar(units2) &&
            !identical(units1, units2)) {
        warning("The first array has y units `", units1,
                "`, but the second array has y units `", units2, "`.")
    }
}

array_y_label <- function(x, default = "Value") {
    value_name <- attr(x, "value_name") %||% default
    units_str <- attr(x, "units")
    if (!is.null(units_str) && nzchar(units_str)) {
        value_name <- paste0(value_name, " [", units_str, "]")
    }
    value_name
}

#' Add lines to an existing plot
#'
#' `r lifecycle::badge("experimental")`
#' `addPlot()` adds another set of values to an existing ggplot. The first
#' method supports adding an `ArraySpeciesBySize` object to a compatible plot,
#' for example to compare the same rate before and after a model change.
#' The method checks whether the existing plot uses a compatible x variable,
#' and warns if the y variable or y-axis units appear to differ.
#'
#' @param plot A ggplot2 object to which the new values should be added.
#' @param x An object containing the values to add.
#' @param species Character vector of species to include. `NULL` (default) means
#'   all species.
#' @param total A boolean value that determines whether the total over all
#'   selected species is plotted as well. Default is `FALSE`.
#' @param background A boolean value that determines whether background species
#'   are included. Ignored if the model does not contain background species.
#'   Default is `TRUE`.
#' @param colour Optional fixed colour for the added lines. If `NULL`, the
#'   species colours from the existing plot are used.
#' @param linetype Optional fixed line type for the added lines. If `NULL`, the
#'   species line types from the existing plot are used.
#' @param linewidth Width of the added lines.
#' @param alpha Transparency of the added lines.
#' @param ... Further arguments used by only some of the methods:
#'
#'   **For `ArraySpeciesBySize` methods:**
#'   \describe{
#'     \item{`all.sizes`}{If `FALSE` (default), values outside a species' size
#'       range (`w_min` to `w_max`) are removed.}
#'     \item{`wlim`}{A numeric vector of length two providing lower and upper
#'       limits for the weight (x) axis. Use `NA` to refer to the existing
#'       minimum or maximum.}
#'     \item{`llim`}{A numeric vector of length two providing lower and upper
#'       limits for the length (x) axis when `size_axis = "l"`. Use `NA` to
#'       refer to the existing minimum or maximum.}
#'     \item{`size_axis`}{Whether to plot size as weight (`"w"`, default) or
#'       length (`"l"`), using the allometric weight-length relationship.}
#'   }
#'
#'   **For `ArrayTimeBySpecies` methods:**
#'   \describe{
#'     \item{`tlim`}{A numeric vector of length two providing lower and upper
#'       limits for the time axis, e.g. `c(1980, 2000)`. Use `NA` to apply no
#'       limit at that end. Default is `c(NA, NA)`.}
#'     \item{`ylim`}{A numeric vector of length two providing lower and upper
#'       limits for the value (y) axis.}
#'   }
#'
#' @return A ggplot2 object.
#' @export
#' @family plotting functions
#'
#' @examples
#' \donttest{
#' p <- plot(getEncounter(NS_params), species = "Cod")
#' addPlot(p, getEncounter(NS_params), species = "Cod")
#' }
addPlot <- function(plot, x, species = NULL, total = FALSE,
                    background = TRUE, colour = NULL, linetype = "dashed",
                    linewidth = 0.8, alpha = 1, ...) {
    UseMethod("addPlot", x)
}

#' @rdname addPlot
#' @usage NULL
#' @export
addPlot.ArraySpeciesBySize <- function(plot, x, species = NULL,
                                       total = FALSE,
                                       background = TRUE,
                                       colour = NULL,
                                       linetype = "dashed",
                                       linewidth = 0.8,
                                       alpha = 1,
                                       all.sizes = FALSE,
                                       wlim = c(NA, NA),
                                       llim = c(NA, NA),
                                       size_axis = c("w", "l"),
                                       ...) {
    if (!inherits(plot, "ggplot")) {
        stop("The `plot` argument must be a ggplot object.")
    }
    assert_that(is.number(linewidth),
                is.number(alpha),
                alpha >= 0,
                alpha <= 1)
    size_axis <- plot_size_axis(size_axis)
    assert_that(length(wlim) == 2,
                length(llim) == 2)

    plot <- deep_copy(plot)
    plot_dat <- prepare_ArraySpeciesBySize_plot_data(
        x, species = species, all.sizes = all.sizes, wlim = wlim,
        total = total, background = background)
    params <- attr(x, "params")
    plot_dat <- convert_plot_size_axis(plot_dat, params, size_axis)
    if (identical(size_axis, "l")) {
        plot_dat <- filter_plot_length_limits(plot_dat, llim)
    }
    x_var <- plot_size_x_var(size_axis)
    y_var <- names(plot_dat)[2]
    check_addPlot_compatible(plot, x_var = x_var, y_var = y_var,
                             units = attr(x, "units"))

    mapping <- aes(x = .data[[x_var]], y = .data[[y_var]],
                   group = .data[["Species"]])
    if (is.null(colour)) {
        mapping$colour <- rlang::quo(.data[["Legend"]])
    }
    if (is.null(linetype)) {
        mapping$linetype <- rlang::quo(.data[["Legend"]])
    }

    layer_args <- list(
        data = plot_dat,
        mapping = mapping,
        linewidth = linewidth,
        alpha = alpha,
        inherit.aes = FALSE
    )
    if (!is.null(colour)) {
        layer_args$colour <- colour
    }
    if (!is.null(linetype)) {
        layer_args$linetype <- linetype
    }

    plot + do.call(geom_line, layer_args)
}

deep_copy <- function(x) {
    unserialize(serialize(x, NULL))
}

check_addPlot_compatible <- function(plot, x_var, y_var, units = NULL) {
    mapping <- plot$mapping
    if (length(plot$layers) > 0) {
        layer_mapping <- plot$layers[[1]]$mapping
        if (!is.null(layer_mapping$x)) {
            mapping$x <- layer_mapping$x
        }
        if (!is.null(layer_mapping$y)) {
            mapping$y <- layer_mapping$y
        }
    }

    plot_x_var <- plot_mapping_var(mapping$x)
    if (!is.null(plot_x_var) && !identical(plot_x_var, x_var)) {
        stop("The data can only be added to a plot with x variable `", x_var,
             "`. The existing plot uses x variable `", plot_x_var, "`.")
    }

    plot_y_var <- plot_mapping_var(mapping$y)
    if (!is.null(plot_y_var) && !identical(plot_y_var, y_var)) {
        warning("The existing plot appears to use y variable `", plot_y_var,
                "`, but the added data uses `", y_var, "`.")
    }

    plot_units <- plot_y_units(plot)
    if (!is.null(plot_units) && !is.null(units) &&
            nzchar(plot_units) && nzchar(units) &&
            !identical(plot_units, units)) {
        warning("The existing plot appears to use y units `", plot_units,
                "`, but the added data uses `", units, "`.")
    }
}

plot_mapping_var <- function(mapping) {
    if (is.null(mapping)) {
        return(NULL)
    }

    label <- rlang::as_label(mapping)
    match <- regmatches(label, regexec("\\.data\\[\\[\"([^\"]+)\"\\]\\]", label))[[1]]
    if (length(match) == 2) {
        return(match[[2]])
    }
    match <- regmatches(label, regexec("\\.data\\$([^ ]+)$", label))[[1]]
    if (length(match) == 2) {
        return(match[[2]])
    }
    if (grepl("^[[:alnum:]_.]+$", label)) {
        return(label)
    }

    NULL
}

plot_y_units <- function(plot) {
    scales <- plot$scales$scales
    for (scale in scales) {
        if ("y" %in% scale$aesthetics && is.character(scale$name)) {
            return(label_units(scale$name))
        }
    }
    label_units(plot$labels$y)
}

label_units <- function(label) {
    if (is.null(label) || !is.character(label) || !nzchar(label)) {
        return(NULL)
    }
    match <- regmatches(label, regexec("\\[([^]]+)\\]\\s*$", label))[[1]]
    if (length(match) == 2) {
        return(match[[2]])
    }
    NULL
}

apply_wlim <- function(data, wlim) {
    if (!is.na(wlim[1])) data <- data[data$w >= wlim[1], ]
    if (!is.na(wlim[2])) data <- data[data$w <= wlim[2], ]
    data
}

prepare_ArraySpeciesBySize_plot_data <- function(x, species = NULL,
                                                 all.sizes = FALSE,
                                                 wlim = c(NA, NA),
                                                 total = FALSE,
                                                 background = TRUE) {
    params <- attr(x, "params")
    w <- get_ArraySpeciesBySize_w(x)

    all_species <- rownames(x)
    if (is.null(species)) {
        species <- all_species
    } else {
        species <- intersect(species, all_species)
        if (length(species) == 0) {
            stop("None of the selected species are in the rate array.")
        }
    }

    value_name <- attr(x, "value_name") %||% "value"
    sel <- all_species %in% species
    mat <- unclass(x)[sel, , drop = FALSE]

    # Compute total across all selected species before size-range trimming
    if (total) {
        total_row <- colSums(mat, na.rm = TRUE)
    }

    plot_dat <- data.frame(
        w = rep(w, each = sum(sel)),
        value = c(mat),
        Species = rownames(mat)
    )

    if (!all.sizes && !is.null(params)) {
        sp_params <- params@species_params
        for (sp in species) {
            if (sp %in% sp_params$species) {
                sp_row <- sp_params[sp_params$species == sp, ]
                plot_dat$value[plot_dat$Species == sp &
                                   (plot_dat$w < sp_row$w_min[1] |
                                        plot_dat$w > sp_row$w_max[1])] <- NA
            }
        }
        plot_dat <- plot_dat[complete.cases(plot_dat), ]
    }

    plot_dat <- apply_wlim(plot_dat, wlim)

    plot_dat$Legend <- plot_dat$Species

    # Handle background species
    if (!is.null(params) && isTRUE(any(params@species_params$is_background))) {
        bkgrd_sp <- params@species_params$species[params@species_params$is_background]
        if (background) {
            plot_dat$Legend[plot_dat$Species %in% bkgrd_sp] <- "Background"
        } else {
            plot_dat <- plot_dat[!plot_dat$Species %in% bkgrd_sp, ]
        }
    }

    # Add total line
    if (total) {
        total_dat <- data.frame(
            w = w,
            value = total_row,
            Species = "Total",
            Legend = "Total"
        )
        total_dat <- apply_wlim(total_dat, wlim)
        plot_dat <- rbind(plot_dat, total_dat)
    }

    names(plot_dat)[2] <- value_name

    plot_dat
}

#' @rdname plotHover
#' @usage NULL
#' @examples
#' \donttest{
#' plotHover(getEncounter(NS_params))
#' }
#' @export
plotHover.ArraySpeciesBySize <- function(x, ...) {
    plotHover(plot(x, ...), ...)
}

#' @export
as.data.frame.ArraySpeciesBySize <- function(x, row.names = NULL,
                                     optional = FALSE, ...) {
    w <- get_ArraySpeciesBySize_w(x)
    sp_names <- rownames(x)
    mat <- unclass(x)
    data.frame(
        w = rep(w, each = nrow(mat)),
        value = c(mat),
        Species = sp_names,
        row.names = row.names,
        check.names = !optional,
        stringsAsFactors = FALSE
    )
}

#' Get the size grid for an ArraySpeciesBySize object
#'
#' Internal helper that returns the consumer size grid `params@w` or the full
#' prey/resource size grid `params@w_full`, depending on the number of columns
#' in the array.
#'
#' @param x An `ArraySpeciesBySize` object.
#'
#' @return A numeric vector giving the size represented by each column.
#' @keywords internal
get_ArraySpeciesBySize_w <- function(x) {
    params <- attr(x, "params")
    if (is.null(params)) {
        w <- as.numeric(colnames(x))
        if (any(is.na(w))) {
            w <- seq_len(ncol(x))
        }
        return(w)
    }
    if (ncol(x) == length(params@w)) {
        return(params@w)
    }
    if (ncol(x) == length(params@w_full)) {
        return(params@w_full)
    }
    stop("Can not determine the size grid for this ArraySpeciesBySize object. ",
         "The number of columns is ", ncol(x), ", but the params object has ",
         length(params@w), " consumer sizes and ", length(params@w_full),
         " full-spectrum sizes.")
}

#' @export
`[.ArraySpeciesBySize` <- function(x, i, j, ..., drop = TRUE) {
    result <- NextMethod()
    # Preserve class only if result is still a 2D matrix
    if (is.matrix(result) && length(dim(result)) == 2) {
        attr(result, "value_name") <- attr(x, "value_name")
        attr(result, "units") <- attr(x, "units")
        attr(result, "params") <- attr(x, "params")
        class(result) <- c("ArraySpeciesBySize", "matrix", "array")
    }
    result
}

#' @export
Ops.ArraySpeciesBySize <- function(e1, e2) {
    # Strip ArraySpeciesBySize class so that arithmetic returns a plain matrix.
    # We unclass both operands and call the generic directly.
    if (is.ArraySpeciesBySize(e1)) e1 <- unclass_rate(e1)
    if (!missing(e2) && is.ArraySpeciesBySize(e2)) e2 <- unclass_rate(e2)
    op <- match.fun(.Generic)
    if (missing(e2)) op(e1) else op(e1, e2)
}

# Helper to strip all ArraySpeciesBySize attributes
unclass_rate <- function(x) {
    x <- unclass(x)
    attr(x, "value_name") <- NULL
    attr(x, "units") <- NULL
    attr(x, "params") <- NULL
    x
}

#' @export
str.ArraySpeciesBySize <- function(object, ...) {
    params <- attr(object, "params")
    attr(object, "params") <- NULL
    class(object) <- c("matrix", "array")
    out <- capture.output(utils::str(object, ...))
    out[1] <- paste0(" 'ArraySpeciesBySize' ", sub("^ ", "", out[1]))
    cat(paste0(out, collapse = "\n"), "\n", sep = "")
    if (!is.null(params)) {
        cat(" - attr(*, \"params\")=Formal class 'MizerParams' [package \"mizer\"] with ",
            length(slotNames(params)), " slots\n", sep = "")
    }
    invisible(NULL)
}
