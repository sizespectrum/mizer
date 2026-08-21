# MizerScan S3 class for the result of a parameter scan
#
# Copyright 2026 Gustav Delius.
# Distributed under the GPL 3 or later.

#' S3 class for the result of a parameter scan
#'
#' `r lifecycle::badge("experimental")`
#' [scanModel()] varies one aspect of a model over a range of values and
#' measures a quantity on the attractor the model settles on at each of them.
#' It returns its result as a `MizerScan` object, which is a data frame
#' carrying, in addition, everything that [plot()] needs to draw it.
#'
#' A `MizerScan` object behaves like an ordinary data frame with one row for
#' each combination of scanned value and series. Its columns are, in order:
#' \describe{
#'   \item{1}{The scanned value. The column is named after the `scan_name`, so
#'     for example a scan over fishing effort has a column called
#'     `"Fishing effort"`. Use `names(scan)[[1]]` rather than hard-coding it.}
#'   \item{2}{The measured quantity, averaged over the attractor. Named after
#'     the `value_name`.}
#'   \item{`Species`}{The series the row refers to. Named `Species` whatever the
#'     series are, because that is the column that mizer's colour and line-type
#'     machinery reads.}
#'   \item{`ymin`, `ymax`}{The smallest and largest value over the sampling
#'     window. On a fixed point these both equal the value; on a limit cycle
#'     they give the range of the oscillation.}
#'   \item{`type`, `settled`}{What kind of attractor was reached, from the
#'     `"convergence"` attribute that [projectToSteady()] attaches to its
#'     result.}
#'   \item{`period`}{The period of the limit cycle in years, or `NA`.}
#'   \item{`residual`}{How far the state still is from a fixed point, as a
#'     per-capita rate in 1/year, see [getSteadyResidual()].}
#' }
#'
#' The first three columns are the x, y and grouping variable in that order,
#' which is the layout [plotDataFrame()] expects.
#'
#' It also carries these attributes:
#' \itemize{
#'   \item `scan_name`, `scan_units` – name and units of the quantity that was
#'     varied, used for the x-axis label.
#'   \item `value_name`, `value_units` – name and units of the quantity that was
#'     measured, used for the y-axis label.
#'   \item `type` – the kind of quantity the values are, see [array_types].
#'   \item `params` – the `MizerParams` object the scan started from, used for
#'     species colours and line types.
#'   \item `reference_lines` – an optional named numeric vector of positions on
#'     the x axis to mark with vertical lines, for example
#'     `c(F_MSY = 0.32)`.
#'   \item `at_max`, `max_value` – for each series, the scanned value at which
#'     the measured quantity is largest, and the value it takes there. See the
#'     section below.
#'   \item `settings` – a list recording the settings the scan was run with.
#' }
#'
#' @section Where the maximum is:
#'
#' The `at_max` attribute holds, for each series, the scanned value at which
#' that series' measured quantity is largest. On a yield-versus-fishing-mortality
#' scan that is \eqn{F_{MSY}}; on a scan over fishing effort it is the effort
#' that maximises the quantity being plotted. `max_value` holds the value
#' attained there.
#'
#' This is the largest value **among those that were scanned**, not the maximum
#' of the underlying curve. It is therefore only as good as the grid you gave in
#' `scan_values`, and the way to sharpen it is to scan a finer grid near the
#' maximum, not to interpolate a coarse one. Subsetting a `MizerScan` with `[`
#' recomputes both attributes from the rows that remain, so they never go stale.
#'
#' @section Limitations:
#'
#' Because the object is a data frame subclass, its attributes survive base R
#' subsetting with `[` but are dropped by functions that rebuild the data frame,
#' including `dplyr::filter()`, `dplyr::mutate()`, `subset()` and `transform()`.
#' A scan that has lost its attributes can no longer be plotted. Subset with `[`,
#' or rebuild the object with `MizerScan()`.
#'
#' @param x A data frame with at least three columns, laid out as described
#'   above. For `is.MizerScan()`, any object to test.
#' @param scan_name A string naming the quantity that was varied.
#' @param scan_units A string giving its units, for example `"1/year"`.
#' @param value_name A string naming the quantity that was measured.
#' @param value_units A string giving its units, for example `"g/year"`.
#' @param type The kind of quantity the measured values are, see [array_types].
#' @param params The `MizerParams` object the scan started from.
#' @param reference_lines An optional named numeric vector of x positions to
#'   mark with vertical lines.
#' @param settings An optional list recording the settings used.
#'
#' @return A `MizerScan` object, which inherits from `data.frame`.
#'
#' @seealso [scanModel()], [plot.MizerScan()]
#' @family scan functions
#' @concept scan
#' @export
#' @examples
#' \donttest{
#' scan <- scanModel(NS_params, scan_values = c(0, 0.5, 1),
#'                   set_func = scanEffort(), species = "Cod")
#' scan
#' summary(scan)
#' attr(scan, "at_max")
#' }
MizerScan <- function(x, scan_name = NULL, scan_units = NULL,
                      value_name = NULL, value_units = NULL,
                      type = NULL, params = NULL,
                      reference_lines = NULL, settings = NULL) {
    if (!is.data.frame(x)) {
        stop("`x` must be a data frame.")
    }
    if (ncol(x) < 3) {
        stop("`x` must have at least three columns: the scanned value, the ",
             "measured value and the series.")
    }
    if (!("Species" %in% names(x))) {
        stop("`x` must have a `Species` column naming the series.")
    }
    # The type is resolved once and stored, so that array_type() and the
    # helpers built on it never need to consult the units attribute, which
    # this class calls `value_units` rather than `units`.
    type <- resolve_array_type(type, value_name, value_units)
    x <- structure(x,
        class = c("MizerScan", "data.frame"),
        scan_name = scan_name,
        scan_units = scan_units,
        value_name = value_name,
        value_units = value_units,
        type = type,
        params = params,
        reference_lines = reference_lines,
        settings = settings
    )
    set_scan_maximum(x)
}

#' @rdname MizerScan
#' @return `is.MizerScan()` returns `TRUE` if `x` is a `MizerScan` object,
#'   `FALSE` otherwise.
#' @export
is.MizerScan <- function(x) {
    inherits(x, "MizerScan")
}

# ---------------------------------------------------------------------------
# Internal accessors
#
# The first two column names are set from `scan_name` and `value_name` and are
# therefore user data. Nothing outside these two functions may hard-code them.
# ---------------------------------------------------------------------------

#' The x and y variables of a MizerScan
#'
#' @param x A MizerScan object.
#' @return The name of the column holding the scanned value / the measured
#'   value.
#' @keywords internal
scan_x_var <- function(x) names(x)[[1]]

#' @rdname scan_x_var
#' @keywords internal
scan_y_var <- function(x) names(x)[[2]]

#' Axis labels for a MizerScan
#'
#' Assembles `"<name> [<units>]"`, the same way `array_y_label()` does for the
#' array classes.
#'
#' @param x A MizerScan object.
#' @param default The label to use when the name is missing.
#' @return A string.
#' @keywords internal
scan_y_label <- function(x, default = "Value") {
    label_with_units(attr(x, "value_name") %||% default,
                     attr(x, "value_units"))
}

#' @rdname scan_y_label
#' @keywords internal
scan_x_label <- function(x, default = "Scan value") {
    label_with_units(attr(x, "scan_name") %||% default,
                     attr(x, "scan_units"))
}

#' @rdname scan_y_label
#' @param name The name of the quantity.
#' @param units The units, possibly NULL.
#' @keywords internal
label_with_units <- function(name, units) {
    if (!is.null(units) && nzchar(units)) {
        paste0(name, " [", units, "]")
    } else {
        name
    }
}

#' Record where each series attains its maximum
#'
#' Sets the `at_max` and `max_value` attributes from the rows currently in the
#' object, so that subsetting cannot leave a stale maximum behind.
#'
#' @param x A MizerScan object.
#' @return `x` with the two attributes set.
#' @keywords internal
set_scan_maximum <- function(x) {
    x_var <- scan_x_var(x)
    y_var <- scan_y_var(x)
    series <- unique(as.character(x[["Species"]]))
    at_max <- rep(NA_real_, length(series))
    max_value <- rep(NA_real_, length(series))
    names(at_max) <- names(max_value) <- series
    for (s in series) {
        rows <- which(as.character(x[["Species"]]) == s)
        y <- x[[y_var]][rows]
        if (all(is.na(y))) next
        best <- rows[[which.max(y)]]
        at_max[[s]] <- x[[x_var]][[best]]
        max_value[[s]] <- x[[y_var]][[best]]
    }
    attr(x, "at_max") <- at_max
    attr(x, "max_value") <- max_value
    x
}

#' Copy the metadata of a MizerScan onto another object
#'
#' @param to The object to receive the attributes.
#' @param from The MizerScan to take them from.
#' @return `to` with the attributes and class of a MizerScan.
#' @keywords internal
copy_scan_attributes <- function(to, from) {
    for (a in c("scan_name", "scan_units", "value_name", "value_units",
                "type", "params", "reference_lines", "settings")) {
        attr(to, a) <- attr(from, a)
    }
    class(to) <- c("MizerScan", "data.frame")
    to
}

# ---------------------------------------------------------------------------
# Methods
# ---------------------------------------------------------------------------

#' @export
print.MizerScan <- function(x, ...) {
    header <- paste0(scan_y_label(x, default = "Value"), " vs ",
                     scan_x_label(x))
    cat(header, "\n")
    n_values <- length(unique(x[[scan_x_var(x)]]))
    n_series <- length(unique(as.character(x[["Species"]])))
    cat(n_values, "scan values x", n_series, "series\n")

    frame <- as.data.frame(x)
    n_show <- min(nrow(frame), mizer_print_defaults$time_max %||% 10)
    print(frame[seq_len(n_show), , drop = FALSE], row.names = FALSE)
    if (n_show < nrow(frame)) {
        cat("# ... ", nrow(frame) - n_show, " more rows\n", sep = "")
    }
    unsettled <- unique(x[[scan_x_var(x)]][!x$settled])
    if (length(unsettled) > 0) {
        cat("Did not settle at ", scan_x_var(x), " = ",
            paste(signif(unsettled, 3), collapse = ", "), "\n", sep = "")
    }
    invisible(x)
}

#' @export
summary.MizerScan <- function(object, ...) {
    x_var <- scan_x_var(object)
    y_var <- scan_y_var(object)
    at_max <- attr(object, "at_max")
    max_value <- attr(object, "max_value")
    series <- names(at_max)

    per_series <- data.frame(
        Species = series,
        Min = vapply(series, function(s) {
            min(object[[y_var]][as.character(object$Species) == s], na.rm = TRUE)
        }, numeric(1)),
        Max = unname(max_value[series]),
        at_max = unname(at_max[series]),
        row.names = NULL,
        stringsAsFactors = FALSE
    )

    result <- list(
        value_label = scan_y_label(object),
        scan_label = scan_x_label(object),
        x_var = x_var,
        n_values = length(unique(object[[x_var]])),
        range = range(object[[x_var]], na.rm = TRUE),
        per_series = per_series,
        attractor = table(object$type[!duplicated(object[[x_var]])]),
        settings = attr(object, "settings")
    )
    class(result) <- "summary.MizerScan"
    result
}

#' @export
print.summary.MizerScan <- function(x, ...) {
    cat(x$value_label, "vs", x$scan_label, "\n")
    cat(x$n_values, "scan values from", signif(x$range[[1]], 3), "to",
        signif(x$range[[2]], 3), "\n\n")
    print(x$per_series, row.names = FALSE)
    cat("\n`at_max` is the scanned value with the largest value, over the",
        "\nvalues that were scanned. Scan a finer grid to sharpen it.\n")
    if (length(x$attractor) > 0) {
        cat("\nAttractors reached:\n")
        print(x$attractor)
    }
    invisible(x)
}

#' @export
as.data.frame.MizerScan <- function(x, ...) {
    for (a in c("scan_name", "scan_units", "value_name", "value_units",
                "type", "params", "reference_lines", "settings",
                "at_max", "max_value")) {
        attr(x, a) <- NULL
    }
    class(x) <- "data.frame"
    x
}

#' @export
str.MizerScan <- function(object, ...) {
    params <- attr(object, "params")
    attr(object, "params") <- NULL
    attr(object, "settings") <- NULL
    class(object) <- "data.frame"
    out <- utils::capture.output(utils::str(object, ...))
    out[1] <- paste0(" 'MizerScan' ", sub("^ ", "", out[1]))
    cat(paste0(out, collapse = "\n"), "\n", sep = "")
    if (!is.null(params)) {
        cat(" - attr(*, \"params\")=Formal class 'MizerParams' [package \"mizer\"] with ",
            length(slotNames(params)), " slots\n", sep = "")
    }
    invisible(NULL)
}

#' @export
`[.MizerScan` <- function(x, i, j, ..., drop = if (missing(i)) TRUE else FALSE) {
    old_names <- names(x)[1:3]
    result <- NextMethod()
    # `[.data.frame` keeps the class of a subclass whatever it does to the
    # contents, so the class has to be taken away again explicitly once the
    # result has stopped being a scan.
    still_a_scan <- is.data.frame(result) && ncol(result) >= 3 &&
        identical(names(result)[1:3], old_names)
    if (still_a_scan) {
        result <- copy_scan_attributes(result, x)
        # Recompute rather than copy, so a subset never carries a maximum
        # taken over rows it no longer contains.
        result <- set_scan_maximum(result)
    } else if (is.data.frame(result)) {
        result <- as.data.frame(result)
    }
    result
}

# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

#' Plot method for `MizerScan` objects
#'
#' `r lifecycle::badge("experimental")`
#' Draws the result of a [scanModel()] run: the measured quantity against the
#' quantity that was scanned, with a band showing the range the quantity takes
#' over the attractor wherever that attractor is not a fixed point.
#'
#' A model that settles on a fixed point contributes a single value, so the band
#' has zero width there. A model that settles on a limit cycle contributes the
#' average over one period as the line and the range over that period as the
#' band, so a Hopf bifurcation shows up as the scan value at which the band
#' opens up.
#'
#' Scan values where the model reached neither a fixed point nor a limit cycle
#' within the time allowed are marked with a cross, because the value plotted
#' there is only an average over the last few years of a run that was still
#' changing.
#'
#' @param x A `MizerScan` object.
#' @param species The species to show. By default all series in the scan.
#' @param style One of `"ribbon"` (default: the average as a line inside the
#'   band), `"envelope"` (lines along the edges of the band, no average) or
#'   `"line"` (no band).
#' @param highlight Name or vector of names of the species to be highlighted
#'   with a thicker line.
#' @param log_x,log_y,log Whether to use logarithmic axes, see [parsePlotLog()].
#' @param xlim,ylim Numeric vectors of length two giving the axis limits. Use
#'   `NA` to refer to the existing minimum or maximum.
#' @param y_ticks The approximate number of ticks desired on the y axis.
#' @param reference_lines Whether to draw the reference lines stored in the
#'   scan, or a named numeric vector of x positions to draw instead.
#' @param mark_max Whether to mark, for each series, the scanned value at which
#'   the measured quantity is largest. See [MizerScan()].
#' @param show_unsettled Whether to mark the scan values where the model did not
#'   settle onto an attractor.
#' @param return_data Whether to return the data frame used for the plot instead
#'   of the plot itself.
#' @param ... Unused.
#'
#' @return A ggplot2 object, unless `return_data = TRUE`, in which case the data
#'   frame used for the plot is returned.
#' @family scan functions
#' @concept scan
#' @seealso [scanModel()], [MizerScan()], [plotting_functions]
#' @export
#' @examples
#' \donttest{
#' scan <- scanModel(NS_params, scan_values = seq(0, 1, 0.25),
#'                   set_func = scanEffort(), species = c("Cod", "Herring"))
#' plot(scan)
#' plot(scan, style = "envelope", mark_max = TRUE)
#' }
plot.MizerScan <- function(x, species = NULL,
                           style = c("ribbon", "envelope", "line"),
                           highlight = NULL,
                           log_x = FALSE, log_y = TRUE, log = NULL,
                           xlim = c(NA, NA), ylim = c(NA, NA), y_ticks = 6,
                           reference_lines = TRUE, mark_max = FALSE,
                           show_unsettled = TRUE, return_data = FALSE, ...) {
    style <- match.arg(style)
    log_y <- array_log_y(x, log_y, log, !missing(log_y))
    log_axes <- parsePlotLog(log, log_x = log_x, log_y = log_y)

    plot_dat <- prepare_MizerScan_plot_data(x, species = species)
    if (return_data) return(plot_dat)

    params <- scan_plot_params(x, plot_dat)
    p <- plotDataFrame(plot_dat, params,
                       style = if (style == "line") "line" else style,
                       xlab = scan_x_label(x), ylab = scan_y_label(x),
                       xtrans = if (log_axes$log_x) "log10" else "identity",
                       ytrans = if (log_axes$log_y) "log10" else "identity",
                       xlim = xlim,
                       ylim = array_ylim(x, ylim, log_axes$log_y, plot_dat[[2]]),
                       y_ticks = y_ticks, highlight = highlight)
    add_scan_annotations(p, x, plot_dat, reference_lines = reference_lines,
                         mark_max = mark_max,
                         show_unsettled = show_unsettled)
}

#' Prepare the data frame for plotting a MizerScan
#'
#' @param x A MizerScan object.
#' @param species The series to keep, or NULL for all of them.
#' @return A data frame with the x, y and grouping variable in the first three
#'   columns, as [plotDataFrame()] requires.
#' @keywords internal
prepare_MizerScan_plot_data <- function(x, species = NULL) {
    frame <- as.data.frame(x)
    if (!is.null(species)) {
        params <- attr(x, "params")
        if (!is.null(params)) {
            species <- valid_species_arg(params, species,
                                         error_on_empty = TRUE)
        }
        keep <- as.character(frame$Species) %in% as.character(species)
        if (!any(keep)) {
            stop("None of the requested species are in this scan.")
        }
        frame <- frame[keep, , drop = FALSE]
    }
    rownames(frame) <- NULL
    frame
}

#' The params object to use when plotting a MizerScan
#'
#' Series that are not species have no colour in the model, and
#' [plotDataFrame()] silently drops any legend level it cannot find a colour
#' for. So any such series is given a colour here, using the ordinary
#' [setColours()] interface, which also leaves the user free to choose a
#' different one.
#'
#' @param x A MizerScan object.
#' @param plot_dat The data frame that will be plotted.
#' @return A MizerParams object with a colour for every series in `plot_dat`.
#' @keywords internal
scan_plot_params <- function(x, plot_dat) {
    params <- attr(x, "params")
    if (is.null(params)) {
        stop("This scan has lost its `params` attribute, so it cannot be ",
             "plotted. Functions that rebuild a data frame, such as ",
             "dplyr::filter() or subset(), drop it. Subset with `[` instead, ",
             "or rebuild the object with MizerScan().")
    }
    series <- unique(as.character(plot_dat$Species))
    missing <- setdiff(series, names(getColours(params)))
    if (length(missing) > 0) {
        colours <- grDevices::hcl.colors(max(3, length(missing)), "Dark 3")
        params <- setColours(params,
                             stats::setNames(as.list(colours[seq_along(missing)]),
                                             missing))
    }
    # A line style is needed too, or the manual linetype scale has no level to
    # match and ggplot2 warns.
    missing <- setdiff(series, names(getLinetypes(params)))
    if (length(missing) > 0) {
        params <- setLinetypes(params,
                               stats::setNames(as.list(rep("solid",
                                                           length(missing))),
                                               missing))
    }
    params
}

#' Add the annotation layers to a MizerScan plot
#'
#' @param p The ggplot object so far.
#' @param x The MizerScan object.
#' @param plot_dat The data frame being plotted.
#' @param reference_lines TRUE to use the stored reference lines, FALSE for
#'   none, or a named numeric vector to use instead.
#' @param mark_max Whether to mark where each series attains its maximum.
#' @param show_unsettled Whether to mark the scan values that did not settle.
#' @return The ggplot object with the extra layers, still a `mizer_plot`.
#' @keywords internal
add_scan_annotations <- function(p, x, plot_dat, reference_lines = TRUE,
                                 mark_max = FALSE, show_unsettled = TRUE) {
    tooltip <- attr(p, "mizer_tooltip")
    x_var <- scan_x_var(x)
    y_var <- scan_y_var(x)

    if (isTRUE(show_unsettled)) {
        # Base subsetting rather than subset() or dplyr::filter(), so that no
        # variable name has to be registered with globalVariables().
        bad <- plot_dat[!plot_dat$settled, , drop = FALSE]
        if (nrow(bad) > 0) {
            p <- p +
                geom_point(data = bad,
                           aes(x = .data[[x_var]], y = .data[[y_var]]),
                           shape = 4, size = 2, colour = "black",
                           inherit.aes = FALSE, show.legend = FALSE) +
                labs(caption = paste("x marks values where the model did not",
                                     "settle onto an attractor"))
        }
    }

    refs <- if (isTRUE(reference_lines)) {
        attr(x, "reference_lines")
    } else if (isFALSE(reference_lines)) {
        NULL
    } else {
        reference_lines
    }
    if (length(refs) > 0) {
        p <- p + geom_vline(xintercept = unname(refs), linetype = "dashed",
                            colour = "grey50")
        if (!is.null(names(refs)) && all(nzchar(names(refs)))) {
            p <- p + annotate("text", x = unname(refs),
                              y = Inf, label = names(refs),
                              hjust = -0.1, vjust = 1.5, size = 3,
                              colour = "grey40")
        }
    }

    if (isTRUE(mark_max)) {
        at_max <- attr(x, "at_max")
        max_value <- attr(x, "max_value")
        keep <- names(at_max) %in% unique(as.character(plot_dat$Species))
        at_max <- at_max[keep]
        max_value <- max_value[keep]
        if (length(at_max) > 0 && any(!is.na(at_max))) {
            marks <- data.frame(xx = unname(at_max), yy = unname(max_value))
            marks <- marks[!is.na(marks$xx), , drop = FALSE]
            p <- p + geom_point(data = marks, aes(x = .data[["xx"]],
                                                  y = .data[["yy"]]),
                                shape = 1, size = 3, colour = "black",
                                inherit.aes = FALSE, show.legend = FALSE)
        }
    }

    make_mizer_plot(p, tooltip)
}

#' @rdname plotHover
#' @export
plotHover.MizerScan <- function(x, ...) {
    plotHover(plot(x, ...))
}
