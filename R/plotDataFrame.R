# Compatibility adapters, shared ggplot metadata and hover conversion.

# Hackiness to get past the 'no visible binding ... ' warning when running check
utils::globalVariables(c("time", "value", "Species", "w", "l", "gear", "Age",
                         "x", "y", "Year", "Yield", "Biomass", "Size",
                         "Proportion", "Prey", "Legend", "Type", "Gear",
                         "Predator", "weight", "a", "b", "age", "w_max",
                         "Model", "rel_diff"))

# A limit of zero or below cannot be shown on a logarithmic axis: it maps to
# -Inf, which ggplot2 drops with a warning about infinite values. Such a limit
# is only ever a request to start the axis at the bottom of the data, so on a
# log10 axis we auto-scale instead of passing the impossible limit on.
drop_nonpositive_limits <- function(lim, trans) {
    if (identical(trans, "log10")) {
        lim[!is.na(lim) & lim <= 0] <- NA
    }
    lim
}

# Values of zero or below have no place on a logarithmic axis either. ggplot2
# discards them anyway, but only after warning about infinite values, which is
# noise for the many mizer rates that are legitimately zero somewhere (fishing
# mortality outside the selected sizes, predation mortality on the largest
# individuals). We discard them ourselves instead.
drop_nonpositive_values <- function(frame, var, trans) {
    if (identical(trans, "log10")) {
        keep <- is.na(frame[[var]]) | frame[[var]] > 0
        # If nothing at all can be shown there is nothing to be gained by
        # emptying the data frame, which would only trade the warning for the
        # complaints that an empty plot draws.
        if (any(keep)) frame <- frame[keep, , drop = FALSE]
    }
    frame
}

#' Make a plot from a data frame
#'
#' This is used internally by most plotting functions.
#'
#' @param frame A data frame with at least three variables.
#'   The first three variables are used, in that order, as:
#'   1. Variable to be plotted on x-axis
#'   2. Variable to be plotted on y-axis
#'   3. Grouping variable
#' @param params A MizerParams object, which is used for the line colours and
#'   line types.
#' @param style The style of the plot. Available options are `"line"` for
#'   `geom_line()`, `"area"` for `geom_area()`, and `"ribbon"` and `"envelope"`
#'   for a band drawn between the `ymin` and `ymax` variables, which must then
#'   be present in `frame`. `"ribbon"` draws the y variable as a line inside the
#'   band; `"envelope"` draws lines along the two edges of the band instead.
#'   Default is `"line"`.
#' @param legend_var The name of the variable that should be used in the legend
#'   and to determine the line style. If NULL then the grouping variable is
#'   used for this purpose.
#' @param wrap_var Optional. The name of the variable that should be used for
#'  creating wrapped facets.
#' @param wrap_scale Optional. Used to pass the scales argument to facet_wrap().
#' @param xlab Label for the x-axis
#' @param ylab Label for the y-axis
#' @param xtrans Transformation for the x-axis. Often "log10" may be useful
#'   instead of the default of "identity".
#' @param ytrans Transformation for the y-axis.
#' @param xlim A numeric vector of length two providing lower and upper limits
#'   for the x-axis viewport. Use `NA` for the existing minimum or maximum.
#'   The limits do not filter data.
#' @param ylim A numeric vector of length two providing lower and upper limits
#'   for the y-axis viewport. Use `NA` for the existing minimum or maximum.
#'   The limits do not filter data.
#' @param y_ticks The approximate number of ticks desired on the y axis
#' @param highlight Name or vector of names of the species to be highlighted.
#' @return A ggplot2 object
#' @keywords internal
#' @export
plotDataFrame <- function(frame, params, style = "line", xlab = waiver(),
                          ylab = waiver(), xtrans = "identity", ytrans = "identity",
                          xlim = c(NA, NA), ylim = c(NA, NA),
                          y_ticks = 6, highlight = NULL, legend_var = NULL,
                          wrap_var = NULL, wrap_scale = NULL) {
    spec <- legacy_plot_spec(
        frame, params, legend_var = legend_var, style = style,
        xlab = xlab, ylab = ylab, xtrans = xtrans, ytrans = ytrans,
        xlim = xlim, ylim = ylim, y_ticks = y_ticks, highlight = highlight,
        wrap_var = wrap_var, wrap_scale = wrap_scale
    )
    render_plot_spec(spec)
}

#' Check that a data frame has the variables a band style needs
#'
#' @param var_names The names of the variables in the data frame.
#' @param style The style that was requested.
#' @return Nothing; called for the error.
#' @keywords internal
check_band_vars <- function(var_names, style) {
    missing <- setdiff(c("ymin", "ymax"), var_names)
    if (length(missing) > 0) {
        stop('style = "', style, '" needs the variables ',
             paste0("`", missing, "`", collapse = " and "),
             ", which the data frame does not have.")
    }
}

#' Construct a named vector of line widths for a plot
#'
#' Helper used by the plotting functions to give highlighted species a thicker
#' line than the rest.
#'
#' @param levels Character vector of the legend levels (usually species names).
#' @param highlight Name or vector of names of the levels to be highlighted with
#'   a thicker line.
#' @return A named numeric vector of line widths, one for each entry in
#'   `levels`, with highlighted entries set to a larger value.
#' @keywords internal
make_linesize <- function(levels, highlight) {
    linesize <- rep(0.8, length(levels))
    names(linesize) <- levels
    linesize[highlight] <- 1.6
    linesize
}

#' Tag a ggplot object as a mizer plot
#'
#' Attaches the tooltip information to a ggplot object and adds the
#' `"mizer_plot"` class so that [plotHover()] knows which aesthetics to show.
#' It also gives the plot the layout that fits an open end of its y limits to
#' the data visible within its x limits when the plot is built.
#'
#' @param plot A ggplot object.
#' @param tooltip Character vector of variable names to include in the plotly
#'   tooltip.
#' @return The `plot` object with the tooltip stored as an attribute and
#'   `"mizer_plot"` prepended to its class.
#' @keywords internal
make_mizer_plot <- function(plot, tooltip) {
    plot <- with_mizer_layout(plot)
    attr(plot, "mizer_tooltip") <- tooltip
    class(plot) <- unique(c("mizer_plot", class(plot)))
    plot
}

#' Print a mizer plot
#'
#' Suppresses the uninformative ggplot2 warning about log transformations
#' introducing infinite values, which occurs when zero values are present
#' on a logged axis.
#'
#' @param x A `mizer_plot` object.
#' @param ... Further arguments passed to the ggplot2 print method.
#' @return The plot object, invisibly.
#' @export
print.mizer_plot <- function(x, ...) {
    class(x) <- setdiff(class(x), "mizer_plot")
    withCallingHandlers(
        print(x, ...),
        warning = function(w) {
            if (grepl("transformation introduced infinite values",
                      conditionMessage(w))) {
                tryInvokeRestart("muffleWarning")
            }
        }
    )
    invisible(x)
}

#' Determine the tooltip variables for a mizer plot
#'
#' Works out which variables should appear in the plotly tooltip, including the
#' legend variable only when it differs from the grouping variable, plus any
#' extra variables.
#'
#' @param frame The data frame underlying the plot.
#' @param group_var Name of the grouping variable.
#' @param x_var Name of the variable on the x-axis.
#' @param y_var Name of the variable on the y-axis.
#' @param legend_var Optional name of the legend variable.
#' @param extra Optional character vector of additional variable names to
#'   include.
#' @return A character vector of variable names for the tooltip.
#' @keywords internal
mizer_tooltip_vars <- function(frame, group_var, x_var, y_var,
                               legend_var = NULL, extra = NULL) {
    tooltip <- c(group_var, x_var, y_var)
    if (!is.null(legend_var) && legend_var %in% names(frame) &&
            !identical(legend_var, group_var) &&
            any(as.character(frame[[legend_var]]) !=
                    as.character(frame[[group_var]]), na.rm = TRUE)) {
        tooltip <- c(tooltip, legend_var)
    }
    unique(c(tooltip, extra))
}

#' Create a hover-enabled plotly plot from a mizer object
#'
#' Creates an interactive plotly version of a mizer plot. Can be called on any
#' mizer array object (such as those returned by [getEncounter()],
#' [getBiomass()], etc.) or on any `mizer_plot` object returned by the named
#' plot functions such as [plotBiomass()], [plotSpectra()], etc.
#'
#' @param x A `mizer_plot`, `ArraySpeciesBySize`, `ArrayTimeBySpecies`, or
#'   `ArrayTimeBySpeciesBySize` object.
#' @param ... Arguments passed to the corresponding `plot()` method for mizer
#'   array objects, or to [plotly::ggplotly()] for `mizer_plot` objects.
#' @return A plotly object.
#' @seealso [plot()], [plotBiomass()], [plotSpectra()], [plotting_functions]
#' @concept plotting functions
#' @export
plotHover <- function(x, ...) UseMethod("plotHover")

# Plotly keys the data held inside each widget with an id returned by
# `basename(tempfile())`. That makes two otherwise identical widgets differ
# across R sessions and causes their pkgdown pages to be rewritten on every
# build. The ids are private to a widget, so replace them with a stable sequence
# while preserving all links between the internal lists.
stable_plotly_ids <- function(widget) {
    keyed_fields <- c("attrs", "visdat", "layoutAttrs")
    old_ids <- unique(unlist(lapply(keyed_fields, function(field) {
        names(widget$x[[field]])
    }), use.names = FALSE))
    old_ids <- old_ids[nzchar(old_ids)]
    if (length(old_ids) == 0) return(widget)

    new_ids <- paste0("mizer-plotly-", seq_along(old_ids))
    for (field in keyed_fields) {
        field_ids <- names(widget$x[[field]])
        if (length(field_ids) > 0) {
            names(widget$x[[field]]) <- new_ids[match(field_ids, old_ids)]
        }
    }
    current <- match(widget$x$cur_data, old_ids)
    if (length(current) == 1 && !is.na(current)) {
        widget$x$cur_data <- new_ids[[current]]
    }
    widget
}

#' @rdname plotHover
#' @usage NULL
#' @export
plotHover.mizer_plot <- function(x = ggplot2::last_plot(), ...,
                                 tooltip = attr(x, "mizer_tooltip") %||%
                                     "all") {
    class(x) <- setdiff(class(x), "mizer_plot")
    widget <- withCallingHandlers(
        plotly::ggplotly(hover_viewport(x), ..., tooltip = tooltip),
        warning = function(w) {
            if (grepl("transformation introduced infinite values",
                      conditionMessage(w))) {
                tryInvokeRestart("muffleWarning")
            }
        }
    )
    stable_plotly_ids(widget)
}

# Every `plotly...()` wrapper ends with a call to `plotly_forward()`, which
# calls the wrapper's plot function with the wrapper's arguments and converts
# the result with `plotHover()`. The wrapper keeps its documented signature
# while the forwarding is implemented only here. It reads the wrapper's formals
# and arguments from the calling frame, so it must be called directly from the
# body of the wrapper.
#
# Each formal of the calling wrapper other than `...` and those in `exclude` is
# forwarded by name with its current value. That value is the wrapper's own
# default if the argument was not supplied, so a wrapper default that differs
# from the plot function's is kept, as is any change the wrapper made to the
# value before calling this. An argument that was not supplied and has no
# default, such as `time_range`, is left missing, so that the plot function's
# own handling of a missing argument applies, and so is an argument that was
# not supplied and whose default is `lifecycle::deprecated()`.
#
# plot_fun: The plot function, given by name so that its body is not spelled
#   out in the call shown by errors and tracebacks.
# tooltip: The tooltip variables passed to `plotHover()`. `NULL` keeps the
#   tooltip that the plot function stored with the plot.
# exclude: Names of formals of the wrapper that are not forwarded.
# extra: A named list of further arguments for the plot function.
# dots: Where the wrapper's `...` go: to the plot function (`"plot"`) or to
#   `plotHover()` (`"hover"`).
plotly_forward <- function(plot_fun, tooltip = NULL, exclude = character(),
                           extra = list(), dots = c("plot", "hover")) {
    dots <- match.arg(dots)
    env <- parent.frame()
    fmls <- formals(sys.function(sys.parent()))
    skip_if_missing <- vapply(fmls, function(default) {
        identical(default, quote(expr = )) ||
            identical(default, quote(lifecycle::deprecated()))
    }, logical(1))
    args <- list()
    for (name in setdiff(names(fmls), c("...", exclude))) {
        if (skip_if_missing[[name]] &&
                eval(call("missing", as.name(name)), env)) {
            next
        }
        args[name] <- list(get(name, envir = env, inherits = FALSE))
    }
    dot_args <- list()
    if ("..." %in% names(fmls)) {
        dot_args <- eval(quote(list(...)), env)
    }
    hover_args <- list()
    if (dots == "hover") {
        hover_args <- dot_args
        dot_args <- list()
    }
    p <- do.call(plot_fun, c(args, extra, dot_args))
    if (!is.null(tooltip)) {
        hover_args$tooltip <- tooltip
    }
    do.call("plotHover", c(list(p), hover_args))
}

#' Make a plot comparing two data frames
#'
#' Used internally by the comparison plotting functions such as [plotSpectra2()]
#' and [plotCDF2()]. The two data frames are combined and drawn with colour
#' identifying the species or group and linetype identifying the object.
#'
#' Both data frames must arrive ready to plot: on the axis they will be drawn
#' against, with any total line already among their rows. Each operand is
#' prepared by whatever produced it, using *its own* model, because a length
#' axis and a density Jacobian both depend on the weight-length relationship of
#' the model the values came from. Doing it here instead would silently impose
#' the first model's parameters on the second.
#'
#' @param frame1,frame2 Data frames whose first three variables are the x, y and
#'   grouping variables. The columns of `frame2` are matched to those of
#'   `frame1` by these roles and by `legend_var`, not by name.
#' @param params A MizerParams object, used for the line colours.
#' @param name1,name2 Labels for the two data frames, used in the linetype
#'   legend.
#' @param xlab,ylab Labels for the x and y axes.
#' @param xtrans,ytrans Transformations for the x and y axes, e.g. `"log10"` or
#'   `"identity"`.
#' @param xlim,ylim Numeric vectors of length two giving the axis limits. Use
#'   `NA` to refer to the existing minimum or maximum.
#' @param y_ticks The approximate number of ticks desired on the y axis.
#' @param highlight Name or vector of names of the species to be highlighted.
#' @param legend_var Name of the variable used in the legend and to determine the
#'   line colour.
#' @return A `mizer_plot` (ggplot2) object.
#' @keywords internal
plotComparisonDataFrame <- function(frame1, frame2, params,
                                    name1 = "First", name2 = "Second",
                                    xlab = waiver(), ylab = waiver(),
                                    xtrans = "identity", ytrans = "identity",
                                    xlim = c(NA, NA), ylim = c(NA, NA),
                                    y_ticks = 6, highlight = NULL,
                                    legend_var = "Legend") {
    assert_that(is.data.frame(frame1), is.data.frame(frame2),
                is(params, "MizerParams"))
    first <- legacy_plot_spec(
        frame1, params, legend_var = legend_var,
        xlab = xlab, ylab = ylab, xtrans = xtrans, ytrans = ytrans,
        xlim = xlim, ylim = ylim, y_ticks = y_ticks, highlight = highlight
    )
    second <- legacy_plot_spec(frame2, params, legend_var = legend_var)
    render_plot_spec(compare_plot_specs(first, second, name1, name2))
}

#' Make a plot of the relative difference between two data frames
#'
#' Used internally by [plotSpectraRelative()] and similar functions. The two
#' data frames are matched up on their shared variables and the relative
#' difference of their y-values is plotted against the x-variable.
#'
#' Both data frames must arrive ready to plot, on the axis they will be drawn
#' against and with any total line already among their rows. See
#' [plotComparisonDataFrame()] for why the preparation belongs to whatever
#' produced each operand rather than here.
#'
#' @param frame1,frame2 Data frames whose first three variables are the x, y and
#'   grouping variables. The columns of `frame2` are matched to those of
#'   `frame1` by these roles and by `legend_var`, not by name.
#' @param params A MizerParams object, used for the line colours.
#' @param xlab Label for the x-axis.
#' @param xtrans Transformation for the x-axis, e.g. `"log10"` or `"identity"`.
#' @param xlim,ylim Numeric vectors of length two giving the axis limits. Use
#'   `NA` to refer to the existing minimum or maximum.
#' @param highlight Name or vector of names of the species to be highlighted.
#' @param legend_var Name of the variable used in the legend and to determine the
#'   line colour.
#' @param interpolate Whether the two series may sit on different x-grids and
#'   should be interpolated onto a common one, see
#'   [interpolate_relative_frames()]. `TRUE` for a size axis, where the two
#'   models can convert weight to length differently; `FALSE` for a time axis,
#'   where they share the saved times or share nothing.
#' @return A `mizer_plot` (ggplot2) object showing the relative difference.
#' @keywords internal
plotRelativeDataFrame <- function(frame1, frame2, params,
                                  xlab = waiver(),
                                  xtrans = "identity",
                                  xlim = c(NA, NA),
                                  ylim = c(NA, NA),
                                  highlight = NULL,
                                  legend_var = "Legend",
                                  interpolate = FALSE) {
    assert_that(is.data.frame(frame1), is.data.frame(frame2),
                is(params, "MizerParams"))
    first <- legacy_plot_spec(
        frame1, params, legend_var = legend_var,
        xlab = xlab, xtrans = xtrans, xlim = xlim, ylim = ylim,
        highlight = highlight
    )
    second <- legacy_plot_spec(frame2, params, legend_var = legend_var)
    render_plot_spec(relative_plot_specs(first, second, interpolate))
}

#' Symmetric relative difference between two values
#'
#' @param first,second Numeric vectors to compare.
#' @return The symmetric relative difference
#'   `2 * (second - first) / (first + second)`.
#' @keywords internal
relative_difference <- function(first, second) {
    2 * (second - first) / (first + second)
}

#' Put two series on a common size grid before comparing them
#'
#' A relative difference can only be formed where both series have a value. On a
#' weight axis they always do: both models share the size grid, so the two
#' frames match up row for row. On a length axis they need not, because each
#' model converts weight to length with its own allometric relationship, so the
#' same weight grid lands on different lengths. Matching the frames by equality
#' of the size coordinate then throws away nearly every point — an inner join on
#' two grids that merely overlap keeps only their exact coincidences.
#'
#' Each series is therefore interpolated, linearly in the logarithm of size
#' since the grid is logarithmic, onto the sorted union of the two sets of
#' coordinates, restricted to the interval both series cover. Outside that
#' interval one of them would have to be extrapolated, which is not a comparison
#' but a guess. When the two grids already coincide the union is that grid, the
#' overlap is all of it, and the interpolation reproduces the values exactly, so
#' the matching case is unchanged.
#'
#' @param frame1,frame2 Data frames of prepared plotting data, sharing their
#'   variable names.
#' @param x_var Name of the size column.
#' @param y_var Name of the value column.
#' @param by_vars Names of the columns identifying a series, typically the
#'   species and the legend group.
#' @return A data frame with the `by_vars`, the size column, and the two value
#'   columns named `<y_var>.x` and `<y_var>.y`, holding only the series present
#'   in both frames.
#' @keywords internal
interpolate_relative_frames <- function(frame1, frame2, x_var, y_var, by_vars) {
    series_key <- function(frame) {
        do.call(paste, c(lapply(frame[, by_vars, drop = FALSE], as.character),
                         sep = "\r"))
    }
    split1 <- split(frame1, series_key(frame1))
    split2 <- split(frame2, series_key(frame2))
    parts <- lapply(intersect(names(split1), names(split2)), function(key) {
        first <- split1[[key]]
        second <- split2[[key]]
        x1 <- first[[x_var]]
        x2 <- second[[x_var]]
        if (length(x1) == 0 || length(x2) == 0) return(NULL)
        xout <- sort(unique(c(x1, x2)))
        xout <- xout[xout >= max(min(x1), min(x2)) &
                         xout <= min(max(x1), max(x2))]
        if (length(xout) == 0) return(NULL)
        out <- first[rep(1, length(xout)), by_vars, drop = FALSE]
        out[[x_var]] <- xout
        out[[paste0(y_var, ".x")]] <-
            interpolate_in_log_size(x1, first[[y_var]], xout)
        out[[paste0(y_var, ".y")]] <-
            interpolate_in_log_size(x2, second[[y_var]], xout)
        out
    })
    parts <- parts[!vapply(parts, is.null, logical(1))]
    if (length(parts) == 0) {
        out <- frame1[0, by_vars, drop = FALSE]
        out[[x_var]] <- numeric(0)
        out[[paste0(y_var, ".x")]] <- numeric(0)
        out[[paste0(y_var, ".y")]] <- numeric(0)
        return(out)
    }
    out <- do.call(rbind, parts)
    rownames(out) <- NULL
    out
}

#' Interpolate a series linearly in the logarithm of size
#'
#' The size grid is logarithmic, so the interpolation is too. A series of a
#' single point can only speak for the coordinate it sits at, and a coordinate
#' that is not positive has no logarithm, so those two cases fall back to
#' matching the coordinates exactly and to a linear interpolation respectively.
#'
#' @param x,y The coordinates and values of the series.
#' @param xout The coordinates to interpolate onto.
#' @return The interpolated values, `NA` where `xout` is outside the range of
#'   `x`.
#' @keywords internal
interpolate_in_log_size <- function(x, y, xout) {
    order_x <- order(x)
    x <- x[order_x]
    y <- y[order_x]
    if (length(x) < 2) {
        out <- rep(NA_real_, length(xout))
        hit <- match(x, xout)
        if (length(hit) == 1 && !is.na(hit)) out[hit] <- y
        return(out)
    }
    if (any(x <= 0) || any(xout <= 0)) {
        return(stats::approx(x, y, xout = xout)$y)
    }
    stats::approx(log(x), y, xout = log(xout))$y
}

#' Helper function to produce nice breaks on logarithmic axes
#'
#' This is needed when the logarithmic y-axis spans less than one order of
#' magnitude, in which case the ggplot2 default produces no ticks.
#'
#' Thanks to Heather Turner at
#' \url{https://stackoverflow.com/questions/14255533/pretty-ticks-for-log-normal-scale-using-ggplot2-dynamic-not-manual}
#'
#' @param n Approximate number of ticks
#'
#' @return A function that can be used as the break argument in calls to
#'   scale_y_continuous() or scale_x_continuous()
#' @export
#' @keywords internal
#' @family helper
log_breaks <- function(n = 6) {
    n <- max(1, n)  # Because n=0 could lead to R crash
    function(x) {
        grDevices::axisTicks(log10(range(x, na.rm = TRUE)),
                             log = TRUE, nint = n)
    }
}
