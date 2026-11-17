# Internal plotting contract. Data arrive fully prepared: selection, totals,
# density conversion and integration belong to the producer, not the renderer.
# This is a plain list, never a new public class or a return_data schema.
# Every renderer reads its columns through the named roles in `variables`.
# Column positions matter only to the data-frame adapters (plotDataFrame() and
# its comparison and relative counterparts) and to the public layout of the
# data that `return_data = TRUE` returns, see plot_data_value_var().
plot_spec_modes <- c("ordinary", "comparison", "relative")
plot_spec_styles <- c("line", "area", "ribbon", "envelope")

new_plot_spec <- function(data, params, x_var, y_var, group_var,
                          legend_var = group_var, frame_var = NULL,
                          xlab = waiver(), ylab = waiver(),
                          xtrans = "identity", ytrans = "identity",
                          xlim = c(NA, NA), ylim = c(NA, NA),
                          style = "line", y_ticks = 6, highlight = NULL,
                          wrap_var = NULL, wrap_scale = NULL,
                          tooltip = NULL, options = list()) {
    variables <- list(x = x_var, value = y_var, group = group_var,
                      legend = legend_var, frame = frame_var)
    missing <- setdiff(unlist(variables, use.names = FALSE), names(data))
    if (length(missing) > 0) {
        stop("Plot variables missing from data: ",
             paste(missing, collapse = ", "))
    }
    if (!style %in% plot_spec_styles) {
        stop("unknown style selected")
    }
    if (!is.null(wrap_var) && !(wrap_var %in% names(data))) {
        stop("The `wrap_var` argument must be the name of a variable ",
             "in the data frame.")
    }
    options <- utils::modifyList(
        list(mode = "ordinary", style = style, y_ticks = y_ticks,
             highlight = highlight, wrap_var = wrap_var,
             wrap_scale = wrap_scale),
        options
    )
    list(data = data, params = params, variables = variables,
         labels = list(x = xlab, y = ylab),
         scales = list(x = xtrans, y = ytrans),
         limits = list(x = xlim, y = ylim),
         tooltip = tooltip, options = options)
}

# Most mizer plots draw one line for each `Species`, coloured by its `Legend`
# entry, so only the size or time coordinate and the value need naming.
species_plot_spec <- function(data, params, x_var, y_var, ...) {
    new_plot_spec(data, params, x_var = x_var, y_var = y_var,
                  group_var = "Species", legend_var = "Legend", ...)
}

# The plotting data that `return_data = TRUE` returns keep their public layout,
# with the value in the second column under its descriptive name. Producers
# that do not set that name themselves look it up here.
plot_data_value_var <- function(plot_dat) names(plot_dat)[[2]]

# The positional contract is confined to the compatibility boundary. In
# particular, keep the public column names and the automatic Legend column.
legacy_plot_spec <- function(frame, params, legend_var = NULL,
                             wrap_var = NULL, ...) {
    assert_that(is.data.frame(frame), is(params, "MizerParams"))
    if (ncol(frame) < 3) {
        stop("The data frame needs to have at least 3 variables.")
    }
    variables <- names(frame)[1:3]
    # A facet must be one of the caller's own columns, not the Legend column
    # that this adapter adds.
    if (!is.null(wrap_var) && !(wrap_var %in% names(frame))) {
        stop("The `wrap_var` argument must be the name of a variable ",
             "in the data frame.")
    }
    if (is.null(legend_var)) {
        frame$Legend <- frame[[variables[[3]]]]
        legend_var <- "Legend"
    } else if (!(legend_var %in% names(frame))) {
        stop("The `legend_var` argument must be the name of a variable ",
             "in the data frame.")
    }
    new_plot_spec(frame, params, x_var = variables[[1]],
                  y_var = variables[[2]], group_var = variables[[3]],
                  legend_var = legend_var, wrap_var = wrap_var, ...)
}

# Align semantic roles, not column positions. Each producer has already used
# its own model to prepare its operand, including its size axis and density.
align_plot_spec_data <- function(first, second) {
    roles <- c("x", "value", "group", "legend")
    from <- unlist(second$variables[roles], use.names = FALSE)
    to <- unlist(first$variables[roles], use.names = FALSE)
    data <- second$data
    positions <- match(names(data), from)
    rename <- !is.na(positions)
    names(data)[rename] <- to[positions[rename]]
    data
}

compare_plot_specs <- function(first, second, name1, name2) {
    data1 <- first$data
    data2 <- align_plot_spec_data(first, second)
    data1$Model <- name1
    data2$Model <- name2
    first$data <- rbind(data1, data2)
    first$data$Model <- factor(first$data$Model, levels = c(name1, name2))
    first$variables$model <- "Model"
    first$options$mode <- "comparison"
    first$tooltip <- NULL
    first
}

relative_plot_specs <- function(first, second, interpolate = FALSE) {
    vars <- first$variables
    by_vars <- unique(c(vars$group, vars$legend))
    data2 <- align_plot_spec_data(first, second)
    if (isTRUE(interpolate)) {
        data <- interpolate_relative_frames(
            first$data, data2, x_var = vars$x, y_var = vars$value,
            by_vars = by_vars
        )
    } else {
        data <- dplyr::inner_join(first$data, data2,
                                  by = c(vars$x, by_vars),
                                  suffix = c(".x", ".y"))
    }
    data$rel_diff <- relative_difference(data[[paste0(vars$value, ".x")]],
                                         data[[paste0(vars$value, ".y")]])
    first$data <- data[is.finite(data$rel_diff), ]
    first$variables$value <- "rel_diff"
    first$labels$y <- "Relative difference"
    first$scales$y <- "identity"
    first$options$mode <- "relative"
    first$tooltip <- NULL
    first
}

plot_spec_legend_levels <- function(spec) {
    intersect(names(spec$params@linecolour),
              spec$data[[spec$variables$legend]])
}

plot_spec_tooltip <- function(spec) {
    if (!is.null(spec$tooltip)) return(spec$tooltip)
    vars <- spec$variables
    if (identical(spec$options$mode, "comparison")) {
        return(mizer_tooltip_vars(spec$data, vars$legend, vars$x, vars$value,
                                  extra = vars$model))
    }
    mizer_tooltip_vars(spec$data, vars$group, vars$x, vars$value, vars$legend)
}

# Axis transforms belong to scales; limits belong to the viewport so that
# clipping cannot remove vertices or change a statistic calculated by a geom.
plot_spec_scales <- function(spec) {
    xbreaks <- if (spec$scales$x == "log10") log_breaks() else waiver()
    ybreaks <- if (spec$scales$y == "log10") {
        log_breaks(n = spec$options$y_ticks)
    } else {
        waiver()
    }
    ylabels <- if (spec$options$mode == "relative") waiver() else prettyNum
    list(scale_y_continuous(trans = spec$scales$y, breaks = ybreaks,
                            labels = ylabels, name = spec$labels$y),
         scale_x_continuous(trans = spec$scales$x, breaks = xbreaks,
                            name = spec$labels$x),
         coord_cartesian(xlim = spec$limits$x, ylim = spec$limits$y))
}

# Limits only move the viewport; they never remove data. An end of the y
# viewport that the caller leaves as NA is fitted to what is visible inside the
# x viewport, as it was when limits filtered the data. Otherwise zooming in on
# x would leave the y axis spanning values that are no longer in view.
#
# The viewport records the requested limits and every data source drawn on the
# plot, so that addPlot() can refit the open ends when it adds a layer. A
# source names its columns by role: `y` lists every column that sets a height,
# `group` every column that separates one line from another, `stack` marks
# stacked areas and `include` holds values that are always shown, such as the
# zero line of a relative plot. `data = NULL` refers to the plot's own data.
plot_viewport <- function(xlim, ylim, xtrans, ytrans, sources) {
    list(x = xlim, y = ylim, xtrans = xtrans, ytrans = ytrans,
         sources = sources)
}

viewport_source <- function(x, y, group = character(), data = NULL,
                            stack = FALSE, include = numeric()) {
    list(data = data, x = x, y = y, group = group, stack = stack,
         include = include)
}

viewport_ylim <- function(viewport, plot_data) {
    ylim <- viewport$y
    if (!anyNA(ylim) || all(is.na(viewport$x))) return(ylim)
    values <- unlist(lapply(viewport$sources, function(source) {
        visible_values(source$data %||% plot_data, source, viewport)
    }))
    if (identical(viewport$ytrans, "log10")) values <- values[values > 0]
    values <- values[is.finite(values)]
    if (length(values) == 0) return(ylim)
    fitted <- ylim
    fitted[is.na(ylim)] <- range(values)[is.na(ylim)]
    # An open end cannot cross a requested one; leave it to ggplot2.
    if (fitted[1] > fitted[2]) return(ylim)
    fitted
}

# The heights a source shows inside the x viewport: its vertices there and the
# points where its lines cross the edges. Lines are straight on the transformed
# axes, so the crossings are interpolated on those axes.
visible_values <- function(data, source, viewport) {
    if (nrow(data) == 0) return(source$include)
    to_axis <- function(v, trans) {
        if (identical(trans, "log10")) suppressWarnings(log10(v)) else v
    }
    edges <- to_axis(viewport$x, viewport$xtrans)
    lower <- if (is.na(edges[1])) -Inf else edges[1]
    upper <- if (is.na(edges[2])) Inf else edges[2]
    crossings <- edges[!is.na(edges)]
    x <- data[[source$x]]
    heights <- lapply(source$y, function(var) data[[var]])
    group <- if (length(source$group) > 0) {
        interaction(data[source$group], drop = TRUE)
    } else {
        factor(rep(1L, nrow(data)))
    }
    if (isTRUE(source$stack)) {
        # Stacked areas reach the sum of their values at each x.
        x_levels <- sort(unique(x))
        heights <- list(vapply(split(heights[[1]],
                                     factor(x, levels = x_levels)),
                               sum, numeric(1), na.rm = TRUE))
        x <- x_levels
        group <- factor(rep(1L, length(x)))
    }
    x <- to_axis(x, viewport$xtrans)
    values <- numeric(0)
    for (height in heights) {
        y <- to_axis(height, viewport$ytrans)
        ok <- is.finite(x) & is.finite(y)
        values <- c(values, y[ok & x >= lower & x <= upper])
        if (length(crossings) == 0) next
        for (rows in split(which(ok), group[ok])) {
            if (length(unique(x[rows])) < 2) next
            values <- c(values, stats::approx(x[rows], y[rows],
                                              xout = crossings,
                                              ties = mean)$y)
        }
    }
    values <- values[is.finite(values)]
    if (identical(viewport$ytrans, "log10")) values <- 10^values
    c(values, source$include)
}

# The source describing a specification's own data.
plot_spec_source <- function(spec) {
    vars <- spec$variables
    mode <- spec$options$mode
    style <- if (mode == "ordinary") spec$options$style else "line"
    y <- switch(style,
                ribbon = c(vars$value, "ymin", "ymax"),
                envelope = c("ymin", "ymax"),
                vars$value)
    group <- unique(c(vars$group, vars$legend, vars$model, vars$frame,
                      spec$options$wrap_var))
    viewport_source(x = vars$x, y = y, group = group,
                    stack = style == "area",
                    include = if (mode == "relative" || style == "area") 0)
}

# Added layers share the existing viewport. Only supplied endpoints change it,
# and the open ends of the y viewport are refitted to include `source`. The
# plot is returned with a new coordinate system; the caller's is not modified.
set_plot_viewport <- function(plot, xlim = c(NA, NA), ylim = c(NA, NA),
                              source = NULL) {
    assert_that(length(xlim) == 2, length(ylim) == 2)
    coord <- plot$coordinates
    merge_limits <- function(existing, requested, axis) {
        existing[!is.na(requested)] <- requested[!is.na(requested)]
        if (plot_axis_is_log(plot, axis)) {
            existing <- drop_nonpositive_limits(existing, "log10")
        }
        existing
    }
    viewport <- attr(plot, "mizer_viewport")
    if (is.null(viewport)) {
        # A plot mizer did not lay out: its other layers are unknown, so only
        # the supplied endpoints change.
        limits <- list(
            x = merge_limits(coord$limits$x %||% c(NA, NA), xlim, "x"),
            y = merge_limits(coord$limits$y %||% c(NA, NA), ylim, "y")
        )
    } else {
        viewport$x <- merge_limits(viewport$x, xlim, "x")
        viewport$y <- merge_limits(viewport$y, ylim, "y")
        if (!is.null(source)) {
            viewport$sources <- c(viewport$sources, list(source))
        }
        limits <- list(x = viewport$x, y = viewport_ylim(viewport, plot$data))
        attr(plot, "mizer_viewport") <- viewport
    }
    plot$coordinates <- ggproto(NULL, coord, limits = limits)
    plot
}

render_plot_spec <- function(spec) {
    vars <- spec$variables
    mode <- spec$options$mode
    if (!mode %in% plot_spec_modes) {
        stop("Unknown plot specification mode `", mode, "`.")
    }
    frame <- spec$data
    legend_levels <- plot_spec_legend_levels(spec)
    frame[[vars$legend]] <- factor(frame[[vars$legend]],
                                   levels = legend_levels)
    if (anyNA(frame[[vars$legend]])) {
        warning("missing legend in params@linecolour, ",
                "some groups won't be displayed")
    }
    spec$limits$x <- drop_nonpositive_limits(spec$limits$x, spec$scales$x)
    spec$limits$y <- drop_nonpositive_limits(spec$limits$y, spec$scales$y)
    if (mode != "relative") {
        frame <- drop_nonpositive_values(frame, vars$x, spec$scales$x)
        frame <- drop_nonpositive_values(frame, vars$value, spec$scales$y)
    }
    if (mode != "ordinary") frame$LineSpec <- frame[[vars$legend]]
    spec$data <- frame
    viewport <- plot_viewport(spec$limits$x, spec$limits$y,
                              spec$scales$x, spec$scales$y,
                              list(plot_spec_source(spec)))
    spec$limits$y <- viewport_ylim(viewport, frame)

    p <- ggplot(frame, plot_spec_mapping(spec)) + plot_spec_scales(spec) +
        plot_spec_layers(spec, legend_levels)
    wrap_var <- spec$options$wrap_var
    if (!is.null(wrap_var)) {
        p <- p + facet_wrap(wrap_var, scales = spec$options$wrap_scale)
    }
    p <- make_mizer_plot(p, plot_spec_tooltip(spec))
    attr(p, "mizer_viewport") <- viewport
    p
}

# Preserve mapping order: ggplotly uses it to order the hover fields.
plot_spec_mapping <- function(spec) {
    vars <- spec$variables
    if (spec$options$mode == "comparison") {
        return(aes(group = interaction(.data[[vars$legend]],
                                       .data[[vars$model]]),
                   colour = .data[[vars$legend]]))
    }
    aes(group = .data[[vars$group]])
}

# The geoms of each mode, followed by the scales that style them. Comparison
# and relative plots distinguish models by linetype and draw the legend
# without line widths; ordinary plots carry the full species styling.
plot_spec_layers <- function(spec, legend_levels) {
    vars <- spec$variables
    x_var <- vars$x
    y_var <- vars$value
    legend_var <- vars$legend
    linecolour <- spec$params@linecolour[legend_levels]
    linesize <- make_linesize(legend_levels, spec$options$highlight)
    hidden_linewidth <- scale_discrete_manual("linewidth", values = linesize,
                                              guide = "none")
    switch(
        spec$options$mode,
        comparison = list(
            geom_line(aes(x = .data[[x_var]], y = .data[[y_var]],
                          linetype = .data[[vars$model]],
                          linewidth = .data[["LineSpec"]])),
            scale_colour_manual(values = linecolour),
            scale_linetype_discrete(drop = FALSE),
            hidden_linewidth
        ),
        relative = list(
            geom_hline(yintercept = 0, linetype = 1,
                       colour = "dark grey", linewidth = 0.75),
            geom_line(aes(x = .data[[x_var]], y = .data[[y_var]],
                          colour = .data[[legend_var]],
                          linewidth = .data[["LineSpec"]])),
            scale_colour_manual(values = linecolour),
            hidden_linewidth
        ),
        ordinary = ordinary_plot_layers(spec, legend_levels, linecolour,
                                        linesize)
    )
}

ordinary_plot_layers <- function(spec, legend_levels, linecolour, linesize) {
    vars <- spec$variables
    x_var <- vars$x
    legend_var <- vars$legend
    style <- spec$options$style
    if (style == "area") {
        return(list(
            geom_area(aes(x = .data[[x_var]], y = .data[[vars$value]],
                          fill = .data[[legend_var]])),
            scale_fill_manual(values = linecolour)
        ))
    }
    band <- list()
    if (style != "line") {
        check_band_vars(names(spec$data), style)
        band <- list(
            geom_ribbon(aes(x = .data[[x_var]], ymin = .data[["ymin"]],
                            ymax = .data[["ymax"]],
                            fill = .data[[legend_var]]),
                        alpha = 0.25, colour = NA),
            scale_fill_manual(values = linecolour)
        )
    }
    # A function rather than an inline aes(), so that each layer keeps its own
    # binding of yv, including both edges of an envelope.
    value_line <- function(yv) {
        geom_line(aes(x = .data[[x_var]], y = .data[[yv]],
                      colour = .data[[legend_var]],
                      linetype = .data[[legend_var]],
                      linewidth = .data[[legend_var]]))
    }
    y_vars <- if (style == "envelope") c("ymax", "ymin") else vars$value
    linetype <- spec$params@linetype[legend_levels]
    c(band, lapply(y_vars, value_line),
      list(scale_colour_manual(values = linecolour),
           scale_linetype_manual(values = linetype),
           scale_discrete_manual("linewidth", values = linesize)))
}

# Animation consumes the same named contract. Static hover plots continue to
# convert the actual ggplot via plotHover(), including user-added layers and
# scales, rather than reconstructing a plot from an earlier specification.
render_plotly_spec <- function(spec) {
    vars <- spec$variables
    df <- spec$data
    legend_levels <- plot_spec_legend_levels(spec)
    group_order <- unlist(lapply(legend_levels, function(legend) {
        intersect(names(spec$params@linecolour),
                  unique(df[[vars$group]][df[[vars$legend]] == legend]))
    }))
    formula_for <- function(variable) {
        rlang::new_formula(NULL, as.name(variable))
    }
    p <- plotly::plot_ly()
    shown_legends <- character(0)
    for (group in group_order) {
        group_data <- df[df[[vars$group]] == group, , drop = FALSE]
        legend <- unique(group_data[[vars$legend]])
        p <- plotly::add_lines(
            p, data = group_data,
            x = formula_for(vars$x), y = formula_for(vars$value),
            frame = formula_for(vars$frame),
            name = legend, legendgroup = legend,
            line = list(color = spec$params@linecolour[[legend]],
                        simplify = FALSE),
            showlegend = !(legend %in% shown_legends)
        )
        shown_legends <- c(shown_legends, legend)
    }
    viewport <- plot_viewport(spec$limits$x, spec$limits$y,
                              spec$scales$x, spec$scales$y,
                              list(plot_spec_source(spec)))
    p <- plotly::layout(
        p,
        xaxis = plotly_axis(df[[vars$x]], spec$limits$x,
                            spec$scales$x == "log10", spec$labels$x),
        yaxis = plotly_axis(df[[vars$value]], viewport_ylim(viewport, df),
                            spec$scales$y == "log10", spec$labels$y),
        legend = list(traceorder = "normal")
    )
    p <- plotly::animation_opts(p, frame = spec$options$frame_duration,
                                transition = spec$options$transition_duration,
                                easing = spec$options$easing)
    stable_plotly_ids(p)
}
