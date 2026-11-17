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
# viewport that the caller leaves open is fitted to what is visible inside the
# x viewport, as it was when limits filtered the data. Otherwise zooming in on
# x would leave the y axis spanning values that are no longer in view.
#
# The fit happens when the plot is built, so it sees what is actually drawn:
# every layer, including those a user adds with `+` after mizer has returned
# the plot, after stats and positions have run (so stacked areas reach the top
# of the stack ggplot2 draws) and separately for each free y scale of a facet.
# ggplot2 builds a plot with the Layout object stored in it, so every mizer plot
# carries a MizerLayout.
mizer_layout_cache <- new.env(parent = emptyenv())

mizer_layout <- function() {
    mizer_layout_cache$layout <- mizer_layout_cache$layout %||%
        new_mizer_layout()
    mizer_layout_cache$layout
}

# Created when first needed rather than when mizer is built, so that it extends
# the Layout of the installed ggplot2.
new_mizer_layout <- function() {
    parent <- get("Layout", envir = asNamespace("ggplot2"))
    ggproto("MizerLayout", parent,
        train_position = function(self, data, x_scale, y_scale) {
            ggproto_parent(parent, self)$train_position(data, x_scale, y_scale)
            # The last training sees the data as drawn.
            self$position_data <- data
        },
        # The first hook that knows the layers, which tell lines from points.
        # Guides are made from the panel parameters, so fit before that.
        setup_panel_guides = function(self, guides, layers) {
            fit_panel_ylim(self, layers)
            ggproto_parent(parent, self)$setup_panel_guides(guides, layers)
        }
    )
}

with_mizer_layout <- function(plot) {
    plot$layout <- mizer_layout()
    plot
}

# Refit the open ends of the y limits of each panel to the heights visible in
# its x range. Panels sharing a y scale share the fit, so a fixed scale stays
# fixed across facets.
fit_panel_ylim <- function(layout, layers) {
    coord <- layout$coord
    if (!inherits(coord, "CoordCartesian") || inherits(coord, "CoordFlip")) {
        return(invisible())
    }
    xlim <- coord$limits$x %||% c(NA, NA)
    ylim <- coord$limits$y %||% c(NA, NA)
    # Without an x limit everything is visible, so ggplot2 chooses.
    if (all(is.na(xlim)) || !anyNA(ylim)) return(invisible())
    scales_y <- layout$panel_scales_y
    if (length(scales_y) == 0 || scales_y[[1]]$is_discrete() ||
            layout$panel_scales_x[[1]]$is_discrete()) {
        return(invisible())
    }
    x_aes <- layout$panel_scales_x[[1]]$aesthetics
    y_aes <- scales_y[[1]]$aesthetics
    path <- vapply(layers, function(layer) {
        inherits(layer$geom, c("GeomPath", "GeomRibbon", "GeomSmooth"))
    }, logical(1))
    panels <- layout$layout
    data <- layout$position_data
    values <- rep(list(numeric()), length(scales_y))
    for (i in seq_len(nrow(panels))) {
        window <- layout$panel_params[[i]]$x.range
        s <- panels$SCALE_Y[[i]]
        for (l in seq_along(data)) {
            d <- data[[l]][data[[l]]$PANEL == panels$PANEL[[i]], ,
                           drop = FALSE]
            visible <- visible_heights(d, window, x_aes, y_aes, path[[l]])
            values[[s]] <- c(values[[s]], visible)
        }
    }
    fitted <- FALSE
    for (s in seq_along(scales_y)) {
        trans <- scales_y[[s]]$get_transformation()
        limits <- fit_open_limits(values[[s]], ylim, trans$transform,
                                  trans$inverse)
        if (is.null(limits)) next
        fitted_coord <- ggproto(NULL, coord,
                                limits = list(x = coord$limits$x, y = limits))
        for (i in which(panels$SCALE_Y == s)) {
            scale_x <- layout$panel_scales_x[[panels$SCALE_X[[i]]]]
            layout$panel_params[[i]] <- fitted_coord$setup_panel_params(
                scale_x, scales_y[[s]], params = layout$coord_params)
        }
        fitted <- TRUE
    }
    if (!fitted) return(invisible())
    if (!is.null(layout$facet$setup_panel_params)) {
        layout$panel_params <-
            layout$facet$setup_panel_params(layout$panel_params, coord)
    }
    # For plotHover(), whose conversion does not use this Layout.
    layout$fitted_panel_params <- layout$panel_params
    invisible()
}

# The y limits in data space with their open ends replaced by the range of
# `values`, which are on the scale that `transform` maps to and `inverse` maps
# back from. NULL when there is nothing to fit or an open end would cross a
# requested one.
fit_open_limits <- function(values, ylim, transform, inverse) {
    if (length(values) == 0) return(NULL)
    requested <- suppressWarnings(transform(ylim))
    fitted <- range(values)
    fitted[!is.na(requested)] <- requested[!is.na(requested)]
    if (fitted[1] > fitted[2]) return(NULL)
    limits <- inverse(fitted)
    limits[!is.na(ylim)] <- ylim[!is.na(ylim)]
    limits
}

# The heights that the layer data `d` shows inside the x range `window`, both
# on the transformed scales: every height of a row reaching into the window,
# such as a vertex, a reference line or a rectangle, and for a layer drawn as
# a line or band the heights where its segments cross the edges of the window.
# Segments join consecutive rows of a group in the order in which they are
# drawn, and are straight on the transformed scales.
visible_heights <- function(d, window, x_aes, y_aes, path) {
    y_cols <- intersect(y_aes, names(d))
    if (nrow(d) == 0 || length(y_cols) == 0) return(numeric())
    x_cols <- intersect(x_aes, names(d))
    if (length(x_cols) == 0) {
        inside <- rep(TRUE, nrow(d))
    } else {
        # An infinite position is drawn at the edge of the panel.
        xs <- lapply(d[x_cols], function(v) {
            v[v == -Inf] <- window[[1]]
            v[v == Inf] <- window[[2]]
            v
        })
        lo <- do.call(pmin, c(unname(xs), na.rm = TRUE))
        hi <- do.call(pmax, c(unname(xs), na.rm = TRUE))
        inside <- !is.na(lo) & lo <= window[[2]] & hi >= window[[1]]
    }
    values <- unlist(d[inside, y_cols], use.names = FALSE)
    if (path && "x" %in% x_cols && nrow(d) > 1) {
        n <- nrow(d)
        x1 <- d$x[-n]
        x2 <- d$x[-1]
        joined <- if ("group" %in% names(d)) {
            d$group[-n] == d$group[-1]
        } else {
            rep(TRUE, n - 1)
        }
        for (edge in window[is.finite(window)]) {
            cross <- joined & x1 != x2 & (x1 - edge) * (x2 - edge) <= 0
            cross[is.na(cross)] <- FALSE
            if (!any(cross)) next
            t <- (edge - x1[cross]) / (x2[cross] - x1[cross])
            for (col in y_cols) {
                y1 <- d[[col]][-n][cross]
                y2 <- d[[col]][-1][cross]
                values <- c(values, y1 + t * (y2 - y1))
            }
        }
    }
    values[is.finite(values)]
}

# plotly converts a ggplot with its own copy of ggplot2's build, which does not
# use the plot's Layout and so would not fit the y limits to the visible data.
# Build the plot here instead and hand plotly the fitted panels through the
# facet, which ggplot2 (>= 4.0.0) lets adjust the panel parameters.
hover_viewport <- function(plot) {
    if (!inherits(plot$layout, "MizerLayout") ||
            is.null(plot$facet$setup_panel_params) ||
            all(is.na(plot$coordinates$limits$x %||% NA))) {
        return(plot)
    }
    panel_params <- ggplot_build(plot)$layout$fitted_panel_params
    if (is.null(panel_params)) return(plot)
    plot$facet <- facet_with_panel_params(plot$facet, panel_params)
    plot
}

# A separate function so that the method's environment holds only the panel
# parameters, not the plot it is attached to.
facet_with_panel_params <- function(facet, panel_params) {
    ggproto(NULL, facet,
            setup_panel_params = function(self, ...) panel_params)
}

# Added layers share the existing viewport. Only supplied endpoints change it;
# the open ends of the y limits are fitted when the plot is built. The plot is
# returned with a new coordinate system; the caller's is not modified.
set_plot_viewport <- function(plot, xlim = c(NA, NA), ylim = c(NA, NA)) {
    assert_that(length(xlim) == 2, length(ylim) == 2)
    coord <- plot$coordinates
    merge_limits <- function(existing, requested, axis) {
        existing <- existing %||% c(NA, NA)
        existing[!is.na(requested)] <- requested[!is.na(requested)]
        if (plot_axis_is_log(plot, axis)) {
            existing <- drop_nonpositive_limits(existing, "log10")
        }
        existing
    }
    limits <- list(x = merge_limits(coord$limits$x, xlim, "x"),
                   y = merge_limits(coord$limits$y, ylim, "y"))
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

    p <- ggplot(frame, plot_spec_mapping(spec)) + plot_spec_scales(spec) +
        plot_spec_layers(spec, legend_levels)
    wrap_var <- spec$options$wrap_var
    if (!is.null(wrap_var)) {
        p <- p + facet_wrap(wrap_var, scales = spec$options$wrap_scale)
    }
    make_mizer_plot(p, plot_spec_tooltip(spec))
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
    p <- plotly::layout(
        p,
        xaxis = plotly_axis(df[[vars$x]], spec$limits$x,
                            spec$scales$x == "log10", spec$labels$x),
        yaxis = plotly_axis(df[[vars$value]], animation_ylim(spec),
                            spec$scales$y == "log10", spec$labels$y),
        legend = list(traceorder = "normal")
    )
    p <- plotly::animation_opts(p, frame = spec$options$frame_duration,
                                transition = spec$options$transition_duration,
                                easing = spec$options$easing)
    stable_plotly_ids(p)
}

# An animation draws only its own lines, without stacks, facets or added
# layers, so its open y limits are fitted to its data directly: to the heights
# that any frame shows between the x limits, which plotly draws without
# expansion. add_lines() joins each trace in the order of x.
animation_ylim <- function(spec) {
    vars <- spec$variables
    xlim <- spec$limits$x
    ylim <- spec$limits$y
    if (all(is.na(xlim)) || !anyNA(ylim)) return(ylim)
    to_axis <- function(v, trans) {
        if (identical(trans, "log10")) suppressWarnings(log10(v)) else v
    }
    df <- spec$data
    d <- data.frame(
        x = to_axis(df[[vars$x]], spec$scales$x),
        y = to_axis(df[[vars$value]], spec$scales$y),
        group = interaction(df[unique(c(vars$group, vars$frame))], drop = TRUE)
    )
    d <- d[order(d$group, d$x), , drop = FALSE]
    window <- to_axis(xlim, spec$scales$x)
    window[is.na(window)] <- c(-Inf, Inf)[is.na(window)]
    values <- visible_heights(d, window, "x", "y", path = TRUE)
    log_y <- identical(spec$scales$y, "log10")
    fit_open_limits(values, ylim,
                    transform = function(v) to_axis(v, spec$scales$y),
                    inverse = if (log_y) function(v) 10^v else identity) %||%
        ylim
}
