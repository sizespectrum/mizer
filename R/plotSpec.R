# Internal plotting contract. Data arrive fully prepared: selection, totals,
# density conversion and integration belong to the producer, not the renderer.
# This is a plain list, never a new public class or a return_data schema.
# density_wrt describes the prepared values, when known; it does not request
# another conversion. Legacy data-frame adapters need not infer it from labels.
new_plot_spec <- function(data, params, x_var, y_var, group_var,
                          legend_var = group_var, frame_var = NULL,
                          xlab = waiver(), ylab = waiver(),
                          xtrans = "identity", ytrans = "identity",
                          xlim = c(NA, NA), ylim = c(NA, NA),
                          style = "line", y_ticks = 6, highlight = NULL,
                          wrap_var = NULL, wrap_scale = NULL,
                          density_wrt = NA_character_, tooltip = NULL,
                          options = list()) {
    variables <- list(x = x_var, value = y_var, group = group_var,
                      legend = legend_var, frame = frame_var)
    missing <- setdiff(unlist(variables, use.names = FALSE), names(data))
    if (length(missing) > 0) {
        stop("Plot variables missing from data: ",
             paste(missing, collapse = ", "))
    }
    options <- utils::modifyList(
        list(mode = "ordinary", style = style, y_ticks = y_ticks,
             highlight = highlight, wrap_var = wrap_var,
             wrap_scale = wrap_scale, legend_warning_var = legend_var),
        options
    )
    list(data = data, params = params, variables = variables,
         labels = list(x = xlab, y = ylab),
         scales = list(x = xtrans, y = ytrans),
         limits = list(x = xlim, y = ylim),
         density_wrt = density_wrt, tooltip = tooltip, options = options)
}

# The positional contract is confined to the compatibility boundary. In
# particular, keep the public column names and the automatic Legend column.
legacy_plot_spec <- function(frame, params, legend_var = NULL, ...) {
    assert_that(is.data.frame(frame), is(params, "MizerParams"))
    if (ncol(frame) < 3) {
        stop("The data frame needs to have at least 3 variables.")
    }
    variables <- names(frame)[1:3]
    if (is.null(legend_var)) {
        frame$Legend <- frame[[variables[[3]]]]
        legend_var <- "Legend"
    } else if (!(legend_var %in% names(frame))) {
        stop("The `legend_var` argument must be the name of a variable ",
             "in the data frame.")
    }
    new_plot_spec(frame, params, x_var = variables[[1]],
                  y_var = variables[[2]], group_var = variables[[3]],
                  legend_var = legend_var, ...)
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
    first$density_wrt <- NA_character_
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

# Keep scale limits (and the relative renderer's historical treatment of
# non-positive coordinates) here. Viewport-only limits are a separate change.
plot_spec_scales <- function(spec) {
    xbreaks <- if (spec$scales$x == "log10") log_breaks() else waiver()
    ybreaks <- if (spec$scales$y == "log10") {
        log_breaks(n = spec$options$y_ticks)
    } else {
        waiver()
    }
    ylabels <- if (spec$options$mode == "relative") waiver() else prettyNum
    list(scale_y_continuous(trans = spec$scales$y, breaks = ybreaks,
                            labels = ylabels, name = spec$labels$y,
                            limits = spec$limits$y),
         scale_x_continuous(trans = spec$scales$x, breaks = xbreaks,
                            name = spec$labels$x, limits = spec$limits$x))
}

render_plot_spec <- function(spec) {
    vars <- spec$variables
    x_var <- vars$x
    y_var <- vars$value
    group_var <- vars$group
    legend_var <- vars$legend
    mode <- spec$options$mode
    style <- spec$options$style
    frame <- spec$data
    legend_levels <- plot_spec_legend_levels(spec)
    frame[[legend_var]] <- factor(frame[[legend_var]], levels = legend_levels)
    if (sum(is.na(frame[[spec$options$legend_warning_var]]))) {
        warning("missing legend in params@linecolour, ",
                "some groups won't be displayed")
    }
    linecolour <- spec$params@linecolour[legend_levels]
    linetype <- spec$params@linetype[legend_levels]
    linesize <- make_linesize(legend_levels, spec$options$highlight)

    if (mode != "relative") {
        spec$limits$x <- drop_nonpositive_limits(spec$limits$x, spec$scales$x)
        spec$limits$y <- drop_nonpositive_limits(spec$limits$y, spec$scales$y)
        frame <- drop_nonpositive_values(frame, x_var, spec$scales$x)
        frame <- drop_nonpositive_values(frame, y_var, spec$scales$y)
    }
    if (mode != "ordinary") frame$LineSpec <- frame[[legend_var]]
    spec$data <- frame

    # Preserve mapping order: ggplotly uses it to order the hover fields.
    mapping <- aes(group = .data[[group_var]])
    if (mode == "comparison") {
        model_var <- vars$model
        mapping <- aes(group = interaction(.data[[legend_var]],
                                           .data[[model_var]]),
                       colour = .data[[legend_var]])
    }
    p <- ggplot(frame, mapping) + plot_spec_scales(spec)

    if (mode == "comparison") {
        p <- p + geom_line(aes(x = .data[[x_var]], y = .data[[y_var]],
                               linetype = .data[[model_var]],
                               linewidth = .data[["LineSpec"]]))
    } else if (mode == "relative") {
        p <- p + geom_hline(yintercept = 0, linetype = 1,
                            colour = "dark grey", linewidth = 0.75) +
            geom_line(aes(x = .data[[x_var]], y = .data[[y_var]],
                          colour = .data[[legend_var]],
                          linewidth = .data[["LineSpec"]]))
    } else if (style == "area") {
        p <- p + geom_area(aes(x = .data[[x_var]], y = .data[[y_var]],
                               fill = .data[[legend_var]])) +
            scale_fill_manual(values = linecolour)
    } else {
        if (!style %in% c("line", "ribbon", "envelope")) {
            stop("unknown style selected")
        }
        if (style != "line") {
            check_band_vars(names(frame), style)
            p <- p + geom_ribbon(aes(x = .data[[x_var]],
                                     ymin = .data[["ymin"]],
                                     ymax = .data[["ymax"]],
                                     fill = .data[[legend_var]]),
                                 alpha = 0.25, colour = NA) +
                scale_fill_manual(values = linecolour)
        }
        # Give each layer its own binding of yv, including both band edges.
        value_line <- function(yv) {
            geom_line(aes(x = .data[[x_var]], y = .data[[yv]],
                          colour = .data[[legend_var]],
                          linetype = .data[[legend_var]],
                          linewidth = .data[[legend_var]]))
        }
        y_vars <- if (style == "envelope") c("ymax", "ymin") else y_var
        for (yv in y_vars) p <- p + value_line(yv)
    }
    if (mode != "ordinary" || style != "area") {
        p <- p + scale_colour_manual(values = linecolour)
        if (mode == "ordinary") {
            p <- p + scale_linetype_manual(values = linetype) +
                scale_discrete_manual("linewidth", values = linesize)
        } else {
            if (mode == "comparison") {
                p <- p + scale_linetype_discrete(drop = FALSE)
            }
            p <- p + scale_discrete_manual("linewidth", values = linesize,
                                           guide = "none")
        }
    }
    wrap_var <- spec$options$wrap_var
    if (!is.null(wrap_var)) {
        if (!(wrap_var %in% names(frame))) {
            stop("The `wrap_var` argument must be the name of a variable ",
                 "in the data frame.")
        }
        p <- p + facet_wrap(wrap_var, scales = spec$options$wrap_scale)
    }
    make_mizer_plot(p, plot_spec_tooltip(spec))
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
        yaxis = plotly_axis(df[[vars$value]], spec$limits$y,
                            spec$scales$y == "log10", spec$labels$y),
        legend = list(traceorder = "normal")
    )
    p <- plotly::animation_opts(p, frame = spec$options$frame_duration,
                                transition = spec$options$transition_duration,
                                easing = spec$options$easing)
    stable_plotly_ids(p)
}
