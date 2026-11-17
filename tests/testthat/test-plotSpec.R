test_that("plot specifications use named roles regardless of column order", {
    frame <- data.frame(unused = 99, fish = rep(c("Cod", "Sprat"), each = 3),
                        response = c(1, 2, 4, 2, 3, 5),
                        mass = rep(c(1, 10, 100), 2))
    spec <- new_plot_spec(frame, NS_params_small,
                          x_var = "mass", y_var = "response",
                          group_var = "fish",
                          xlab = "Mass", ylab = "Response", xtrans = "log10",
                          highlight = "Cod")
    expect_identical(class(spec), "list")
    before <- spec
    p <- render_plot_spec(spec)
    expect_identical(spec, before)
    expect_s3_class(p, "mizer_plot")
    expect_equal(p$data$response, frame$response)
    expect_equal(p$data$mass, frame$mass)
    expect_identical(attr(p, "mizer_tooltip"), c("fish", "mass", "response"))
    layer <- ggplot2::ggplot_build(p)$data[[1]]
    expect_equal(sort(layer$x), sort(log10(frame$mass)))
    expect_equal(sort(layer$y), sort(frame$response))
    cod <- layer$colour == NS_params_small@linecolour[["Cod"]]
    expect_equal(unique(layer$linewidth[cod]), 1.6)
    reordered <- spec
    reordered$data <- frame[, rev(names(frame))]
    expect_equal(ggplot2::ggplot_build(render_plot_spec(reordered))$data,
                 ggplot2::ggplot_build(p)$data)
    custom_tooltip <- spec
    custom_tooltip$tooltip <- c("response", "fish")
    expect_identical(attr(render_plot_spec(custom_tooltip), "mizer_tooltip"),
                     c("response", "fish"))
    skip_if_not_installed("vdiffr")
    vdiffr::expect_doppelganger("Named plot specification", p)
})

test_that("comparison and relative specifications align named roles", {
    first <- new_plot_spec(
        data.frame(fish = "Cod", response = c(1, 2, 4), mass = c(1, 10, 100)),
        NS_params_small, x_var = "mass", y_var = "response", group_var = "fish",
        xlab = "Mass", ylab = "Response", xtrans = "log10", highlight = "Cod"
    )
    second <- new_plot_spec(
        data.frame(size = c(1, 10, 100), animal = "Cod", rate = c(2, 4, 8)),
        NS_params_small, x_var = "size", y_var = "rate", group_var = "animal"
    )
    comparison <- compare_plot_specs(first, second, "Before", "After")
    p_comparison <- render_plot_spec(comparison)
    expect_equal(p_comparison$data$response, c(1, 2, 4, 2, 4, 8))
    expect_equal(p_comparison$data$mass, rep(c(1, 10, 100), 2))
    expect_identical(levels(p_comparison$data$Model), c("Before", "After"))
    expect_identical(attr(p_comparison, "mizer_tooltip"),
                     c("fish", "mass", "response", "Model"))

    relative <- relative_plot_specs(first, second)
    expect_equal(relative$data$rel_diff, rep(2 / 3, 3))
    expect_identical(relative$variables$value, "rel_diff")
    p_relative <- render_plot_spec(relative)
    expect_identical(attr(p_relative, "mizer_tooltip"),
                     c("fish", "mass", "rel_diff"))

    # Interpolation uses the named size and value fields too.
    second$data <- data.frame(size = c(10, 1000), animal = "Cod",
                              rate = c(2, 6))
    interpolated <- relative_plot_specs(first, second, interpolate = TRUE)
    expect_equal(interpolated$data$mass, c(10, 100))
    expect_equal(interpolated$data$rel_diff, c(0, 0))

    skip_if_not_installed("vdiffr")
    vdiffr::expect_doppelganger("Named comparison specification", p_comparison)
    vdiffr::expect_doppelganger("Named relative specification", p_relative)
})

test_that("comparison data frames are aligned by role, not by position", {
    frame1 <- data.frame(w = c(1, 10), value = c(1, 2), Species = "Cod",
                         Legend = "Cod")
    frame2 <- data.frame(length = c(1, 10), biomass = c(3, 4),
                         Species = "Cod", Legend = "Cod")
    p <- plotComparisonDataFrame(frame1, frame2, NS_params_small)
    expect_named(p$data, c("w", "value", "Species", "Legend", "Model",
                           "LineSpec"))
    expect_equal(p$data$value, c(1, 2, 3, 4))
    rel <- plotRelativeDataFrame(frame1, frame2, NS_params_small)
    expect_equal(rel$data$rel_diff, relative_difference(c(1, 2), c(3, 4)))
})

test_that("specifications are validated when they are built", {
    frame <- data.frame(mass = 1:2, response = 1:2, fish = "Cod")
    expect_error(new_plot_spec(frame, NS_params_small, "mass", "response",
                               "fish", style = "bogus"),
                 "unknown style selected")
    expect_error(new_plot_spec(frame, NS_params_small, "mass", "response",
                               "fish", wrap_var = "missing"),
                 "wrap_var")
    expect_error(new_plot_spec(frame, NS_params_small, "mass", "missing",
                               "fish"),
                 "Plot variables missing from data: missing")
    spec <- new_plot_spec(frame, NS_params_small, "mass", "response", "fish",
                          options = list(mode = "bogus"))
    expect_error(render_plot_spec(spec), "Unknown plot specification mode")
})

test_that("the legend warning checks the legend variable that is used", {
    frame <- data.frame(mass = 1:2, response = 1:2, fish = "Cod",
                        key = "Not a colour")
    expect_warning(plotDataFrame(frame, NS_params_small, legend_var = "key"),
                   "missing legend in params@linecolour")
    expect_no_warning(plotDataFrame(frame, NS_params_small))
})

test_that("the legacy adapter preserves data columns and log filtering", {
    frame <- data.frame(w = c(1, 10, 100), value = c(-1, 0, 2),
                        Species = "Cod")
    spec <- legacy_plot_spec(frame, NS_params_small, ytrans = "log10",
                             ylim = c(0, 4))
    p <- render_plot_spec(spec)
    expect_named(p$data, c("w", "value", "Species", "Legend"))
    expect_equal(p$data$value, 2)
    expect_equal(p$coordinates$limits$y, c(NA, 4))
    # Rendering does not mutate the producer's data or introduce a viewport
    # policy into preparation. Linear plots retain signed values as before.
    expect_equal(spec$data$value, frame$value)
    spec$scales$y <- "identity"
    p <- render_plot_spec(spec)
    expect_equal(p$data$value, frame$value)
    expect_equal(p$coordinates$limits$y, c(0, 4))
})

test_that("animation specifications use named fields and preserve frame data", {
    frame <- data.frame(step = rep(0:1, each = 6),
                        fish = rep(rep(c("Cod", "Sprat"), each = 3), 2),
                        response = rep(c(1, 2, 4, 2, 3, 5), 2),
                        mass = rep(c(1, 10, 100), 4),
                        key = "Background")
    spec <- new_plot_spec(
        frame, NS_params_small, x_var = "mass", y_var = "response",
        group_var = "fish", legend_var = "key", frame_var = "step",
        xlab = "Mass", ylab = "Response", xtrans = "log10",
        xlim = c(2, 90), ylim = c(0, 6),
        options = list(frame_duration = 800, transition_duration = 0,
                       easing = "linear")
    )
    animation <- render_plotly_spec(spec)
    built <- plotly::plotly_build(animation)
    expect_equal(built$x$layout$xaxis$range, log10(c(2, 90)))
    expect_equal(built$x$layout$yaxis$range, c(0, 6))
    expect_identical(animation$animation$frame$duration, 800)
    expect_identical(animation$animation$transition$duration, 0)
    for (frame_data in built$x$frames) {
        expect_length(frame_data$data, 2)
        expect_equal(vapply(frame_data$data, function(trace) length(trace$x),
                            integer(1)), c(3L, 3L))
        expect_equal(sum(vapply(frame_data$data, `[[`, logical(1),
                                "showlegend")), 1)
    }
    reordered <- spec
    reordered$data <- reordered$data[, rev(names(reordered$data))]
    expect_equal(plotly::plotly_build(render_plotly_spec(reordered))$x$frames,
                 built$x$frames)
})

test_that("viewport limits preserve vertices and agree with hover ranges", {
    data <- data.frame(x = c(1, 10, 100), y = c(1, 10, 100), fish = "Cod")
    spec <- new_plot_spec(data, NS_params_small, x_var = "x", y_var = "y",
                          group_var = "fish", xtrans = "log10",
                          ytrans = "log10", xlim = c(2, 50), ylim = c(3, 60))
    p <- render_plot_spec(spec)
    expect_null(p$scales$get_scales("x")$limits)
    expect_null(p$scales$get_scales("y")$limits)
    expect_equal(p$data$x, data$x)
    built <- ggplot2::ggplot_build(p)
    expect_equal(built$data[[1]]$x, log10(data$x))
    expect_equal(built$data[[1]]$y, log10(data$y))
    hover <- plotly::plotly_build(plotHover(p))
    panel <- built$layout$panel_params[[1]]
    expect_equal(hover$x$layout$xaxis$range, panel$x.range)
    expect_equal(hover$x$layout$yaxis$range, panel$y.range)
    expect_equal(as.numeric(hover$x$data[[1]]$x), log10(data$x))
    skip_if_not_installed("vdiffr")
    vdiffr::expect_doppelganger("Viewport keeps crossing lines", p)
})

# The y range of each panel of a built plot.
built_panels <- function(p) ggplot2::ggplot_build(p)$layout$panel_params
# ggplot2's default expansion of a continuous range, 5% at each end.
expand5 <- function(r) r + c(-1, 1) * 0.05 * diff(r)
# Two facet panels that are flat at the edges of the window x in [3, 5], at
# heights 1 and 3 and at 100 and 300, and rise steeply outside it.
two_panel_data <- function() {
    shape <- c(1000, 1, 1, 2, 3, 3, 1000)
    data.frame(x = rep(1:7, 2), y = c(shape, 100 * shape),
               fish = rep(c("Cod", "Sprat"), each = 7),
               panel = rep(c("a", "b"), each = 7))
}

test_that("an open end of the y viewport is fitted to the visible data", {
    # y = x on log-log axes, so the line crosses the edges of the visible x
    # range at the heights of those edges.
    data <- data.frame(x = c(1, 10, 100), y = c(1, 10, 100), fish = "Cod")
    spec <- new_plot_spec(data, NS_params_small, x_var = "x", y_var = "y",
                          group_var = "fish", xtrans = "log10",
                          ytrans = "log10", xlim = c(2, 50))
    p <- render_plot_spec(spec)
    # The plot keeps its limits and data; the fit is made when it is built.
    expect_equal(p$coordinates$limits$y, c(NA, NA))
    expect_equal(p$data$y, data$y)
    panel <- built_panels(p)[[1]]
    expect_equal(panel$y.range, expand5(panel$x.range))
    # A requested end is kept; only the open one is fitted.
    spec$limits$y <- c(1, NA)
    panel <- built_panels(render_plot_spec(spec))[[1]]
    expect_equal(panel$y.range, expand5(c(0, panel$x.range[2])))
    # Without an x limit everything is visible, so ggplot2 chooses.
    spec$limits <- list(x = c(NA, NA), y = c(NA, NA))
    built <- ggplot2::ggplot_build(render_plot_spec(spec))
    expect_null(built$layout$fitted_panel_params)
    # An open end that would cross a requested one is left to ggplot2.
    spec$limits <- list(x = c(2, 50), y = c(NA, 1))
    built <- ggplot2::ggplot_build(render_plot_spec(spec))
    expect_null(built$layout$fitted_panel_params)
})

test_that("fitted viewports keep reference lines, bands and stacks in view", {
    first <- new_plot_spec(
        data.frame(x = 1:4, y = c(1, 2, 3, 4), fish = "Cod"),
        NS_params_small, x_var = "x", y_var = "y", group_var = "fish",
        xlim = c(2, 3))
    second <- first
    second$data$y <- c(1.1, 2.2, 3.3, 4.4)
    # The zero line of a relative plot stays visible.
    rel <- render_plot_spec(relative_plot_specs(first, second))
    expect_equal(built_panels(rel)[[1]]$y.range,
                 expand5(c(0, max(rel$data$rel_diff))))
    # Stacked areas reach the top of the stack, here y = 2x.
    stacked <- new_plot_spec(
        data.frame(x = rep(1:4, 2), y = rep(1:4, 2),
                   fish = rep(c("Cod", "Sprat"), each = 4)),
        NS_params_small, x_var = "x", y_var = "y", group_var = "fish",
        xlim = c(2, 3), style = "area")
    panel <- built_panels(render_plot_spec(stacked))[[1]]
    expect_equal(panel$y.range, expand5(c(0, 2 * panel$x.range[2])))
    # Both edges of a band count, here at y = x - 1 and y = x + 1.
    band <- new_plot_spec(
        data.frame(x = 1:4, y = 1:4, ymin = 0:3, ymax = 2:5, fish = "Cod"),
        NS_params_small, x_var = "x", y_var = "y", group_var = "fish",
        xlim = c(2, 3), style = "ribbon")
    panel <- built_panels(render_plot_spec(band))[[1]]
    expect_equal(panel$y.range, expand5(panel$x.range + c(-1, 1)))
})

test_that("stacked areas on different size grids are fitted as drawn", {
    # ggplot2 aligns the two flat series before stacking them, so the stack
    # reaches 200 although no x value is shared.
    stacked <- new_plot_spec(
        data.frame(x = c(1:10, 1:10 + 0.5), y = 100,
                   fish = rep(c("Cod", "Sprat"), each = 10)),
        NS_params_small, x_var = "x", y_var = "y", group_var = "fish",
        xlim = c(3, 8), style = "area")
    expect_equal(built_panels(render_plot_spec(stacked))[[1]]$y.range,
                 expand5(c(0, 200)))
})

test_that("layers added with ggplot2 are fitted when the plot is built", {
    base <- render_plot_spec(new_plot_spec(
        data.frame(x = 1:4, y = c(1, 2, 3, 4), fish = "Cod"),
        NS_params_small, x_var = "x", y_var = "y", group_var = "fish",
        xlim = c(2, 3)))
    window <- built_panels(base)[[1]]$x.range
    # A reference line stays in view.
    panel <- built_panels(base + ggplot2::geom_hline(yintercept = 10))[[1]]
    expect_equal(panel$y.range, expand5(c(window[1], 10)))
    # Only the visible part of an added line counts, up to where it crosses
    # the edge of the window.
    added <- data.frame(x = 1:4, y = c(100, 5, 6, 100))
    line <- base + ggplot2::geom_line(ggplot2::aes(x, y), data = added,
                                      inherit.aes = FALSE)
    crossing <- 6 + (window[2] - 3) * (100 - 6)
    expect_equal(built_panels(line)[[1]]$y.range,
                 expand5(c(window[1], crossing)))
    # Points are not joined across the edge.
    points <- base + ggplot2::geom_point(ggplot2::aes(x, y), data = added,
                                         inherit.aes = FALSE)
    expect_equal(built_panels(points)[[1]]$y.range, expand5(c(window[1], 6)))
})

test_that("each free y scale of a facet is fitted on its own", {
    spec <- new_plot_spec(two_panel_data(), NS_params_small, x_var = "x", y_var = "y",
                          group_var = "fish", xlim = c(3, 5),
                          wrap_var = "panel", wrap_scale = "free_y")
    panels <- built_panels(render_plot_spec(spec))
    expect_equal(panels[[1]]$y.range, expand5(c(1, 3)))
    expect_equal(panels[[2]]$y.range, expand5(c(100, 300)))
    # A fixed scale stays shared by the panels.
    spec$options$wrap_scale <- "fixed"
    panels <- built_panels(render_plot_spec(spec))
    expect_equal(panels[[1]]$y.range, expand5(c(1, 300)))
    expect_equal(panels[[2]]$y.range, expand5(c(1, 300)))
})

test_that("plotHover() shows the ranges fitted when the plot is built", {
    skip_if(utils::packageVersion("ggplot2") < "4.0.0",
            "ggplot2 < 4.0.0 cannot hand panel parameters to plotly")
    p <- render_plot_spec(new_plot_spec(
        two_panel_data(), NS_params_small, x_var = "x", y_var = "y", group_var = "fish",
        xlim = c(3, 5), wrap_var = "panel", wrap_scale = "free_y")) +
        ggplot2::geom_hline(yintercept = 10)
    panels <- built_panels(p)
    hover <- plotly::plotly_build(plotHover(p))
    expect_equal(hover$x$layout$yaxis$range, panels[[1]]$y.range)
    expect_equal(hover$x$layout$yaxis2$range, panels[[2]]$y.range)
    # The plot's own facet still leaves the panel parameters to the build.
    expect_identical(p$facet$setup_panel_params(list("unchanged"), NULL),
                     list("unchanged"))
})

test_that("set_plot_viewport() only changes the supplied endpoints", {
    base <- render_plot_spec(new_plot_spec(
        data.frame(x = 1:4, y = c(1, 2, 3, 4), fish = "Cod"),
        NS_params_small, x_var = "x", y_var = "y", group_var = "fish",
        xlim = c(2, 3)))
    p <- set_plot_viewport(base, xlim = c(1, NA), ylim = c(NA, 5))
    expect_equal(p$coordinates$limits, list(x = c(1, 3), y = c(NA, 5)))
    # The input plot is not modified.
    expect_equal(base$coordinates$limits$x, c(2, 3))
    # Plots that mizer did not lay out are not fitted.
    plain <- ggplot2::ggplot(data.frame(x = 1:4, y = c(100, 5, 6, 100)),
                             ggplot2::aes(x, y)) +
        ggplot2::geom_line()
    plain <- set_plot_viewport(plain, xlim = c(2, 3))
    expect_equal(plain$coordinates$limits$y, c(NA, NA))
    expect_null(ggplot2::ggplot_build(plain)$layout$fitted_panel_params)
})

test_that("animations fit an open y range to the visible sizes", {
    frame <- data.frame(step = rep(0:1, each = 3), fish = "Cod",
                        response = c(1, 10, 1000, 2, 20, 2000),
                        mass = rep(c(1, 10, 100), 2))
    spec <- new_plot_spec(
        frame, NS_params_small, x_var = "mass", y_var = "response",
        group_var = "fish", frame_var = "step", xtrans = "log10",
        ytrans = "log10", xlim = c(1, 10),
        options = list(frame_duration = 1, transition_duration = 0,
                       easing = "linear"))
    built <- plotly::plotly_build(render_plotly_spec(spec))
    expect_equal(built$x$layout$yaxis$range, log10(c(1, 20)))
})
