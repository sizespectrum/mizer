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
