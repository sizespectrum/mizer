# Plot fixtures use edition 1, as do the existing visual baselines.
withr::local_options(mizer_defaults_edition = 1)
params <- plot_params_small
sim <- plot_sim_small
species <- c(1, 2)

test_that("plotHover on array objects has correct tooltip fields", {
    expect_equal(tooltip_fields(plotHover(getEncounter(params),
                                          species = species)),
                 c("Species", "w", "Encounter rate"))
    expect_equal(tooltip_fields(plotHover(getBiomass(sim), species = species)),
                 c("Species", "Year", "Biomass"))
    expect_equal(tooltip_fields(plotHover(getEncounter(sim),
                                          species = species[[1]])),
                 c("Species", "w", "Encounter rate"))
})

test_that("plotHover(plot(...)) uses concise mizer tooltips", {
    p <- plot(getEncounter(NS_params_small), species = "Cod")
    expect_s3_class(p, "mizer_plot")
    gp <- plotHover(p)
    expect_equal(tooltip_fields(gp), c("Species", "w", "Encounter rate"))
})

test_that("plotHover uses reproducible internal plotly ids", {
    p1 <- plotHover(getBiomass(NS_sim_small))
    p2 <- plotHover(getBiomass(NS_sim_small))

    expect_identical(names(p1$x$attrs), "mizer-plotly-1")
    expect_identical(names(p1$x$attrs), names(p2$x$attrs))
    expect_identical(names(p1$x$visdat), names(p2$x$visdat))
    expect_identical(p1$x$cur_data, p2$x$cur_data)
    expect_error(plotly::plotly_build(p1), NA)
})

test_that("comparison helpers preserve non-size x variables", {
    p <- plot2(getBiomass(NS_sim_small), getBiomass(NS_sim_small), species = "Cod")
    expect_s3_class(p, "ggplot")
    expect_true("Year" %in% names(p$data))
    expect_error(ggplot2::ggplot_build(p), NA)

    p_rel <- plotRelative(getBiomass(NS_sim_small), getBiomass(NS_sim_small),
                          species = "Cod")
    expect_s3_class(p_rel, "ggplot")
    expect_true("Year" %in% names(p_rel$data))
    expect_error(ggplot2::ggplot_build(p_rel), NA)
})

test_that("log_breaks returns monotonic breaks", {
    f <- log_breaks(5)
    b <- f(c(1e-3, 1e3))
    expect_true(is.numeric(b))
    expect_true(all(diff(b) > 0))
})

test_that("plotDataFrame builds ggplot with species legend", {
    params <- NS_params_small
    sp <- species_params(params)$species[1:2]
    df <- data.frame(
        x = c(1, 2, 3, 1, 2, 3),
        y = c(1, 2, 1, 2, 3, 2),
        Species = rep(sp, each = 3)
    )
    p <- plotDataFrame(
        df, params,
        style = "line",
        xlab = "x", ylab = "y",
        legend_var = "Species"
    )
    expect_true(is_ggplot(p))
})

test_that("plotDataFrame validates helper arguments", {
    params <- NS_params_small
    df <- data.frame(x = 1:3, y = 1:3, Species = species_params(params)$species[1])
    expect_error(plotDataFrame(df[, 1:2], params),
                 "at least 3 variables")
    expect_error(plotDataFrame(df, params, legend_var = "missing"),
                 "legend_var")
    expect_error(plotDataFrame(df, params, wrap_var = "missing"),
                 "wrap_var")
    expect_error(plotDataFrame(df, params, wrap_var = "Legend"),
                 "wrap_var")
    expect_error(plotDataFrame(df, params, style = "bogus"),
                 "unknown style selected")
})

test_that("plotDataFrame supports area plots, wrapping and log x breaks", {
    params <- NS_params_small
    sp <- species_params(params)$species[1:2]
    df <- data.frame(
        x = c(1, 10, 100, 1, 10, 100),
        y = c(1, 2, 3, 2, 3, 4),
        Species = rep(sp, each = 3),
        Panel = rep(c("A", "B"), each = 3)
    )
    p <- plotDataFrame(df, params,
                       style = "area",
                       xtrans = "log10",
                       legend_var = "Species",
                       wrap_var = "Panel",
                       wrap_scale = "free_y")
    expect_true(is_ggplot(p))
    expect_true(p$facet$params$free$y)
    expect_false(p$facet$params$free$x)
    expect_false(inherits(p$scales$scales[[2]]$breaks, "waiver"))
})

test_that("a log axis drops the values and limits it cannot show", {
    params <- NS_params_small
    sp <- species_params(params)$species[1:2]
    df <- data.frame(
        x = c(1, 10, 100, 1, 10, 100),
        y = c(1, 2, 0, 2, 3, 4),
        Species = rep(sp, each = 3)
    )
    # A y-limit of zero is impossible on a log axis and is auto-scaled instead,
    # and the zero value is dropped rather than warned about.
    p <- expect_no_warning(
        ggplot2::ggplot_build(
            plotDataFrame(df, params, ytrans = "log10", ylim = c(0, NA),
                          legend_var = "Species")))
    expect_equal(nrow(p$plot$data), 5)
    expect_true(all(is.na(p$plot$coordinates$limits$y)))
    # On a linear axis nothing is dropped and the limit is kept.
    p2 <- plotDataFrame(df, params, ylim = c(0, NA), legend_var = "Species")
    expect_equal(nrow(p2$data), 6)
    expect_equal(p2$coordinates$limits$y, c(0, NA))
})

test_that("a log axis that would hide everything is left alone", {
    params <- NS_params_small
    df <- data.frame(x = 1:3, y = c(0, 0, 0),
                     Species = species_params(params)$species[1])
    # All values are unplottable, so dropping them would only replace the
    # warning with the complaints of an empty plot.
    p <- plotDataFrame(df, params, ytrans = "log10", legend_var = "Species")
    expect_equal(nrow(p$data), 3)
})

test_that("plotDataFrame() needs ymin and ymax for the band styles", {
    frame <- data.frame(x = 1:2, y = 1:2, Species = "Cod")
    expect_error(plotDataFrame(frame, NS_params_small, style = "ribbon"),
                 "needs the variables")
})

test_that("interpolation reproduces the values exactly on a shared grid", {
    grid <- c(1, 2, 4, 8)
    frame1 <- data.frame(w = grid, value = c(1, 2, 3, 4),
                         Species = "Cod", Legend = "Cod")
    frame2 <- data.frame(w = grid, value = c(2, 4, 6, 8),
                         Species = "Cod", Legend = "Cod")
    out <- interpolate_relative_frames(frame1, frame2, "w", "value",
                                       c("Species", "Legend"))
    expect_equal(out$w, grid)
    expect_equal(out$value.x, frame1$value)
    expect_equal(out$value.y, frame2$value)
})

test_that("interpolation keeps only the domain both series cover", {
    frame1 <- data.frame(w = c(1, 2, 4), value = c(1, 2, 4),
                         Species = "Cod", Legend = "Cod")
    frame2 <- data.frame(w = c(2, 4, 8), value = c(2, 4, 8),
                         Species = "Cod", Legend = "Cod")
    out <- interpolate_relative_frames(frame1, frame2, "w", "value",
                                       c("Species", "Legend"))
    # Outside [2, 4] one of the two would have to be extrapolated.
    expect_equal(out$w, c(2, 4))
})

test_that("interpolation puts two offset grids on their union", {
    frame1 <- data.frame(w = c(1, 4, 16), value = c(0, 2, 4),
                         Species = "Cod", Legend = "Cod")
    frame2 <- data.frame(w = c(2, 8), value = c(1, 3),
                         Species = "Cod", Legend = "Cod")
    out <- interpolate_relative_frames(frame1, frame2, "w", "value",
                                       c("Species", "Legend"))
    expect_equal(out$w, c(2, 4, 8))
    # Both series are linear in log2 of the size, so the interpolation is exact.
    expect_equal(out$value.x, c(1, 2, 3))
    expect_equal(out$value.y, c(1, 2, 3))
})

test_that("interpolation drops a series that only one frame has", {
    frame1 <- data.frame(w = c(1, 2), value = c(1, 2),
                         Species = c("Cod", "Cod"), Legend = "Cod")
    frame2 <- data.frame(w = c(1, 2), value = c(1, 2),
                         Species = c("Sprat", "Sprat"), Legend = "Sprat")
    out <- interpolate_relative_frames(frame1, frame2, "w", "value",
                                       c("Species", "Legend"))
    expect_equal(nrow(out), 0)
    expect_setequal(names(out), c("Species", "Legend", "w", "value.x",
                                  "value.y"))
})

test_that("interpolate_in_log_size handles the awkward series", {
    # A single point speaks only for the coordinate it sits at.
    expect_equal(interpolate_in_log_size(2, 5, c(1, 2, 4)), c(NA, 5, NA))
    # A non-positive coordinate has no logarithm, so the fallback is linear.
    expect_equal(interpolate_in_log_size(c(0, 2), c(0, 2), c(0, 1, 2)),
                 c(0, 1, 2))
    # Unsorted input is sorted first.
    expect_equal(interpolate_in_log_size(c(4, 1), c(4, 1), c(1, 2, 4)),
                 c(1, 2.5, 4))
})

test_that("plotComparisonDataFrame draws the frames it is given", {
    frame1 <- data.frame(l = c(1, 2), value = c(1, 2),
                         Species = "Cod", Legend = "Cod")
    frame2 <- data.frame(l = c(3, 4), value = c(3, 4),
                         Species = "Cod", Legend = "Cod")
    p <- plotComparisonDataFrame(frame1, frame2, NS_params_small,
                                 name1 = "A", name2 = "B")
    # No conversion of its own: the coordinates are exactly what it was handed.
    expect_equal(p$data$l, c(1, 2, 3, 4))
    expect_equal(as.character(p$data$Model), c("A", "A", "B", "B"))
})

test_that("plotRelativeDataFrame joins unless told to interpolate", {
    frame1 <- data.frame(w = c(1, 2, 4), value = c(1, 1, 1),
                         Species = "Cod", Legend = "Cod")
    frame2 <- data.frame(w = c(2, 4, 8), value = c(3, 3, 3),
                         Species = "Cod", Legend = "Cod")
    joined <- plotRelativeDataFrame(frame1, frame2, NS_params_small)
    expect_equal(joined$data$w, c(2, 4))
    interpolated <- plotRelativeDataFrame(frame1, frame2, NS_params_small,
                                          interpolate = TRUE)
    expect_equal(interpolated$data$w, c(2, 4))
    expect_equal(unique(interpolated$data$rel_diff), 1)
})

test_that("explicit tooltips and user ggplot edits survive hover conversion", {
    frame <- data.frame(mass = c(1, 10, 100), response = c(1, 2, 4),
                        fish = "Cod")
    p <- plotDataFrame(frame, NS_params_small,
                       xlab = "Mass", ylab = "Response") +
        ggplot2::labs(title = "Custom title") +
        ggplot2::geom_hline(yintercept = 3, colour = "red")
    hover <- plotly::plotly_build(plotHover(p, tooltip = "response"))
    expect_identical(hover$x$layout$title$text, "Custom title")
    expect_true(any(vapply(hover$x$data, function(trace) {
        length(trace$y) > 0 && all(trace$y == 3)
    }, logical(1))))
    text <- unlist(lapply(hover$x$data, `[[`, "text"))
    expect_true(any(grepl("response", text)))
    expect_false(any(grepl("fish:|mass:", text)))
})

test_that("comparison renderers hide the LineSpec legend", {
    p_rel <- plotRelative(getBiomass(NS_sim_small), getBiomass(NS_sim_small), species = "Cod")
    expect_equal(p_rel$scales$get_scales("linewidth")$guide, "none")

    p2 <- plot2(getBiomass(NS_sim_small), getBiomass(NS_sim_small), species = "Cod")
    expect_equal(p2$scales$get_scales("linewidth")$guide, "none")
})
