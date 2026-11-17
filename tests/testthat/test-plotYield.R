# Plot fixtures use edition 1, as do the existing visual baselines.
withr::local_options(mizer_defaults_edition = 1)
params <- plot_params_small
sim <- plot_sim_small
sim0 <- plot_sim0_small
species <- c(1, 2)

test_that("yield plotting helpers validate comparison and gear selection", {
    sim_shifted <- sim
    dimnames(sim_shifted@n)$time <- as.character(10:13)
    lifecycle::expect_deprecated(
        expect_error(plotYield(sim, sim_shifted),
                     "do not have the same times")
    )

    y <- plotYieldGear(sim, species = species,
                       gears = dimnames(sim@params@selectivity)$gear[[1]],
                       return_data = TRUE)
    expect_true(all(y$Gear == dimnames(sim@params@selectivity)$gear[[1]]))
})

test_that("yield plotly wrappers return plotly objects", {
    expect_s3_class(plotlyYield(sim, species = species, ylim = c(1e5, 1e12)), "plotly")
    expect_s3_class(plotlyYieldGear(sim, species = species, ylim = c(1e5, 1e12)), "plotly")
})

test_that("yield plotting helpers accept ylim", {
    p <- plotYield(sim, species = species, ylim = c(1e-5, 1))
    expect_equal(p$coordinates$limits$y, c(1e-5, 1))

    p_gear <- plotYieldGear(sim, species = species, ylim = c(1e-5, 1))
    expect_equal(p_gear$coordinates$limits$y, c(1e-5, 1))
})

test_that("plotYieldGear supports log_x, log_y, and log arguments", {
    p_default <- plotYieldGear(sim, species = species)
    expect_identical(p_default$scales$get_scales("x")$trans$name, "identity")
    expect_identical(p_default$scales$get_scales("y")$trans$name, "log-10")

    p_log_x <- plotYieldGear(sim, species = species, log_x = TRUE)
    expect_identical(p_log_x$scales$get_scales("x")$trans$name, "log-10")

    p_log_y_false <- plotYieldGear(sim, species = species, log_y = FALSE)
    expect_identical(p_log_y_false$scales$get_scales("y")$trans$name, "identity")

    p_log_xy <- plotYieldGear(sim, species = species, log = "xy")
    expect_identical(p_log_xy$scales$get_scales("x")$trans$name, "log-10")
    expect_identical(p_log_xy$scales$get_scales("y")$trans$name, "log-10")
})

test_that("plotYield delegates and preserves its public return-data schemas", {
    sim <- NS_sim_small

    # Single sim plot
    p <- plotYield(sim)
    expect_s3_class(p, "ggplot")
    expect_identical(p$scales$get_scales("x")$name, "Year")
    expect_identical(p$scales$get_scales("y")$name, "Yield [g/year]")

    y <- plotYield(sim, return_data = TRUE)
    expect_named(y, c("Year", "Yield", "Species"))

    lifecycle::expect_deprecated(
        y2 <- plotYield(sim, sim, return_data = TRUE)
    )
    expect_named(y2, c("Year", "Yield", "Species", "Simulation"))
})

test_that("yield plots return ggplots", {
    p <- plotYield(sim, species = species, total = TRUE)
    expect_ggplot(p)

    p <- plotYieldGear(sim, species = species)
    expect_ggplot(p)
})

test_that("interactive wrappers preserve tooltip fields", {
    gp <- plotlyYield(sim, species = species)
    expect_s3_class(gp, "plotly")
    expect_equal(tooltip_fields(gp), c("Species", "Year", "Yield"))

    lifecycle::expect_deprecated(
        expect_s3_class(plotlyYield(sim, sim), "plotly")
    )

    gp <- plotlyYieldGear(sim, species = species)
    expect_s3_class(gp, "plotly")
    expect_equal(tooltip_fields(gp), c("Year", "Yield", "Species", "Gear"))
})

test_that("public return-data shape is preserved", {
    expect_warning(p <- plotYield(sim, sim0, species = species, return_data = TRUE))
    expect_equal(dim(p), c(8,4))

    expect_equal(dim(plotYieldGear(sim, species = species, return_data = TRUE)), c(8,4))
})

test_that("time limits only set yield viewports", {
    all_rows <- nrow(plotYield(sim, species = species, return_data = TRUE))
    limited_rows <- nrow(plotYield(sim, species = species,
                                   tlim = c(1, 3), return_data = TRUE))
    expect_equal(limited_rows, all_rows)

    # plotYieldGear: tlim retains all rows
    all_rows <- nrow(plotYieldGear(sim, species = species, return_data = TRUE))
    limited_rows <- nrow(plotYieldGear(sim, species = species,
                                       tlim = c(1, 3), return_data = TRUE))
    expect_equal(limited_rows, all_rows)

    # The deprecation of animate(time_range) is tested in
    # test-animateSpectra.R, where animate() is defined. It cannot be tested
})

test_that("time viewports preserve returned data and totals", {
    for (fun in list(plotYield, plotYieldGear)) {
        full <- fun(NS_sim_small, total = TRUE, return_data = TRUE)
        zoom <- fun(NS_sim_small, total = TRUE, tlim = c(0.5, 2.5),
                    ylim = c(1, 10), return_data = TRUE)
        expect_equal(zoom, full)
    }
})
