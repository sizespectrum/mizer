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
