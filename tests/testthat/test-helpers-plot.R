test_that("log_breaks returns monotonic breaks", {
    f <- log_breaks(5)
    b <- f(c(1e-3, 1e3))
    expect_true(is.numeric(b))
    expect_true(all(diff(b) > 0))
})

test_that("plotDataFrame builds ggplot with species legend", {
    params <- NS_params
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


