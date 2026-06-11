# second_order_w ----
test_that("second_order_w slot defaults to c(flux_limiter=FALSE, bin_average=FALSE)", {
    params <- NS_params_small
    sow <- second_order_w(params)
    expect_identical(sow, c(flux_limiter = FALSE, bin_average = FALSE))
})

test_that("second_order_w getter works", {
    params <- NS_params_small
    params@second_order_w <- c(flux_limiter = TRUE, bin_average = TRUE)
    sow <- second_order_w(params)
    expect_identical(sow, c(flux_limiter = TRUE, bin_average = TRUE))
})

test_that("second_order_w setter with single TRUE sets both entries", {
    params <- NS_params_small
    second_order_w(params) <- TRUE
    expect_identical(second_order_w(params),
                     c(flux_limiter = TRUE, bin_average = TRUE))
})

test_that("second_order_w setter with single FALSE sets both entries", {
    params <- NS_params_small
    second_order_w(params) <- TRUE
    second_order_w(params) <- FALSE
    expect_identical(second_order_w(params),
                     c(flux_limiter = FALSE, bin_average = FALSE))
})

test_that("second_order_w setter with named vector sets individual entries", {
    params <- NS_params_small
    second_order_w(params) <- c(flux_limiter = TRUE)
    expect_true(second_order_w(params)[["flux_limiter"]])
    expect_false(second_order_w(params)[["bin_average"]])

    second_order_w(params) <- c(bin_average = TRUE)
    expect_true(second_order_w(params)[["flux_limiter"]])
    expect_true(second_order_w(params)[["bin_average"]])
})

test_that("second_order_w setter validates input", {
    params <- NS_params_small
    expect_error(second_order_w(params) <- "yes")
    expect_error(second_order_w(params) <- NA)
    expect_error(second_order_w(params) <- c(unknown = TRUE))
})

test_that("second_order_w setter re-runs setParams", {
    params <- NS_params_small
    second_order_w(params) <- TRUE
    expect_identical(second_order_w(params),
                     c(flux_limiter = TRUE, bin_average = TRUE))
    expect_s4_class(params, "MizerParams")
})

test_that("second_order_w slot is preserved by upgradeParams", {
    params <- NS_params_small
    params2 <- upgradeParams(params)
    expect_identical(second_order_w(params2),
                     c(flux_limiter = FALSE, bin_average = FALSE))
})

test_that("validObject accepts second_order_w = c(TRUE, TRUE)", {
    params <- NS_params_small
    params@second_order_w <- c(flux_limiter = TRUE, bin_average = TRUE)
    expect_silent(validObject(params))
})

test_that("validObject rejects invalid second_order_w", {
    params <- NS_params_small
    params@second_order_w <- logical(0)
    expect_error(validObject(params))

    params@second_order_w <- TRUE  # unnamed single value
    expect_error(validObject(params))
})
