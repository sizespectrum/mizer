# second_order_w ----
test_that("second_order_w slot defaults to FALSE", {
    params <- NS_params_small
    expect_false(second_order_w(params))
})

test_that("second_order_w getter works", {
    params <- NS_params_small
    params@second_order_w <- TRUE
    expect_true(second_order_w(params))
})

test_that("second_order_w setter validates input", {
    params <- NS_params_small
    expect_error(second_order_w(params) <- "yes")
    expect_error(second_order_w(params) <- NA)
})

test_that("second_order_w setter re-runs setParams", {
    params <- NS_params_small
    # Changing second_order_w should return a valid params object
    second_order_w(params) <- TRUE
    expect_true(second_order_w(params))
    expect_s4_class(params, "MizerParams")
})

test_that("second_order_w slot is preserved by upgradeParams", {
    params <- NS_params_small
    # upgradeParams should not alter the slot
    params2 <- upgradeParams(params)
    expect_false(second_order_w(params2))
})

test_that("validObject accepts second_order_w = TRUE", {
    params <- NS_params_small
    params@second_order_w <- TRUE
    expect_silent(validObject(params))
})

test_that("validObject rejects invalid second_order_w", {
    params <- NS_params_small
    params@second_order_w <- logical(0)
    expect_error(validObject(params))
})
