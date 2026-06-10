# high_order ----
test_that("high_order slot defaults to FALSE", {
    params <- NS_params_small
    expect_false(high_order(params))
})

test_that("high_order getter works", {
    params <- NS_params_small
    params@high_order <- TRUE
    expect_true(high_order(params))
})

test_that("high_order setter validates input", {
    params <- NS_params_small
    expect_error(high_order(params) <- "yes")
    expect_error(high_order(params) <- NA)
})

test_that("high_order setter re-runs setParams", {
    params <- NS_params_small
    # Changing high_order should return a valid params object
    high_order(params) <- TRUE
    expect_true(high_order(params))
    expect_s4_class(params, "MizerParams")
})

test_that("high_order slot is preserved by upgradeParams", {
    params <- NS_params_small
    # upgradeParams should not alter the slot
    params2 <- upgradeParams(params)
    expect_false(high_order(params2))
})

test_that("validObject accepts high_order = TRUE", {
    params <- NS_params_small
    params@high_order <- TRUE
    expect_silent(validObject(params))
})

test_that("validObject rejects invalid high_order", {
    params <- NS_params_small
    params@high_order <- logical(0)
    expect_error(validObject(params))
})
