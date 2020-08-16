params <- NS_params
no_sp <- nrow(params@species_params)

## setMaxIntakeRate ----
test_that("ssetMaxIntakeRate works", {
    expect_identical(setMaxIntakeRate(params, params@intake_max), params)
    params@species_params$h <- 2 * params@species_params$h
    p2 <- setMaxIntakeRate(params)
    expect_identical(2 * params@intake_max, p2@intake_max)
})
test_that("Comment works on intake_max", {
    comment(params@intake_max) <- "test"
    params <- setMaxIntakeRate(params, intake_max = params@intake_max)
    expect_identical(comment(params@intake_max), "test")
    expect_message(setMaxIntakeRate(params), NA)
    params@species_params$h <- 1
    expect_message(setMaxIntakeRate(params),
                   "has been commented")
})

# getMaxIntakeRate ----
test_that("getMaxIntakeRate works", {
    p <- setMaxIntakeRate(params, intake_max = getMaxIntakeRate(params))
    expect_identical(params, p)
})