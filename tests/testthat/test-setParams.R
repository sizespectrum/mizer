params <- NS_params
no_sp <- nrow(params@species_params)

## setParams ----
test_that("setParams can leave params unchanged", {
    expect_equal(setParams(params), params)
})
