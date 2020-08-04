params <- NS_params
no_sp <- nrow(params@species_params)

## setMetabolicRate ----
test_that("setMetabolicRate works", {
    expect_identical(setMetabolicRate(params, params@metab), params)
    params@species_params$ks <- 2 * params@species_params$ks
    p2 <- setMetabolicRate(params)
    expect_identical(2 * params@metab, p2@metab)
})
test_that("Comment works on metab", {
    comment(params@metab) <- "test"
    params <- setMetabolicRate(params, metab = params@metab)
    expect_identical(comment(params@metab), "test")
    expect_message(setMetabolicRate(params), NA)
    params@species_params$k <- 1
    expect_message(setMetabolicRate(params),
                   "has been commented")
})