params <- NS_params
no_sp <- nrow(params@species_params)

## setExtMort ----
test_that("setExtMort works", {
    expect_identical(setExtMort(params, params@mu_b), params)
    params@species_params$z0 <- 2 * params@species_params$z0
    p2 <- setExtMort(params)
    expect_identical(2 * params@mu_b, p2@mu_b)
})
test_that("Comment works on mu_b", {
    comment(params@mu_b) <- "test"
    params <- setExtMort(params, z0 = params@mu_b)
    expect_identical(comment(params@mu_b), "test")
    expect_message(setExtMort(params), NA)
    params@species_params$z0 <- 1
    expect_message(setExtMort(params),
                   "has been commented")
})