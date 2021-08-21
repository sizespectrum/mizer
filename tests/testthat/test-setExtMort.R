params <- NS_params

## setExtMort ----
test_that("setExtMort works", {
    params@species_params$z0 <- 2 * params@species_params$z0
    p2 <- setExtMort(params)
    expect_identical(2 * params@mu_b, p2@mu_b)
    # only mu_b changed
    p2@mu_b <- params@mu_b
    expect_identical(p2, params)
})
test_that("Comment works on mu_b", {
    params <- NS_params
    # if no comment, it is set automatically
    z0 <- params@mu_b
    params <- setExtMort(params, z0 = z0)
    expect_identical(comment(params@mu_b), "set manually")
    
    # comment is stored
    comment(z0) <- "test"
    params <- setExtMort(params, z0 = z0)
    expect_identical(comment(params@mu_b), "test")
    
    # if no comment, previous comment is kept
    comment(z0) <- NULL
    params <- setExtMort(params, z0 = z0)
    expect_identical(comment(params@mu_b), "test")
    
    # no message when nothing changes
    expect_message(setExtMort(params), NA)
    # but message when a change is not stored due to comment
    params@species_params$z0 <- 1
    expect_message(setExtMort(params),  "has been commented")
    # Can reset
    p <- setExtMort(params, reset = TRUE)
    expect_equal(p@mu_b[1, 1], 1,
                 check.attributes = FALSE)
    expect_warning(setExtMort(params, z0 = z0,
                                    reset = TRUE),
                   "Because you set `reset = TRUE`, the")
})

# getExtMort ----
test_that("getExtMort works", {
    expect_identical(getExtMort(NS_params),
                     NS_params@mu_b)
})