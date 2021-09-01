params <- NS_params

## setExtMort ----
test_that("setExtMort works", {
    params <- NS_params
    params@species_params$z0 <- 2 * params@species_params$z0
    p2 <- setExtMort(params)
    expect_identical(2 * params@mu_b, p2@mu_b)
    # only mu_b changed
    p2@mu_b <- params@mu_b
    expect_unchanged(p2, params)
    
    # supplying ext_mort
    p2 <- setExtMort(params, 3 * params@mu_b)
    comment(p2@mu_b) <- NULL
    expect_equal(p2@mu_b, 3 * params@mu_b)
})
test_that("Comment works on mu_b", {
    params <- NS_params
    # if no comment, it is set automatically
    ext_mort <- params@mu_b
    params <- setExtMort(params, ext_mort = ext_mort)
    expect_identical(comment(params@mu_b), "set manually")
    
    # comment is stored
    comment(ext_mort) <- "test"
    params <- setExtMort(params, ext_mort = ext_mort)
    expect_identical(comment(params@mu_b), "test")
    
    # if no comment, previous comment is kept
    comment(ext_mort) <- NULL
    params <- setExtMort(params, ext_mort = ext_mort)
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
    expect_warning(setExtMort(params, ext_mort = ext_mort,
                                    reset = TRUE),
                   "Because you set `reset = TRUE`, the")
})

# getExtMort ----
test_that("getExtMort works", {
    expect_identical(getExtMort(NS_params),
                     NS_params@mu_b)
})

test_that("Can get and set slot", {
    params <- NS_params
    ext_mort <- getExtMort(params)
    expect_identical(ext_mort(params), ext_mort)
    new <- 2 * ext_mort
    comment(new) <- "test"
    ext_mort(params) <- new
    expect_identical(ext_mort(params), new)
})