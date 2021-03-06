params <- NS_params
no_sp <- nrow(params@species_params)

## setExtMort ----
test_that("setExtMort works", {
    expect_identical(setExtMort(params, params@mu_b, comment_z0 = NULL), params)
    params@species_params$z0 <- 2 * params@species_params$z0
    p2 <- setExtMort(params)
    expect_identical(2 * params@mu_b, p2@mu_b)
})
test_that("Comment works on mu_b", {
    mu_b <- params@mu_b
    comment(mu_b) <- "test"
    params <- setExtMort(params, z0 = mu_b)
    expect_identical(comment(params@mu_b), "test")
    
    # no message when nothing changes
    expect_message(setExtMort(params), NA)
    # but message when a change is not stored due to comment
    params@species_params$z0 <- 1
    expect_message(setExtMort(params),
                   "has been commented")
    
    params <- setExtMort(params, z0 = mu_b,
                         comment_z0 = "overwrite")
    expect_identical(comment(params@mu_b), "test")
    # but it is used otherwise
    comment(mu_b) <- NULL
    params <- setExtMort(params, z0 = mu_b,
                         comment_z0 = "overwrite")
    expect_identical(comment(params@mu_b), "overwrite")
})

# getExtMort ----
test_that("getExtMort works", {
    p <- setExtMort(params, z0 = getExtMort(params), comment_z0 = NULL)
    expect_identical(params, p)
})