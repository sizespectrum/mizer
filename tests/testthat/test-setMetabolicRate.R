params <- NS_params
no_sp <- nrow(params@species_params)

## setMetabolicRate ----
test_that("setMetabolicRate works", {
    expect_identical(setMetabolicRate(params, params@metab), params)
    params@species_params$ks <- 2 * params@species_params$ks
    p2 <- setMetabolicRate(params)
    expect_identical(2 * params@metab, p2@metab)
})
test_that("setMetabolicRate can set exponent p", {
    # no change where p is already set in species_params
    params <- setMetabolicRate(params, p = 1)
    expect_identical(params@species_params$p, rep(0.7, 12))
    # but change where it is not
    params@species_params$p[[1]] <- NA
    params <- setMetabolicRate(params, p = 1)
    expect_identical(params@species_params$p, c(1, rep(0.7, 11)))
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

# getMetabolicRate ----
test_that("getMetabolicRate works", {
    p <- setMetabolicRate(params, metab = getMetabolicRate(params))
    expect_identical(params, p)
})