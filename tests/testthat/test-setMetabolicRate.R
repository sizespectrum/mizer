params <- NS_params

# setMetabolicRate ----
test_that("setMetabolicRate works", {
    params@species_params$ks <- 2 * params@species_params$ks
    p2 <- setMetabolicRate(params)
    expect_identical(2 * params@metab, p2@metab)
    # only metab changed
    p2@metab <- params@metab
    expect_unchanged(p2, params)
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
    params <- NS_params
    # if no comment, it is set automatically
    metab <- params@metab
    params <- setMetabolicRate(params, metab = metab)
    expect_identical(comment(params@metab), "set manually")
    
    # comment is stored
    comment(metab) <- "test"
    params <- setMetabolicRate(params, metab = metab)
    expect_identical(comment(params@metab), "test")
    
    # if no comment, previous comment is kept
    comment(metab) <- NULL
    params <- setMetabolicRate(params, metab = metab)
    expect_identical(comment(params@metab), "test")
    
    # no message when nothing changes
    expect_message(setMetabolicRate(params), NA)
    # but message when a change is not stored due to comment
    params@species_params$ks <- 1
    expect_message(setMetabolicRate(params),  "has been commented")
    # Can reset
    p <- setMetabolicRate(params, reset = TRUE)
    expect_equal(p@metab[, 1], params@w[1]^params@species_params$p,
                 ignore_attr = TRUE)
    expect_warning(setMetabolicRate(params, metab = metab,
                                    reset = TRUE),
                   "Because you set `reset = TRUE`, the")
})

# getMetabolicRate ----
test_that("getMetabolicRate works", {
    expect_identical(getMetabolicRate(NS_params),
                     NS_params@metab)
})

test_that("Can get and set metab slot", {
    params <- NS_params
    new <- 2 * metab(params)
    comment(new) <- "test"
    metab(params) <- new
    expect_identical(metab(params), new)
})
