params <- NS_params_small

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
    expect_identical(params@species_params$p, rep(0.7, 3))
    # but change where it is not
    params@species_params$p[[1]] <- NA
    params <- setMetabolicRate(params, p = 1)
    expect_identical(params@species_params$p, c(1, rep(0.7, 2)))
})
test_that("Comment works on metab", {
    params <- NS_params_small
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
    expect_true(is.ArraySpeciesBySize(getMetabolicRate(NS_params_small)))
    expect_equal(getMetabolicRate(NS_params_small), NS_params_small@metab,
                 ignore_attr = TRUE)
})

test_that("Can get and set metab slot", {
    params <- NS_params_small
    new <- 2 * metab(params)
    comment(new) <- "test"
    metab(params) <- new
    expect_equal(metab(params), new, ignore_attr = TRUE)
    expect_identical(comment(params@metab), "test")
})

test_that("setMetabolicRate uses the documented defaults and validates inputs", {
    params <- NS_params_small
    params@species_params$p[] <- NA
    params <- setMetabolicRate(params)
    expect_identical(species_params(params)$p, rep(3 / 4, nrow(species_params(params))))

    expect_error(setMetabolicRate(NS_params_small, p = "x"), "p must be numeric")

    new <- metab(NS_params_small)
    bad_names <- new
    dimnames(bad_names)[[1]] <- rev(dimnames(bad_names)[[1]])
    expect_error(setMetabolicRate(NS_params_small, metab = bad_names),
                 "same ordering of species")

    bad_values <- new
    bad_values[1, 1] <- -1
    expect_error(setMetabolicRate(NS_params_small, metab = bad_values))
})
