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
    expect_identical(params@species_params$p, rep(0.7, 3), ignore_attr = TRUE)
    # but change where it is not
    params@species_params$p[[1]] <- NA
    params <- setMetabolicRate(params, p = 1)
    expect_identical(params@species_params$p, c(1, rep(0.7, 2)), ignore_attr = TRUE)
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
    # `p` defaults to `n`. NS_params_small has n = 2/3, so this also pins that
    # the default tracks `n` rather than a hardcoded 3/4.
    params <- NS_params_small
    params@species_params$p[] <- NA
    params <- setMetabolicRate(params)
    expect_equal(species_params(params)$p, species_params(params)$n,
                 ignore_attr = TRUE)

    expect_error(setMetabolicRate(NS_params_small, p = "x"), "p must be numeric")

    # The default must not depend on how the model was reached: filling `p`
    # via a standalone setter call must agree with the construction path.
    sp <- data.frame(species = c("a", "b"), w_inf = c(100, 200), n = 0.7)
    params <- suppressMessages(newMultispeciesParams(sp))
    expect_equal(species_params(params)$p, c(0.7, 0.7), ignore_attr = TRUE)

    params@species_params$p[] <- NA
    params <- suppressMessages(setMetabolicRate(params))
    expect_equal(species_params(params)$p, c(0.7, 0.7), ignore_attr = TRUE)

    new <- metab(NS_params_small)
    bad_names <- new
    dimnames(bad_names)[[1]] <- rev(dimnames(bad_names)[[1]])
    expect_error(setMetabolicRate(NS_params_small, metab = bad_names),
                 "same ordering of species")

    bad_values <- new
    bad_values[1, 1] <- -1
    expect_error(setMetabolicRate(NS_params_small, metab = bad_values))
})
