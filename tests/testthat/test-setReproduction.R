params <- setReproduction(NS_params_small)
no_sp <- nrow(params@species_params)

# setReproduction ----
test_that("setReproduction works", {
    expect_unchanged(setReproduction(params), params)
    maturity <- array(1, dim = c(no_sp, length(params@w)))
    p2 <- setReproduction(params, maturity = maturity)
    p3 <- setReproduction(p2, maturity = maturity,
                          repro_prop = p2@psi)
    comment(p3@psi) <- comment(p2@psi)
    expect_unchanged(p2, p3)
    p3 <- setReproduction(params, repro_prop = repro_prop(params))
    comment(p3@psi) <- comment(params@psi)
    expect_unchanged(params, p3)
    expect_error(setReproduction(params, RDD = "str"),
                 "Arguments of RDD function can only contain 'rdi', 'species_params' and `t`.")
    params@species_params$erepro[1] <- NA
    p2 <- setReproduction(params, RDD = "SheperdRDD")
    expect_equal(p2@species_params$erepro[1], 1, ignore_attr = TRUE)
    p2@species_params$sheperd_b <- 0
    expect_error(getRDD(p2),
                 "The species_params dataframe must contain columns sheperd_b and sheperd_c.")
    p2 <- setReproduction(params, RDD = "RickerRDD")
    expect_error(getRDD(p2),
                 "The ricker_b column is missing in species_params")
    p2@species_params$ricker_b <- 0
    expect_equal(getRDI(p2), getRDD(p2))
})
test_that("setReproduction checks arguments", {
    params <- NS_params_small
    params@species_params$w_max[[2]] <- NA
    params@species_params$w_max[[2]] <- NA
    sp_name <- params@species_params$species[2]
    expect_error(setReproduction(params),
                 paste0("The following species are missing their upper size-grid boundary `w_max`: ", sp_name))
    params@species_params$w_max[[2]] <- 1e-5
    expect_error(setReproduction(params),
                 "Some of the upper size-grid boundaries \\(`w_max`\\) are smaller than the egg sizes.")
    params@species_params$w_max <- NULL
    expect_error(setReproduction(params),
                 "The upper size-grid boundary must be specified in the `w_max` column of the species parameter data frame.")

    params <- NS_params_small
    params@species_params$w_mat[[3]] <- NA
    sp_name <- params@species_params$species[3]
    expect_message(pa <- setReproduction(params),
                 paste0("Note: The following species were missing data for their maturity size w_mat: ", sp_name, "."))
})

# * Comments ----
test_that("Comment works on maturity", {
    params <- NS_params_small
    # if no comment, it is set automatically
    maturity <- params@maturity
    params <- setReproduction(params, maturity = maturity)
    expect_identical(comment(params@maturity), "set manually")

    # comment is stored
    comment(maturity) <- "test"
    params <- setReproduction(params, maturity = maturity)
    expect_identical(comment(params@maturity), "test")

    # if no comment, previous comment is kept
    comment(maturity) <- NULL
    params <- setReproduction(params, maturity = maturity)
    expect_identical(comment(params@maturity), "test")

    # no message when nothing changes
    expect_message(setReproduction(params), NA)
    # but message when a change is not stored due to comment
    params@species_params$w_mat <- params@species_params$w_mat * 1.1
    expect_message(setReproduction(params),  "has been commented")
    # Can reset
    p <- setReproduction(params, reset = TRUE)
    # The increase in w_mat should lower the maturity curve
    expect_gt(max(params@maturity - p@maturity), 0)
    expect_warning(setReproduction(params, maturity = maturity,
                                    reset = TRUE),
                   "Because you set `reset = TRUE`, the")
})

test_that("Comment works on psi", {
    params <- NS_params_small
    # if no comment, it is set automatically
    repro_prop <- getReproductionProportion(params)
    params <- setReproduction(params, repro_prop = repro_prop)
    expect_identical(comment(params@psi), "set manually")

    # comment is stored
    comment(repro_prop) <- "test"
    params <- setReproduction(params, repro_prop = repro_prop)
    expect_identical(comment(params@psi), "test")

    # if no comment, previous comment is kept
    comment(repro_prop) <- NULL
    params <- setReproduction(params, repro_prop = repro_prop)
    expect_identical(comment(params@psi), "test")

    # no message when nothing changes
    expect_message(setReproduction(params), NA)
    # but message when a change is not stored due to comment
    params@species_params$w_max <- params@species_params$w_max / 1.1
    params@species_params$w_repro_max <- params@species_params$w_repro_max / 1.1
    expect_message(setReproduction(params),  "has been commented")
    # Can reset
    p <- setReproduction(params, reset = TRUE)
    # The decrease in w_repro_max should increase the psi curve
    expect_lt(min(params@psi - p@psi), 0)
    expect_warning(setReproduction(params, repro_prop = repro_prop,
                                    reset = TRUE),
                   "Because you set `reset = TRUE`, the")
})

# getMaturityProportion ----
test_that("getMaturityProportion works", {
    params <- setReproduction(NS_params_small)
    maturity <- getMaturityProportion(params)
    expect_identical(maturity(params), maturity)
    params2 <- setReproduction(params, maturity =  maturity)
    comment(params2@maturity) <- NULL
    expect_unchanged(params, params2)
})

# getReproductionProportion ----
test_that("getReproductionProportion works", {
    params <- setReproduction(NS_params_small)
    repro_prop <- getReproductionProportion(params)
    expect_identical(repro_prop(params), repro_prop)
    params2 <- setReproduction(params, repro_prop = repro_prop)
    comment(params2@psi) <- NULL
    expect_unchanged(params, params2)
})

test_that("getReproductionProportion returns a proportion",{
    params <- NS_params_small
    # Make extremely wide maturity ogive
    species_params(params)$w_mat25 <- 1
    repro_prop <- getReproductionProportion(params)
    expect_true(all(repro_prop >= 0))
    expect_true(all(repro_prop <= 1))
})

test_that("Can get and set repro_prop", {
    params <- NS_params_small
    new <- repro_prop(params) ^ 2
    comment(new) <- "test"
    repro_prop(params) <- new
    expect_equal(repro_prop(params)[2, 10], new[2, 10])
    expect_equal(getReproductionProportion(params)[2, 10], new[2, 10])
})

test_that("Can get and set maturity", {
    params <- NS_params_small
    new <- 1/2 * maturity(params)
    comment(new) <- "test"
    maturity(params) <- new
    expect_equal(maturity(params), new, ignore_attr = TRUE)
    expect_equal(getMaturityProportion(params), new, ignore_attr = TRUE)
    expect_identical(comment(params@maturity), "test")
})
