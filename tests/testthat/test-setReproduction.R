params <- NS_params
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
    p3 <- setReproduction(params, repro_prop = p2@psi)
    comment(p3@psi) <- comment(params@psi)
    expect_unchanged(params, p3)
    expect_error(setReproduction(params, RDD = "str"),
                 "Arguments of RDD function can only contain 'rdi', 'species_params' and `t`.")
    expect_error(setReproduction(params, RDD = "sum"),
                 "The RDD function needs to have at least arguments `rdi` and `...`.")
    params@species_params$erepro[1] <- NA
    p2 <- setReproduction(params, RDD = "SheperdRDD")
    expect_equal(p2@species_params$erepro[1], 1)
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
    params <- NS_params
    params@species_params$w_inf[[2]] <- NA
    expect_error(setReproduction(params),
                 "The following species are missing data for their maximum size w_inf: Sandeel")
    params@species_params$w_inf[[2]] <- 1e-5
    expect_error(setReproduction(params),
                 "Some of the asymptotic sizes are smaller than the egg sizes.")
    params@species_params$w_inf <- NULL
    expect_error(setReproduction(params),
                 "The maximum sizes of the species must be specified in the w_inf column of the species parameter data frame.")

    params <- NS_params
    params@species_params$w_mat[[2]] <- NA
    expect_message(pa <- setReproduction(params),
                 "Note: The following species were missing data for their maturity size w_mat: Sandeel.")
})

# * Comments ----
test_that("Comment works on maturity", {
    params <- NS_params
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
    params <- NS_params
    # if no comment, it is set automatically
    repro_prop <- getReproductionProportion(params)
    params <- setReproduction(params, repro_prop = repro_prop)
    expect_identical(comment(getReproductionProportion(params)), "set manually")
    
    # comment is stored
    comment(repro_prop) <- "test"
    params <- setReproduction(params, repro_prop = repro_prop)
    expect_identical(comment(getReproductionProportion(params)), "test")
    
    # if no comment, previous comment is kept
    comment(repro_prop) <- NULL
    params <- setReproduction(params, repro_prop = repro_prop)
    expect_identical(comment(getReproductionProportion(params)), "test")
    
    # no message when nothing changes
    expect_message(setReproduction(params), NA)
    # but message when a change is not stored due to comment
    params@species_params$w_inf <- params@species_params$w_inf * 1.1
    expect_message(setReproduction(params),  "has been commented")
    # Can reset
    p <- setReproduction(params, reset = TRUE)
    # The increase in w_inf should lower the psi curve
    expect_gt(max(params@psi - p@psi), 0)
    expect_warning(setReproduction(params, repro_prop = repro_prop,
                                    reset = TRUE),
                   "Because you set `reset = TRUE`, the")
})

# getMaturityProportion ----
test_that("getMaturityProportion works", {
    params <- setReproduction(NS_params)
    maturity <- getMaturityProportion(params)
    params2 <- setReproduction(params, maturity =  maturity)
    comment(params2@maturity) <- NULL
    expect_unchanged(params, params2)
})

# getReproductionProportion ----
test_that("getReproductionProportion works", {
    params <- setReproduction(NS_params)
    repro_prop <- getReproductionProportion(params)
    params2 <- setReproduction(params, repro_prop = repro_prop)
    comment(params2@psi) <- NULL
    expect_unchanged(params, params2)
})

test_that("Can get and set repro_prop", {
    params <- NS_params
    new <- repro_prop(params) ^ 2 
    comment(new) <- "test"
    repro_prop(params) <- new
    expect_equal(repro_prop(params)[2,50], new[2, 50])
})
test_that("Can get and set maturity", {
    params <- NS_params
    new <- 1/2 * maturity(params)
    comment(new) <- "test"
    maturity(params) <- new
    expect_identical(maturity(params), new)
})