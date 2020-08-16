params <- NS_params
no_sp <- nrow(params@species_params)

# setReproduction ----
test_that("setReproduction works", {
    expect_equal(setReproduction(params), params)
    maturity <- array(1, dim = c(no_sp, length(params@w)))
    p2 <- setReproduction(params, maturity = maturity)
    expect_equal(p2, setReproduction(p2, maturity = maturity,
                                     repro_prop = p2@psi))
    expect_equal(params, setReproduction(params, repro_prop = p2@psi))
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
    params <- setReproduction(NS_params)
    comment(params@maturity) <- "test"
    params <- setReproduction(params, maturity = params@maturity)
    expect_identical(comment(params@maturity), "test")
    expect_message(setReproduction(params), NA)
    params@species_params$w_mat <- params@species_params$w_mat + 1
    expect_message(setReproduction(params),
                   "maturity ogive has been commented")
})
test_that("Comment works on psi", {
    params <- NS_params
    repro_prop <- params@psi
    repro_prop[] <- 1
    comment(repro_prop) <- "test"
    params <- setReproduction(params, repro_prop = repro_prop)
    expect_identical(comment(params@psi), "test")
    # We use entry 50 below to be above the low-value cutoff in psi
    expect_equal(params@psi[1, 50], params@maturity[1, 50])
    expect_message(setReproduction(params),
                   "has been commented")
    # Test that we can still update maturity
    maturity <- getMaturityProportion(params)
    maturity[] <- 1
    params <- setReproduction(params, maturity = maturity)
    expect_equal(params@psi[1, 50], 1)
})

# getMaturityProportion ----
test_that("getMaturityProportion works", {
    params <- setReproduction(NS_params)
    maturity <- getMaturityProportion(params)
    params2 <- setReproduction(params, maturity =  maturity)
    expect_identical(params, params2)
})

# getReproductionProportion ----
test_that("getReproductionProportion works", {
    params <- setReproduction(NS_params)
    repro_prop <- getReproductionProportion(params)
    params2 <- setReproduction(params, repro_prop = repro_prop)
    expect_identical(params, params2)
})