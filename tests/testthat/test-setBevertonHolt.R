# setBevertonHolt ----
test_that("setBevertonHolt sets erepro correctly when setting all values", {
    params <- setBevertonHolt(NS_params,
                              erepro = 10 * NS_params@species_params$erepro)
    expect_identical(params@species_params$erepro,
                     10 * NS_params@species_params$erepro)
    expect_equal(getRequiredRDD(params), getRDD(params))
})

test_that("setBevertonHolt sets erepro correctly when setting same value for all species", {
    expect_warning(params <- setBevertonHolt(NS_params, erepro = 0.1),
                   "For the following species `erepro` has been")
    expect_identical(params@species_params$R_max[params@species_params$species == "Gurnard"],
                     Inf)
    expect_equal(getRequiredRDD(NS_params), getRDD(params))
})

test_that("setBevertonHolt sets erepro correctly when setting same value for single species by using named vectors", {
    params <- NS_params
    erepro_old <- params@species_params$erepro
    erepro_new <- 0.1
    names(erepro_new) <- params@species_params$species[1]
    params <- setBevertonHolt(NS_params, erepro = erepro_new)
    expect_identical(params@species_params$erepro[2:12], erepro_old[2:12])
    expect_identical(params@species_params$erepro[1], 0.1)
})

test_that("setBevertonHolt sets erepro correctly when setting same value for some species by using named vector", {
    params <- NS_params
    erepro_old <- params@species_params$erepro
    erepro_new <- c(0.1, 0.2, NA)
    names(erepro_new) <- params@species_params$species[1:3]
    params <- setBevertonHolt(NS_params, erepro = erepro_new)
    expect_identical(params@species_params$erepro[3:12], erepro_old[3:12])
    expect_identical(params@species_params$erepro[1:2], c(0.1, 0.2))
    expect_equal(getRequiredRDD(NS_params)[1:2], getRDD(params)[1:2])
})

test_that("setBevertonHolt sets erepro correctly when setting same value for some species some species by NAs", {
    params <- NS_params
    erepro_old <- params@species_params$erepro
    erepro_new <- erepro_old
    erepro_new[3:12] <- NA
    erepro_new[1:2] <- c(0.1, 0.2)
    params <- setBevertonHolt(NS_params, erepro = erepro_new)
    expect_identical(params@species_params$erepro[3:12], erepro_old[3:12])
    expect_identical(params@species_params$erepro[1:2], c(0.1, 0.2))
})

# R_max ----
test_that("setBevertonHolt sets R_max correctly when setting all values", {
    params <- setBevertonHolt(NS_params,
                              R_max = 10 * NS_params@species_params$R_max)
    expect_identical(params@species_params$R_max,
                     10 * NS_params@species_params$R_max)
    expect_equal(getRequiredRDD(params), getRDD(params))
})

test_that("setBevertonHolt sets R_max correctly when setting values for some species by using named vector", {
    params <- NS_params
    R_max_old <- params@species_params$R_max
    R_max_new <- c(1e12, 1e13, NA)
    names(R_max_new) <- params@species_params$species[1:3]
    params <- setBevertonHolt(NS_params, R_max = R_max_new)
    expect_identical(params@species_params$R_max[3:12], R_max_old[3:12])
    expect_identical(params@species_params$R_max[1:2], c(1e12, 1e13))
})

test_that("setBevertonHolt issues warning when an R_max leads to an erepro > 1", {
    R_max_new <- NS_params@species_params$R_max * 1.02
    expect_warning(params <- setBevertonHolt(NS_params, R_max = R_max_new),
                   "The following species require an unrealistic reproductive efficiency greater than 1: Plaice")
    expect_gt(params@species_params$erepro[params@species_params$species == "Plaice"], 1)
    expect_identical(params@species_params$R_max, R_max_new)
})

# reproduction_level ----
test_that("setBevertonHolt sets reproduction_level correctly", {
    expect_warning(params <- setBevertonHolt(NS_params, reproduction_level = 0.4),
                   "The following species require an unrealistic reproductive efficiency greater than 1: Plaice")
    expect_equal(getRDD(params), params@species_params$R_max * 0.4, ignore_attr = TRUE)
    expect_equal(getRequiredRDD(params), getRDD(params))
    expect_equal(getReproductionLevel(params)[[1]], 0.4)
})

# R_factor ----
test_that("setBevertonHolt sets R_factor correctly", {
    expect_warning(params <- setBevertonHolt(NS_params, R_factor = 4),
                   "The following species require an unrealistic reproductive efficiency greater than 1: Plaice")
    expect_equal(getRDD(params), params@species_params$R_max / 4, ignore_attr = TRUE)
    expect_equal(getRequiredRDD(params), getRDD(params))
})

# special values ----
test_that("setBevertonHolt does nothing when called with only NA values", {
    erepro_new <- rep(NA, 12)
    params <- setBevertonHolt(NS_params, erepro = erepro_new)
    expect_identical(params, NS_params)
    erepro_new <- NA
    names(erepro_new) <- "Cod"
    params <- setBevertonHolt(NS_params, erepro = erepro_new)
    expect_identical(params, NS_params)
})

test_that("setBevertonHolt throws error when called with invalid values", {
    expect_error(setBevertonHolt(NS_params, erepro = rep("x", 12)),
                 "You provided invalid non-numeric values.")
    expect_error(setBevertonHolt(NS_params, erepro = 1, R_max = 2),
                 "You should only provide `params` and one other argument.")
    expect_error(setBevertonHolt(NS_params, R_max = 1:4),
                 "You need to supply a vector of length 12 or a single number or a named vector.")
    expect_error(setBevertonHolt(NS_params, reproduction_level = 2),
                 "The reproduction level must be smaller than 1 and non-negative.")
    expect_error(setBevertonHolt(NS_params, R_factor = 1/2),
                 "The R_factor must be greater than 1.")
})

test_that("reproduction_level of 0 works", {
    params <- setBevertonHolt(NS_params, reproduction_level = 0)
    expect_equal(getRDD(params), mizer:::getRequiredRDD(params))
    expect_equal(params@species_params$R_max[1], Inf)
    
    params <- setBevertonHolt(NS_params, R_factor = Inf)
    expect_equal(getRDD(params), mizer:::getRequiredRDD(params))
    expect_equal(params@species_params$R_max[1], Inf)
})

test_that("R_max is increased when needed", {
    expect_warning(p <- setBevertonHolt(NS_params, R_max = c(1, 2, rep(NA, 10))),
                 "has been increased to give a reproduction level of 0.99: Sprat, Sandeel")
    expect_gt(p@species_params$R_max[1], NS_params@species_params$R_max[1])
})
