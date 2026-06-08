# setBevertonHolt ----
NS_required_rdd <- getRequiredRDD(NS_params)
NS_rdd <- getRDD(NS_params)

test_that("setBevertonHolt sets erepro correctly when setting all values", {
    params <- setBevertonHolt(NS_params,
                              erepro = 10 * NS_params@species_params$erepro)
    expect_identical(params@species_params$erepro,
                     10 * NS_params@species_params$erepro)
    rdd <- getRDD(params)
    expect_equal(getRequiredRDD(params), rdd)
})

test_that("setBevertonHolt sets erepro correctly when setting same value for all species", {
    sp_name <- NS_params@species_params$species[3]
    expect_warning(params <- setBevertonHolt(NS_params, erepro = 0.1),
                   "For the following species `erepro` has been")
    expect_identical(params@species_params$R_max[params@species_params$species == sp_name],
                     Inf)
    expect_equal(NS_required_rdd, getRDD(params))
})

test_that("setBevertonHolt sets erepro correctly when setting same value for single species by using named vectors", {
    params <- NS_params
    erepro_old <- params@species_params$erepro
    erepro_new <- 0.1
    names(erepro_new) <- params@species_params$species[1]
    params <- setBevertonHolt(NS_params, erepro = erepro_new)
    expect_identical(params@species_params$erepro[2:3], erepro_old[2:3])
    expect_identical(params@species_params$erepro[1], 0.1)
})

test_that("setBevertonHolt sets erepro correctly when setting same value for some species by using named vector", {
    params <- NS_params
    erepro_old <- params@species_params$erepro
    erepro_new <- c(0.1, 0.2, NA)
    names(erepro_new) <- params@species_params$species[1:3]
    params <- setBevertonHolt(NS_params, erepro = erepro_new)
    expect_identical(params@species_params$erepro[3], erepro_old[3])
    expect_identical(params@species_params$erepro[1:2], c(0.1, 0.2))
    expect_equal(NS_required_rdd[1:2], getRDD(params)[1:2])
})

test_that("setBevertonHolt sets erepro correctly when setting same value for some species some species by NAs", {
    params <- NS_params
    erepro_old <- params@species_params$erepro
    erepro_new <- erepro_old
    erepro_new[3] <- NA
    erepro_new[1:2] <- c(0.1, 0.2)
    params <- setBevertonHolt(NS_params, erepro = erepro_new)
    expect_identical(params@species_params$erepro[3], erepro_old[3])
    expect_identical(params@species_params$erepro[1:2], c(0.1, 0.2))
})

test_that("setBevertonHolt updates `time_modified`", {
    expect_false(identical(setBevertonHolt(NS_params, erepro = 1)@time_modified, 
                           NS_params@time_modified))
})

test_that("setBevertonHolt without extra arguments keeps previous erepro", {
    params <- setBevertonHolt(NS_params)
    expect_identical(params@species_params$erepro,
                     NS_params@species_params$erepro)
})

# R_max ----
test_that("setBevertonHolt sets R_max correctly when setting all values", {
    params <- setBevertonHolt(NS_params,
                              R_max = 10 * NS_params@species_params$R_max)
    expect_identical(params@species_params$R_max,
                     10 * NS_params@species_params$R_max)
    rdd <- getRDD(params)
    expect_equal(getRequiredRDD(params), rdd)
})

test_that("setBevertonHolt sets R_max correctly when setting values for some species by using named vector", {
    params <- NS_params
    R_max_old <- params@species_params$R_max
    R_max_new <- c(1e12, 1e13, NA)
    names(R_max_new) <- params@species_params$species[1:3]
    params <- setBevertonHolt(NS_params, R_max = R_max_new)
    expect_identical(params@species_params$R_max[3], R_max_old[3])
    expect_identical(params@species_params$R_max[1:2], c(1e12, 1e13))
})

test_that("setBevertonHolt issues warning when an R_max leads to an erepro > 1", {
    sp_name <- NS_params@species_params$species[3]
    rdi_sp <- getRDI(NS_params)[sp_name]
    req_rdd_sp <- getRequiredRDD(NS_params)[sp_name]
    # Set R_max just below the threshold where erepro > 1 is needed
    # Use a named single-species vector to avoid floating-point issues with other species
    threshold <- req_rdd_sp * rdi_sp / (rdi_sp - req_rdd_sp)
    R_max_new <- setNames(threshold * 0.99, sp_name)
    expect_warning(params <- setBevertonHolt(NS_params, R_max = R_max_new),
                   paste0("The following species require an unrealistic value greater than 1 for `erepro`: ", sp_name))
    expect_gt(params@species_params$erepro[params@species_params$species == sp_name], 1)
    expect_equal(params@species_params$R_max[3], unname(R_max_new))
})

# reproduction_level ----
test_that("setBevertonHolt sets reproduction_level correctly", {
    sp_name <- NS_params@species_params$species[3]
    expect_warning(params <- setBevertonHolt(NS_params, reproduction_level = 0.8),
                   paste0("The following species require an unrealistic value greater than 1 for `erepro`: ", sp_name))
    rdd <- getRDD(params)
    expect_equal(rdd, params@species_params$R_max * 0.8, ignore_attr = TRUE)
    expect_equal(getRequiredRDD(params), rdd)
    expect_equal(getReproductionLevel(params)[[1]], 0.8)
})

test_that("getReproductionLevel is getRDD divided by R_max", {
    expect_equal(getReproductionLevel(NS_params),
                 NS_rdd / NS_params@species_params$R_max)
})

# R_factor ----
test_that("setBevertonHolt sets R_factor correctly", {
    sp_name <- NS_params@species_params$species[3]
    expect_warning(params <- setBevertonHolt(NS_params, R_factor = 1.5),
                   paste0("The following species require an unrealistic value greater than 1 for `erepro`: ", sp_name))
    rdd <- getRDD(params)
    expect_equal(rdd, params@species_params$R_max / 1.5, ignore_attr = TRUE)
    expect_equal(getRequiredRDD(params), rdd)
})

test_that("setRmax is an exact alias for setBevertonHolt", {
    p1 <- suppressWarnings(setRmax(NS_params, reproduction_level = 0.3))
    p2 <- suppressWarnings(setBevertonHolt(NS_params, reproduction_level = 0.3))
    expect_equal(species_params(p1)$erepro, species_params(p2)$erepro)
    expect_equal(species_params(p1)$R_max, species_params(p2)$R_max)
    expect_identical(getReproductionLevel(p1), getReproductionLevel(p2))
})

# special values ----
test_that("setBevertonHolt does nothing when called with only NA values", {
    erepro_new <- rep(NA, 3)
    params <- setBevertonHolt(NS_params, erepro = erepro_new)
    expect_identical(params, NS_params)
    erepro_new <- NA
    sp_name <- NS_params@species_params$species[3]
    names(erepro_new) <- sp_name
    params <- setBevertonHolt(NS_params, erepro = erepro_new)
    expect_identical(params, NS_params)
})

test_that("setBevertonHolt throws error when called with invalid values", {
    expect_error(setBevertonHolt(NS_params, erepro = rep("x", 3)),
                 "You provided invalid non-numeric values.")
    expect_error(setBevertonHolt(NS_params, erepro = 1, R_max = 2),
                 "You should only provide `params` and one other argument.")
    expect_error(setBevertonHolt(NS_params, R_max = 1:4),
                 "You need to supply a vector of length 3 or a single number or a named vector.")
    expect_error(setBevertonHolt(NS_params, reproduction_level = 2),
                 "The reproduction level must be smaller than 1 and non-negative.")
    expect_error(setBevertonHolt(NS_params, R_factor = 1/2),
                 "The R_factor must be greater than 1.")
})

test_that("reproduction_level of 0 works", {
    params <- setBevertonHolt(NS_params, reproduction_level = 0)
    expect_equal(getRDD(params), getRequiredRDD(params))
    expect_equal(params@species_params$R_max[1], Inf)
    
    params <- setBevertonHolt(NS_params, R_factor = Inf)
    expect_equal(getRDD(params), getRequiredRDD(params))
    expect_equal(params@species_params$R_max[1], Inf)
})

test_that("R_max is increased when needed", {
    sp_name1 <- NS_params@species_params$species[1]
    sp_name2 <- NS_params@species_params$species[2]
    expect_warning(p <- setBevertonHolt(NS_params, R_max = c(1, 2, NA)),
                 paste0("has been increased to give a reproduction level of 0.99: ", sp_name1, ", ", sp_name2))
    # R_max was increased from the requested values (1 and 2) to rdd/0.99
    expect_gt(p@species_params$R_max[1], 1)
    expect_gt(p@species_params$R_max[2], 2)
})


