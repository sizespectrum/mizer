local_edition(3)

test_that("calibrateBiomass works", {
    params <- NS_params
    # Does nothing when no observed biomass
    expect_identical(calibrateBiomass(params), params)
    species_params(params)$biomass_observed <- NA
    expect_identical(calibrateBiomass(params), params)
    # Does nothing if observed already equals model
    species_params(params)$biomass_cutoff <- 1e-4
    species_params(params)$biomass_observed <- 
        rowSums(sweep(params@initial_n, 2, params@w * params@dw, "*"))
    expect_unchanged(calibrateBiomass(params), params)
    # Even if only partially observed
    species_params(params)$biomass_observed[1:5] <- NA
    expect_unchanged(calibrateBiomass(params), params)
    # If we double the observations, we get twice the abundance
    species_params(params)$biomass_observed <- 
        species_params(params)$biomass_observed * 2
    params2 <- calibrateBiomass(params)
    expect_equal(params2@initial_n, params@initial_n * 2)
    # We don't need to check other slots because this function uses
    # `scaleModel()` which is unit-tested separately.
})


test_that("calibrateNumber works", {
    params <- NS_params
    # Does nothing when no observed Number
    expect_identical(calibrateNumber(params), params)
    species_params(params)$number_observed <- NA
    expect_identical(calibrateNumber(params), params)
    # Does nothing if observed already equals model
    species_params(params)$number_cutoff <- 1e-4
    species_params(params)$number_observed <- 
        rowSums(sweep(params@initial_n, 2, params@dw, "*"))
    expect_unchanged(calibrateNumber(params), params)
    # Even if only partially observed
    species_params(params)$number_observed[1:5] <- NA
    expect_unchanged(calibrateNumber(params), params)
    # If we double the observations, we get twice the abundance
    species_params(params)$number_observed <- 
        species_params(params)$number_observed * 2
    params2 <- calibrateNumber(params)
    expect_equal(params2@initial_n, params@initial_n * 2)
    # We don't need to check other slots because this function uses
    # `scaleModel()` which is unit-tested separately.
})


test_that("calibrateYield works", {
    params <- NS_params
    # Does nothing when no observed yield
    expect_identical(calibrateYield(params), params)
    species_params(params)$yield_observed <- NA
    expect_identical(calibrateYield(params), params)
    # Does nothing if observed already equals model
    biomass <- sweep(params@initial_n, 2, params@w * params@dw, "*")
    yield_model <- rowSums(biomass * getFMort(params))
    species_params(params)$yield_observed <- yield_model
    expect_unchanged(calibrateYield(params), params)
    # Even if only partially observed
    species_params(params)$yield_observed[1:5] <- NA
    expect_unchanged(calibrateYield(params), params)
    # If we double the observations, we get twice the abundance
    species_params(params)$yield_observed <- 
        species_params(params)$yield_observed * 2
    params2 <- calibrateYield(params)
    expect_equal(params2@initial_n, params@initial_n * 2)
    # If we double the catchability as well, we get the original abundance
    gear_params(params)$catchability <- 
        gear_params(params)$catchability * 2
    params2 <- calibrateYield(params)
    expect_equal(params2@initial_n, params@initial_n)
    # We don't need to check other slots because this function uses
    # `scaleModel()` which is unit-tested separately.
})


test_that("scaleModel does not change dynamics.", {
    factor <- 10
    sim <- project(NS_params, t_max = 1)
    params2 <- scaleModel(NS_params, factor)
    sim2 <- project(params2, t_max = 1)
    expect_equal(sim2@n[1, , ], sim@n[1, , ] * factor)
    expect_equal(sim2@n[2, , ], sim@n[2, , ] * factor)
})
