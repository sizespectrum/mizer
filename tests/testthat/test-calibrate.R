# calibrateBiomass ----
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
    expect_equal(calibrateBiomass(params), params)
    # Even if only partially observed
    species_params(params)$biomass_observed[1:5] <- NA
    expect_equal(calibrateBiomass(params), params)
    # If we double the observations, we get twice the abundance
    species_params(params)$biomass_observed <- 
        species_params(params)$biomass_observed * 2
    params2 <- calibrateBiomass(params)
    expect_equal(params2@initial_n, params@initial_n * 2)
    # We don't need to check other slots because this function uses
    # `scaleModel()` which is unit-tested separately.
})

# scaleModel ----
test_that("scaleModel does not change dynamics.", {
    factor <- 10
    sim <- project(NS_params, t_max = 1)
    params2 <- scaleModel(NS_params, factor)
    sim2 <- project(params2, t_max = 1)
    expect_equal(sim2@n[1, , ], sim@n[1, , ] * factor)
    expect_equal(sim2@n[2, , ], sim@n[2, , ] * factor)
})
