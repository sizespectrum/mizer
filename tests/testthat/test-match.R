# matchBiomasses ----
test_that("matchBiomasses works", {
    params <- NS_params
    
    # Does nothing when no observed biomass
    expect_identical(matchBiomasses(params), params)
    species_params(params)$biomass_observed <- NA
    expect_identical(matchBiomasses(params), params)
    # Does nothing if observed already equals model
    species_params(params)$biomass_cutoff <- 1e-4
    biomass_actual <-
        rowSums(sweep(params@initial_n, 2, params@w * params@dw, "*"))
    species_params(params)$biomass_observed <- biomass_actual
    expect_equal(matchBiomasses(params), params)
    # Even if only partially observed
    species_params(params)$biomass_observed[1:5] <- NA
    expect_equal(matchBiomasses(params), params)
    
    # If we double the observations, we get twice the abundance
    species <- 1:9
    species_params(params)$biomass_observed <- 
        species_params(params)$biomass_observed * 2
    params2 <- matchBiomasses(params, species)
    expect_equal(params2@initial_n[6:9, ], params@initial_n[6:9, ] * 2)
    # but unobserved species don't change
    expect_equal(params2@initial_n[1:5, ], params@initial_n[1:5, ])
    # and unselected species don't change
    expect_equal(params2@initial_n[10:12, ], params@initial_n[10:12, ])
    
    # Actually matches biomasses
    species_params(params)$biomass_observed <- 1:12
    params2 <- matchBiomasses(params)
    biomasses_new <-
        rowSums(sweep(params2@initial_n, 2, params@w * params@dw, "*"))
    expect_equal(biomasses_new, 1:12, check.attributes = FALSE)
    
})

# matchYields ----
test_that("matchYields works", {
    params <- NS_params
    
    # Does nothing when no observed yield
    expect_identical(matchYields(params), params)
    species_params(params)$yield_observed <- NA
    expect_identical(matchYields(params), params)
    # Does nothing if observed already equals model
    yield_actual <-
        rowSums(sweep(params@initial_n, 2, params@w * params@dw, "*"))
    species_params(params)$yield_observed <- yield_actual
    expect_equal(matchYields(params), params)
    # Even if only partially observed
    species_params(params)$yield_observed[1:5] <- NA
    expect_equal(matchYields(params), params)
    
    # If we double the observations, we get twice the abundance
    species <- 1:9
    species_params(params)$yield_observed <- 
        species_params(params)$yield_observed * 2
    params2 <- matchYields(params, species)
    expect_equal(params2@initial_n[6:9, ], params@initial_n[6:9, ] * 2)
    # but unobserved species don't change
    expect_equal(params2@initial_n[1:5, ], params@initial_n[1:5, ])
    # and unselected species don't change
    expect_equal(params2@initial_n[10:12, ], params@initial_n[10:12, ])
    
    # Actually matches yieldes
    species_params(params)$yield_observed <- 1:12
    params2 <- matchYields(params)
    yieldes_new <-
        rowSums(sweep(params2@initial_n, 2, params@w * params@dw, "*"))
    expect_equal(yieldes_new, 1:12, check.attributes = FALSE)
    
})
