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
    expect_message(
        expect_identical(matchYields(params), params),
        "not be changed: Sprat, Sandeel, N.pout, Herring, Dab, Whiting, Sole, Gurnard, Plaice, Haddock, Cod, Saithe.")
    # Does nothing if observed already equals model
    biomass <- sweep(params@initial_n, 2, params@w * params@dw, "*")
    yield_model <- rowSums(biomass * getFMort(params))
    species_params(params)$yield_observed <- yield_model
    expect_message(
        expect_equal(matchYields(params), params),
        "not be changed: Sprat, Sandeel, N.pout.")
    # Even if only partially observed
    species_params(params)$yield_observed[1:5] <- NA
    expect_message(
        expect_equal(matchYields(params), params),
        "not be changed: Sprat, Sandeel, N.pout, Herring, Dab.")
    
    # If we double the observations, we get twice the abundance
    species <- 1:9
    species_params(params)$yield_observed <- 
        species_params(params)$yield_observed * 2
    expect_message(
        params2 <- matchYields(params, species),
        "not be changed: Sprat, Sandeel, N.pout, Herring, Dab.")
    expect_equal(params2@initial_n[6:9, ], params@initial_n[6:9, ] * 2)
    # but unobserved species don't change
    expect_equal(params2@initial_n[1:5, ], params@initial_n[1:5, ])
    # and unselected species don't change
    expect_equal(params2@initial_n[10:12, ], params@initial_n[10:12, ])
    
    # Actually matches yieldes
    species_params(params)$yield_observed <- 1:12
    expect_message(params2 <- matchYields(params))
    biomass2 <- sweep(params2@initial_n, 2, params@w * params@dw, "*")
    yield_model2 <- rowSums(biomass2 * getFMort(params2))
    expect_equal(yield_model2[4:12], 4:12, check.attributes = FALSE)
    
})
