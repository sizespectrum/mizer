devtools::load_all()
test_that("plotBiomass works with use_cutoff", {
    params <- NS_params
    species_params(params)$biomass_cutoff <- 10
    sim <- project(params, t_max = 1, effort = 1)
    
    # Test with return_data = TRUE to check values
    # Default behavior (use_cutoff = FALSE)
    p_default <- plotBiomass(sim, return_data = TRUE)
    bm_default <- getBiomass(sim)
    # Check total for a species matches
    expect_equal(p_default$Biomass[p_default$Species == "Cod" & p_default$Year == 1],
                 bm_default["1", "Cod"], ignore_attr = TRUE)
    
    # With use_cutoff = TRUE
    p_cutoff <- plotBiomass(sim, use_cutoff = TRUE, return_data = TRUE)
    bm_cutoff <- getBiomass(sim, use_cutoff = TRUE)
    expect_equal(p_cutoff$Biomass[p_cutoff$Species == "Cod" & p_cutoff$Year == 1],
                 bm_cutoff["1", "Cod"], ignore_attr = TRUE)
    
    # Check that values are different (since cutoff is 10g)
    expect_true(p_default$Biomass[p_default$Species == "Cod" & p_default$Year == 1] >
                p_cutoff$Biomass[p_cutoff$Species == "Cod" & p_cutoff$Year == 1])
    
    # Test plotlyBiomass accepts the argument
    expect_error(plotlyBiomass(sim, use_cutoff = TRUE), NA)
})
