test_that("plotBiomass works with use_cutoff", {
    params <- NS_params
    species_params(params)$biomass_cutoff <- 10
    sim <- project(params, t_max = 1, effort = 1)
    
    sp_name <- species_params(params)$species[11]
    # Test with return_data = TRUE to check values
    # Default behavior (use_cutoff = FALSE)
    p_default <- plotBiomass(sim, return_data = TRUE)
    bm_default <- getBiomass(sim)
    # Check total for a species matches
    expect_equal(p_default$Biomass[p_default$Species == sp_name & p_default$Year == 1],
                 bm_default["1", sp_name], ignore_attr = TRUE)
    
    # With use_cutoff = TRUE
    p_cutoff <- plotBiomass(sim, use_cutoff = TRUE, return_data = TRUE)
    bm_cutoff <- getBiomass(sim, use_cutoff = TRUE)
    expect_equal(p_cutoff$Biomass[p_cutoff$Species == sp_name & p_cutoff$Year == 1],
                 bm_cutoff["1", sp_name], ignore_attr = TRUE)
    
    # Check that values are different (since cutoff is 10g)
    expect_true(p_default$Biomass[p_default$Species == sp_name & p_default$Year == 1] >
                p_cutoff$Biomass[p_cutoff$Species == sp_name & p_cutoff$Year == 1])
    
    # Test plotlyBiomass accepts the argument
    expect_error(plotlyBiomass(sim, use_cutoff = TRUE), NA)
})

test_that("plotBiomass validates time range and can include total", {
    sim <- project(NS_params, t_max = 1, effort = 1, progress_bar = FALSE)
    expect_error(plotBiomass(sim, start_time = 1, end_time = 1),
                 "must be less than")

    d <- plotBiomass(sim, species = "Cod", total = TRUE, return_data = TRUE)
    expect_true("Total" %in% d$Species)
})
