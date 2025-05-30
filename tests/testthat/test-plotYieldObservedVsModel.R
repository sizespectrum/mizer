test_that("plotYieldObservedVsModel works", {
    
  # Set up parameters
    params <- NS_params
    
    # check you get error without yield_observed column
    expect_message(
        expect_error(plotYieldObservedVsModel(params),
                     "You have not provided values"))
    
    # pull out data frame for yield comparison
    species_params(params)$yield_observed <-
        c(0.8, 61, 12, 35, 1.6, 20, 10, 7.6, 135, 60, 30, 78)
    expect_warning(params <- calibrateYield(params))
    expect_message(dummy <- plotYieldObservedVsModel(params, return_data = T))
    
    # Check yields equal those put in
    expect_equal(dummy$observed, species_params(params)$yield_observed[4:12])
    
    # check that you get error with no species
    expect_error(plotYieldObservedVsModel(params, species = rep(F, 12)),
                 "No species selected, please fix.")
    
    # Try removing observed yields.
    params2 <- params
    species_params(params2)$yield_observed[c(1, 7, 10)] <- NA
    # plot without unobserved species
    expect_message(dummy <- plotYieldObservedVsModel(params2, return_data = T))
    expect_equal(as.character(dummy$species), 
                 species_params(params)$species[c(4, 5, 6, 8, 9, 11, 12)])
    expect_equal(dummy$observed, 
                 species_params(params2)$yield_observed[c(4, 5, 6, 8, 9, 11, 12)])
    # plot with unobserved species
    expect_message(dummy <- plotYieldObservedVsModel(params2, return_data = T,
                                                    show_unobserved = TRUE))
    expect_equal(as.character(dummy$species), 
                 species_params(params)$species[4:12])
    
    # # Try removing species, check it still checks out
    sp_select <- c(4, 7, 10, 11, 12) # choose some species
    dummy <- plotYieldObservedVsModel(params, species = sp_select, 
                                      return_data = T)
    expect_equal(nrow(dummy), length(sp_select))
    expect_equal(dummy$observed, 
                 species_params(params)$yield_observed[sp_select])
    
    # Finally, look at plot
    expect_message(p <- plotYieldObservedVsModel(params))
    expect_true(is_ggplot(p))
    expect_identical(p$labels$x, "observed yield [g/year]")
    expect_identical(p$labels$y, "model yield [g/year]")
    
    vdiffr::expect_doppelganger("plotYieldObservedVsModel", p)
})
