test_that("plotYieldObservedVsModel works", {

  # Set up parameters
    params <- NS_params_small
    # Note: NS_params_small has Industrial gear effort=0, so Sprat is not fished.
    # plotYieldObservedVsModel will exclude Sprat with a message.

    # check you get error without yield_observed column;
    # "Sprat not fished" message is emitted before the error
    expect_message(
        expect_error(plotYieldObservedVsModel(params),
                     "You have not provided values"))

    # pull out data frame for yield comparison
    # Sprat (species 1) is not fished so will be excluded; set observed to NA
    species_params(params)$yield_observed <-
        c(NA, 61, 12)
    expect_warning(params <- calibrateYield(params))
    expect_message(dummy <- plotYieldObservedVsModel(params, return_data = T))

    # Check yields equal those put in for fished species (Herring and Cod)
    expect_equal(dummy$observed, species_params(params)$yield_observed[2:3])

    # check that you get error with no species
    expect_error(plotYieldObservedVsModel(params, species = rep(F, 3)),
                 "No species selected, please fix.")

    # Try removing observed yield for Herring (species 2)
    params2 <- params
    species_params(params2)$yield_observed[2] <- NA
    # plot without unobserved species (only Cod shown)
    expect_message(dummy <- plotYieldObservedVsModel(params2, return_data = T))
    expect_equal(as.character(dummy$species),
                 species_params(params)$species[3])
    expect_equal(dummy$observed,
                 species_params(params2)$yield_observed[3])
    # plot with unobserved species (Herring + Cod shown; Sprat still excluded as unfished)
    expect_message(dummy <- plotYieldObservedVsModel(params2, return_data = T,
                                                    show_unobserved = TRUE))
    expect_equal(as.character(dummy$species),
                 species_params(params)$species[2:3])

    # Try selecting species, check it still checks out
    sp_select <- c(2, 3) # Herring and Cod (fished, have yield_observed)
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

test_that("plotYieldObservedVsModel methods for MizerSim and plotly work", {
    params <- NS_params_small
    species_params(params)$yield_observed <-
        c(NA, 61, 12)
    expect_warning(params <- calibrateYield(params))
    sim <- project(params, t_max = 0.1, progress_bar = FALSE)

    expect_message(dummy_sim <- plotYieldObservedVsModel(sim, return_data = TRUE))
    expect_message(dummy_params <- plotYieldObservedVsModel(finalParams(sim),
                                                            return_data = TRUE))
    expect_equal(dummy_sim, dummy_params, ignore_attr = TRUE)

    expect_message(p <- plotYieldObservedVsModel(sim, ratio = TRUE))
    expect_true(is_ggplot(p))
    expect_identical(p$labels$y, "model yield / observed yield")

    expect_message(pp <- plotlyYieldObservedVsModel(params))
    expect_s3_class(pp, "plotly")
})
