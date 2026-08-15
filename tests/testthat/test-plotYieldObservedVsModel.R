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

    # pull out data frame for yield comparison. The observed yields are set
    # close to the model yields so that the plot is legible.
    # Sprat (species 1) is not fished so will be excluded; set observed to NA
    gear_params(params)$yield_observed <- c(NA, 1.2, 0.8) * getYield(params)
    expect_message(dummy <- plotYieldObservedVsModel(params, return_data = T))

    # Check yields equal those put in for fished species (Herring and Cod)
    expect_equal(dummy$observed, gear_params(params)$yield_observed[2:3],
                 ignore_attr = TRUE)

    # check that you get error with no species
    expect_error(plotYieldObservedVsModel(params, species = rep(F, 3)),
                 "No species selected, please fix.")

    # Try removing observed yield for Herring (species 2)
    params2 <- params
    gear_params(params2)$yield_observed[2] <- NA
    # plot without unobserved species (only Cod shown)
    expect_message(dummy <- plotYieldObservedVsModel(params2, return_data = T))
    expect_equal(as.character(dummy$species),
                 species_params(params)$species[3])
    expect_equal(dummy$observed,
                 gear_params(params2)$yield_observed[3],
                 ignore_attr = TRUE)
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
                 gear_params(params)$yield_observed[sp_select],
                 ignore_attr = TRUE)

    # The observed yields are also accepted in the species parameters
    params3 <- NS_params_small
    species_params(params3)$yield_observed <-
        gear_params(params)$yield_observed
    expect_message(dummy3 <- plotYieldObservedVsModel(params3, return_data = T))
    expect_equal(dummy3, dummy)

    # Finally, look at plot
    expect_message(p <- plotYieldObservedVsModel(params))
    expect_true(is_ggplot(p))
    expect_identical(p$labels$x, "observed yield [g/year]")
    expect_identical(p$labels$y, "model yield [g/year]")

    vdiffr::expect_doppelganger("plotYieldObservedVsModel", p)
})

test_that("plotYieldObservedVsModel methods for MizerSim and plotly work", {
    params <- NS_params_small
    gear_params(params)$yield_observed <- c(NA, 1.2, 0.8) * getYield(params)
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

test_that("get_yield_observed sums over gears", {
    params <- NS_params_small
    # No observations anywhere
    expect_null(get_yield_observed(params))

    # Give Cod a second gear
    gp <- gear_params(params)
    gp2 <- rbind(gp, gp["Cod, Otter", ])
    gp2$gear[[4]] <- "Extra"
    gear_params(params) <- gp2

    gear_params(params)$yield_observed <- c(NA, 5, 3, 4)
    expect_equal(get_yield_observed(params),
                 c(Sprat = NA, Herring = 5, Cod = 7))

    # A species with observations in the gear parameters ignores the species
    # parameter, a species without one uses it.
    species_params(params)$yield_observed <- c(1, 2, 3)
    expect_equal(get_yield_observed(params),
                 c(Sprat = 1, Herring = 5, Cod = 7))
})

test_that("get_yield_observed reads the species parameters", {
    params <- NS_params_small
    species_params(params)$yield_observed <- c(1, NA, 3)
    expect_equal(get_yield_observed(params),
                 c(Sprat = 1, Herring = NA, Cod = 3))

    species_params(params)$yield_observed <- c("a", "b", "c")
    expect_error(get_yield_observed(params),
                 "The column 'yield_observed' in the species parameter data frame is not numeric")
})

test_that("get_yield_observed rejects non-numeric gear observations", {
    params <- NS_params_small
    gear_params(params)$yield_observed <- c("a", "b", "c")
    expect_error(get_yield_observed(params),
                 "The column 'yield_observed' in the gear parameter data frame is not numeric")
})

test_that("plotYieldObservedVsModel agrees with getYield in both modes", {
    check <- function(params) {
        # Sprat is not fished in this fixture, so it is dropped from the plot
        gear_params(params)$yield_observed <- c(NA, 1.2, 0.8) * getYield(params)
        expect_message(dummy <- plotYieldObservedVsModel(params,
                                                        return_data = TRUE))
        expect_equal(dummy$model, unname(getYield(params)[c(2, 3)]))
    }
    p <- NS_params_small
    check(p)
    second_order_w(p) <- c(bin_average = TRUE)
    check(p)
})
