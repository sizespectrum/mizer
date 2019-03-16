context("Plotting methods")

# Initialisation ----------------
data(NS_species_params_gears)
data(inter)
params <- MizerParams(NS_species_params_gears, inter, no_w = 30)
sim <- project(params, effort = 1, t_max = 3, dt = 1, t_save = 1)
species <- c("Cod", "Haddock")

test_that("plots have not changed", {
expect_known_value(plotBiomass(sim, species = species, total = TRUE,
                               start_time = 0, end_time = 2.8,
                               y_ticks = 4),
                   "values/plotBiomass")
expect_known_value(plotYield(sim, species = species, total = TRUE),
                   "values/plotYield")
expect_known_value(plotYieldGear(sim, species = species), "values/plotYieldGear")
expect_known_value(plotSpectra(sim, species = species, total = TRUE,
                               time_range = 1:3, power = 2,
                               ylim = c(1e6, NA)),
                   "values/plotSpectra")
expect_known_value(plotFeedingLevel(sim, species = species, time_range = 1:3),
                   "values/plotFeedingLevel")
expect_known_value(plotM2(sim, species = species, time_range = 1:3),
                   "values/plotM2")
expect_known_value(plotFMort(sim, species = species, time_range = 1:3),
                   "values/plotFMort")
expect_known_value(plotGrowthCurves(sim, species = species, percentage = TRUE,
                                    max_age = 50),
                   "values/plotGrowthCurves")
})

test_that("plotly function do not throw error", {
    expect_error(plotlyBiomass(sim, species = species), NA)
    expect_error(plotlyFeedingLevel(sim, species = species), NA)
    expect_error(plotlyYield(sim, species = species), NA)
    expect_error(plotlyYield(sim, sim), NA)
    expect_error(plotlyYieldGear(sim, species = species), NA)
    expect_error(plotlySpectra(params, species = species), NA)
    expect_error(plotlyM2(sim, species = species), NA)
    expect_error(plotlyFMort(sim, species = species), NA)
    expect_error(plotlyGrowthCurves(sim, species = species), NA)
    expect_error(plotlyGrowthCurves(params, species = species), NA)
})

