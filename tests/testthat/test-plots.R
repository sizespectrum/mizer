context("Plotting methods")

# Initialisation ----------------
data(NS_species_params_gears)
data(inter)
params <- newMultispeciesParams(NS_species_params_gears, inter, no_w = 30)
sim <- project(params, effort = 1, t_max = 3, dt = 1, t_save = 1)
sim0 <- project(params, effort = 0, t_max = 3, dt = 1, t_save = 1)
species <- c("Cod", "Haddock")

# plots have not changed ----
test_that("plots have not changed", {
p <- plotBiomass(sim, species = species, total = TRUE,
                 start_time = 0, end_time = 2.8,
                 y_ticks = 4)
p$plot_env <- NULL
expect_known_value(p, "values/plotBiomass")

p <- plotYield(sim, species = species, total = TRUE)
p$plot_env <- NULL
expect_known_value(p, "values/plotYield")

p <- plotYieldGear(sim, species = species)
p$plot_env <- NULL
expect_known_value(p, "values/plotYieldGear")

p <- plotSpectra(sim, species = species, total = TRUE,
                 time_range = 1:3, power = 2,
                 ylim = c(1e6, NA))
p$plot_env <- NULL
expect_known_value(p, "values/plotSpectra")

p <- plotFeedingLevel(sim, species = species, time_range = 1:3)
p$plot_env <- NULL
expect_known_value(p, "values/plotFeedingLevel")

p <- plotPredMort(sim, species = species, time_range = 1:3)
p$plot_env <- NULL
expect_known_value(p, "values/plotPredMort")

p <- plotFMort(sim, species = species, time_range = 1:3)
p$plot_env <- NULL
expect_known_value(p, "values/plotFMort")

p <- plotGrowthCurves(sim, species = species, percentage = TRUE,
                      max_age = 50)
p$plot_env <- NULL
expect_known_value(p, "values/plotGrowthCurves")

sim@params@species_params$a <- 0.0058
sim@params@species_params$b <- 3.13
p <- plotGrowthCurves(sim, species = "Haddock", max_age = 50)
p$plot_env <- NULL
expect_known_value(p, "values/plotGrowthCurvesSingle")

p <- displayFrames(getBiomassFrame(sim0, 
                                   species = species,
                                   start_time = 1,
                                   total = TRUE), 
                   getBiomassFrame(sim,
                                   end_time = 3,
                                   ylim = c(1e12)),
                   params)
p$plot_env <- NULL
expect_known_value(p, "values/displayFrames")

expect_known_value(getSSBFrame(sim, species = "Cod", total = TRUE),
                   "values/getSSBFrame")
})

test_that("plot function do not throw error", {
    expect_error(plot(sim, species = species, wlim = c(10, 100), w_min = 10), NA)
})

# plotly functions do not throw error
test_that("plotly functions do not throw error", {
    expect_error(plotlyBiomass(sim, species = species), NA)
    expect_error(plotlyFeedingLevel(sim, species = species), NA)
    expect_error(plotlyYield(sim, species = species), NA)
    expect_error(plotlyYield(sim, sim), NA)
    expect_error(plotlyYieldGear(sim, species = species), NA)
    expect_error(plotlySpectra(params, species = species), NA)
    expect_error(plotlyPredMort(sim, species = species), NA)
    expect_error(plotlyFMort(sim, species = species), NA)
    expect_error(plotlyGrowthCurves(sim, species = species), NA)
    expect_error(plotlyGrowthCurves(params, species = species), NA)
})

