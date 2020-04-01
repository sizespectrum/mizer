context("Plotting methods")

# Initialisation ----------------
species_params <- NS_species_params_gears
species_params$pred_kernel_type <- "truncated_lognormal"
params <- newMultispeciesParams(species_params, inter, no_w = 30,
                                n = 2/3, p = 0.7, lambda = 2.8 - 2/3)
sim <- project(params, effort = 1, t_max = 3, dt = 1, t_save = 1)
sim0 <- project(params, effort = 0, t_max = 3, dt = 1, t_save = 1)
species <- c("Cod", "Haddock")

# plots have not changed ----
test_that("plots have not changed", {
p <- plotBiomass(sim, species = species, total = TRUE,
                 start_time = 0, end_time = 2.8,
                 y_ticks = 4)
vdiffr::expect_doppelganger("Plot Biomass", p)

p <- plotYield(sim, species = species, total = TRUE)
vdiffr::expect_doppelganger("Plot Yield", p)

p <- plotYieldGear(sim, species = species)
vdiffr::expect_doppelganger("Plot Yield by Gear", p)

p <- plotSpectra(sim, species = species, total = TRUE,
                 time_range = 1:3, power = 2,
                 ylim = c(1e6, NA))
vdiffr::expect_doppelganger("Plot Spectra", p)

p <- plotFeedingLevel(sim, species = species, time_range = 1:3)
vdiffr::expect_doppelganger("Plot Feeding Level", p)

p <- plotPredMort(sim, species = species, time_range = 1:3)
vdiffr::expect_doppelganger("PlotPredation Mortality", p)

p <- plotFMort(sim, species = species, time_range = 1:3)
vdiffr::expect_doppelganger("PlotFishing Mortality", p)

p <- plotGrowthCurves(sim, species = species, percentage = TRUE,
                      max_age = 50)
vdiffr::expect_doppelganger("Plot Growth Curves", p)

sim@params@species_params[["a"]] <- 0.0058
sim@params@species_params[["b"]] <- 3.13
p <- plotGrowthCurves(sim, species = "Haddock", max_age = 50)
vdiffr::expect_doppelganger("Plot Single Growth Curve", p)
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

