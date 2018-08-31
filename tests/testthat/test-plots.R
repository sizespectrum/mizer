context("Plotting methods")

# Initialisation ----------------
data(NS_species_params_gears)
data(inter)
params <- MizerParams(NS_species_params_gears, inter, no_w = 30)
sim <- project(params, effort = 1, t_max = 2, dt = 1, t_save = 1)

test_that("plotting function do not throw error", {
    expect_error(plot(sim), NA)
    expect_error(plotYield(sim), NA)
    expect_error(plotYield(sim, sim), NA)
    expect_error(plotYieldGear(sim), NA)
    expect_error(plotSpectra(params), NA)
    expect_error(plotFMort(sim), NA)
    expect_error(plotGrowthCurves(sim), NA)
    expect_error(plotGrowthCurves(params), NA)
})

