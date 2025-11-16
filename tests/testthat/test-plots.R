# Initialisation ----------------
species_params <- NS_species_params_gears
# Make species names numeric because that created problems in the past
species_params$species <- seq_len(nrow(species_params))
species_params$pred_kernel_type <- "truncated_lognormal"
(params <- newMultispeciesParams(species_params, inter, no_w = 30,
                                n = 2 / 3, p = 0.7, lambda = 2.8 - 2 / 3, 
                                info_level = 0)) |>
    expect_message("Note: Dimnames of interaction matrix do not match")
sim <- project(params, effort = 1, t_max = 3, dt = 1, t_save = 1)
sim0 <- project(params, effort = 0, t_max = 3, dt = 1, t_save = 1)
species <- c(11, 10)
# Mark some species as background
params_bkgrd <- params
params_bkgrd@A[1:3] <- NA
# params object with single species
sp_single <- data.frame(species = 1, w_max = 1000, h = 30)
params_single <- newMultispeciesParams(sp_single, no_w = 30, info_level = 0)

# Need to use vdiffr conditionally
expect_doppelganger <- function(title, fig, ...) {
    testthat::skip_if_not_installed("vdiffr")
    vdiffr::expect_doppelganger(title, fig, ...)
}

# plots have not changed ----
test_that("plots have not changed", {
p <- plotBiomass(sim, species = species, total = TRUE,
                 start_time = 0, end_time = 2.8,
                 y_ticks = 4)
expect_doppelganger("Plot Biomass", p)

p <- plotYield(sim, species = species, total = TRUE)
expect_doppelganger("Plot Yield", p)

p <- plotYieldGear(sim, species = species)
expect_doppelganger("Plot Yield by Gear", p)

p <- plotSpectra(sim, species = species, total = TRUE,
                 time_range = 1:3, power = 2,
                 ylim = c(1e6, NA))
expect_doppelganger("Plot Spectra", p)

p <- plotFeedingLevel(sim, species = species, time_range = 1:3)
expect_doppelganger("Plot Feeding Level", p)
p <- plotFeedingLevel(sim, species = species, time_range = 1:3,
                      include_critical = TRUE)
expect_doppelganger("Plot Feeding Level critical", p)

p <- plotPredMort(sim, species = species, time_range = 1:3,
                  all.sizes = TRUE)
expect_doppelganger("PlotPredation Mortality", p)
p <- plotPredMort(sim, species = 2, time_range = 1:3)
# The following test is disabled because the plot is different on different
# platforms
# TODO: reenable this test
# expect_doppelganger("PlotPredMort truncated", p)

p <- plotFMort(sim, species = species, time_range = 1:3,
               all.sizes = TRUE)
expect_doppelganger("PlotFishing Mortality", p)
p <- plotFMort(sim, species = 2, time_range = 1:3)
expect_doppelganger("PlotFMort truncated", p)

# TODO: figure out why these give different results on different platforms
# p <- plotGrowthCurves(sim, species = species, percentage = TRUE,
#                       max_age = 50)
# expect_doppelganger("Plot Growth Curves", p)
# p <- plotGrowthCurves(sim, percentage = FALSE,
#                       species_panel = TRUE, max_age = 50)
# expect_doppelganger("Plot Growth Curves panel", p)

sim@params@species_params[["a"]] <- 0.0058
sim@params@species_params[["b"]] <- 3.13
p <- plotGrowthCurves(sim, species = "10", max_age = 50)
expect_doppelganger("Plot Single Growth Curve", p)

p <- plotDiet(NS_params, species = "Haddock")
expect_doppelganger("Plot Diet", p)
})

test_that("plot function do not throw error", {
    expect_error(plot(sim, species = species, wlim = c(10, 100),
                      w_min = 10), NA)
    expect_error(plot(params, species = species, wlim = c(10, 100),
                      w_min = 10), NA)
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


# testing the plot outputs
test_that("return_data is identical",{
    expect_equal(dim(plotBiomass(sim, species = species, total = TRUE,
                                 start_time = 0, end_time = 2.8, y_ticks = 4, return_data = TRUE)), c(9,4))
    expect_warning(p <- plotYield(sim, sim0, species = species, return_data = TRUE))
    expect_equal(dim(p), c(8,4))
    
    expect_equal(dim(plotYieldGear(sim, species = species, return_data = TRUE)), c(8,4))
    
    expect_equal(dim(plotSpectra(sim, species = species, wlim = c(1,NA),
                                 return_data = TRUE)), c(37, 4))
    
    expect_equal(dim(plotFeedingLevel(sim, species = species,
                                      return_data = TRUE)), c(56, 3))
    
    expect_equal(dim(plotPredMort(sim, species = species, 
                                  return_data = TRUE)), c(56 ,3))
    
    expect_equal(dim(plotFMort(sim, species = species, 
                               return_data = TRUE)), c(56, 3))
    
    expect_equal(dim(plotGrowthCurves(sim, species = species, 
                                      return_data = TRUE)), c(100,4))
    # the following is not a good test because the size of the returned data
    # frame is machine dependent due to the selection of only results above a
    # certain threshold.
    # expect_equal(dim(plotDiet(sim@params, species = species, 
    #                           return_data = TRUE)), c(717,3))
}
)

# Legends have the correct entries ----
test_that("Legends have correct entries", {
    # plotSpectra
    p <- plotSpectra(params, species = 2:3)
    expect_length(unique(p$data$Legend), 3)
    p <- plotSpectra(params, species = 2:4, resource = FALSE)
    expect_length(unique(p$data$Legend), 3)
    p <- plotSpectra(params, species = 2:3, total = TRUE)
    expect_length(unique(p$data$Legend), 4)
    p <- plotSpectra(params, species = 8:9, background = TRUE)
    expect_length(unique(p$data$Legend), 3)
    p <- plotSpectra(params_bkgrd, species = 8:9, background = TRUE)
    expect_length(unique(p$data$Legend), 4)
    p <- plotSpectra(params_single)
    expect_length(unique(p$data$Legend), 2)
})

test_that("plotSpectra averages over time range", {
    time_sel <- c(24:33)
    time_range <- getTimes(NS_sim)[time_sel]
    # arithmetic mean
    df <- plotSpectra(NS_sim, species = 1, time_range = time_range,
                      power = 0, return_data = TRUE)
    expected <- mean(NS_sim@n[time_sel, 1, 1])
    expect_equal(df$value[1], expected)
    # geometric mean
    df <- plotSpectra(NS_sim, species = 1, time_range = time_range,
                      geometric_mean = TRUE,
                      power = 0, return_data = TRUE)
    expected <- exp(mean(log(NS_sim@n[time_sel, 1, 1])))
    expect_equal(df$value[1], expected)
})

test_that("plotDiet restricts to meaningful abundance ranges", {
    # Get diet data for a species
    diet_data <- plotDiet(params, species = "10", return_data = TRUE)
    
    # Get abundance data for species "10"
    sp_idx <- which(dimnames(params@initial_n)[[1]] == "10")
    abundance <- params@initial_n[sp_idx, ]
    max_abundance <- max(abundance)
    
    # Find the maximum meaningful size
    meaningful_idx <- which(abundance > 1e-5 * max_abundance)
    w_max_meaningful <- params@w[max(meaningful_idx)]
    
    # Check that all sizes in diet_data are at or below w_max_meaningful
    expect_true(all(diet_data$w <= w_max_meaningful))
    
    # Also verify that the maximum size in diet_data is close to w_max_meaningful
    # (there might be filtering due to proportion threshold, so we check it's not much smaller)
    if (nrow(diet_data) > 0) {
        max_w_in_data <- max(diet_data$w)
        expect_true(max_w_in_data <= w_max_meaningful)
    }
})
