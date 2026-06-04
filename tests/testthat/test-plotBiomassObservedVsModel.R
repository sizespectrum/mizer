test_that("plotBiomassObservedVsModel works", {
# Set up parameters
params <- NS_params

# check you get error without biomass_observed column
expect_error(plotBiomassObservedVsModel(params))

# pull out data frame for biomass comparison
species_params(params)$biomass_observed <-
  c(0.8, 61, 12, 35, 1.6, 20, 10, 7.6, 135, 60, 30, 78)
species_params(params)$biomass_cutoff <- 10
params <- calibrateBiomass(params)
dummy <- plotBiomassObservedVsModel(params, return_data = TRUE)

# Check biomasses equal those put in
expect_equal(dummy$observed, species_params(params)$biomass_observed)

# check that you get error with no species
expect_error(plotBiomassObservedVsModel(params, species = rep(FALSE, 12)))

# Try removing observed biomasses.
params2 <- params # copy over
species_params(params2)$biomass_observed[c(1, 7, 10)] <- NA
# plot without unobserved species
dummy <- plotBiomassObservedVsModel(params2, return_data = TRUE)
expect_equal(as.character(dummy$species), 
             species_params(params)$species[!is.na(species_params(params2)$biomass_observed)])
expect_equal(dummy$observed, 
             species_params(params2)$biomass_observed
             [!is.na(species_params(params2)$biomass_observed)])
# plot with unobserved species
dummy <- plotBiomassObservedVsModel(params2, return_data = TRUE,
                                   show_unobserved = TRUE)
expect_equal(as.character(dummy$species), 
             species_params(params)$species)

# # Try removing species, check it still checks out
sp_select <- c(1, 4, 7, 10, 11, 12) # choose some species
dummy <- plotBiomassObservedVsModel(params, species = sp_select, return_data = T)
expect_equal(nrow(dummy), length(sp_select))
expect_equal(dummy$observed, 
             species_params(params)$biomass_observed[sp_select])

# Finally, look at default plot (model biomass vs observed biomass)
dummy <- plotBiomassObservedVsModel(params, return_data = TRUE)
p <- plotBiomassObservedVsModel(params)
expect_true(is_ggplot(p))
expect_identical(p$labels$x, "observed biomass [g]")
expect_identical(p$labels$y, "model biomass [g]")
expect_identical(p$data, dummy)

# Look at plot of ratio
dummy <- plotBiomassObservedVsModel(params, ratio = TRUE, return_data = TRUE)
p <- plotBiomassObservedVsModel(params, ratio = TRUE)
expect_identical(p$labels$y, "model biomass / observed biomass")
expect_identical(p$data, dummy)

vdiffr::expect_doppelganger("plotBiomassObservedVsModel", p)

})

test_that("plotBiomassObservedVsModel methods for MizerSim and plotly work", {
params <- NS_params
species_params(params)$biomass_observed <-
  c(0.8, 61, 12, 35, 1.6, 20, 10, 7.6, 135, 60, 30, 78)
species_params(params)$biomass_cutoff <- 10
params <- calibrateBiomass(params)
sim <- project(params, t_max = 0.1, progress_bar = FALSE)

dummy_sim <- plotBiomassObservedVsModel(sim, return_data = TRUE)
dummy_params <- plotBiomassObservedVsModel(setInitialValues(sim@params, sim),
                                           return_data = TRUE,
                                           ratio = FALSE)
expect_identical(dummy_sim, dummy_params)

    p <- plotBiomassObservedVsModel(sim)
    expect_true(is_ggplot(p))
    expect_identical(p$labels$y, "model biomass [g]")

    pp <- plotlyBiomassObservedVsModel(params)
    expect_s3_class(pp, "plotly")
})
