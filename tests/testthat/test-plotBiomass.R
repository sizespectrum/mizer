# Plot fixtures use edition 1, as do the existing visual baselines.
withr::local_options(mizer_defaults_edition = 1)
params <- plot_params_small
sim <- plot_sim_small
species <- c(1, 2)

test_that("plotBiomass works with use_cutoff", {
    params <- NS_params_small
    species_params(params)$biomass_cutoff <- 10
    sim <- project(params, t_max = 1, effort = 1)

    sp_name <- species_params(params)$species[3]
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
    sim <- project(NS_params_small, t_max = 1, effort = 1, progress_bar = FALSE)
    expect_error(plotBiomass(sim, tlim = c(1,1)),
                 "must be less than")

    d <- plotBiomass(sim, species = "Cod", total = TRUE, return_data = TRUE)
    expect_true("Total" %in% d$Species)
})

test_that("biomass plot appearance is unchanged", {
    p <- plotBiomass(sim, species = species, total = TRUE,
                     tlim = c(0, 2.8),
                     y_ticks = 4)
    expect_doppelganger("Plot Biomass", p)
})

test_that("interactive wrappers preserve tooltip fields", {
    gp <- plotlyBiomass(sim, species = species)
    expect_s3_class(gp, "plotly")
    expect_equal(tooltip_fields(gp), c("Species", "Year", "Biomass"))
})

test_that("public return-data shape is preserved", {
    expect_equal(dim(plotBiomass(sim, species = species, total = TRUE,
                                 tlim = c(0, 2.8), y_ticks = 4, return_data = TRUE)), c(12,4))
})

test_that("time limits only set the biomass viewport", {
    # plotBiomass: deprecated start_time/end_time triggers warning
    expect_warning(plotBiomass(sim, start_time = 0), "deprecated")
    expect_warning(plotBiomass(sim, end_time = 2), "deprecated")

    # plotBiomass: tlim retains all rows
    all_rows <- nrow(plotBiomass(sim, species = species, return_data = TRUE))
    limited_rows <- nrow(plotBiomass(sim, species = species,
                                     tlim = c(1, 3), return_data = TRUE))
    expect_equal(limited_rows, all_rows)
})

test_that("biomass viewport limits are set correctly", {
    # Test ylim in plotBiomass
    p <- plotBiomass(sim, species = species, ylim = c(1e3, 1e6))
    expect_equal(p$coordinates$limits$y, c(1e3, 1e6))

    # Test with NA values
    p <- plotBiomass(sim, species = species, ylim = c(NA, 1e6))
    expect_equal(p$coordinates$limits$y[2], 1e6)
    expect_true(is.na(p$coordinates$limits$y[1]))
})

test_that("time viewports preserve returned data and totals", {
    for (fun in list(plotBiomass)) {
        full <- fun(NS_sim_small, total = TRUE, return_data = TRUE)
        zoom <- fun(NS_sim_small, total = TRUE, tlim = c(0.5, 2.5),
                    ylim = c(1, 10), return_data = TRUE)
        expect_equal(zoom, full)
    }
})

test_that("zooming the time axis fits the y axis to the visible years", {
    # The y axis is logarithmic.
    y_range <- function(p) {
        10^ggplot2::ggplot_build(p)$layout$panel_params[[1]]$y.range
    }
    bm <- plotBiomass(sim, tlim = c(2, 3))
    years <- bm$data$Year >= 2 & bm$data$Year <= 3
    expect_lt(diff(log(y_range(bm))), diff(log(y_range(plotBiomass(sim)))))
    expect_lte(y_range(bm)[1], min(bm$data[years, 2]))
    expect_gte(y_range(bm)[2], max(bm$data[years, 2]))
})
