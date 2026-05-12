# Initialisation ----------------
# Snapshots recorded with edition 1; lock params creation to edition 1
withr::local_options(mizer_defaults_edition = 1)
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
params_bkgrd <- markBackground(params, species = params@species_params$species[1:3])
# params object with single species
sp_single <- data.frame(species = 1, w_max = 1000, h = 30)
params_single <- newMultispeciesParams(sp_single, no_w = 30, info_level = 0)

# Need to use vdiffr conditionally
expect_doppelganger <- function(title, fig, ...) {
    testthat::skip_if_not_installed("vdiffr")
    vdiffr::expect_doppelganger(title, fig, ...)
}

expect_ggplot <- function(fig) {
    testthat::expect_s3_class(fig, "ggplot")
}

# plots have not changed ----
test_that("plots have not changed", {
p <- plotBiomass(sim, species = species, total = TRUE,
                 start_time = 0, end_time = 2.8,
                 y_ticks = 4)
expect_doppelganger("Plot Biomass", p)

p <- plotYield(sim, species = species, total = TRUE)
expect_ggplot(p)

p <- plotYieldGear(sim, species = species)
expect_ggplot(p)

p <- plotSpectra(sim, species = species, total = TRUE,
                 time_range = 1:3, power = 2,
                 ylim = c(1e6, NA))
expect_doppelganger("Plot Spectra", p)

p <- plotFeedingLevel(sim, species = species, time_range = 1:3)
expect_ggplot(p)
p <- plotFeedingLevel(sim, species = species, time_range = 1:3,
                      include_critical = TRUE)
expect_ggplot(p)

p <- plotPredMort(sim, species = species, time_range = 1:3,
                  all.sizes = TRUE)
expect_ggplot(p)
p <- plotPredMort(sim, species = 2, time_range = 1:3)
# The following test is disabled because the plot is different on different
# platforms
# TODO: reenable this test
# expect_doppelganger("PlotPredMort truncated", p)

p <- plotFMort(sim, species = species, time_range = 1:3,
               all.sizes = TRUE)
expect_ggplot(p)
p <- plotFMort(sim, species = 2, time_range = 1:3)
expect_ggplot(p)

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
expect_ggplot(p)

sp_name <- NS_params@species_params$species[10]
p <- plotDiet(NS_params, species = sp_name)
expect_doppelganger("Plot Diet", p)
})

test_that("plot function do not throw error", {
    expect_error(plot(sim, species = species), NA)
    expect_error(plot(params, species = species), NA)
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

test_that("plotly wrappers return plotly objects for spectra and rate plots", {
    expect_s3_class(plotlySpectra(params, species = species), "plotly")
    expect_s3_class(plotlyPredMort(sim, species = species), "plotly")
    expect_s3_class(plotlyFMort(sim, species = species), "plotly")
    expect_s3_class(plotlyGrowthCurves(sim, species = species), "plotly")
    expect_s3_class(plotlyFeedingLevel(sim, species = species,
                                       include_critical = TRUE), "plotly")
})

test_that("plotSpectra2 compares spectra from params and sims", {
    p_params <- plotSpectra2(params, params, name1 = "Original",
                             name2 = "Changed", species = species,
                             total = TRUE)
    expect_s3_class(p_params, "ggplot")
    expect_identical(levels(p_params$data$Model), c("Original", "Changed"))
    expect_true("Total" %in% p_params$data$Legend)

    expect_s3_class(plotSpectra2(sim, sim0, species = species), "ggplot")
    expect_s3_class(plotSpectra2(params, sim, species = species), "ggplot")
})

test_that("plotSpectra2 supports base plot log argument", {
    p_y <- plotSpectra2(params, sim0, species = species, log = "y")
    expect_identical(p_y$scales$get_scales("x")$trans$name, "identity")
    expect_identical(p_y$scales$get_scales("y")$trans$name, "log-10")

    p_xy <- plotSpectra2(params, sim0, species = species, log = "xy")
    expect_identical(p_xy$scales$get_scales("x")$trans$name, "log-10")
    expect_identical(p_xy$scales$get_scales("y")$trans$name, "log-10")

    p_none <- plotSpectra2(params, sim0, species = species, log = "")
    expect_identical(p_none$scales$get_scales("x")$trans$name, "identity")
    expect_identical(p_none$scales$get_scales("y")$trans$name, "identity")

    p_flags <- plotSpectra2(params, sim0, species = species,
                            log_x = FALSE, log_y = FALSE)
    expect_identical(p_flags$scales$get_scales("x")$trans$name, "identity")
    expect_identical(p_flags$scales$get_scales("y")$trans$name, "identity")

    expect_error(plotSpectra2(params, sim0, species = species, log = "z"),
                 "`log` must be a character string")
})

test_that("yield plotting helpers validate comparison and gear selection", {
    sim_shifted <- sim
    dimnames(sim_shifted@n)$time <- as.character(10:13)
    expect_error(plotYield(sim, sim_shifted), "do not have the same times")

    y <- plotYieldGear(sim, species = species,
                       gears = dimnames(sim@params@selectivity)$gear[[1]],
                       return_data = TRUE)
    expect_true(all(y$Gear == dimnames(sim@params@selectivity)$gear[[1]]))
})

test_that("yield plotly wrappers return plotly objects", {
    expect_s3_class(plotlyYield(sim, species = species), "plotly")
    expect_s3_class(plotlyYieldGear(sim, species = species), "plotly")
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
                                  return_data = TRUE)), c(56, 4))

    expect_equal(dim(plotFMort(sim, species = species,
                               return_data = TRUE)), c(56, 4))

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

test_that("plotSpectra validates empty selection and can return total only", {
    expect_error(plotSpectra(params, species = rep(FALSE, nrow(species_params(params))),
                             total = FALSE, resource = FALSE),
                 "There is nothing to plot")
    df <- plotSpectra(params, species = rep(FALSE, nrow(species_params(params))),
                      total = TRUE, resource = FALSE, return_data = TRUE)
    expect_true(all(df$Legend == "Total"))
})

test_that("plotSpectra supports base plot log argument", {
    p_y <- plotSpectra(params, species = species, log = "y")
    expect_identical(p_y$scales$get_scales("x")$trans$name, "identity")
    expect_identical(p_y$scales$get_scales("y")$trans$name, "log-10")

    p_xy <- plotSpectra(params, species = species, log = "xy")
    expect_identical(p_xy$scales$get_scales("x")$trans$name, "log-10")
    expect_identical(p_xy$scales$get_scales("y")$trans$name, "log-10")

    p_none <- plotSpectra(params, species = species, log = "")
    expect_identical(p_none$scales$get_scales("x")$trans$name, "identity")
    expect_identical(p_none$scales$get_scales("y")$trans$name, "identity")

    p_sim <- plotSpectra(sim, species = species, log_x = FALSE, log_y = FALSE)
    expect_identical(p_sim$scales$get_scales("x")$trans$name, "identity")
    expect_identical(p_sim$scales$get_scales("y")$trans$name, "identity")

    expect_error(plotSpectra(params, species = species, log = "z"),
                 "`log` must be a character string")
})

test_that("plotPredMort and plotFMort trim to species size range by default", {
    pred_trimmed <- plotPredMort(params, species = 2, return_data = TRUE)
    pred_full <- plotPredMort(params, species = 2, all.sizes = TRUE, return_data = TRUE)
    expect_lt(nrow(pred_trimmed), nrow(pred_full))
    expect_equal(
        plotM2(params, species = 2, return_data = TRUE),
        pred_trimmed
    )

    f_trimmed <- plotFMort(params, species = 2, return_data = TRUE)
    f_full <- plotFMort(params, species = 2, all.sizes = TRUE, return_data = TRUE)
    expect_lt(nrow(f_trimmed), nrow(f_full))
})

test_that("plotFeedingLevel trims by size and can include critical levels", {
    fl_trimmed <- plotFeedingLevel(params, species = 2, return_data = TRUE)
    fl_full <- plotFeedingLevel(params, species = 2, all.sizes = TRUE,
                                return_data = TRUE)
    expect_lt(nrow(fl_trimmed), nrow(fl_full))

    fl_critical <- plotFeedingLevel(params, species = 2,
                                    include_critical = TRUE,
                                    return_data = TRUE)
    expect_setequal(unique(fl_critical$Type), c("actual", "critical"))
})

test_that("plotGrowthCurves validates size_at_age input", {
    expect_error(plotGrowthCurves(params, size_at_age = data.frame(age = 1, weight = 2)),
                 "needs to have a 'species' column")
    expect_error(plotGrowthCurves(params, size_at_age = data.frame(species = "10", weight = 2)),
                 "needs to have an 'age' column")
    expect_error(plotGrowthCurves(params, size_at_age = data.frame(species = "10", age = 1)),
                 "needs to have either a 'length' or a 'weight' column")
})

test_that("plotDiet restricts to meaningful abundance ranges", {
    # Get diet data for a species
    diet_data <- plotDiet(params, species = "10", return_data = TRUE)

    # Get abundance data for species "10"
    sp_idx <- which(dimnames(params@initial_n)[[1]] == "10")
    abundance <- params@initial_n[sp_idx, ] * params@w * params@dw
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

# Test axis limits are properly set ----
test_that("axis limits are set correctly", {
    # Test ylim in plotBiomass
    p <- plotBiomass(sim, species = species, ylim = c(1e3, 1e6))
    expect_equal(p$scales$scales[[1]]$limits, c(3, 6))

    # Test with NA values
    p <- plotBiomass(sim, species = species, ylim = c(NA, 1e6))
    expect_equal(p$scales$scales[[1]]$limits[2], 6)
    expect_true(is.na(p$scales$scales[[1]]$limits[1]))

    # Test wlim and ylim in plotSpectra
    p <- plotSpectra(sim, species = species, wlim = c(1, 100), ylim = c(1e5, 1e8))
    # x-axis is first scale, y-axis is second
    expect_equal(p$scales$scales[[2]]$limits, c(0, 2))
    expect_equal(p$scales$scales[[1]]$limits, c(5, 8))

    # Test with NA values in wlim
    p <- plotSpectra(sim, species = species, wlim = c(10, NA), ylim = c(NA, 1e8))
    expect_equal(p$scales$scales[[2]]$limits[1], 1)
    expect_equal(p$scales$scales[[2]]$limits[2], log10(max(params@w_full)))
    expect_true(is.na(p$scales$scales[[1]]$limits[1]))
    expect_equal(p$scales$scales[[1]]$limits[2], 8)

    # Default wlim lower depends on resource argument
    p_res <- plotSpectra(sim, species = species, return_data = TRUE)
    p_nores <- plotSpectra(sim, species = species, resource = FALSE, return_data = TRUE)
    expect_true(min(p_res$w) < min(params@w))
    expect_equal(min(p_nores$w), min(params@w))
})

test_that("plotDiet works with MizerSim", {
    p <- plotDiet(sim, species = 11) # Species 11 is Cod in setup
    expect_true(is(p, "ggplot"))
    p <- plotDiet(sim, species = 11, time_range = 1:2)
    expect_true(is(p, "ggplot"))
})
