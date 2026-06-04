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
    expect_error(plotlySpectra2(params, sim, species = species), NA)
    expect_error(plotlyPredMort(sim, species = species), NA)
    expect_error(plotlyFMort(sim, species = species), NA)
    expect_error(plotlyGrowthCurves(sim, species = species), NA)
    expect_error(plotlyGrowthCurves(params, species = species), NA)
    expect_error(plotlyDiet(params, species = species[[1]]), NA)
})

test_that("plotly wrappers return plotly objects for spectra and rate plots", {
    expect_s3_class(plotlySpectra(params, species = species), "plotly")
    expect_s3_class(plotlySpectra2(params, sim, species = species), "plotly")
    expect_s3_class(plotlyCDF(params, species = species,
                              resource = FALSE), "plotly")
    expect_s3_class(plotlyCDF2(params, params, species = species,
                               resource = FALSE), "plotly")
    expect_s3_class(plotlySpectraRelative(params, params, species = species,
                                          resource = FALSE), "plotly")
    expect_s3_class(plotlyPredMort(sim, species = species), "plotly")
    expect_s3_class(plotlyFMort(sim, species = species), "plotly")
    expect_s3_class(plotlyGrowthCurves(sim, species = species), "plotly")
    expect_s3_class(plotlyFeedingLevel(sim, species = species,
                                       include_critical = TRUE), "plotly")
    expect_s3_class(plotlyDiet(params, species = species[[1]]), "plotly")
})

test_that("ggplotly(plot(...)) uses concise mizer tooltips", {
    p <- plot(getEncounter(NS_params), species = "Cod")
    expect_s3_class(p, "mizer_plot")
    gp <- ggplotly(p)
    first_tip <- gp$x$data[[1]]$text[[1]]
    expect_true(grepl("Species: Cod", first_tip, fixed = TRUE))
    expect_true(grepl("w:", first_tip, fixed = TRUE))
    expect_true(grepl("value:", first_tip, fixed = TRUE))
    legend_matches <- gregexpr("Legend:", first_tip, fixed = TRUE)[[1]]
    expect_lte(sum(legend_matches > 0), 1)

    ggplot2::set_last_plot(p)
    gp_last <- ggplotly()
    expect_identical(gp_last$x$data[[1]]$text[[1]], first_tip)
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

test_that("comparison helpers preserve non-size x variables", {
    p <- plot2(getBiomass(NS_sim), getBiomass(NS_sim), species = "Cod")
    expect_s3_class(p, "ggplot")
    expect_true("Year" %in% names(p$data))
    expect_error(ggplot2::ggplot_build(p), NA)

    p_rel <- plotRelative(getBiomass(NS_sim), getBiomass(NS_sim),
                          species = "Cod")
    expect_s3_class(p_rel, "ggplot")
    expect_true("Year" %in% names(p_rel$data))
    expect_error(ggplot2::ggplot_build(p_rel), NA)
})

test_that("plotSpectraRelative plots symmetric relative difference", {
    params2 <- params
    params2@initial_n[] <- params@initial_n * 2

    p <- plotSpectraRelative(params, params2, species = species,
                             resource = FALSE)
    expect_s3_class(p, "ggplot")
    expect_true(all(abs(p$data$rel_diff - 2 / 3) < 1e-12))
    expect_identical(p$scales$get_scales("x")$trans$name, "log-10")

    p_linear <- plotSpectraRelative(params, params2, species = species,
                                    resource = FALSE, log_x = FALSE)
    expect_identical(p_linear$scales$get_scales("x")$trans$name, "identity")

    expect_s3_class(plotSpectraRelative(sim, sim0, species = species,
                                        resource = FALSE), "ggplot")
    expect_s3_class(plotSpectraRelative(params, sim, species = species,
                                        resource = FALSE), "ggplot")
})

test_that("plotCDF plots cumulative spectra from small to large sizes", {
    p <- plotCDF(params, species = species, resource = FALSE, power = 0,
                 wlim = c(1, NA), return_data = TRUE)
    expect_true(all(p$w >= 1))
    expect_true(all(p$value >= 0))
    for (sp in unique(p$Species)) {
        sp_dat <- p[p$Species == sp, ]
        expect_equal(max(sp_dat$value), 1)
        expect_true(all(diff(sp_dat$value) >= -1e-12))
    }

    spectra <- plotSpectra(params, species = species, resource = FALSE,
                           power = 1, wlim = c(1, NA), return_data = TRUE)
    cdf <- plotCDF(params, species = species, resource = FALSE,
                   power = 1, wlim = c(1, NA), normalise = FALSE,
                   return_data = TRUE)
    widths <- params@dw_full[match(spectra$w, params@w_full)]
    expected <- sum(spectra$value[spectra$Species == species[[1]]] *
                        widths[spectra$Species == species[[1]]])
    observed <- max(cdf$value[cdf$Species == species[[1]]])
    expect_equal(observed, expected)

    p_plot <- plotCDF(params, species = species, resource = FALSE)
    expect_s3_class(p_plot, "ggplot")
    expect_identical(p_plot$scales$get_scales("x")$trans$name, "log-10")
    expect_match(p_plot$scales$get_scales("y")$name, "biomass",
                 ignore.case = TRUE)

    p_linear <- plotCDF(params, species = species, resource = FALSE,
                        log_x = FALSE)
    expect_identical(p_linear$scales$get_scales("x")$trans$name, "identity")

    p_log_none <- plotCDF(params, species = species, resource = FALSE,
                          log = "")
    expect_identical(p_log_none$scales$get_scales("x")$trans$name,
                     "identity")

    p_log_x <- plotCDF(params, species = species, resource = FALSE,
                       log = "x")
    expect_identical(p_log_x$scales$get_scales("x")$trans$name, "log-10")

    p_log_y <- plotCDF(params, species = species, resource = FALSE,
                       log_y = TRUE)
    expect_identical(p_log_y$scales$get_scales("y")$trans$name, "log-10")

    p_log_xy <- plotCDF(params, species = species, resource = FALSE,
                        log = "xy")
    expect_identical(p_log_xy$scales$get_scales("x")$trans$name, "log-10")
    expect_identical(p_log_xy$scales$get_scales("y")$trans$name, "log-10")

    p_abundance <- plotCDF(params, species = species, resource = FALSE,
                           power = 0)
    expect_match(p_abundance$scales$get_scales("y")$name, "abundance",
                 ignore.case = TRUE)
})

test_that("plotCDF supports simulations, resource, total and unnormalised output", {
    p <- plotCDF(sim, species = species, time_range = 1:3,
                 total = TRUE, resource = TRUE, normalise = FALSE,
                 return_data = TRUE)
    expect_true(all(c("Resource", "Total") %in% p$Legend))
    expect_true(max(p$value) > 1)

    p_plot <- plotCDF(sim, species = species, time_range = 1:3,
                      total = TRUE, resource = TRUE, normalise = FALSE)
    expect_s3_class(p_plot, "ggplot")
})

test_that("plotCDF2 compares cumulative distributions", {
    p <- plotCDF2(params, params, name1 = "Original", name2 = "Changed",
                  species = species, total = TRUE, resource = FALSE,
                  wlim = c(1, NA), normalise = FALSE, log = "")
    expect_s3_class(p, "ggplot")
    expect_identical(levels(p$data$Model), c("Original", "Changed"))
    expect_true("Total" %in% p$data$Legend)
    expect_true(all(p$data$w >= 1))
    expect_identical(p$scales$get_scales("x")$trans$name, "identity")

    p_log <- plotCDF2(params, sim, species = species, resource = FALSE,
                      log = "x")
    expect_s3_class(p_log, "ggplot")
    expect_identical(p_log$scales$get_scales("x")$trans$name, "log-10")

    expect_s3_class(plotCDF2(sim, sim0, species = species,
                             time_range = 1:3, resource = FALSE),
                    "ggplot")
    p_log_y2 <- plotCDF2(params, sim, species = species, resource = FALSE,
                         log = "y")
    expect_identical(p_log_y2$scales$get_scales("y")$trans$name, "log-10")
})

test_that("size-based plots support length axes", {
    params_len <- params
    params_len@species_params$a <- 0.01
    params_len@species_params$b <- 3
    sim_len <- sim
    sim_len@params <- params_len

    spectra_w <- plotSpectra(params_len, species = species, resource = FALSE,
                             total = FALSE, return_data = TRUE)
    spectra_l <- plotSpectra(params_len, species = species, resource = FALSE,
                             total = FALSE, size_axis = "l",
                             return_data = TRUE)
    expect_true("l" %in% names(spectra_l))
    expect_false("w" %in% names(spectra_l))
    sp_idx <- match(as.character(spectra_w$Species),
                    as.character(params_len@species_params$species))
    expected_l <- w2l(spectra_w$w,
                      params_len@species_params[sp_idx, , drop = FALSE])
    expect_equal(spectra_l$l, expected_l)

    llim <- stats::quantile(spectra_l$l, c(0.25, 0.75), names = FALSE)
    spectra_l_limited <- plotSpectra(params_len, species = species,
                                     resource = FALSE, total = FALSE,
                                     size_axis = "l", llim = llim,
                                     return_data = TRUE)
    expect_true(all(spectra_l_limited$l >= llim[1]))
    expect_true(all(spectra_l_limited$l <= llim[2]))

    spectra_hidden <- plotSpectra(params_len, species = species,
                                  resource = TRUE, total = TRUE,
                                  size_axis = "l", return_data = TRUE)
    expect_false(any(spectra_hidden$Legend %in% c("Resource", "Total")))

    p <- plotSpectra(params_len, species = species, resource = FALSE,
                     size_axis = "l")
    expect_identical(p$scales$get_scales("x")$name, "Length [cm]")
    expect_identical(p$scales$get_scales("x")$trans$name, "log-10")

    expect_true("l" %in% names(plotCDF(params_len, species = species,
                                       resource = FALSE, size_axis = "l",
                                       return_data = TRUE)))
    cdf_l_limited <- plotCDF(params_len, species = species, resource = FALSE,
                             size_axis = "l", llim = llim,
                             return_data = TRUE)
    expect_true(all(cdf_l_limited$l >= llim[1]))
    expect_true(all(cdf_l_limited$l <= llim[2]))
    expect_equal(stats::aggregate(value ~ Species, cdf_l_limited, max)$value,
                 rep(1, length(unique(cdf_l_limited$Species))))
    expect_true("l" %in% names(plotSpectra2(params_len, params_len,
                                            species = species,
                                            resource = FALSE,
                                            size_axis = "l")$data))
    expect_true("l" %in% names(plotCDF2(params_len, params_len,
                                        species = species,
                                        resource = FALSE,
                                        size_axis = "l")$data))
    expect_true("l" %in% names(plotSpectraRelative(params_len, params_len,
                                                   species = species,
                                                   resource = FALSE,
                                                   size_axis = "l")$data))

    expect_true("l" %in% names(plotFeedingLevel(params_len, species = species,
                                                size_axis = "l",
                                                return_data = TRUE)))
    expect_true("l" %in% names(plotPredMort(params_len, species = species,
                                            size_axis = "l",
                                            return_data = TRUE)))
    expect_true("l" %in% names(plotFMort(params_len, species = species,
                                         size_axis = "l",
                                         return_data = TRUE)))
    expect_true("l" %in% names(plotDiet(params_len, species = species[[1]],
                                        size_axis = "l",
                                        return_data = TRUE)))
    expect_true(all(plotDiet(params_len, species = species[[1]],
                             size_axis = "l", llim = llim,
                             return_data = TRUE)$l >= llim[1]))
    expect_true("l" %in% names(plot(getPredMort(params_len),
                                    species = species, size_axis = "l",
                                    return_data = TRUE)))
    rate_l_limited <- plot(getPredMort(params_len), species = species,
                           size_axis = "l", llim = llim,
                           return_data = TRUE)
    expect_true(all(rate_l_limited$l >= llim[1]))
    expect_true(all(rate_l_limited$l <= llim[2]))
    expect_true("l" %in% names(plot(getFMort(sim_len), species = species,
                                    size_axis = "l", return_data = TRUE)))
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

    p <- plotSpectra(sim, species = species, resource = FALSE,
                     size_axis = "l", llim = c(10, 100))
    expect_equal(p$scales$scales[[2]]$limits, c(1, 2))

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
