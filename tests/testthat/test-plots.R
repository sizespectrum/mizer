# Initialisation ----------------
# Snapshots recorded with edition 1; lock params creation to edition 1
withr::local_options(mizer_defaults_edition = 1)
species_params <- NS_species_params_gears_small
# Make species names numeric because that created problems in the past
species_params$species <- seq_len(nrow(species_params))
species_params$pred_kernel_type <- "truncated_lognormal"
(params <- newMultispeciesParams(species_params, inter_small, no_w = 30,
                                n = 2 / 3, p = 0.7, lambda = 2.8 - 2 / 3,
                                info_level = 0)) |>
    expect_message("Note: Dimnames of interaction matrix do not match")
sim <- project(params, effort = 1, t_max = 3, dt = 1, t_save = 1)
sim0 <- project(params, effort = 0, t_max = 3, dt = 1, t_save = 1)
species <- c(1, 2)
# Mark some species as background
params_bkgrd <- markBackground(params, species = params@species_params$species[1])
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
                 tlim = c(0, 2.8),
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
p <- plotGrowthCurves(sim, species = 3, max_age = 50)
expect_ggplot(p)

sp_name <- NS_params_small@species_params$species[2]
p <- plotDiet(NS_params_small, species = sp_name)
expect_doppelganger("Plot Diet", p)
})

tooltip_fields <- function(gp, trace = 1) {
    tip <- gp$x$data[[trace]]$text[[1]]
    trimws(sub(":.*", "", strsplit(tip, "<br />", fixed = TRUE)[[1]]))
}

test_that("plotly wrappers return plotly objects with correct tooltips", {
    gp <- plotlyBiomass(sim, species = species)
    expect_s3_class(gp, "plotly")
    expect_equal(tooltip_fields(gp), c("Species", "Year", "Biomass"))

    gp <- plotlyYield(sim, species = species)
    expect_s3_class(gp, "plotly")
    expect_equal(tooltip_fields(gp), c("Species", "Year", "Yield"))

    expect_s3_class(plotlyYield(sim, sim), "plotly")

    gp <- plotlyYieldGear(sim, species = species)
    expect_s3_class(gp, "plotly")
    expect_equal(tooltip_fields(gp), c("Year", "Yield", "Species", "Gear"))

    gp <- plotlySpectra(params, species = species)
    expect_s3_class(gp, "plotly")
    expect_equal(tooltip_fields(gp), c("Species", "w", "Biomass density"))

    gp <- plotlySpectra2(params, sim, species = species)
    expect_s3_class(gp, "plotly")
    expect_equal(tooltip_fields(gp), c("Legend", "w", "Biomass density", "Model"))

    gp <- plotlyCDF(params, species = species, resource = FALSE)
    expect_s3_class(gp, "plotly")
    expect_equal(tooltip_fields(gp), c("Species", "w", "Cumulative proportion of biomass"))

    gp <- plotlyCDF2(params, params, species = species, resource = FALSE)
    expect_s3_class(gp, "plotly")
    expect_equal(tooltip_fields(gp), c("Legend", "w", "Cumulative proportion of biomass", "Model"))

    gp <- plotlySpectraRelative(params, params, species = species, resource = FALSE)
    expect_s3_class(gp, "plotly")
    expect_equal(tooltip_fields(gp, trace = 2), c("w", "rel_diff", "Legend"))

    gp <- plotlyPredMort(sim, species = species)
    expect_s3_class(gp, "plotly")
    expect_equal(tooltip_fields(gp), c("Species", "w", "Predation mortality"))

    gp <- plotlyFMort(sim, species = species)
    expect_s3_class(gp, "plotly")
    expect_equal(tooltip_fields(gp), c("Species", "w", "Fishing mortality"))

    gp <- plotlyFeedingLevel(sim, species = species, include_critical = TRUE)
    expect_s3_class(gp, "plotly")
    expect_equal(tooltip_fields(gp), c("Species", "w", "Feeding level"))

    gp <- plotlyGrowthCurves(sim, species = species)
    expect_s3_class(gp, "plotly")
    expect_equal(tooltip_fields(gp), c("Species", "Age", "Size [g]"))

    expect_s3_class(plotlyGrowthCurves(params, species = species), "plotly")

    gp <- plotlyDiet(params, species = species[[1]])
    expect_s3_class(gp, "plotly")
    expect_equal(tooltip_fields(gp), c("w", "Proportion", "Prey"))
})

test_that("plotHover on array objects has correct tooltip fields", {
    expect_equal(tooltip_fields(plotHover(getEncounter(params),
                                          species = species)),
                 c("Species", "w", "Encounter rate"))
    expect_equal(tooltip_fields(plotHover(getBiomass(sim), species = species)),
                 c("Species", "Year", "Biomass"))
    expect_equal(tooltip_fields(plotHover(getEncounter(sim),
                                          species = species[[1]])),
                 c("Species", "w", "Encounter rate"))
})

test_that("plotHover(plot(...)) uses concise mizer tooltips", {
    p <- plot(getEncounter(NS_params_small), species = "Cod")
    expect_s3_class(p, "mizer_plot")
    gp <- plotHover(p)
    expect_equal(tooltip_fields(gp), c("Species", "w", "Encounter rate"))
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
                 "containing only")
})

test_that("comparison helpers preserve non-size x variables", {
    p <- plot2(getBiomass(NS_sim_small), getBiomass(NS_sim_small), species = "Cod")
    expect_s3_class(p, "ggplot")
    expect_true("Year" %in% names(p$data))
    expect_error(ggplot2::ggplot_build(p), NA)

    p_rel <- plotRelative(getBiomass(NS_sim_small), getBiomass(NS_sim_small),
                          species = "Cod")
    expect_s3_class(p_rel, "ggplot")
    expect_true("Year" %in% names(p_rel$data))
    expect_error(ggplot2::ggplot_build(p_rel), NA)
})

test_that("plotRelative, plot2, and plotGrowthCurves do not have LineSpec legend", {
    p_rel <- plotRelative(getBiomass(NS_sim_small), getBiomass(NS_sim_small), species = "Cod")
    expect_equal(p_rel$scales$get_scales("linewidth")$guide, "none")

    p2 <- plot2(getBiomass(NS_sim_small), getBiomass(NS_sim_small), species = "Cod")
    expect_equal(p2$scales$get_scales("linewidth")$guide, "none")

    p_growth <- plotGrowthCurves(NS_params_small, species = c("Cod", "Herring"))
    expect_equal(p_growth$scales$get_scales("linewidth")$guide, "none")
    expect_equal(rlang::as_label(p_growth$layers[[1]]$mapping$colour), "LineSpec")
    expect_equal(rlang::as_label(p_growth$layers[[1]]$mapping$linetype), "LineSpec")
    expect_equal(rlang::as_label(p_growth$layers[[1]]$mapping$linewidth), "LineSpec")
    expect_equal(p_growth$scales$get_scales("colour")$name, "Species")
    expect_equal(p_growth$scales$get_scales("linetype")$name, "Species")
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
    expect_true(all(p[[2]] >= 0))
    for (sp in unique(p$Species)) {
        sp_dat <- p[p$Species == sp, ]
        expect_equal(max(sp_dat[[2]]), 1)
        expect_true(all(diff(sp_dat[[2]]) >= -1e-12))
    }

    spectra <- plotSpectra(params, species = species, resource = FALSE,
                           power = 1, wlim = c(1, NA), return_data = TRUE)
    cdf <- plotCDF(params, species = species, resource = FALSE,
                   power = 1, wlim = c(1, NA), normalise = FALSE,
                   return_data = TRUE)
    widths <- params@dw_full[match(spectra$w, params@w_full)]
    expected <- sum(spectra[[2]][spectra$Species == species[[1]]] *
                        widths[spectra$Species == species[[1]]])
    observed <- max(cdf[[2]][cdf$Species == species[[1]]])
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
    expect_true(max(p[[2]]) > 1)

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
    expect_equal(spectra_l$l, expected_l, ignore_attr = TRUE)

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
    # The cumulative now sits on each bin's upper edge (#383), so the top
    # in-window point — where the normalised CDF reaches 1 — may exceed llim[2]
    # by up to one bin, but no more (one weight-bin ratio is a safe length bound
    # since the length ratio beta^(1/b) < beta).
    beta_w <- params_len@w_full[2] / params_len@w_full[1]
    expect_true(all(cdf_l_limited$l <= llim[2] * beta_w))
    y_var <- names(cdf_l_limited)[2]
    expect_equal(stats::aggregate(cdf_l_limited[[y_var]] ~ cdf_l_limited$Species,
                                  FUN = max)[[2]],
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
    expect_s3_class(plotlyYield(sim, species = species, ylim = c(1e-5, 1)), "plotly")
    expect_s3_class(plotlyYieldGear(sim, species = species, ylim = c(1e-5, 1)), "plotly")
})

test_that("yield plotting helpers accept ylim", {
    p <- plotYield(sim, species = species, ylim = c(1e-5, 1))
    expect_equal(p$scales$get_scales("y")$limits, log10(c(1e-5, 1)))

    p_gear <- plotYieldGear(sim, species = species, ylim = c(1e-5, 1))
    expect_equal(p_gear$scales$get_scales("y")$limits, log10(c(1e-5, 1)))
})

test_that("plotYieldGear supports log_x, log_y, and log arguments", {
    p_default <- plotYieldGear(sim, species = species)
    expect_identical(p_default$scales$get_scales("x")$trans$name, "identity")
    expect_identical(p_default$scales$get_scales("y")$trans$name, "log-10")

    p_log_x <- plotYieldGear(sim, species = species, log_x = TRUE)
    expect_identical(p_log_x$scales$get_scales("x")$trans$name, "log-10")

    p_log_y_false <- plotYieldGear(sim, species = species, log_y = FALSE)
    expect_identical(p_log_y_false$scales$get_scales("y")$trans$name, "identity")

    p_log_xy <- plotYieldGear(sim, species = species, log = "xy")
    expect_identical(p_log_xy$scales$get_scales("x")$trans$name, "log-10")
    expect_identical(p_log_xy$scales$get_scales("y")$trans$name, "log-10")
})



# testing the plot outputs
test_that("return_data is identical",{
    expect_equal(dim(plotBiomass(sim, species = species, total = TRUE,
                                 tlim = c(0, 2.8), y_ticks = 4, return_data = TRUE)), c(9,4))
    expect_warning(p <- plotYield(sim, sim0, species = species, return_data = TRUE))
    expect_equal(dim(p), c(8,4))

    expect_equal(dim(plotYieldGear(sim, species = species, return_data = TRUE)), c(8,4))

    expect_equal(dim(plotSpectra(sim, species = species, wlim = c(1,NA),
                                 return_data = TRUE)), c(22, 4))

    expect_equal(dim(plotFeedingLevel(sim, species = species,
                                      return_data = TRUE)), c(40, 3))

    expect_equal(dim(plotPredMort(sim, species = species,
                                  return_data = TRUE)), c(40, 4))

    expect_equal(dim(plotFMort(sim, species = species,
                               return_data = TRUE)), c(40, 4))

    # 200 rows = 50 ages x 2 species x 2 legends (model + von Bertalanffy);
    # the von Bertalanffy curve is added because a, b, k_vb and w_inf are all
    # present (a, b and w_inf are now supplied by default).
    expect_equal(dim(plotGrowthCurves(sim, species = species,
                                      return_data = TRUE)), c(200,4))
    # the following is not a good test because the size of the returned data
    # frame is machine dependent due to the selection of only results above a
    # certain threshold.
    # expect_equal(dim(plotDiet(sim@params, species = species,
    #                           return_data = TRUE)), c(717,3))
}
)

# tlim parameter ----
test_that("tlim filters time axis correctly", {
    # plotBiomass: deprecated start_time/end_time triggers warning
    expect_warning(plotBiomass(sim, start_time = 0), "deprecated")
    expect_warning(plotBiomass(sim, end_time = 2), "deprecated")

    # plotBiomass: tlim reduces number of rows
    all_rows <- nrow(plotBiomass(sim, species = species, return_data = TRUE))
    limited_rows <- nrow(plotBiomass(sim, species = species,
                                     tlim = c(1, 3), return_data = TRUE))
    expect_lt(limited_rows, all_rows)

    # plotYield: tlim reduces number of rows
    all_rows <- nrow(plotYield(sim, species = species, return_data = TRUE))
    limited_rows <- nrow(plotYield(sim, species = species,
                                   tlim = c(1, 3), return_data = TRUE))
    expect_lt(limited_rows, all_rows)

    # plotYieldGear: tlim reduces number of rows
    all_rows <- nrow(plotYieldGear(sim, species = species, return_data = TRUE))
    limited_rows <- nrow(plotYieldGear(sim, species = species,
                                       tlim = c(1, 3), return_data = TRUE))
    expect_lt(limited_rows, all_rows)

    # animate: deprecated time_range triggers warning
    expect_warning(animate(sim, time_range = 1:3), "deprecated")
})

# Legends have the correct entries ----
test_that("Legends have correct entries", {
    # plotSpectra
    p <- plotSpectra(params, species = 2:3)
    expect_length(unique(p$data$Legend), 3)
    p <- plotSpectra(params, species = 1:3, resource = FALSE)
    expect_length(unique(p$data$Legend), 3)
    p <- plotSpectra(params, species = 2:3, total = TRUE)
    expect_length(unique(p$data$Legend), 4)
    p <- plotSpectra(params, species = 2:3, background = TRUE)
    expect_length(unique(p$data$Legend), 3)
    p <- plotSpectra(params_bkgrd, species = 2:3, background = TRUE)
    expect_length(unique(p$data$Legend), 4)
    p <- plotSpectra(params_single)
    expect_length(unique(p$data$Legend), 2)
})

test_that("plotSpectra averages over time range", {
    time_sel <- c(2:4)
    time_range <- getTimes(NS_sim_small)[time_sel]
    # arithmetic mean
    df <- plotSpectra(NS_sim_small, species = 1, time_range = time_range,
                      power = 0, return_data = TRUE)
    expected <- mean(NS_sim_small@n[time_sel, 1, 1])
    expect_equal(df[[2]][1], expected)
    # geometric mean
    df <- plotSpectra(NS_sim_small, species = 1, time_range = time_range,
                      geometric_mean = TRUE,
                      power = 0, return_data = TRUE)
    expected <- exp(mean(log(NS_sim_small@n[time_sel, 1, 1])))
    expect_equal(df[[2]][1], expected)
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
                 "containing only")
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
    diet_data <- plotDiet(params, species = "3", return_data = TRUE)

    # Get abundance data for species "3"
    sp_idx <- which(dimnames(params@initial_n)[[1]] == "3")
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
    p <- plotDiet(sim, species = 2) # Species 11 is Cod in setup
    expect_true(is(p, "ggplot"))
    p <- plotDiet(sim, species = 2, time_range = 1:2)
    expect_true(is(p, "ggplot"))
})

test_that("plotDiet.MizerSim uses the simulated abundance, not the initial", {
    # sim is projected with fishing effort = 1, so its diet at later times
    # differs from the diet computed from the initial (params) abundances.
    times <- as.numeric(dimnames(sim@n)$time)
    diet_params <- getDiet(sim@params)
    diet_final <- getDiet(sim, time_range = max(times), drop = TRUE)
    expect_false(isTRUE(all.equal(diet_params, diet_final,
                                  check.attributes = FALSE)))
    # By default the MizerSim plot uses the final time step, so its data matches
    # the MizerParams plot built from the diet and abundance at that time step.
    last <- dim(sim@n)[1]
    n_last <- apply(sim@n[last, , , drop = FALSE], c(2, 3), mean)
    expect_equal(
        plotDiet(sim, species = 2, return_data = TRUE),
        plot_diet(sim@params, n = n_last, diet = diet_final, species = 2,
                  log_x = TRUE, log_y = FALSE, wlim = c(NA, NA),
                  llim = c(NA, NA), size_axis = "w", return_data = TRUE))
})

test_that("plotDiet.MizerSim averages rates, not proportions, over a range", {
    tr <- range(as.numeric(dimnames(sim@n)$time))
    # Correct aggregation: average the consumption rates then normalise
    rates <- getDiet(sim, time_range = tr, drop = FALSE, proportion = FALSE)
    rates_mean <- apply(rates, c(2, 3, 4), mean)
    total <- rowSums(rates_mean, dims = 2)
    diet_expected <- sweep(rates_mean, c(1, 2), total, "/")
    diet_expected[is.nan(diet_expected)] <- 0
    # Naive (wrong) aggregation: average the per-step proportions
    props <- getDiet(sim, time_range = tr, drop = FALSE, proportion = TRUE)
    diet_naive <- apply(props, c(2, 3, 4), mean)
    # The two genuinely differ for this projection
    expect_false(isTRUE(all.equal(diet_expected, diet_naive)))
    # plotDiet uses the correct (rate-averaged) proportions
    time_elements <- get_time_elements(sim, tr)
    n_avg <- apply(sim@n[time_elements, , , drop = FALSE], c(2, 3), mean)
    expect_equal(
        plotDiet(sim, species = 2, time_range = tr, return_data = TRUE),
        plot_diet(sim@params, n = n_avg, diet = diet_expected, species = 2,
                  log_x = TRUE, log_y = FALSE, wlim = c(NA, NA),
                  llim = c(NA, NA), size_axis = "w", return_data = TRUE))
})

# Second-order power weighting in plotSpectra / plotCDF (#383) --------------

test_that("plotSpectra draws the spectrum at bin centres with the w^power weight there", {
    p0 <- params
    p1 <- params
    second_order_w(p1) <- c(bin_average = TRUE)
    beta <- p0@w_full[2] / p0@w_full[1]
    for (pw in c(0, 1, 2)) {
        d0 <- plotSpectra(p0, power = pw, return_data = TRUE)
        d1 <- plotSpectra(p1, power = pw, return_data = TRUE)
        d0 <- d0[order(d0$Species, d0$w), ]
        d1 <- d1[order(d1$Species, d1$w), ]
        expect_equal(nrow(d0), nrow(d1))
        # x moves to the geometric bin centre (a uniform sqrt(beta) shift) ...
        expect_equal(unname(d1$w / d0$w), rep(sqrt(beta), nrow(d1)))
        # ... and the w^power weight is evaluated there, scaling the value
        # (column 2, named by the y-label) by (w*/w)^power = beta^(power/2).
        expect_equal(unname(d1[[2]] / d0[[2]]),
                     rep(beta^(pw / 2), nrow(d1)))
    }
})

test_that("plotSpectra default (first order) is unchanged", {
    expect_identical(plotSpectra(params, power = 2, return_data = TRUE),
                     plotSpectra(params, power = 2, return_data = TRUE))
    # The default model never shifts: x stays on the model grid nodes.
    d <- plotSpectra(params, power = 2, resource = FALSE, total = FALSE,
                     background = FALSE, return_data = TRUE)
    expect_true(all(d$w %in% params@w))
})

test_that("plotCDF places the cumulative on upper bin edges (inclusive convention)", {
    upper_edges <- round(params@w_full + params@dw_full, 6)
    # The inclusive cumulative sum belongs on each bin's *upper* edge w_k + dw_k,
    # in both the default and the second-order schemes.
    cdf0 <- plotCDF(params, power = 2, return_data = TRUE)
    expect_true(all(round(cdf0$w, 6) %in% upper_edges))
    # The left bin edges (nodes) are no longer used for the placement.
    nodes <- round(params@w, 6)
    expect_false(all(round(cdf0$w, 6) %in% nodes))

    p1 <- params
    second_order_w(p1) <- c(bin_average = TRUE)
    cdf <- plotCDF(p1, power = 2, return_data = TRUE)
    centres <- round(c(bin_midpoints(p1), bin_midpoints(p1, w_full = TRUE)), 6)
    expect_true(all(round(cdf$w, 6) %in% upper_edges))
    # never on the geometric bin centres
    expect_false(any(round(cdf$w, 6) %in% setdiff(centres, upper_edges)))
    # Cumulative is monotonic increasing per species.
    sp1 <- cdf[cdf$Species == p1@species_params$species[1], ]
    sp1 <- sp1[order(sp1$w), ]
    expect_true(all(diff(sp1[[2]]) >= -1e-12))
})

test_that("plotCDF cumulative is second-order and free of the one-bin offset", {
    # For a community spectrum N = C w^-2 with power = 2 the integrand N w^power
    # is constant, so the exact CDF is linear in w: F(w) = V (w - w_min). The
    # bin-average of C w^-2 is exact at the geometric centre, so the plotted
    # density value is the constant V at every bin. We can therefore check the
    # cumulative analytically, which pins down the edge placement: a one-bin
    # offset (plotting F at the *left* edge) would break F(w_min-bin) != 0 and
    # the slope, so this is a sharp regression test for issue #383.
    p <- params
    second_order_w(p) <- c(bin_average = TRUE)
    sp <- 1
    p@initial_n[sp, ] <- p@w^(-2)            # community spectrum N proportional to w^-2

    # The weighted density is the constant V = beta at every bin centre.
    dens <- plotSpectra(p, species = sp, power = 2, resource = FALSE,
                        total = FALSE, background = FALSE, return_data = TRUE)
    V <- mean(dens[[2]])
    expect_equal(unname(dens[[2]]), rep(V, nrow(dens)))

    cdf <- plotCDF(p, species = sp, power = 2, resource = FALSE, total = FALSE,
                   background = FALSE, normalise = FALSE, return_data = TRUE)
    cdf <- cdf[order(cdf$w), ]
    w_min <- min(p@w)
    # F plotted at the upper edge equals the exact integral V (w_upper - w_min).
    expect_equal(unname(cdf[[2]]), V * (cdf$w - w_min), tolerance = 1e-8)
    # The smallest plotted x is the upper edge of the first bin, and its value is
    # one full bin integral (not zero, and not two bins) — i.e. no offset.
    expect_equal(cdf$w[1], unname(p@w[1] + p@dw_full[match(p@w[1], p@w_full)]))
    expect_equal(cdf[[2]][1], V * (cdf$w[1] - w_min), tolerance = 1e-8)
})

test_that("plotSpectraRelative shifts x to centres", {
    p1a <- params
    p1b <- params
    p1b@initial_n <- p1b@initial_n * 1.5
    second_order_w(p1a) <- c(bin_average = TRUE)
    second_order_w(p1b) <- c(bin_average = TRUE)
    p_centre <- plotSpectraRelative(p1a, p1b, species = species, resource = FALSE)
    d_centre <- p_centre$data[order(p_centre$data$Species, p_centre$data$w), ]
    # The x-location is the geometric bin centre, not the node.
    p_node <- plotSpectraRelative(params, params, species = species,
                                  resource = FALSE)
    expect_true(all(p_node$data$w %in% params@w))
    expect_false(all(d_centre$w %in% params@w))
})
