# Plot fixtures use edition 1, as do the existing visual baselines.
withr::local_options(mizer_defaults_edition = 1)
params <- plot_params_small
sim <- plot_sim_small
species <- c(1, 2)

test_that("plotGrowthCurves validates size_at_age input", {
    expect_error(plotGrowthCurves(params, size_at_age = data.frame(age = 1, weight = 2)),
                 "needs to have a 'species' column")
    expect_error(plotGrowthCurves(params, size_at_age = data.frame(species = "10", weight = 2)),
                 "needs to have an 'age' column")
    expect_error(plotGrowthCurves(params, size_at_age = data.frame(species = "10", age = 1)),
                 "needs to have either a 'length' or a 'weight' column")
})

test_that("growth curves return a ggplot", {
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
})

test_that("interactive wrappers preserve tooltip fields", {
    gp <- plotlyGrowthCurves(sim, species = species)
    expect_s3_class(gp, "plotly")
    expect_equal(tooltip_fields(gp), c("Species", "Age", "Size [g]"))

    expect_s3_class(plotlyGrowthCurves(params, species = species), "plotly")
})

test_that("growth curves hide the LineSpec legend", {
    p_growth <- plotGrowthCurves(NS_params_small, species = c("Cod", "Herring"))
    expect_equal(p_growth$scales$get_scales("linewidth")$guide, "none")
    expect_equal(rlang::as_label(p_growth$layers[[1]]$mapping$colour), "LineSpec")
    expect_equal(rlang::as_label(p_growth$layers[[1]]$mapping$linetype), "LineSpec")
    expect_equal(rlang::as_label(p_growth$layers[[1]]$mapping$linewidth), "LineSpec")
    expect_equal(p_growth$scales$get_scales("colour")$name, "Species")
    expect_equal(p_growth$scales$get_scales("linetype")$name, "Species")
})

test_that("public return-data shape is preserved", {
    # 200 rows = 50 ages x 2 species x 2 legends (model + von Bertalanffy);
    # the von Bertalanffy curve is added because a, b, k_vb and w_inf are all
    # present (a, b and w_inf are now supplied by default).
    expect_equal(dim(plotGrowthCurves(sim, species = species,
                                      return_data = TRUE)), c(200,4))
})
