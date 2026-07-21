# plotBifurcation ------------------------------------------------------------

test_that("plotBifurcation returns a tidy data frame with return_data", {
    d <- plotBifurcation(NS_params_small, effort = c(0, 0.5, 1),
                         t_max = 6, t_sample_default = 3, t_save = 0.5,
                         progress_bar = FALSE, return_data = TRUE)
    expect_s3_class(d, "data.frame")
    expect_named(d, c("Effort", "Species", "ymin", "ymax", "type"))
    # one row per (effort, species)
    expect_equal(nrow(d), 3 * length(valid_species_arg(NS_params_small)))
    expect_setequal(unique(d$Effort), c(0, 0.5, 1))
    # envelope is well ordered and positive for biomass
    expect_true(all(d$ymax >= d$ymin))
    expect_true(all(d$ymin > 0))
})

test_that("plotBifurcation returns a mizer_plot object", {
    p <- plotBifurcation(NS_params_small, effort = c(0, 1),
                         t_max = 6, t_sample_default = 3, t_save = 0.5,
                         progress_bar = FALSE)
    expect_s3_class(p, "mizer_plot")
    expect_s3_class(p, "ggplot")
})

test_that("plotBifurcation respects the species and value arguments", {
    d <- plotBifurcation(NS_params_small, effort = c(0, 1), species = "Cod",
                         value = "yield", t_max = 6, t_sample_default = 3,
                         t_save = 0.5, progress_bar = FALSE, return_data = TRUE)
    expect_setequal(unique(d$Species), "Cod")
})

test_that("plotBifurcation validates its arguments", {
    expect_error(
        plotBifurcation(NS_params_small, effort = 1),
        "length\\(effort\\) not greater than or equal to 2"
    )
    expect_error(
        plotBifurcation(NS_params_small, t_max = -1),
        "t_max not greater than 0"
    )
})
