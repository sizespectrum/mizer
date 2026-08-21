# plotBifurcation ------------------------------------------------------------
#
# plotBifurcation() is now a thin wrapper over scanModel(); the machinery itself
# is tested in test-scanModel.R.

fast <- list(t_max = 6, t_per = 1.5, dt = 0.5, t_save = 0.5, t_sample = 3,
             progress_bar = FALSE)

bifurcation <- function(...) {
    suppressMessages(suppressWarnings(
        do.call(plotBifurcation, c(list(...), fast))))
}

test_that("plotBifurcation returns a MizerScan with return_data", {
    d <- bifurcation(NS_params_small, effort = c(0, 0.5, 1), return_data = TRUE)
    expect_true(is.MizerScan(d))
    expect_s3_class(d, "data.frame")
    expect_named(d, c("Effort", "Biomass", "Species", "ymin", "ymax", "type",
                      "settled", "period", "residual"))
    # one row per (effort, species)
    expect_equal(nrow(d), 3 * length(valid_species_arg(NS_params_small)))
    expect_setequal(unique(d$Effort), c(0, 0.5, 1))
    # envelope is well ordered and positive for biomass
    expect_true(all(d$ymax >= d$ymin))
    expect_true(all(d$ymin > 0))
})

test_that("plotBifurcation returns a mizer_plot object", {
    p <- bifurcation(NS_params_small, effort = c(0, 1))
    expect_s3_class(p, "mizer_plot")
    expect_s3_class(p, "ggplot")
})

test_that("plotBifurcation respects the species and value arguments", {
    d <- bifurcation(NS_params_small, effort = c(0, 1), species = "Cod",
                     value = "yield", return_data = TRUE)
    expect_setequal(unique(d$Species), "Cod")
    # The y-axis label comes from the measured quantity itself
    expect_identical(attr(d, "value_name"), "Yield rate")
    expect_identical(names(d)[[2]], "Yield rate")
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

test_that("the arguments replaced by scanModel() are soft-deprecated", {
    # Called directly rather than through the helper above, which suppresses
    # warnings and already supplies t_save.
    base <- list(NS_params_small, effort = c(0, 1), t_max = 6, t_per = 1.5,
                 dt = 0.5, progress_bar = FALSE)

    expect_warning(
        d <- suppressMessages(do.call(plotBifurcation,
                                      c(base, list(return_data = TRUE,
                                                   t_sample_default = 3)))),
        class = "lifecycle_warning_deprecated")
    # and it still does what t_sample now does
    expect_identical(attr(d, "settings")$t_sample, 3)

    expect_warning(
        suppressMessages(do.call(plotBifurcation,
                                 c(base, list(t_save = 0.25)))),
        class = "lifecycle_warning_deprecated")

    expect_warning(
        p <- suppressMessages(do.call(plotBifurcation,
                                      c(base, list(ytrans = "identity")))),
        class = "lifecycle_warning_deprecated")
    expect_s3_class(p, "ggplot")
})

test_that("plotBifurcation() resolves `species` the way mizer does elsewhere", {
    # Background species are not target species, so they are left out by
    # default, and a logical or numeric argument indexes the model's species.
    p <- NS_params_small
    sp <- species_params(p)
    sp$is_background <- c(FALSE, TRUE, FALSE)
    species_params(p) <- sp
    d <- suppressMessages(suppressWarnings(
        plotBifurcation(p, effort = c(0, 1), t_max = 4, t_sample = 2,
                        return_data = TRUE, progress_bar = FALSE)))
    expect_identical(unique(as.character(d$Species)), c("Sprat", "Cod"))

    d <- suppressMessages(suppressWarnings(
        plotBifurcation(NS_params_small, effort = c(0, 1),
                        species = c(TRUE, FALSE, TRUE), t_max = 4,
                        t_sample = 2, return_data = TRUE,
                        progress_bar = FALSE)))
    expect_identical(unique(as.character(d$Species)), c("Sprat", "Cod"))
})
