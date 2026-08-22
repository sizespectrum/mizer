# plotYieldVsF ----------------------------------------------------------------
#
# plotYieldVsF() is a thin wrapper over scanModel() with scanFishingMortality()
# as its setter; the machinery itself is tested in test-scanModel.R and
# test-scan_setters.R. What is tested here is the wrapper: that it measures the
# yield of the one species asked for, that it hands scanModel() the setter and
# the sweep order it is supposed to, and that its own arguments are checked.

fast <- list(t_max = 6, t_per = 1.5, dt = 0.5, t_save = 0.5, t_sample = 3,
             progress_bar = FALSE)

yield_vs_f <- function(...) {
    suppressMessages(suppressWarnings(
        do.call(plotYieldVsF, c(list(...), fast))))
}

test_that("plotYieldVsF returns a MizerScan with return_data", {
    d <- yield_vs_f(NS_params_small, "Cod", F_range = c(0, 0.4, 0.8),
                    return_data = TRUE)
    expect_true(is.MizerScan(d))
    expect_s3_class(d, "data.frame")
    expect_named(d, c("Fishing mortality on Cod", "Yield rate", "Species",
                      "ymin", "ymax", "termination", "converged", "attractor",
                      "period", "residual"))
    # One row per fishing mortality, for the one species asked for
    expect_setequal(unique(d$Species), "Cod")
    expect_equal(nrow(d), 3)
    expect_setequal(unique(d[[1]]), c(0, 0.4, 0.8))
    # The envelope is well ordered, and the yield vanishes at zero fishing
    expect_true(all(d$ymax >= d$ymin))
    expect_equal(d[[2]][d[[1]] == 0], 0)
    expect_true(all(d[[2]][d[[1]] > 0] > 0))
    # The x axis carries the units the setter supplies
    expect_identical(attr(d, "scan_units"), "1/year")
})

test_that("plotYieldVsF returns a mizer_plot object", {
    p <- yield_vs_f(NS_params_small, "Cod", F_range = c(0.2, 0.6))
    expect_s3_class(p, "mizer_plot")
    expect_s3_class(p, "ggplot")
})

test_that("plotYieldVsF leaves the fishing on the other species alone", {
    # The whole point of scanning the fishing mortality rather than the effort:
    # only Cod's yield responds, because only Cod's fishing is varied.
    d <- yield_vs_f(NS_params_small, "Cod", F_range = c(0.2, 1.2),
                    return_data = TRUE)
    yield <- d[[2]][order(d[[1]])]
    expect_length(yield, 2)
    expect_false(isTRUE(all.equal(yield[[1]], yield[[2]])))

    # The effort of the gears catching the other species is untouched
    p <- attr(d, "params")
    expect_identical(initial_effort(p)[c("Industrial", "Pelagic")],
                     initial_effort(NS_params_small)[c("Industrial", "Pelagic")])
})

test_that("F_range is built from F_min, F_max and no_steps when missing", {
    d <- yield_vs_f(NS_params_small, "Cod", F_max = 0.9, no_steps = 4,
                    return_data = TRUE)
    expect_setequal(unique(d[[1]]), seq(0, 0.9, length.out = 4))
})

test_that("an F_MSY species parameter is carried through as a reference line", {
    p <- NS_params_small
    p@species_params$F_MSY <- c(NA, NA, 0.35)
    d <- yield_vs_f(p, "Cod", F_range = c(0.2, 0.6), return_data = TRUE)
    expect_identical(attr(d, "reference_lines"), c(F_MSY = 0.35))
})

test_that("plotYieldVsF validates its arguments", {
    expect_error(
        plotYieldVsF(NS_params_small, c("Cod", "Herring")),
        "one species at a time"
    )
    expect_error(
        suppressWarnings(plotYieldVsF(NS_params_small, "Nonexistent")),
        "No species have been selected"
    )
    expect_error(
        plotYieldVsF(NS_params_small, "Cod", F_range = 0.5),
        "length\\(F_range\\) not greater than or equal to 2"
    )
    expect_error(
        plotYieldVsF(NS_params_small, "Cod", F_min = 1, F_max = 0),
        "F_max not greater than F_min"
    )
})
