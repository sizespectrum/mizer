# Plot fixtures use edition 1, as do the existing visual baselines.
withr::local_options(mizer_defaults_edition = 1)
params <- plot_params_small
sim <- plot_sim_small
species <- c(1, 2)

test_that("rate preparation averages saved times and keeps plotting metadata", {
    for (second in c(FALSE, TRUE)) {
        sim <- NS_sim_small
        second_order_w(sim@params) <- c(bin_average = second)
        getter <- function(object, time_range, drop) {
            expect_equal(time_range, c(1, 3))
            expect_false(drop)
            # Each species and size has 0.2 at the first time and 0.6 at the
            # second. Averaging any dimension except time gives a different
            # shape or leaves those two values in the result.
            array(rep(c(0.2, 0.6), nrow(object@params@initial_n) *
                          length(object@params@w)),
                  dim = c(2, dim(object@params@initial_n)),
                  dimnames = c(list(time = c("1", "3")),
                               dimnames(object@params@initial_n)))
        }
        rate <- prepare_plot_rate(sim, getter, c(1, 3),
                                  value_name = "Feeding level", units = "",
                                  type = "proportion")
        expect_s3_class(rate, "ArraySpeciesBySize")
        expect_equal(as.numeric(rate), rep(0.4, length(sim@params@initial_n)))
        expect_equal(dimnames(rate), dimnames(sim@params@initial_n))
        expect_identical(array_type(rate), "proportion")
        expect_identical(attr(rate, "value_name"), "Feeding level")
        expect_identical(attr(rate, "units"), "")
        expected_w <- if (second) bin_midpoints(sim@params) else sim@params@w
        expect_equal(get_ArraySpeciesBySize_w(rate), expected_w)
    }
})

test_that("rate preparation accepts plain matrices at the final time", {
    sim <- NS_sim_small
    values <- sim@params@initial_n["Cod", , drop = FALSE]
    getter <- function(object, time_range, drop) {
        expect_equal(time_range, max(getTimes(object)))
        expect_false(drop)
        values
    }
    rate <- prepare_plot_rate(sim, getter, value_name = "Custom rate",
                              units = "1/year")
    expect_s3_class(rate, "ArraySpeciesBySize")
    expect_equal(as.numeric(rate), as.numeric(values))
    expect_equal(dimnames(rate), dimnames(values))
    expect_identical(attr(rate, "value_name"), "Custom rate")
    expect_identical(attr(rate, "representation"), "average")
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

test_that("plotFeedingLevel shows the whole of [0, 1] without hiding data", {
    expect_identical(plotFeedingLevel(params)$coordinates$limits$y, c(0, 1))

    # A critical feeding level above 1 stays visible, where the fixed [0, 1]
    # window used to draw it off the top of the plot
    params_hungry <- params
    params_hungry@metab[] <- params@metab * 5
    plotted <- plotFeedingLevel(params_hungry, include_critical = TRUE,
                                return_data = TRUE)[["Feeding level"]]
    expect_gt(max(plotted), 1)
    limits <- plotFeedingLevel(params_hungry,
                               include_critical = TRUE)$coordinates$limits$y
    expect_identical(limits[[1]], 0)
    expect_gte(limits[[2]], max(plotted))

    # A logarithmic axis is left to the data, as before
    expect_true(all(is.na(
        plotFeedingLevel(params, log_y = TRUE)$coordinates$limits$y)))
})

test_that("plotPredMort tolerates a log y axis", {
    params <- NS_params_small
    sim <- project(params, t_max = 1, effort = 1)
    # Predation mortality is zero on the largest individuals and the plot asks
    # for a y-axis starting at zero, so this used to warn twice.
    expect_no_warning(ggplot2::ggplot_build(plotPredMort(sim, log_y = TRUE)))
    expect_no_warning(ggplot2::ggplot_build(plotFMort(sim, log_y = TRUE)))
})

test_that("plotFeedingLevel and plotlyFeedingLevel with include_critical", {
    params <- NS_params_small

    p <- plotFeedingLevel(params, include_critical = TRUE)
    expect_s3_class(p, "ggplot")
    expect_equal(length(p$layers), 2)

    df <- plotFeedingLevel(params, include_critical = TRUE, return_data = TRUE)
    expect_true(all(c("actual", "critical") %in% df$Type))

    pl <- plotlyFeedingLevel(params, include_critical = TRUE)
    expect_s3_class(pl, "plotly")
})

test_that("rate plots return ggplots", {
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
})

test_that("interactive wrappers preserve tooltip fields", {
    gp <- plotlyPredMort(sim, species = species)
    expect_s3_class(gp, "plotly")
    expect_equal(tooltip_fields(gp), c("Species", "w", "Predation mortality"))

    gp <- plotlyFMort(sim, species = species)
    expect_s3_class(gp, "plotly")
    expect_equal(tooltip_fields(gp), c("Species", "w", "Fishing mortality"))

    gp <- plotlyFeedingLevel(sim, species = species, include_critical = TRUE)
    expect_s3_class(gp, "plotly")
    expect_equal(tooltip_fields(gp), c("Species", "w", "Feeding level"))
})

test_that("rate plots support length axes", {
    params_len <- params
    params_len@species_params$a <- 0.01
    params_len@species_params$b <- 3
    sim_len <- sim
    sim_len@params <- params_len

    expect_true("l" %in% names(plotFeedingLevel(params_len, species = species,
                                                size_axis = "l",
                                                return_data = TRUE)))
    expect_true("l" %in% names(plotPredMort(params_len, species = species,
                                            size_axis = "l",
                                            return_data = TRUE)))
    expect_true("l" %in% names(plotFMort(params_len, species = species,
                                         size_axis = "l",
                                         return_data = TRUE)))
})

test_that("public return-data shape is preserved", {
    expect_equal(dim(plotFeedingLevel(sim, species = species,
                                      return_data = TRUE)), c(40, 4))

    expect_equal(dim(plotPredMort(sim, species = species,
                                  return_data = TRUE)), c(40, 4))

    expect_equal(dim(plotFMort(sim, species = species,
                               return_data = TRUE)), c(40, 4))
})

test_that("size viewports preserve returned data and totals", {
    for (second in c(FALSE, TRUE)) {
        p <- NS_params_small
        second_order_w(p) <- c(bin_average = second)
        for (axis in c("w", "l")) {
            for (fun in list(plotFeedingLevel, plotPredMort, plotFMort)) {
                full <- fun(p, species = "Cod", size_axis = axis,
                            return_data = TRUE)
                zoom <- fun(p, species = "Cod", size_axis = axis,
                            wlim = c(1, 100), llim = c(5, 20),
                            ylim = c(0.2, 0.8), return_data = TRUE)
                expect_equal(zoom, full)
            }
        }
    }
})

test_that("critical feeding levels keep the requested viewport", {
    feeding <- plotFeedingLevel(NS_params_small, wlim = c(1, 100),
                                ylim = c(0.2, 0.8), include_critical = TRUE)
    expect_equal(feeding$coordinates$limits$x, c(1, 100))
    expect_equal(feeding$coordinates$limits$y, c(0.2, 0.8))
})
