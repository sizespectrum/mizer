# Cache resource arrays once to avoid recomputing them in every test block.
mort_small <- getResourceMort(NS_params_small)
n_resource_small <- NResource(NS_sim_small)
# A second array to compare against, with all values doubled.
mort_small2 <- ArrayResourceBySize(unclass_resource(mort_small) * 2,
                                   value_name = attr(mort_small, "value_name"),
                                   units = attr(mort_small, "units"),
                                   params = NS_params_small)

test_that("ArrayResourceBySize constructor works", {
    vec <- 1:10
    names(vec) <- signif(NS_params_small@w_full[1:10], 3)

    rate <- ArrayResourceBySize(vec, value_name = "Test rate", units = "1/year")

    expect_s3_class(rate, "ArrayResourceBySize")
    expect_true(is.ArrayResourceBySize(rate))
    expect_false(is.ArrayResourceBySize(vec))
    expect_identical(length(rate), length(vec))
    expect_identical(attr(rate, "value_name"), "Test rate")
    expect_identical(attr(rate, "units"), "1/year")
})

test_that("ArrayResourceBySize constructor validates input", {
    expect_error(ArrayResourceBySize(matrix(1:4, 2, 2)), "`x` must be a numeric vector.")
})

test_that("Rate functions return ArrayResourceBySize", {
    expect_true(is.ArrayResourceBySize(mort_small))
    expect_true(is.ArrayResourceBySize(initialNResource(NS_params_small)))
})

test_that("print.ArrayResourceBySize works", {
    expect_output(print(mort_small), "Resource mortality")
    expect_output(print(mort_small), "sizes")
    expect_output(print(mort_small), "1/year")
})

test_that("print.ArrayResourceBySize truncates a wide size grid", {
    mort <- getResourceMort(NS_params)
    expect_output(print(mort), "showing .* of \\d+ sizes")
    expect_output(print(mort), "log-spaced")
    expect_output(print(mort), "as.data.frame")
})

test_that("summary.ArrayResourceBySize works", {
    s <- summary(mort_small)
    expect_s3_class(s, "summary.ArrayResourceBySize")
    expect_identical(s$value_name, "Resource mortality")
    expect_output(print(s), "Resource mortality")
})

test_that("str.ArrayResourceBySize works", {
    expect_output(str(mort_small), "ArrayResourceBySize")
    expect_output(str(mort_small), "Resource mortality")
})

test_that("plot2.ArrayResourceBySize compares compatible arrays", {
    p <- plot2(mort_small, mort_small2, name1 = "Original", name2 = "Changed",
               wlim = c(1, NA), log = "xy")
    expect_s3_class(p, "ggplot")
    expect_identical(levels(p$data$Model), c("Original", "Changed"))
    expect_identical(unique(as.character(p$data$Legend)), "Resource")
    expect_true(all(p$data$w >= 1))
    expect_identical(p$scales$get_scales("x")$trans$name, "log-10")
    expect_identical(p$scales$get_scales("y")$trans$name, "log-10")

    p_none <- plot2(mort_small, mort_small2, log = "")
    expect_identical(p_none$scales$get_scales("x")$trans$name, "identity")
    expect_identical(p_none$scales$get_scales("y")$trans$name, "identity")

    expect_error(plot2(mort_small, getEncounter(NS_params_small)),
                 "Both objects must be")
    expect_error(plot2(mort_small, mort_small2, log = "z"), "containing only")
})

test_that("plotRelative.ArrayResourceBySize plots symmetric relative difference", {
    p <- plotRelative(mort_small, mort_small2, wlim = c(1, NA))
    expect_s3_class(p, "ggplot")
    expect_true(all(p$data$w >= 1))
    expect_true(all(abs(p$data$rel_diff - 2 / 3) < 1e-12))
    expect_identical(p$scales$get_scales("x")$trans$name, "log-10")

    p_linear <- plotRelative(mort_small, mort_small2, log_x = FALSE)
    expect_identical(p_linear$scales$get_scales("x")$trans$name, "identity")

    expect_error(plotRelative(mort_small, getEncounter(NS_params_small)),
                 "Both objects must be")
})

test_that("addPlot.ArrayResourceBySize adds lines to an existing ggplot", {
    p <- plot(mort_small)
    p2 <- addPlot(p, mort_small2, linetype = "dashed", alpha = 0.5)

    expect_s3_class(p2, "ggplot")
    expect_equal(length(p2$layers), length(p$layers) + 1)
    added <- p2$layers[[length(p2$layers)]]
    expect_identical(added$aes_params$linetype, "dashed")
    expect_identical(added$aes_params$alpha, 0.5)
    # Colour is fixed rather than mapped, because the existing plot's colour
    # scale need not contain a "Resource" level.
    expect_identical(added$aes_params$colour,
                     NS_params_small@linecolour[["Resource"]])
    expect_null(added$mapping$colour)

    # The original plot is not modified in place
    p2$layers[[1]]$aes_params$alpha <- 0.25
    expect_null(p$layers[[1]]$aes_params$alpha)

    expect_error(addPlot("not a plot", mort_small), "ggplot")
    expect_error(addPlot(plot(getBiomass(NS_sim_small)), mort_small),
                 "x variable `w`")
})

test_that("resource plotting methods warn about species/total/background", {
    expect_warning(plot2(mort_small, mort_small2, species = "Cod"),
                   "`species` is not used for resource arrays")
    expect_warning(plotRelative(mort_small, mort_small2, total = TRUE),
                   "`total` is not used for resource arrays")
    expect_warning(addPlot(plot(mort_small), mort_small, background = FALSE),
                   "`background` is not used for resource arrays")
    expect_warning(animate(n_resource_small, species = "Cod", total = TRUE),
                   "`species`, `total` are not used for resource arrays")
    # The defaults must not warn
    expect_silent(warn_unused_resource_args())
})

test_that("ArrayTimeByResourceBySize constructor works", {
    mat <- matrix(1:20, nrow = 2, ncol = 10)
    rownames(mat) <- as.character(2000:2001)
    colnames(mat) <- signif(NS_params_small@w_full[1:10], 3)

    out <- ArrayTimeByResourceBySize(mat, value_name = "Test rate", units = "1/g")

    expect_s3_class(out, "ArrayTimeByResourceBySize")
    expect_true(is.ArrayTimeByResourceBySize(out))
    expect_false(is.ArrayTimeByResourceBySize(mat))
    expect_identical(dim(out), dim(mat))
    expect_identical(attr(out, "value_name"), "Test rate")
    expect_identical(attr(out, "units"), "1/g")
})

test_that("ArrayTimeByResourceBySize constructor validates input", {
    expect_error(ArrayTimeByResourceBySize(1:10), "`x` must be a matrix.")
})

test_that("NResource returns ArrayTimeByResourceBySize", {
    expect_true(is.ArrayTimeByResourceBySize(n_resource_small))
})

test_that("print.ArrayTimeByResourceBySize works", {
    expect_output(print(n_resource_small), "Number density")
    expect_output(print(n_resource_small), "times x")
    expect_output(print(n_resource_small), "1/g")
})

test_that("print.ArrayTimeByResourceBySize truncates a long time series", {
    sim_long <- suppressMessages(
        project(NS_params_small, t_max = 60, dt = 0.5, t_save = 0.5,
               progress_bar = FALSE))
    n_resource <- NResource(sim_long)
    expect_identical(nrow(n_resource), 121L)
    out <- paste(capture.output(print(n_resource)), collapse = "\n")
    expect_match(out, "\\b0\\b")
    expect_match(out, "\\b60\\b")
    expect_match(out, "showing 8 of 121 times")
    expect_match(out, "evenly spaced")
})

test_that("summary.ArrayTimeByResourceBySize works", {
    s <- summary(n_resource_small)
    expect_s3_class(s, "summary.ArrayTimeByResourceBySize")
    expect_identical(s$value_name, "Number density")
    expect_output(print(s), "Number density")
})

test_that("str.ArrayTimeByResourceBySize works", {
    expect_output(str(n_resource_small), "ArrayTimeByResourceBySize")
    expect_output(str(n_resource_small), "Number density")
})

test_that("ArrayTimeByResourceBySize comparison methods slice and delegate", {
    times <- as.numeric(dimnames(n_resource_small)[[1]])

    p <- plot2(n_resource_small, n_resource_small,
               name1 = "First", name2 = "Second")
    expect_s3_class(p, "ggplot")
    expect_identical(levels(p$data$Model), c("First", "Second"))

    # By default the final time step is shown, and `time` selects another
    expect_equal(p$data[p$data$Model == "First", ][[2]],
                 plot(n_resource_small, return_data = TRUE)[[2]])
    p_first <- plot2(n_resource_small, n_resource_small, time = times[1])
    expect_equal(p_first$data[p_first$data$Model == "First", ][[2]],
                 plot(n_resource_small, time = times[1],
                      return_data = TRUE)[[2]])

    pr <- plotRelative(n_resource_small, n_resource_small)
    expect_s3_class(pr, "ggplot")
    expect_true(all(abs(pr$data$rel_diff) < 1e-12))

    pa <- addPlot(plot(n_resource_small), n_resource_small, time = times[1])
    expect_s3_class(pa, "ggplot")
    expect_equal(length(pa$layers), length(plot(n_resource_small)$layers) + 1)

    expect_error(plot2(n_resource_small, mort_small), "Both objects must be")
})

test_that("animate dispatches on ArrayTimeByResourceBySize", {
    times <- as.numeric(dimnames(n_resource_small)[[1]])

    p <- animate(n_resource_small)
    expect_s3_class(p, "plotly")
    traces <- plotly::plotly_build(p)$x$data
    expect_identical(unique(vapply(traces, function(tr) tr$name, character(1))),
                     "Resource")

    frames <- function(x) {
        sort(unique(vapply(plotly::plotly_build(x)$x$frames,
                           function(f) f$name, character(1))))
    }
    expect_length(frames(p), length(times))
    expect_length(frames(animate(n_resource_small,
                                 tlim = c(times[2], times[3]))), 2)

    expect_error(animate(n_resource_small, size_axis = "l"),
                 "length axis is not available")
})
