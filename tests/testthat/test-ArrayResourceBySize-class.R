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

    # By default the final time step is shown, and `time` selects another.
    # The plots use a logarithmic y axis, which cannot show the sizes where the
    # resource density is zero, so those are dropped before the data reaches
    # the plot but not from `return_data`.
    plottable <- function(x) x[x > 0]
    expect_equal(p$data[p$data$Model == "First", ][[2]],
                 plottable(plot(n_resource_small, return_data = TRUE)[[2]]))
    p_first <- plot2(n_resource_small, n_resource_small, time = times[1])
    expect_equal(p_first$data[p_first$data$Model == "First", ][[2]],
                 plottable(plot(n_resource_small, time = times[1],
                                return_data = TRUE)[[2]]))

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

    # A length axis is now available, using the resource weight-length
    # parameters rather than a species'
    expect_s3_class(animate(n_resource_small, size_axis = "l"), "plotly")
})

test_that("plot() can express a resource density per logarithmic size", {
    resource <- initialNResource(NS_params_small)
    expect_identical(array_type(resource), "density")
    by_weight <- plot(resource, return_data = TRUE)

    # The resource has no length axis, but per log weight needs no
    # weight-length relationship and so is available
    per_log <- plot(resource, per_log_size = TRUE, return_data = TRUE)
    expect_equal(per_log[[2]], by_weight[[2]] * by_weight$w)
    expect_identical(plot(resource, per_log_size = TRUE)$scales$
                         get_scales("y")$name,
                     "Number density in log weight")
    expect_identical(plot(resource)$scales$get_scales("y")$name,
                     "Number density [1/g]")

    # A length axis is still refused, and so is per_log_size on a non-density
    expect_error(plot(resource_level(NS_params_small), per_log_size = TRUE),
                 "only applies to an array that holds a density")
})

test_that("the resource can be plotted against length", {
    params <- setResource(NS_params_small)
    resource <- initialNResource(params)
    rp <- resource_params(params)

    by_weight <- plot(resource, return_data = TRUE)
    by_length <- plot(resource, size_axis = "l", return_data = TRUE)
    expect_identical(names(by_length)[[1]], "l")
    expect_equal(by_length$l, (by_weight$w / rp$a)^(1 / rp$b))

    # It is a density, so the values are converted as well as the axis
    expect_equal(by_length[[2]],
                 by_weight[[2]] * rp$b * by_weight$w / by_length$l)
    p <- plot(resource, size_axis = "l")
    expect_identical(p$scales$get_scales("x")$name, "Length [cm]")
    expect_identical(p$scales$get_scales("y")$name, "Number density [1/cm]")

    # llim trims the length axis
    expect_lt(nrow(plot(resource, size_axis = "l", llim = c(0.1, 1),
                        return_data = TRUE)),
              nrow(by_length))

    # A quantity that is not a density keeps its values
    mort_w <- plot(getResourceMort(params), return_data = TRUE)
    mort_l <- plot(getResourceMort(params), size_axis = "l",
                   return_data = TRUE)
    expect_equal(mort_l[[2]], mort_w[[2]])

    # and animate() no longer refuses a length axis
    expect_s3_class(animate(NResource(NS_sim_small), size_axis = "l"), "plotly")
})

test_that("the resource joins the species on a length-axis spectrum", {
    params <- setResource(NS_params_small)
    on_l <- plotSpectra(params, size_axis = "l", return_data = TRUE)
    expect_true("Resource" %in% on_l$Legend)
    # The resource sits at its own convention, shorter than a fish of the same
    # weight by the ratio of the two `a` values
    on_w <- plotSpectra(params, return_data = TRUE)
    res_w <- on_w[on_w$Legend == "Resource", ]
    res_l <- on_l[on_l$Legend == "Resource", ]
    rp <- resource_params(params)
    expect_equal(res_l$l, (res_w$w / rp$a)^(1 / rp$b))
})

test_that("plot2, plotRelative and addPlot support length axis and per_log_size for resource arrays", {
    params <- setResource(NS_params_small)
    r1 <- initialNResource(params)
    r2 <- r1
    r2[] <- r1 * 1.5

    # ArrayResourceBySize methods
    p2_l <- plot2(r1, r2, size_axis = "l")
    expect_s3_class(p2_l, "ggplot")
    expect_identical(p2_l$scales$get_scales("x")$name, "Length [cm]")
    expect_identical(p2_l$scales$get_scales("y")$name, "Number density [1/cm]")

    p2_log <- plot2(r1, r2, per_log_size = TRUE)
    expect_s3_class(p2_log, "ggplot")
    expect_identical(p2_log$scales$get_scales("y")$name, "Number density in log weight")

    pr_l <- plotRelative(r1, r2, size_axis = "l")
    expect_s3_class(pr_l, "ggplot")
    expect_identical(pr_l$scales$get_scales("x")$name, "Length [cm]")

    pa_l <- addPlot(plot(r1, size_axis = "l"), r2, size_axis = "l")
    expect_s3_class(pa_l, "ggplot")
    expect_equal(length(pa_l$layers), 2)

    # ArrayTimeByResourceBySize methods
    nr1 <- NResource(NS_sim_small)
    nr2 <- nr1
    nr2[] <- nr1 * 1.2
    t2_l <- plot2(nr1, nr2, size_axis = "l")
    expect_s3_class(t2_l, "ggplot")
    expect_identical(t2_l$scales$get_scales("x")$name, "Length [cm]")

    tr_l <- plotRelative(nr1, nr2, size_axis = "l")
    expect_s3_class(tr_l, "ggplot")
    expect_identical(tr_l$scales$get_scales("x")$name, "Length [cm]")

    ta_l <- addPlot(plot(nr1, size_axis = "l"), nr2, size_axis = "l")
    expect_s3_class(ta_l, "ggplot")
    expect_equal(length(ta_l$layers), 2)
})


# Comparing two resource spectra across models ----

test_that("resource comparisons use each model's own length relationship", {
    params2 <- NS_params_small
    params2@resource_params$a <- resource_length_params(NS_params_small)$a / 8
    mort1 <- getResourceMort(NS_params_small)
    mort2 <- getResourceMort(params2)

    dat1 <- ArrayResourceBySize_plot_data(mort1, size_axis = "l")
    dat2 <- ArrayResourceBySize_plot_data(mort2, size_axis = "l")
    b <- resource_length_params(NS_params_small)$b
    expect_equal(dat2$l, dat1$l * 8^(1 / b))

    p2 <- plot2(mort1, mort2, size_axis = "l")
    expect_equal(sort(unique(p2$data$l[p2$data$Model == "Second"])),
                 sort(unique(dat2$l[dat2[[2]] > 0])))
})

test_that("plotRelative interpolates two resource spectra onto a common grid", {
    params2 <- NS_params_small
    params2@resource_params$a <- resource_length_params(NS_params_small)$a / 8
    mort1 <- getResourceMort(NS_params_small)
    mort2 <- getResourceMort(params2)

    rel <- plotRelative(mort1, mort2, size_axis = "l")
    expect_gt(nrow(rel$data), 10)
    dat1 <- ArrayResourceBySize_plot_data(mort1, size_axis = "l")
    dat2 <- ArrayResourceBySize_plot_data(mort2, size_axis = "l")
    expect_true(all(rel$data$l >= max(min(dat1$l), min(dat2$l))))
    expect_true(all(rel$data$l <= min(max(dat1$l), max(dat2$l))))
})

test_that("resource comparisons honour highlight", {
    mort <- getResourceMort(NS_params_small)
    doubled <- ArrayResourceBySize(unclass_resource(mort) * 2,
                                   value_name = attr(mort, "value_name"),
                                   units = attr(mort, "units"),
                                   params = NS_params_small)
    # A linear y axis, because the resource mortality vanishes at the top of
    # the grid and a logarithmic one would only warn about it.
    p2 <- plot2(mort, doubled, highlight = "Resource", log_y = FALSE)
    expect_equal(drawn_linewidth(p2, "Resource", NS_params_small), 1.6)
    plain <- plot2(mort, doubled, log_y = FALSE)
    expect_equal(drawn_linewidth(plain, "Resource", NS_params_small), 0.8)

    pr <- plotRelative(mort, doubled, highlight = "Resource")
    expect_equal(drawn_linewidth(pr, "Resource", NS_params_small, layer = 2),
                 1.6)
})

test_that("plot() and the comparisons share one resource preparation", {
    mort <- getResourceMort(NS_params_small)
    prepared <- ArrayResourceBySize_plot_data(mort, wlim = c(1e-5, NA),
                                              size_axis = "l")
    from_plot <- plot(mort, wlim = c(1e-5, NA), size_axis = "l",
                      return_data = TRUE)
    expect_equal(from_plot, prepared)
})
