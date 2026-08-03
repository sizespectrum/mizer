# Cache the most expensive repeated computation once at file load.
fmort_small <- getFMort(NS_sim_small)

test_that("ArrayTimeBySpeciesBySize constructor works", {
    arr <- array(1:200, dim = c(4, 5, 10))
    dimnames(arr) <- list(
        time = as.character(2000:2003),
        sp = letters[1:5],
        w = as.character(seq(1, 100, by = 10))
    )

    out <- ArrayTimeBySpeciesBySize(arr, value_name = "Test rate",
                                    units = "1/year")

    expect_s3_class(out, "ArrayTimeBySpeciesBySize")
    expect_true(is.ArrayTimeBySpeciesBySize(out))
    expect_false(is.ArrayTimeBySpeciesBySize(arr))
    expect_true(is.array(out))
    expect_identical(dim(out), dim(arr))
    expect_identical(attr(out, "value_name"), "Test rate")
    expect_identical(attr(out, "units"), "1/year")
})

test_that("ArrayTimeBySpeciesBySize constructor validates input", {
    expect_error(ArrayTimeBySpeciesBySize(matrix(1:6, 2, 3)), "`x` must be a 3D array.")
    expect_error(ArrayTimeBySpeciesBySize(1:10), "`x` must be a 3D array.")
})

test_that("Rate functions return ArrayTimeBySpeciesBySize from MizerSim", {
    expect_true(is.ArrayTimeBySpeciesBySize(fmort_small))
    expect_true(is.ArrayTimeBySpeciesBySize(getFeedingLevel(NS_sim_small)))
    expect_true(is.ArrayTimeBySpeciesBySize(getPredMort(NS_sim_small)))
})

test_that("Rate functions from MizerSim have correct dimnames", {
    expect_identical(dimnames(fmort_small)[[2]], NS_params_small@species_params$species)
    expect_identical(dimnames(fmort_small)[[3]], dimnames(NS_params_small@metab)[[2]])
    expect_identical(dim(fmort_small)[1], dim(NS_sim_small@n)[1])
})

test_that("getFMort with drop=TRUE and single time returns ArraySpeciesBySize", {
    times <- dimnames(NS_sim_small@n)$time
    f1 <- getFMort(NS_sim_small, time_range = times[3], drop = TRUE)
    expect_true(is.ArraySpeciesBySize(f1))
    expect_identical(attr(f1, "value_name"), "Fishing mortality")
})

test_that("print.ArrayTimeBySpeciesBySize works", {
    expect_output(print(fmort_small), "Fishing mortality")
    expect_output(print(fmort_small), "times x")
    expect_output(print(fmort_small), "1/year")
})

test_that("print.ArrayTimeBySpeciesBySize shows only the final time slice", {
    expect_output(print(fmort_small), "Showing final time step")
    final_time <- dimnames(fmort_small)[[1]][dim(fmort_small)[1]]
    expect_output(print(fmort_small), final_time)
})

test_that("summary.ArrayTimeBySpeciesBySize works", {
    s <- summary(fmort_small)
    expect_s3_class(s, "summary.ArrayTimeBySpeciesBySize")
    expect_identical(s$value_name, "Fishing mortality")
    expect_identical(nrow(s$per_species), nrow(NS_params_small@species_params))
    expect_output(print(s), "Fishing mortality")
    expect_output(print(s), "1/year")
    expect_output(print(s), "times x")
})

test_that("str.ArrayTimeBySpeciesBySize works", {
    expect_output(str(fmort_small), "ArrayTimeBySpeciesBySize")
    expect_output(str(fmort_small), "Fishing mortality")
    expect_output(str(fmort_small), "1/year")
    expect_output(str(fmort_small), "params")

    out <- capture.output(str(fmort_small))
    expect_false(any(grepl("intake_max", out)))
})

test_that("plot.ArrayTimeBySpeciesBySize supports base plot log argument", {
    p_y <- plot(fmort_small, log = "y")
    expect_s3_class(p_y, "ggplot")
    expect_identical(p_y$scales$get_scales("x")$trans$name, "identity")
    expect_identical(p_y$scales$get_scales("y")$trans$name, "log-10")

    p_xy <- plot(fmort_small, log = "xy")
    expect_identical(p_xy$scales$get_scales("x")$trans$name, "log-10")
    expect_identical(p_xy$scales$get_scales("y")$trans$name, "log-10")
})

test_that("plot.ArrayTimeBySpeciesBySize time argument selects correct slice", {
    times <- as.numeric(dimnames(fmort_small)[[1]])
    p <- plot(fmort_small, time = times[3], return_data = TRUE)
    expect_true(is.data.frame(p))
})

test_that("plot2.ArrayTimeBySpeciesBySize compares selected time slices", {
    times <- as.numeric(dimnames(fmort_small)[[1]])

    p <- plot2(fmort_small, fmort_small, name1 = "Original", name2 = "Changed",
               species = "Cod", time = times[3], total = TRUE,
               wlim = c(1, NA), log = "xy")
    expect_s3_class(p, "ggplot")
    expect_identical(levels(p$data$Model), c("Original", "Changed"))
    expect_true(all(p$data$Species %in% c("Cod", "Total")))
    expect_true(all(p$data$w >= 1))
    expect_identical(p$scales$get_scales("x")$trans$name, "log-10")
    expect_identical(p$scales$get_scales("y")$trans$name, "log-10")

    expect_error(plot2(fmort_small, getBiomass(NS_sim_small)), "Both objects must be")
})

test_that("plotRelative.ArrayTimeBySpeciesBySize compares selected time slices", {
    fmort2 <- fmort_small
    fmort2[] <- unclass(fmort_small) * 2
    times <- as.numeric(dimnames(fmort_small)[[1]])

    p <- plotRelative(fmort_small, fmort2, species = "Cod", time = times[3],
                      total = TRUE, wlim = c(1, NA))
    expect_s3_class(p, "ggplot")
    expect_true(all(p$data$Species %in% c("Cod", "Total")))
    expect_true(all(p$data$w >= 1))
    expect_true(all(abs(p$data$rel_diff - 2 / 3) < 1e-12))
    expect_identical(p$scales$get_scales("x")$trans$name, "log-10")

    expect_error(plotRelative(fmort_small, getBiomass(NS_sim_small)), "Both objects must be")
})

test_that("plot.ArrayTimeBySpeciesBySize preserves single species dimension", {
    arr <- array(seq_len(6), dim = c(2, 1, 3),
                 dimnames = list(time = c("2000", "2001"),
                                 sp = "Cod",
                                 w = c("1", "10", "100")))
    rate <- ArrayTimeBySpeciesBySize(arr, value_name = "Test rate",
                                     units = "1/year")

    p <- plot(rate, time = 2001, return_data = TRUE)

    expect_true(is.data.frame(p))
    expect_identical(p$Species, rep("Cod", 3))
    expect_equal(p[[2]], unname(arr["2001", "Cod", ]))
})

test_that("as.data.frame.ArrayTimeBySpeciesBySize returns correct structure", {
    df <- as.data.frame(fmort_small)
    expect_true(is.data.frame(df))
    expect_named(df, c("time", "Species", "w", "value"))
    expected_rows <- prod(dim(fmort_small))
    expect_equal(nrow(df), expected_rows)
})

test_that("[.ArrayTimeBySpeciesBySize preserves class for 3D result", {
    sub <- fmort_small[1:3, , ]
    expect_true(is.ArrayTimeBySpeciesBySize(sub))
    expect_identical(attr(sub, "value_name"), "Fishing mortality")
})

test_that("[.ArrayTimeBySpeciesBySize returns ArraySpeciesBySize when time is dropped", {
    slice <- fmort_small[1, , ]
    expect_false(is.ArrayTimeBySpeciesBySize(slice))
    expect_true(is.ArraySpeciesBySize(slice))
    expect_identical(attr(slice, "value_name"), "Fishing mortality")
})

test_that("[.ArrayTimeBySpeciesBySize returns ArrayTimeBySpecies when size is dropped", {
    slice <- fmort_small[, , 1]
    expect_false(is.ArrayTimeBySpeciesBySize(slice))
    expect_true(is.ArrayTimeBySpecies(slice))
    expect_identical(attr(slice, "value_name"), "Fishing mortality")
})

test_that("[.ArrayTimeBySpeciesBySize leaves time by size matrices plain", {
    slice <- fmort_small[, 1, ]
    expect_false(is.ArrayTimeBySpeciesBySize(slice))
    expect_false(is.ArraySpeciesBySize(slice))
    expect_false(is.ArrayTimeBySpecies(slice))
    expect_true(is.matrix(slice))
})

test_that("Ops.ArrayTimeBySpeciesBySize strips class for arithmetic", {
    result <- fmort_small * 2
    expect_false(is.ArrayTimeBySpeciesBySize(result))
    expect_true(is.array(result))
})

test_that("animate dispatches on ArrayTimeBySpeciesBySize", {
    p <- animate(fmort_small)
    expect_s3_class(p, "plotly")
})

test_that("animate.ArrayTimeBySpeciesBySize respects species argument", {
    p <- animate(fmort_small, species = c("Cod", "Herring"))
    expect_s3_class(p, "plotly")
    expect_error(animate(fmort_small, species = "NotASpecies"),
                 "None of the selected species are in the array.")
})

test_that("animate.ArrayTimeBySpeciesBySize respects time_range argument", {
    times <- as.numeric(dimnames(fmort_small)[[1]])
    p <- animate(fmort_small, time_range = c(times[1], times[4]))
    expect_s3_class(p, "plotly")
})

test_that("animate.ArrayTimeBySpeciesBySize sets axis ranges without dropping vertices", {
    p <- animate(fmort_small, species = "Cod", time_range = c(2000, 2001),
                 wlim = c(1, 1000), ylim = c(1e-3, 1))
    built_plot <- plotly::plotly_build(p)
    frame_lengths <- lengths(lapply(built_plot$x$frames, function(frame) {
        frame$data[[1]]$x
    }))

    expect_equal(built_plot$x$layout$xaxis$range, log10(c(1, 1000)))
    expect_equal(built_plot$x$layout$yaxis$range, log10(c(1e-3, 1)))
    expect_equal(frame_lengths, rep(dim(fmort_small)[[3]], length(frame_lengths)))
})

test_that("animate.ArrayTimeBySpeciesBySize can disable interpolation between frames", {
    p <- animate(fmort_small, time_range = c(2000, 2001), transition_duration = 0)
    expect_identical(p$animation$transition$duration, 0)
})

test_that("animate.ArrayTimeBySpeciesBySize exposes plotly animation timing controls", {
    p <- animate(fmort_small, time_range = c(2000, 2001),
                 frame_duration = 800,
                 transition_duration = 120,
                 easing = "cubic-in-out")

    expect_identical(p$animation$frame$duration, 800)
    expect_identical(p$animation$transition$duration, 120)
    expect_identical(p$animation$transition$easing, "cubic-in-out")
})

test_that("animate.ArrayTimeBySpeciesBySize respects total argument", {
    trace_names <- function(p) {
        vapply(plotly::plotly_build(p)$x$data, function(t) t$name, character(1))
    }
    expect_true("Total" %in% trace_names(animate(fmort_small, total = TRUE,
                                                  time_range = c(2000, 2001))))
    expect_false("Total" %in% trace_names(animate(fmort_small, total = FALSE,
                                                   time_range = c(2000, 2001))))
})

test_that("animate.ArrayTimeBySpeciesBySize respects background argument", {
    # Use Herring so the background species has non-zero fishing mortality
    params_bkgrd <- markBackground(NS_params_small, species = "Herring")
    sim_bkgrd <- project(params_bkgrd, t_max = 5, t_save = 1)
    fmort_bkgrd <- getFMort(sim_bkgrd)
    trace_names <- function(p) {
        vapply(plotly::plotly_build(p)$x$data, function(t) t$name, character(1))
    }
    expect_true("Background" %in% trace_names(animate(fmort_bkgrd, background = TRUE)))
    expect_false("Background" %in% trace_names(animate(fmort_bkgrd, background = FALSE)))
})

test_that("animateSpectra is a backward-compatible alias for animate", {
    expect_s3_class(animateSpectra(fmort_small), "plotly")
    expect_s3_class(animateSpectra(NS_sim_small), "plotly")
})

# Second-order plotting placement (#382) -----------------------------------

test_that("the time-series size axis honours representation under second_order_w", {
    # getFMort.MizerSim tags its result "average".
    expect_identical(attr(fmort_small, "representation"), "average")
    p <- NS_sim_small@params
    # Default first-order: size axis stays on the nodes (df unchanged). The
    # stored dimnames are rounded to 3 sig figs, so compare with tolerance.
    df_off <- as.data.frame(fmort_small)
    expect_equal(sort(unique(df_off$w)), sort(unname(p@w)),
                 tolerance = 1e-2)

    # Switch the model to second order and recompute on the sim.
    sim2 <- NS_sim_small
    second_order_w(sim2@params) <- c(bin_average = TRUE)
    fmort2 <- getFMort(sim2)
    expect_equal(get_ArrayTimeBySpeciesBySize_w(fmort2), bin_midpoints(p))
    df_on <- as.data.frame(fmort2)
    expect_equal(sort(unique(df_on$w)), sort(unname(bin_midpoints(p))))
})

test_that("a time slice keeps the representation tag", {
    sim2 <- NS_sim_small
    second_order_w(sim2@params) <- c(bin_average = TRUE)
    fmort2 <- getFMort(sim2)
    slice <- ArrayTimeBySpeciesBySize_slice(fmort2)
    expect_s3_class(slice, "ArraySpeciesBySize")
    expect_identical(attr(slice, "representation"), "average")
    expect_equal(get_ArraySpeciesBySize_w(slice), bin_midpoints(sim2@params))
})
