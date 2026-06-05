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
    expect_true(is.ArrayTimeBySpeciesBySize(getFMort(NS_sim)))
    expect_true(is.ArrayTimeBySpeciesBySize(getFeedingLevel(NS_sim)))
    expect_true(is.ArrayTimeBySpeciesBySize(getPredMort(NS_sim)))
})

test_that("Rate functions from MizerSim have correct dimnames", {
    fmort <- getFMort(NS_sim)
    expect_identical(dimnames(fmort)[[2]], NS_params@species_params$species)
    expect_identical(dimnames(fmort)[[3]], dimnames(NS_params@metab)[[2]])
    expect_identical(dim(fmort)[1], dim(NS_sim@n)[1])
})

test_that("getFMort with drop=TRUE and single time returns ArraySpeciesBySize", {
    times <- dimnames(NS_sim@n)$time
    f1 <- getFMort(NS_sim, time_range = times[10], drop = TRUE)
    expect_true(is.ArraySpeciesBySize(f1))
    expect_identical(attr(f1, "value_name"), "Fishing mortality")
})

test_that("print.ArrayTimeBySpeciesBySize works", {
    fmort <- getFMort(NS_sim)
    expect_output(print(fmort), "Fishing mortality")
    expect_output(print(fmort), "times x")
    expect_output(print(fmort), "1/year")
})

test_that("summary.ArrayTimeBySpeciesBySize works", {
    fmort <- getFMort(NS_sim)
    s <- summary(fmort)
    expect_s3_class(s, "summary.ArrayTimeBySpeciesBySize")
    expect_identical(s$value_name, "Fishing mortality")
    expect_identical(nrow(s$per_species), nrow(NS_params@species_params))
    expect_output(print(s), "Fishing mortality")
    expect_output(print(s), "1/year")
    expect_output(print(s), "times x")
})

test_that("str.ArrayTimeBySpeciesBySize works", {
    fmort <- getFMort(NS_sim)
    expect_output(str(fmort), "ArrayTimeBySpeciesBySize")
    expect_output(str(fmort), "Fishing mortality")
    expect_output(str(fmort), "1/year")
    expect_output(str(fmort), "params")
    
    out <- capture.output(str(fmort))
    expect_false(any(grepl("intake_max", out)))
})

test_that("plot.ArrayTimeBySpeciesBySize returns a ggplot", {
    fmort <- getFMort(NS_sim)
    p <- plot(fmort)
    expect_s3_class(p, "ggplot")
})

test_that("plot.ArrayTimeBySpeciesBySize supports base plot log argument", {
    fmort <- getFMort(NS_sim)

    p_y <- plot(fmort, log = "y")
    expect_identical(p_y$scales$get_scales("x")$trans$name, "identity")
    expect_identical(p_y$scales$get_scales("y")$trans$name, "log-10")

    p_xy <- plot(fmort, log = "xy")
    expect_identical(p_xy$scales$get_scales("x")$trans$name, "log-10")
    expect_identical(p_xy$scales$get_scales("y")$trans$name, "log-10")
})

test_that("plot.ArrayTimeBySpeciesBySize time argument selects correct slice", {
    fmort <- getFMort(NS_sim)
    times <- as.numeric(dimnames(fmort)[[1]])
    p <- plot(fmort, time = times[5], return_data = TRUE)
    expect_true(is.data.frame(p))
})

test_that("plot2.ArrayTimeBySpeciesBySize compares selected time slices", {
    fmort <- getFMort(NS_sim)
    times <- as.numeric(dimnames(fmort)[[1]])

    p <- plot2(fmort, fmort, name1 = "Original", name2 = "Changed",
               species = "Cod", time = times[5], total = TRUE,
               wlim = c(1, NA), log = "xy")
    expect_s3_class(p, "ggplot")
    expect_identical(levels(p$data$Model), c("Original", "Changed"))
    expect_true(all(p$data$Species %in% c("Cod", "Total")))
    expect_true(all(p$data$w >= 1))
    expect_identical(p$scales$get_scales("x")$trans$name, "log-10")
    expect_identical(p$scales$get_scales("y")$trans$name, "log-10")

    expect_error(plot2(fmort, getBiomass(NS_sim)), "Both objects must be")
})

test_that("plotRelative.ArrayTimeBySpeciesBySize compares selected time slices", {
    fmort <- getFMort(NS_sim)
    fmort2 <- fmort
    fmort2[] <- unclass(fmort) * 2
    times <- as.numeric(dimnames(fmort)[[1]])

    p <- plotRelative(fmort, fmort2, species = "Cod", time = times[5],
                      total = TRUE, wlim = c(1, NA))
    expect_s3_class(p, "ggplot")
    expect_true(all(p$data$Species %in% c("Cod", "Total")))
    expect_true(all(p$data$w >= 1))
    expect_true(all(abs(p$data$rel_diff - 2 / 3) < 1e-12))
    expect_identical(p$scales$get_scales("x")$trans$name, "log-10")

    expect_error(plotRelative(fmort, getBiomass(NS_sim)), "Both objects must be")
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
    fmort <- getFMort(NS_sim)
    df <- as.data.frame(fmort)
    expect_true(is.data.frame(df))
    expect_named(df, c("time", "Species", "w", "value"))
    expected_rows <- prod(dim(fmort))
    expect_equal(nrow(df), expected_rows)
})

test_that("[.ArrayTimeBySpeciesBySize preserves class for 3D result", {
    fmort <- getFMort(NS_sim)
    sub <- fmort[1:3, , ]
    expect_true(is.ArrayTimeBySpeciesBySize(sub))
    expect_identical(attr(sub, "value_name"), "Fishing mortality")
})

test_that("[.ArrayTimeBySpeciesBySize returns ArraySpeciesBySize when time is dropped", {
    fmort <- getFMort(NS_sim)
    slice <- fmort[1, , ]
    expect_false(is.ArrayTimeBySpeciesBySize(slice))
    expect_true(is.ArraySpeciesBySize(slice))
    expect_identical(attr(slice, "value_name"), "Fishing mortality")
})

test_that("[.ArrayTimeBySpeciesBySize returns ArrayTimeBySpecies when size is dropped", {
    fmort <- getFMort(NS_sim)
    slice <- fmort[, , 1]
    expect_false(is.ArrayTimeBySpeciesBySize(slice))
    expect_true(is.ArrayTimeBySpecies(slice))
    expect_identical(attr(slice, "value_name"), "Fishing mortality")
})

test_that("[.ArrayTimeBySpeciesBySize leaves time by size matrices plain", {
    fmort <- getFMort(NS_sim)
    slice <- fmort[, 1, ]
    expect_false(is.ArrayTimeBySpeciesBySize(slice))
    expect_false(is.ArraySpeciesBySize(slice))
    expect_false(is.ArrayTimeBySpecies(slice))
    expect_true(is.matrix(slice))
})

test_that("Ops.ArrayTimeBySpeciesBySize strips class for arithmetic", {
    fmort <- getFMort(NS_sim)
    result <- fmort * 2
    expect_false(is.ArrayTimeBySpeciesBySize(result))
    expect_true(is.array(result))
})

test_that("animate dispatches on ArrayTimeBySpeciesBySize", {
    fmort <- getFMort(NS_sim)
    p <- animate(fmort)
    expect_s3_class(p, "plotly")
})

test_that("animate.ArrayTimeBySpeciesBySize respects species argument", {
    fmort <- getFMort(NS_sim)
    p <- animate(fmort, species = c("Cod", "Herring"))
    expect_s3_class(p, "plotly")
    expect_error(animate(fmort, species = "NotASpecies"),
                 "None of the selected species are in the array.")
})

test_that("animate.ArrayTimeBySpeciesBySize respects time_range argument", {
    fmort <- getFMort(NS_sim)
    times <- as.numeric(dimnames(fmort)[[1]])
    p <- animate(fmort, time_range = c(times[1], times[10]))
    expect_s3_class(p, "plotly")
})

test_that("animate.ArrayTimeBySpeciesBySize sets axis ranges without dropping vertices", {
    fmort <- getFMort(NS_sim)
    p <- animate(fmort, species = "Cod", time_range = c(2000, 2001),
                 wlim = c(1, 1000), ylim = c(1e-3, 1))
    built_plot <- plotly::plotly_build(p)
    frame_lengths <- lengths(lapply(built_plot$x$frames, function(frame) {
        frame$data[[1]]$x
    }))

    expect_equal(built_plot$x$layout$xaxis$range, log10(c(1, 1000)))
    expect_equal(built_plot$x$layout$yaxis$range, log10(c(1e-3, 1)))
    expect_equal(frame_lengths, rep(dim(fmort)[[3]], length(frame_lengths)))
})

test_that("animate.ArrayTimeBySpeciesBySize can disable interpolation between frames", {
    fmort <- getFMort(NS_sim)
    p <- animate(fmort, time_range = c(2000, 2001), transition_duration = 0)
    expect_identical(p$animation$transition$duration, 0)
})

test_that("animate.ArrayTimeBySpeciesBySize exposes plotly animation timing controls", {
    fmort <- getFMort(NS_sim)
    p <- animate(fmort, time_range = c(2000, 2001),
                 frame_duration = 800,
                 transition_duration = 120,
                 easing = "cubic-in-out")

    expect_identical(p$animation$frame$duration, 800)
    expect_identical(p$animation$transition$duration, 120)
    expect_identical(p$animation$transition$easing, "cubic-in-out")
})

test_that("animate.ArrayTimeBySpeciesBySize respects total argument", {
    fmort <- getFMort(NS_sim)
    trace_names <- function(p) {
        vapply(plotly::plotly_build(p)$x$data, function(t) t$name, character(1))
    }
    expect_true("Total" %in% trace_names(animate(fmort, total = TRUE,
                                                  time_range = c(2000, 2001))))
    expect_false("Total" %in% trace_names(animate(fmort, total = FALSE,
                                                   time_range = c(2000, 2001))))
})

test_that("animate.ArrayTimeBySpeciesBySize respects background argument", {
    # Use Herring so the background species has non-zero fishing mortality
    params_bkgrd <- markBackground(NS_params, species = "Herring")
    sim_bkgrd <- project(params_bkgrd, t_max = 5, t_save = 1)
    fmort <- getFMort(sim_bkgrd)
    trace_names <- function(p) {
        vapply(plotly::plotly_build(p)$x$data, function(t) t$name, character(1))
    }
    expect_true("Background" %in% trace_names(animate(fmort, background = TRUE)))
    expect_false("Background" %in% trace_names(animate(fmort, background = FALSE)))
})

test_that("animateSpectra is a backward-compatible alias for animate", {
    fmort <- getFMort(NS_sim)
    expect_s3_class(animateSpectra(fmort), "plotly")
    expect_s3_class(animateSpectra(NS_sim), "plotly")
})
