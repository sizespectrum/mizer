# Plot fixtures use edition 1, as do the existing visual baselines.
withr::local_options(mizer_defaults_edition = 1)
params <- plot_params_small
species <- c(1, 2)

test_that("proportion_ylim widens to [0, 1] without ever hiding data", {
    # The whole of [0, 1] is shown even for values that use little of it
    expect_identical(proportion_ylim(c(NA, NA), FALSE, c(0.2, 0.3)), c(0, 1))
    # and the range is widened, never narrowed, for values outside it
    expect_identical(proportion_ylim(c(NA, NA), FALSE, c(0.2, 2.5)), c(0, 2.5))
    expect_identical(proportion_ylim(c(NA, NA), FALSE, c(-0.5, 0.3)), c(-0.5, 1))
    # An explicit limit always wins, end by end
    expect_identical(proportion_ylim(c(0.1, 0.5), FALSE, c(0.2, 0.3)),
                     c(0.1, 0.5))
    expect_identical(proportion_ylim(c(NA, 0.5), FALSE, c(0.2, 0.3)), c(0, 0.5))
    expect_identical(proportion_ylim(c(0.1, NA), FALSE, c(0.2, 0.3)), c(0.1, 1))
    # A logarithmic axis has no place for the 0, so it is left alone
    expect_identical(proportion_ylim(c(NA, NA), TRUE, c(0.2, 0.3)), c(NA, NA))
    # Non-finite values are ignored rather than swallowing the range
    expect_identical(proportion_ylim(c(NA, NA), FALSE, c(0.2, NA, Inf)),
                     c(0, 1))
})

test_that("add_total_line sums series on a shared grid exactly", {
    plot_dat <- data.frame(w = rep(c(1, 2, 4), 2), value = 1:6,
                           Species = rep(c("a", "b"), each = 3))
    out <- add_total_line(plot_dat)
    total <- out[out$Species == "Total", ]
    expect_identical(total$w, c(1, 2, 4))
    expect_equal(total$value, c(1 + 4, 2 + 5, 3 + 6))
    # The original rows are untouched
    expect_equal(out[out$Species != "Total", names(plot_dat)], plot_dat,
                 ignore_attr = TRUE)
})

test_that("add_total_line interpolates series that do not share a grid", {
    # `b` lives on a shifted grid, so the total is evaluated on the union and
    # each series is interpolated onto it, contributing nothing outside its own
    # range.
    plot_dat <- rbind(
        data.frame(w = c(1, 4), value = c(1, 4), Species = "a"),
        data.frame(w = c(2, 8), value = c(10, 40), Species = "b"))
    total <- add_total_line(plot_dat)
    total <- total[total$Species == "Total", ]
    expect_identical(total$w, c(1, 2, 4, 8))
    # At w = 1 only `a`; at w = 8 only `b`
    expect_equal(total$value[[1]], 1)
    expect_equal(total$value[[4]], 40)
    # At w = 2, `a` interpolated linearly in log w between (1, 1) and (4, 4)
    a_at_2 <- stats::approx(log(c(1, 4)), c(1, 4), xout = log(2))$y
    expect_equal(total$value[[2]], a_at_2 + 10)
    # At w = 4, `b` interpolated between (2, 10) and (8, 40)
    b_at_4 <- stats::approx(log(c(2, 8)), c(10, 40), xout = log(4))$y
    expect_equal(total$value[[3]], 4 + b_at_4)
})

test_that("add_total_line groups by the columns it is told to", {
    plot_dat <- data.frame(w = rep(c(1, 2), 4), value = 1:8,
                           Species = rep(c("a", "b"), each = 2),
                           time = rep(c(1, 1, 2, 2), each = 2))
    out <- add_total_line(plot_dat, by = "time")
    totals <- out[out$Species == "Total", ]
    expect_identical(nrow(totals), 4L)
    # time 1 holds a = (1, 2) and b = (3, 4); time 2 holds a = (5, 6) and
    # b = (7, 8)
    expect_equal(totals$value[totals$time == 1], c(1 + 3, 2 + 4))
    expect_equal(totals$value[totals$time == 2], c(5 + 7, 6 + 8))
})

test_that("validate_density_wrt accepts only the known measures", {
    expect_identical(validate_density_wrt(NULL), NA_character_)
    expect_identical(validate_density_wrt(NA), NA_character_)
    expect_identical(validate_density_wrt(NA_character_), NA_character_)
    for (measure in density_measures) {
        expect_identical(validate_density_wrt(measure), measure)
    }
    expect_error(validate_density_wrt("weight"), "must be NA or one of")
    expect_error(validate_density_wrt(c("w", "l")), "must be NA or one of")
    expect_error(validate_density_wrt(1), "must be NA or one of")
})

test_that("the target density measure follows the size axis", {
    # The size variable follows the axis, but whether the values are per size
    # or per logarithmic size is a property of the values and is preserved.
    expect_identical(density_target_measure("w", "l"), "l")
    expect_identical(density_target_measure("l", "w"), "w")
    expect_identical(density_target_measure("log_w", "l"), "log_l")
    expect_identical(density_target_measure("log_l", "w"), "log_w")
    expect_identical(density_target_measure("w", "w"), "w")
    expect_identical(density_target_measure("log_l", "l"), "log_l")
    # Values that are not a density have no target measure
    expect_identical(density_target_measure(NA, "l"), NA_character_)
    expect_identical(density_target_measure(NULL, "l"), NA_character_)
})

test_that("density Jacobians are correct, reciprocal and transitive", {
    w <- c(1, 10, 100)
    l <- c(2, 5, 20)
    b <- c(3, 3.1, 2.9)
    expect_equal(density_measure_jacobian("w", "l", w, l, b), b * w / l)
    expect_equal(density_measure_jacobian("l", "w", w, l, b), l / (b * w))
    expect_equal(density_measure_jacobian("log_w", "log_l", w, l, b), b)
    expect_equal(density_measure_jacobian("log_l", "log_w", w, l, b), 1 / b)
    for (measure in density_measures) {
        expect_equal(density_measure_jacobian(measure, measure, w, l, b),
                     rep(1, 3))
    }
    # Converting there and back is the identity
    for (from in density_measures) {
        for (to in density_measures) {
            expect_equal(density_measure_jacobian(from, to, w, l, b) *
                             density_measure_jacobian(to, from, w, l, b),
                         rep(1, 3))
        }
    }
    # A density per log size is the size times the density per size
    expect_equal(density_measure_jacobian("w", "log_w", w, l, b), w)
    expect_equal(density_measure_jacobian("l", "log_l", w, l, b), l)
})

test_that("convert_density_units swaps the size unit", {
    expect_identical(convert_density_units("1/g", "w", "l"), "1/cm")
    expect_identical(convert_density_units("g^-1/year", "w", "l"),
                     "cm^-1/year")
    expect_identical(convert_density_units("1/cm", "l", "w"), "1/g")
    expect_identical(convert_density_units("cm^-1/year", "l", "w"),
                     "g^-1/year")
    # No conversion, so no change
    expect_identical(convert_density_units("1/g", "w", "w"), "1/g")
    expect_identical(convert_density_units("g/year", NA_character_,
                                           NA_character_), "g/year")
    # Nothing declared falls back to the canonical units of the target
    expect_identical(convert_density_units(NULL, "w", "l"), "1/cm")
    expect_identical(convert_density_units("", "l", "w"), "1/g")
    expect_null(convert_density_units(NULL, "log_w", "log_l"))
    # A per log size density carries no size unit, so nothing is swapped
    expect_identical(convert_density_units("1/year", "log_w", "log_l"),
                     "1/year")
    # Units with no per-size factor to identify are left alone
    expect_identical(convert_density_units("g/year", "w", "l"), "g/year")
})

test_that("convert_plot_density_axis converts in both directions", {
    sp <- params@species_params
    w <- params@w[1:5]
    plot_dat <- data.frame(w = rep(w, 2), value = seq_len(10),
                           Species = rep(sp$species[1:2], each = 5))
    sp_idx <- match(as.character(plot_dat$Species), as.character(sp$species))
    l <- w2l(plot_dat$w, sp[sp_idx, ])
    b <- unname(sp$b[sp_idx])
    jacobian <- b * plot_dat$w / l

    # A density with respect to weight, plotted against length
    converted <- convert_plot_density_axis(plot_dat, params, "l",
                                           density_wrt = "w")
    expect_identical(names(converted), c("l", "value", "Species"))
    expect_equal(converted$l, l)
    expect_equal(converted$value, plot_dat$value * jacobian)

    # A density with respect to length, plotted against weight, is the inverse
    length_dat <- plot_dat
    length_dat$value <- converted$value
    back <- convert_plot_density_axis(length_dat, params, "w",
                                      density_wrt = "l")
    expect_identical(names(back), c("w", "value", "Species"))
    expect_equal(back$value, plot_dat$value)

    # Densities per logarithmic size use the b Jacobian
    log_converted <- convert_plot_density_axis(plot_dat, params, "l",
                                               density_wrt = "log_w")
    expect_equal(log_converted$value, plot_dat$value * b)
    log_dat <- plot_dat
    log_dat$value <- log_converted$value
    log_back <- convert_plot_density_axis(log_dat, params, "w",
                                          density_wrt = "log_l")
    expect_equal(log_back$value, plot_dat$value)

    # Values that are not a density only have their size coordinate converted
    not_density <- convert_plot_density_axis(plot_dat, params, "l")
    expect_equal(not_density$value, plot_dat$value)
    expect_equal(not_density$l, l)
    # and are untouched on a weight axis
    expect_identical(convert_plot_density_axis(plot_dat, params, "w"),
                     plot_dat)
    expect_identical(convert_plot_density_axis(plot_dat, params, "w",
                                               density_wrt = "log_w"),
                     plot_dat)

    # The value column can be named or indexed
    expect_equal(convert_plot_density_axis(plot_dat, params, "l",
                                           density_wrt = "w",
                                           value_col = "value"),
                 converted)

    # The Jacobian is a species property, so rows with no species are dropped
    with_total <- rbind(plot_dat,
                        data.frame(w = w, value = 1, Species = "Total"))
    expect_identical(nrow(convert_plot_density_axis(with_total, params, "w",
                                                    density_wrt = "l")),
                     nrow(plot_dat))
})

test_that("per_log_size chooses the measure independently of the size axis", {
    # The size variable follows the axis, the log-ness follows per_log_size
    expect_identical(density_target_measure("w", "w", TRUE), "log_w")
    expect_identical(density_target_measure("w", "l", TRUE), "log_l")
    expect_identical(density_target_measure("log_w", "w", FALSE), "w")
    expect_identical(density_target_measure("log_w", "l", FALSE), "l")
    # NULL keeps whichever the values already are
    expect_identical(density_target_measure("w", "l", NULL), "l")
    expect_identical(density_target_measure("log_w", "l", NULL), "log_l")
    # Still nothing to do for values that are not a density
    expect_identical(density_target_measure(NA, "w", TRUE), NA_character_)
})

test_that("convert_plot_density_axis expresses a density per logarithmic size", {
    sp <- params@species_params
    w <- params@w[1:5]
    plot_dat <- data.frame(w = rep(w, 2), value = seq_len(10),
                           Species = rep(sp$species[1:2], each = 5))
    sp_idx <- match(as.character(plot_dat$Species), as.character(sp$species))
    b <- unname(sp$b[sp_idx])
    l <- w2l(plot_dat$w, sp[sp_idx, ])

    # Per log weight is the density times the weight
    per_log_w <- convert_plot_density_axis(plot_dat, params, "w",
                                           density_wrt = "w",
                                           per_log_size = TRUE)
    expect_identical(names(per_log_w), c("w", "value", "Species"))
    expect_equal(per_log_w$value, plot_dat$value * plot_dat$w)

    # Per log length is the density per length times the length
    per_log_l <- convert_plot_density_axis(plot_dat, params, "l",
                                           density_wrt = "w",
                                           per_log_size = TRUE)
    expect_identical(names(per_log_l), c("l", "value", "Species"))
    expect_equal(per_log_l$value, plot_dat$value * b * plot_dat$w)
    per_l <- convert_plot_density_axis(plot_dat, params, "l",
                                       density_wrt = "w")
    expect_equal(per_log_l$value, per_l$value * l)

    # per_log_size = FALSE undoes it again
    log_dat <- plot_dat
    log_dat$value <- per_log_w$value
    expect_equal(convert_plot_density_axis(log_dat, params, "w",
                                           density_wrt = "log_w",
                                           per_log_size = FALSE)$value,
                 plot_dat$value)

    # Going to a logarithmic measure on a weight axis needs no length, so it
    # keeps rows whose species has no weight-length relationship
    with_total <- rbind(plot_dat,
                        data.frame(w = w, value = 1, Species = "Total"))
    kept <- convert_plot_density_axis(with_total, params, "w",
                                      density_wrt = "w", per_log_size = TRUE)
    expect_identical(nrow(kept), nrow(with_total))
    expect_equal(kept$value, with_total$value * with_total$w)
})

test_that("convert_density_units drops the size unit for a logarithmic measure", {
    # A number per log weight interval is just a number
    expect_identical(convert_density_units("1/g", "w", "log_w"), "")
    expect_identical(convert_density_units("1/cm", "l", "log_l"), "")
    # and the rest of a compound unit survives
    expect_identical(convert_density_units("g^-1/year", "w", "log_w"),
                     "1/year")
    # Units with no per-size factor to remove are left alone
    expect_identical(convert_density_units("g/year", "w", "log_w"), "g/year")
})

test_that("plot_y_is_log reads the transformation off the plot", {
    bio <- getBiomass(NS_sim_small)
    expect_true(plot_y_is_log(plot(bio, log_y = TRUE)))
    expect_false(plot_y_is_log(plot(bio, log_y = FALSE)))
    # A plot with no y scale of its own is linear, ggplot2's default.
    bare <- ggplot2::ggplot(data.frame(x = 1, y = 1),
                            ggplot2::aes(x = x, y = y)) + ggplot2::geom_point()
    expect_false(plot_y_is_log(bare))
})
