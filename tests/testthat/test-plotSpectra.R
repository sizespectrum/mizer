# Plot fixtures use edition 1, as do the existing visual baselines.
withr::local_options(mizer_defaults_edition = 1)
params <- plot_params_small
sim <- plot_sim_small
sim0 <- plot_sim0_small
species <- c(1, 2)
params_bkgrd <- plot_params_background_small
params_single <- plot_params_single_small

test_that("plotSpectra2 compares spectra from params and sims", {
    p_params <- plotSpectra2(params, params, name1 = "Original",
                             name2 = "Changed", species = species,
                             total = TRUE)
    expect_s3_class(p_params, "ggplot")
    expect_identical(levels(p_params$data$Model), c("Original", "Changed"))
    expect_true("Total" %in% p_params$data$Legend)

    expect_s3_class(plotSpectra2(sim, sim0, species = species), "ggplot")
    expect_s3_class(plotSpectra2(params, sim, species = species), "ggplot")
})

test_that("plotSpectra2 supports base plot log argument", {
    p_y <- plotSpectra2(params, sim0, species = species, log = "y")
    expect_identical(p_y$scales$get_scales("x")$trans$name, "identity")
    expect_identical(p_y$scales$get_scales("y")$trans$name, "log-10")

    p_xy <- plotSpectra2(params, sim0, species = species, log = "xy")
    expect_identical(p_xy$scales$get_scales("x")$trans$name, "log-10")
    expect_identical(p_xy$scales$get_scales("y")$trans$name, "log-10")

    p_none <- plotSpectra2(params, sim0, species = species, log = "")
    expect_identical(p_none$scales$get_scales("x")$trans$name, "identity")
    expect_identical(p_none$scales$get_scales("y")$trans$name, "identity")

    p_flags <- plotSpectra2(params, sim0, species = species,
                            log_x = FALSE, log_y = FALSE)
    expect_identical(p_flags$scales$get_scales("x")$trans$name, "identity")
    expect_identical(p_flags$scales$get_scales("y")$trans$name, "identity")

    expect_error(plotSpectra2(params, sim0, species = species, log = "z"),
                 "containing only")
})

test_that("plotSpectraRelative plots symmetric relative difference", {
    params2 <- params
    params2@initial_n[] <- params@initial_n * 2

    p <- plotSpectraRelative(params, params2, species = species,
                             resource = FALSE)
    expect_s3_class(p, "ggplot")
    expect_true(all(abs(p$data$rel_diff - 2 / 3) < 1e-12))
    expect_identical(p$scales$get_scales("x")$trans$name, "log-10")

    p_linear <- plotSpectraRelative(params, params2, species = species,
                                    resource = FALSE, log_x = FALSE)
    expect_identical(p_linear$scales$get_scales("x")$trans$name, "identity")

    expect_s3_class(plotSpectraRelative(sim, sim0, species = species,
                                        resource = FALSE), "ggplot")
    expect_s3_class(plotSpectraRelative(params, sim, species = species,
                                        resource = FALSE), "ggplot")
})

test_that("resolve_spectrum_power splits the power into its two factors", {
    expect_identical(resolve_spectrum_power(),
                     list(power = 1, biomass = TRUE, per_log_size = FALSE))
    expect_identical(resolve_spectrum_power(biomass = FALSE),
                     list(power = 0, biomass = FALSE, per_log_size = FALSE))
    expect_identical(resolve_spectrum_power(per_log_size = TRUE),
                     list(power = 2, biomass = TRUE, per_log_size = TRUE))
    expect_identical(resolve_spectrum_power(biomass = FALSE,
                                            per_log_size = TRUE),
                     list(power = 1, biomass = FALSE, per_log_size = TRUE))

    # `power` on its own keeps the interpretation mizer has always given it
    expect_identical(resolve_spectrum_power(power = 0),
                     list(power = 0, biomass = FALSE, per_log_size = FALSE))
    expect_identical(resolve_spectrum_power(power = 1),
                     list(power = 1, biomass = TRUE, per_log_size = FALSE))
    expect_identical(resolve_spectrum_power(power = 2),
                     list(power = 2, biomass = TRUE, per_log_size = TRUE))
    # ... including for powers that are not the sum of the two flags
    expect_identical(resolve_spectrum_power(power = 1.5),
                     list(power = 1.5, biomass = TRUE, per_log_size = FALSE))

    # `power` together with one flag determines the other
    expect_identical(resolve_spectrum_power(power = 1, biomass = FALSE),
                     list(power = 1, biomass = FALSE, per_log_size = TRUE))
    expect_identical(resolve_spectrum_power(power = 1, per_log_size = TRUE),
                     list(power = 1, biomass = FALSE, per_log_size = TRUE))
    expect_identical(resolve_spectrum_power(power = 2, biomass = TRUE),
                     list(power = 2, biomass = TRUE, per_log_size = TRUE))

    # Contradictions are errors rather than silently ignored arguments
    expect_error(resolve_spectrum_power(power = 0, biomass = TRUE),
                 "but `biomass = TRUE` and `per_log_size = FALSE` imply")
    expect_error(resolve_spectrum_power(power = 1, biomass = TRUE,
                                        per_log_size = TRUE),
                 "imply `power = 2`")
    expect_error(resolve_spectrum_power(power = 1.5, biomass = TRUE),
                 "imply `power = 1`")
})

test_that("plotSpectra flags agree with the corresponding powers", {
    for (flags in list(list(biomass = FALSE, per_log_size = FALSE, power = 0),
                       list(biomass = TRUE, per_log_size = FALSE, power = 1),
                       list(biomass = FALSE, per_log_size = TRUE, power = 1),
                       list(biomass = TRUE, per_log_size = TRUE, power = 2))) {
        by_flags <- plotSpectra(params, species = species, resource = FALSE,
                                biomass = flags$biomass,
                                per_log_size = flags$per_log_size,
                                return_data = TRUE)
        by_power <- plotSpectra(params, species = species, resource = FALSE,
                                power = flags$power, return_data = TRUE)
        expect_equal(by_flags[[2]], by_power[[2]], ignore_attr = TRUE)
        expect_equal(by_flags$w, by_power$w)
    }
})

test_that("the y-axis label distinguishes the two spectra with power 1", {
    biomass_density <- plotSpectra(params, species = species,
                                   resource = FALSE, biomass = TRUE)
    expect_identical(biomass_density$scales$get_scales("y")$name,
                     "Biomass density")
    per_log <- plotSpectra(params, species = species, resource = FALSE,
                           biomass = FALSE, per_log_size = TRUE)
    expect_identical(per_log$scales$get_scales("y")$name,
                     "Number density in log weight")
})

test_that("a density in log size uses the logarithmic Jacobian on a length axis", {
    params_len <- params
    params_len@species_params$a <- 0.01
    params_len@species_params$b <- 3

    by_weight <- plotSpectra(params_len, species = species, resource = FALSE,
                             biomass = FALSE, per_log_size = TRUE,
                             return_data = TRUE)
    by_length <- plotSpectra(params_len, species = species, resource = FALSE,
                             biomass = FALSE, per_log_size = TRUE,
                             size_axis = "l", return_data = TRUE)
    sp_idx <- match(by_length$Species, params_len@species_params$species)
    b <- params_len@species_params$b[sp_idx]
    # d log(w) / d log(l) = b, not the density Jacobian b * w / l that the
    # same power = 1 spectrum gets when it is a biomass density
    expect_equal(by_length[[2]], unname(by_weight[[2]] * b))
    expect_identical(names(by_length)[[2]], "Number density in log length")

    biomass_l <- plotSpectra(params_len, species = species, resource = FALSE,
                             biomass = TRUE, size_axis = "l",
                             return_data = TRUE)
    expect_false(isTRUE(all.equal(biomass_l[[2]], by_length[[2]],
                                  check.attributes = FALSE)))
})

test_that("plotSpectra2 honours the biomass and per_log_size flags", {
    compared <- plotSpectra2(params, params, species = species,
                             resource = FALSE, biomass = FALSE,
                             per_log_size = TRUE)
    expect_identical(compared$scales$get_scales("y")$name,
                     "Number density in log weight")
    single <- plotSpectra(params, species = species, resource = FALSE,
                          biomass = FALSE, per_log_size = TRUE,
                          return_data = TRUE)
    first <- compared$data[compared$data$Model == "First", ]
    expect_equal(first[[2]], single[[2]], ignore_attr = TRUE)
})

test_that("length-axis spectra transform their densities", {
    for (bin_average in c(FALSE, TRUE)) {
        p <- params
        second_order_w(p) <- c(bin_average = bin_average)
        for (power in 0:2) {
            by_weight <- plotSpectra(
                p, species = species, resource = FALSE, power = power,
                size_axis = "w", return_data = TRUE)
            by_length <- plotSpectra(
                p, species = species, resource = FALSE, power = power,
                size_axis = "l", return_data = TRUE)
            sp_idx <- match(by_length$Species, p@species_params$species)
            b <- p@species_params$b[sp_idx]
            jacobian <- if (power == 2) {
                b
            } else {
                b * by_weight$w / by_length$l
            }
            expect_equal(by_length[[2]], unname(by_weight[[2]] * jacobian))
        }
    }

    p <- plotSpectra(params, species = species, resource = FALSE,
                     power = 1, size_axis = "l")
    expect_identical(p$scales$get_scales("y")$name,
                     "Biomass density [g/cm]")

    compared <- plotSpectra2(params, params, species = species,
                             resource = FALSE, power = 0, size_axis = "l")
    compared_data <- compared$data[compared$data$Model == "First", ]
    by_weight <- plotSpectra(params, species = species, resource = FALSE,
                             power = 0, return_data = TRUE)
    sp_idx <- match(compared_data$Species, params@species_params$species)
    jacobian <- params@species_params$b[sp_idx] * by_weight$w /
        compared_data$l
    expect_equal(compared_data[[2]], unname(by_weight[[2]] * jacobian))

    # The viewport uses converted units without filtering either operand.
    limit <- stats::median(compared_data[[2]])
    limited <- plotSpectra2(params, params, species = species,
                            resource = FALSE, power = 0, size_axis = "l",
                            ylim = c(limit, NA))
    expect_equal(limited$data, compared$data)
    expect_equal(limited$coordinates$limits$y, c(limit, NA))
    single <- plotSpectra(params, species = species, resource = FALSE,
                          power = 0, size_axis = "l", ylim = c(limit, NA),
                          return_data = TRUE)
    expect_equal(sort(limited$data[[2]][limited$data$Model == "First"]),
                 sort(single[[2]]))
})

test_that("the spectrum total is unchanged on a weight axis", {
    # The total is now summed after the size axis is converted rather than
    # before, which on a weight axis must make no difference at all.
    for (power in 0:2) {
        d <- plotSpectra(params, total = TRUE, power = power,
                         return_data = TRUE)
        total <- d[d$Legend == "Total", ]
        fish_idx <- (length(params@w_full) - length(params@w) +
                         1):length(params@w_full)
        expected <- params@initial_n_pp
        expected[fish_idx] <- expected[fish_idx] + colSums(params@initial_n)
        expected <- expected * params@w_full^power
        keep <- expected > 1e-20
        expect_equal(total[[2]], unname(expected[keep]))
    }
})

test_that("the spectrum total appears on a length axis", {
    on_l <- plotSpectra(params, total = TRUE, size_axis = "l",
                        return_data = TRUE)
    total <- on_l[on_l$Legend == "Total", ]
    expect_gt(nrow(total), 0)

    # All the species here share one weight-length relationship, so at each
    # length the total is exactly the sum of the series drawn there
    contributors <- on_l[on_l$Legend != "Total", ]
    at <- total$l[[which.max(total$l)]]
    expect_equal(total[[2]][total$l == at],
                 sum(contributors[[2]][contributors$l == at]))

    # The total is the total of everything the model holds, so it does not
    # depend on what is drawn: neither hiding the resource nor selecting
    # species changes it
    small <- min(params@w)
    total_at <- function(...) {
        d <- plotSpectra(params, total = TRUE, return_data = TRUE, ...)
        d[[2]][d$Legend == "Total" & d$w == small]
    }
    expect_equal(total_at(resource = FALSE), total_at(resource = TRUE))
    expect_equal(total_at(species = species), total_at())
})

test_that("the comparison spectra keep their total on a length axis", {
    params2 <- params
    species_params(params2)$gamma <- species_params(params2)$gamma * 1.1

    for (size_axis in c("w", "l")) {
        compared <- plotSpectra2(params, params2, total = TRUE,
                                 size_axis = size_axis)$data
        expect_true("Total" %in% compared$Legend)
        relative <- plotSpectraRelative(params, params2, total = TRUE,
                                        size_axis = size_axis)$data
        expect_true("Total" %in% relative$Legend)
    }

    # Each model in the comparison holds exactly what plotSpectra() draws for
    # it on its own, total included, on either axis
    for (size_axis in c("w", "l")) {
        for (power in 0:2) {
            compared <- plotSpectra2(params, params, total = TRUE,
                                     power = power,
                                     size_axis = size_axis)$data
            single <- plotSpectra(params, total = TRUE, power = power,
                                  size_axis = size_axis, return_data = TRUE)
            expect_equal(sort(compared[[2]][compared$Model == "First"]),
                         sort(single[[2]]))
            expect_equal(sort(compared[[2]][compared$Legend == "Total" &
                                                compared$Model == "First"]),
                         sort(single[[2]][single$Legend == "Total"]))
        }
    }
})

test_that("convert_plot_spectrum_axis agrees with the density conversion", {
    plot_dat <- data.frame(w = params@w[1:5], value = seq_len(5),
                           Species = params@species_params$species[[1]])
    expect_identical(convert_plot_spectrum_axis(plot_dat, params, "l",
                                                power = 0),
                     convert_plot_density_axis(plot_dat, params, "l",
                                               density_wrt = "w"))
    expect_identical(convert_plot_spectrum_axis(plot_dat, params, "l",
                                                power = 2),
                     convert_plot_density_axis(plot_dat, params, "l",
                                               density_wrt = "log_w"))
    expect_identical(spectrum_density_wrt(TRUE), "log_w")
    expect_identical(spectrum_density_wrt(FALSE), "w")
    expect_identical(spectrum_density_wrt(NULL), "w")
})

test_that("Legends have correct entries", {
    # plotSpectra
    p <- plotSpectra(params, species = 2:3)
    expect_length(unique(p$data$Legend), 3)
    p <- plotSpectra(params, species = 1:3, resource = FALSE)
    expect_length(unique(p$data$Legend), 3)
    p <- plotSpectra(params, species = 2:3, total = TRUE)
    expect_length(unique(p$data$Legend), 4)
    p <- plotSpectra(params, species = 2:3, background = TRUE)
    expect_length(unique(p$data$Legend), 3)
    p <- plotSpectra(params_bkgrd, species = 2:3, background = TRUE)
    expect_length(unique(p$data$Legend), 4)
    p <- plotSpectra(params_single)
    expect_length(unique(p$data$Legend), 2)
})

test_that("plotSpectra averages over time range", {
    time_sel <- c(2:4)
    time_range <- getTimes(NS_sim_small)[time_sel]
    # arithmetic mean
    df <- plotSpectra(NS_sim_small, species = 1, time_range = time_range,
                      power = 0, return_data = TRUE)
    expected <- mean(NS_sim_small@n[time_sel, 1, 1])
    expect_equal(df[[2]][1], expected)
    # geometric mean
    df <- plotSpectra(NS_sim_small, species = 1, time_range = time_range,
                      geometric_mean = TRUE,
                      power = 0, return_data = TRUE)
    expected <- exp(mean(log(NS_sim_small@n[time_sel, 1, 1])))
    expect_equal(df[[2]][1], expected)
})

test_that("plotSpectra validates empty selection and can return total only", {
    expect_error(plotSpectra(params, species = rep(FALSE, nrow(species_params(params))),
                             total = FALSE, resource = FALSE),
                 "There is nothing to plot")
    df <- plotSpectra(params, species = rep(FALSE, nrow(species_params(params))),
                      total = TRUE, resource = FALSE, return_data = TRUE)
    expect_true(all(df$Legend == "Total"))
})

test_that("plotSpectra supports base plot log argument", {
    p_y <- plotSpectra(params, species = species, log = "y")
    expect_identical(p_y$scales$get_scales("x")$trans$name, "identity")
    expect_identical(p_y$scales$get_scales("y")$trans$name, "log-10")

    p_xy <- plotSpectra(params, species = species, log = "xy")
    expect_identical(p_xy$scales$get_scales("x")$trans$name, "log-10")
    expect_identical(p_xy$scales$get_scales("y")$trans$name, "log-10")

    p_none <- plotSpectra(params, species = species, log = "")
    expect_identical(p_none$scales$get_scales("x")$trans$name, "identity")
    expect_identical(p_none$scales$get_scales("y")$trans$name, "identity")

    p_sim <- plotSpectra(sim, species = species, log_x = FALSE, log_y = FALSE)
    expect_identical(p_sim$scales$get_scales("x")$trans$name, "identity")
    expect_identical(p_sim$scales$get_scales("y")$trans$name, "identity")

    expect_error(plotSpectra(params, species = species, log = "z"),
                 "containing only")
})

test_that("plotSpectra draws the spectrum at bin centres with the w^power weight there", {
    p0 <- params
    p1 <- params
    second_order_w(p1) <- c(bin_average = TRUE)
    beta <- p0@w_full[2] / p0@w_full[1]
    for (pw in c(0, 1, 2)) {
        d0 <- plotSpectra(p0, power = pw, return_data = TRUE)
        d1 <- plotSpectra(p1, power = pw, return_data = TRUE)
        d0 <- d0[order(d0$Species, d0$w), ]
        d1 <- d1[order(d1$Species, d1$w), ]
        expect_equal(nrow(d0), nrow(d1))
        # x moves to the geometric bin centre (a uniform sqrt(beta) shift) ...
        expect_equal(unname(d1$w / d0$w), rep(sqrt(beta), nrow(d1)))
        # ... and the w^power weight is evaluated there, scaling the value
        # (column 2, named by the y-label) by (w*/w)^power = beta^(power/2).
        expect_equal(unname(d1[[2]] / d0[[2]]),
                     rep(beta^(pw / 2), nrow(d1)))
    }
})

test_that("plotSpectra default (first order) is unchanged", {
    expect_identical(plotSpectra(params, power = 2, return_data = TRUE),
                     plotSpectra(params, power = 2, return_data = TRUE))
    # The default model never shifts: x stays on the model grid nodes.
    d <- plotSpectra(params, power = 2, resource = FALSE, total = FALSE,
                     background = FALSE, return_data = TRUE)
    expect_true(all(d$w %in% params@w))
})

test_that("plotSpectraRelative shifts x to centres", {
    p1a <- params
    p1b <- params
    p1b@initial_n <- p1b@initial_n * 1.5
    second_order_w(p1a) <- c(bin_average = TRUE)
    second_order_w(p1b) <- c(bin_average = TRUE)
    p_centre <- plotSpectraRelative(p1a, p1b, species = species, resource = FALSE)
    d_centre <- p_centre$data[order(p_centre$data$Species, p_centre$data$w), ]
    # The x-location is the geometric bin centre, not the node.
    p_node <- plotSpectraRelative(params, params, species = species,
                                  resource = FALSE)
    expect_true(all(p_node$data$w %in% params@w))
    expect_false(all(d_centre$w %in% params@w))
})

test_that("plotSpectraRelative compares two models with different lengths", {
    params_a <- params
    params_a@species_params$a <- 0.01
    params_a@species_params$b <- 3
    params_b <- params_a
    params_b@species_params$a <- 0.005

    # The two models put the same weights at different lengths, so their length
    # grids only overlap; matching by equality would keep almost nothing.
    rel <- plotSpectraRelative(params_a, params_b, species = species,
                               resource = FALSE, size_axis = "l")
    expect_gt(nrow(rel$data), 10)
    expect_true(all(is.finite(rel$data$rel_diff)))

    # Two identical models are still compared exactly, on either axis.
    for (axis in c("w", "l")) {
        same <- plotSpectraRelative(params_a, params_a, species = species,
                                    resource = FALSE, size_axis = axis)
        expect_true(all(abs(same$data$rel_diff) < 1e-12))
    }
})

test_that("spectrum plot appearance is unchanged", {
    p <- plotSpectra(sim, species = species, total = TRUE,
                     time_range = 1:3, power = 2,
                     ylim = c(1e6, NA))
    expect_doppelganger("Plot Spectra", p)
})

test_that("interactive wrappers preserve tooltip fields", {
    gp <- plotlySpectra(params, species = species)
    expect_s3_class(gp, "plotly")
    expect_equal(tooltip_fields(gp), c("Species", "w", "Biomass density"))

    gp <- plotlySpectra2(params, sim, species = species)
    expect_s3_class(gp, "plotly")
    expect_equal(tooltip_fields(gp), c("Legend", "w", "Biomass density", "Model"))

    gp <- plotlySpectraRelative(params, params, species = species, resource = FALSE)
    expect_s3_class(gp, "plotly")
    expect_equal(tooltip_fields(gp, trace = 2), c("w", "rel_diff", "Legend"))
})

test_that("plotSpectra honours biomass with power (#501)", {
    # Previously `biomass` was silently dropped whenever `power` was supplied
    numbers_in_log_w <- plotSpectra(params, species = species,
                                    resource = FALSE, power = 1,
                                    biomass = FALSE, return_data = TRUE)
    expect_identical(names(numbers_in_log_w)[[2]],
                     "Number density in log weight")
    expect_error(plotSpectra(params, power = 2, biomass = FALSE),
                 "imply `power = 0`")

    # The plotly wrapper passes the flag on rather than swallowing it
    plotly_numbers <- plotlySpectra(params, species = species,
                                    resource = FALSE, biomass = FALSE)
    expect_identical(plotly_numbers$x$layout$yaxis$title$text,
                     "Number density [1/g]")
})

test_that("spectra and comparisons support length axes", {
    params_len <- params
    params_len@species_params$a <- 0.01
    params_len@species_params$b <- 3
    sim_len <- sim
    sim_len@params <- params_len

    spectra_w <- plotSpectra(params_len, species = species, resource = FALSE,
                             total = FALSE, return_data = TRUE)
    spectra_l <- plotSpectra(params_len, species = species, resource = FALSE,
                             total = FALSE, size_axis = "l",
                             return_data = TRUE)
    expect_true("l" %in% names(spectra_l))
    expect_false("w" %in% names(spectra_l))
    sp_idx <- match(as.character(spectra_w$Species),
                    as.character(params_len@species_params$species))
    expected_l <- w2l(spectra_w$w,
                      params_len@species_params[sp_idx, , drop = FALSE])
    expect_equal(spectra_l$l, expected_l, ignore_attr = TRUE)

    llim <- stats::quantile(spectra_l$l, c(0.25, 0.75), names = FALSE)
    spectra_l_limited <- plotSpectra(params_len, species = species,
                                     resource = FALSE, total = FALSE,
                                     size_axis = "l", llim = llim,
                                     return_data = TRUE)
    expect_equal(spectra_l_limited, spectra_l)

    # The resource and the total both appear on a length axis: the resource
    # has its own weight-length relationship, and the total is summed after
    # the conversion, at equal length.
    spectra_all <- plotSpectra(params_len, species = species,
                               resource = TRUE, total = TRUE,
                               size_axis = "l", return_data = TRUE)
    expect_true("Resource" %in% spectra_all$Legend)
    expect_true("Total" %in% spectra_all$Legend)

    p <- plotSpectra(params_len, species = species, resource = FALSE,
                     size_axis = "l")
    expect_identical(p$scales$get_scales("x")$name, "Length [cm]")
    expect_identical(p$scales$get_scales("x")$trans$name, "log-10")

    expect_true("l" %in% names(plotSpectra2(params_len, params_len,
                                            species = species,
                                            resource = FALSE,
                                            size_axis = "l")$data))

    expect_true("l" %in% names(plotSpectraRelative(params_len, params_len,
                                                   species = species,
                                                   resource = FALSE,
                                                   size_axis = "l")$data))
})

test_that("public return-data shape is preserved", {
    expect_equal(plotSpectra(sim, species = species, wlim = c(1, NA),
                             return_data = TRUE),
                 plotSpectra(sim, species = species, return_data = TRUE))
})

test_that("spectrum viewport limits are set correctly", {
    # Test wlim and ylim in plotSpectra
    p <- plotSpectra(sim, species = species, wlim = c(1, 100), ylim = c(1e5, 1e8))
    # x-axis is first scale, y-axis is second
    expect_equal(p$coordinates$limits$x, c(1, 100))
    expect_equal(p$coordinates$limits$y, c(1e5, 1e8))

    # Test with NA values in wlim
    p <- plotSpectra(sim, species = species, wlim = c(10, NA), ylim = c(NA, 1e8))
    expect_equal(p$coordinates$limits$x[1], 10)
    expect_equal(p$coordinates$limits$x[2], max(params@w_full))
    expect_true(is.na(p$coordinates$limits$y[1]))
    expect_equal(p$coordinates$limits$y[2], 1e8)

    p <- plotSpectra(sim, species = species, resource = FALSE,
                     size_axis = "l", llim = c(10, 100))
    expect_equal(p$coordinates$limits$x, c(10, 100))

    # Default wlim lower depends on resource argument
    p_res <- plotSpectra(sim, species = species, return_data = TRUE)
    p_nores <- plotSpectra(sim, species = species, resource = FALSE, return_data = TRUE)
    expect_true(min(p_res$w) < min(params@w))
    expect_equal(min(p_nores$w), min(params@w))
})

test_that("spectrum comparison length-axis totals are retained", {
    params <- NS_params_small

    s2 <- plotSpectra2(params, params, size_axis = "l", total = TRUE)
    expect_s3_class(s2, "ggplot")
    expect_identical(s2$scales$get_scales("x")$name, "Length [cm]")
    expect_true("Total" %in% levels(s2$data$Legend))

    # plotSpectraRelative with size_axis = "l" and total = TRUE
    sr <- plotSpectraRelative(params, params, size_axis = "l", total = TRUE)
    expect_s3_class(sr, "ggplot")
    expect_identical(sr$scales$get_scales("x")$name, "Length [cm]")
    expect_true("Total" %in% levels(sr$data$Legend))
})

test_that("size viewports preserve returned data and totals", {
    for (second in c(FALSE, TRUE)) {
        p <- NS_params_small
        second_order_w(p) <- c(bin_average = second)
        for (axis in c("w", "l")) {
            for (fun in list(plotSpectra)) {
                full <- fun(p, species = "Cod", size_axis = axis,
                            return_data = TRUE)
                zoom <- fun(p, species = "Cod", size_axis = axis,
                            wlim = c(1, 100), llim = c(5, 20),
                            ylim = c(0.2, 0.8), return_data = TRUE)
                expect_equal(zoom, full)
            }
            for (fun in list(plotSpectra2, plotSpectraRelative)) {
                full <- fun(p, p, total = TRUE, size_axis = axis)
                zoom <- fun(p, p, total = TRUE, size_axis = axis,
                            wlim = c(1, 100), llim = c(5, 20),
                            ylim = c(0.2, 0.8))
                expect_equal(zoom$data, full$data)
            }
        }
    }
})
