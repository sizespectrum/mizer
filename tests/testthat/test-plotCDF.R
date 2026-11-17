# Plot fixtures use edition 1, as do the existing visual baselines.
withr::local_options(mizer_defaults_edition = 1)
params <- plot_params_small
sim <- plot_sim_small
sim0 <- plot_sim0_small
species <- c(1, 2)

test_that("plotCDF rejects per_log_size and accepts biomass", {
    expect_error(plotCDF(params, per_log_size = TRUE),
                 "not allowed here")
    expect_error(plotCDF2(params, params, per_log_size = TRUE),
                 "not allowed here")
    by_flag <- plotCDF(params, species = species, resource = FALSE,
                       biomass = FALSE, return_data = TRUE)
    by_power <- plotCDF(params, species = species, resource = FALSE,
                        power = 0, return_data = TRUE)
    expect_equal(by_flag[[2]], by_power[[2]], ignore_attr = TRUE)
    expect_identical(names(by_flag)[[2]], names(by_power)[[2]])
})

test_that("plotCDF plots cumulative spectra from small to large sizes", {
    p <- plotCDF(params, species = species, resource = FALSE, power = 0,
                 wlim = c(1, NA), return_data = TRUE)
    expect_equal(p, plotCDF(params, species = species, resource = FALSE,
                            power = 0, return_data = TRUE))
    expect_true(all(p[[2]] >= 0))
    for (sp in unique(p$Species)) {
        sp_dat <- p[p$Species == sp, ]
        expect_equal(max(sp_dat[[2]]), 1)
        expect_true(all(diff(sp_dat[[2]]) >= -1e-12))
    }

    spectra <- plotSpectra(params, species = species, resource = FALSE,
                           power = 1, wlim = c(1, NA), return_data = TRUE)
    cdf <- plotCDF(params, species = species, resource = FALSE,
                   power = 1, wlim = c(1, NA), normalise = FALSE,
                   return_data = TRUE)
    widths <- params@dw_full[match(spectra$w, params@w_full)]
    expected <- sum(spectra[[2]][spectra$Species == species[[1]]] *
                        widths[spectra$Species == species[[1]]])
    observed <- max(cdf[[2]][cdf$Species == species[[1]]])
    expect_equal(observed, expected)

    p_plot <- plotCDF(params, species = species, resource = FALSE)
    expect_s3_class(p_plot, "ggplot")
    expect_identical(p_plot$scales$get_scales("x")$trans$name, "log-10")
    expect_match(p_plot$scales$get_scales("y")$name, "biomass",
                 ignore.case = TRUE)

    p_linear <- plotCDF(params, species = species, resource = FALSE,
                        log_x = FALSE)
    expect_identical(p_linear$scales$get_scales("x")$trans$name, "identity")

    p_log_none <- plotCDF(params, species = species, resource = FALSE,
                          log = "")
    expect_identical(p_log_none$scales$get_scales("x")$trans$name,
                     "identity")

    p_log_x <- plotCDF(params, species = species, resource = FALSE,
                       log = "x")
    expect_identical(p_log_x$scales$get_scales("x")$trans$name, "log-10")

    p_log_y <- plotCDF(params, species = species, resource = FALSE,
                       log_y = TRUE)
    expect_identical(p_log_y$scales$get_scales("y")$trans$name, "log-10")

    p_log_xy <- plotCDF(params, species = species, resource = FALSE,
                        log = "xy")
    expect_identical(p_log_xy$scales$get_scales("x")$trans$name, "log-10")
    expect_identical(p_log_xy$scales$get_scales("y")$trans$name, "log-10")

    p_abundance <- plotCDF(params, species = species, resource = FALSE,
                           power = 0)
    expect_match(p_abundance$scales$get_scales("y")$name, "abundance",
                 ignore.case = TRUE)
})

test_that("plotCDF supports simulations, resource, total and unnormalised output", {
    p <- plotCDF(sim, species = species, time_range = 1:3,
                 total = TRUE, resource = TRUE, normalise = FALSE,
                 return_data = TRUE)
    expect_true(all(c("Resource", "Total") %in% p$Legend))
    expect_true(max(p[[2]]) > 1)

    p_plot <- plotCDF(sim, species = species, time_range = 1:3,
                      total = TRUE, resource = TRUE, normalise = FALSE)
    expect_s3_class(p_plot, "ggplot")
})

test_that("plotCDF2 compares cumulative distributions", {
    p <- plotCDF2(params, params, name1 = "Original", name2 = "Changed",
                  species = species, total = TRUE, resource = FALSE,
                  wlim = c(1, NA), normalise = FALSE, log = "")
    expect_s3_class(p, "ggplot")
    expect_identical(levels(p$data$Model), c("Original", "Changed"))
    expect_true("Total" %in% p$data$Legend)
    expect_equal(p$coordinates$limits$x, c(1, NA))
    expect_true(any(p$data$w < 1))
    expect_identical(p$scales$get_scales("x")$trans$name, "identity")

    p_log <- plotCDF2(params, sim, species = species, resource = FALSE,
                      log = "x")
    expect_s3_class(p_log, "ggplot")
    expect_identical(p_log$scales$get_scales("x")$trans$name, "log-10")

    expect_s3_class(plotCDF2(sim, sim0, species = species,
                             time_range = 1:3, resource = FALSE),
                    "ggplot")
    p_log_y2 <- plotCDF2(params, sim, species = species, resource = FALSE,
                         log = "y")
    expect_identical(p_log_y2$scales$get_scales("y")$trans$name, "log-10")
})

test_that("plotCDF places the cumulative on upper bin edges (inclusive convention)", {
    upper_edges <- round(params@w_full + params@dw_full, 6)
    # The inclusive cumulative sum belongs on each bin's *upper* edge w_k + dw_k,
    # in both the default and the second-order schemes.
    cdf0 <- plotCDF(params, power = 2, return_data = TRUE)
    expect_true(all(round(cdf0$w, 6) %in% upper_edges))
    # The left bin edges (nodes) are no longer used for the placement.
    nodes <- round(params@w, 6)
    expect_false(all(round(cdf0$w, 6) %in% nodes))

    p1 <- params
    second_order_w(p1) <- c(bin_average = TRUE)
    cdf <- plotCDF(p1, power = 2, return_data = TRUE)
    centres <- round(c(bin_midpoints(p1), bin_midpoints(p1, w_full = TRUE)), 6)
    expect_true(all(round(cdf$w, 6) %in% upper_edges))
    # never on the geometric bin centres
    expect_false(any(round(cdf$w, 6) %in% setdiff(centres, upper_edges)))
    # Cumulative is monotonic increasing per species.
    sp1 <- cdf[cdf$Species == p1@species_params$species[1], ]
    sp1 <- sp1[order(sp1$w), ]
    expect_true(all(diff(sp1[[2]]) >= -1e-12))
})

test_that("plotCDF accumulates the getter weights without a one-bin offset", {
    # With N = w^-2 and power = 2, the trapezoidal bin-average weight that
    # size_integral_weights() gives each bin but the last is a constant
    # multiple, (1 + beta^2) / 2, of its width. The cumulative curve is then
    # linear in w through the lower edge of the grid, which pins down the edge
    # placement: plotting each value at the left edge of its bin (#383) would
    # break it.
    p <- params
    second_order_w(p) <- c(bin_average = TRUE)
    sp <- 1
    p@initial_n[sp, ] <- p@w^(-2)
    cdf <- plotCDF(p, species = sp, power = 2, resource = FALSE, total = FALSE,
                   background = FALSE, normalise = FALSE, return_data = TRUE)
    cdf <- cdf[order(cdf$w), ]
    beta <- p@w[2] / p@w[1]
    inner <- seq_len(nrow(cdf) - 1)
    expect_equal(unname(cdf[[2]][inner]),
                 (1 + beta^2) / 2 * (cdf$w[inner] - min(p@w)),
                 tolerance = 1e-8)
    # The smallest plotted x is the upper edge of the first bin.
    expect_equal(cdf$w[1], unname(p@w[1] + p@dw[1]))
    # The end is the integral that the summary functions report.
    expect_equal(tail(cdf[[2]], 1),
                 unname(sizeIntegral(p, weighting = p@w^2,
                                     exact_bounds = TRUE)[sp]))
})

test_that("CDF bounds include partial bins exactly once in both schemes", {
    for (second in c(FALSE, TRUE)) {
        p <- NS_params_small
        second_order_w(p) <- c(bin_average = second)
        p@initial_n[] <- 0
        bin <- 8
        p@initial_n["Cod", bin] <- 3
        lower <- p@w[bin] + 0.2 * p@dw[bin]
        upper <- p@w[bin] + 0.7 * p@dw[bin]
        for (power in c(0, 1, 2, 1.5)) {
            # The rule of sizeIntegral(exact_bounds = TRUE): the weight at the
            # left edge in a first-order model; otherwise w^power interpolated
            # linearly across the bin and averaged over the covered part.
            edge_weight <- p@w[bin + 0:1]^power
            full_weight <- if (second) mean(edge_weight) else edge_weight[1]
            cut_weight <- edge_weight[1] +
                if (second) (0.2 + 0.7) / 2 * diff(edge_weight) else 0
            full <- plotCDF(p, species = "Cod", power = power,
                            normalise = FALSE, return_data = TRUE)
            expect_equal(max(full[[2]]), 3 * full_weight * p@dw[bin])
            cut <- plotCDF(p, species = "Cod", power = power,
                           min_w = lower, max_w = upper,
                           normalise = FALSE, return_data = TRUE)
            expect_equal(cut$w, upper)
            expect_equal(cut[[2]], 3 * cut_weight * (upper - lower))
            normalised <- plotCDF(p, species = "Cod", power = power,
                                   min_w = lower, max_w = upper,
                                   return_data = TRUE)
            expect_equal(normalised[[2]], 1)
        }
        # Several full bins between two partial ones; no neighbour may leak
        # into the mask, and a viewport may not mask the result a second time.
        p@initial_n["Cod", ] <- 2
        upper <- p@w[11] + 0.6 * p@dw[11]
        cut <- plotCDF(p, species = "Cod", biomass = FALSE,
                       min_w = lower, max_w = upper,
                       normalise = FALSE, return_data = TRUE)
        expect_equal(tail(cut[[2]], 1), 2 * (upper - lower))
        expect_equal(cut[[2]], 2 * (cut$w - lower))
        zoom <- plotCDF(p, species = "Cod", biomass = FALSE,
                        min_w = lower, max_w = upper,
                        wlim = c(lower + 1, upper - 1), ylim = c(0.1, 0.9),
                        normalise = FALSE, return_data = TRUE)
        expect_equal(zoom, cut)
    }
})

test_that("CDF length bounds use each species and the resource's conversion", {
    for (second in c(FALSE, TRUE)) {
        p <- NS_params_small
        second_order_w(p) <- c(bin_average = second)
        p@species_params$a <- c(0.01, 0.03, 0.02)
        p@species_params$b <- c(3, 2.8, 3.2)
        p@initial_n[] <- 1
        p@initial_n_pp[] <- 1
        result <- plotCDF(p, min_l = 1, max_l = 2, size_axis = "l",
                          resource = TRUE, total = TRUE, biomass = FALSE,
                          normalise = FALSE, return_data = TRUE)
        sp <- p@species_params
        rp <- resource_length_params(p)
        amounts <- c(setNames(sp$a * (2^sp$b - 1), sp$species),
                     Resource = rp$a * (2^rp$b - 1))
        amounts <- c(amounts, Total = sum(amounts))
        for (name in names(amounts)) {
            curve <- result[result$Species == name, ]
            expect_equal(tail(curve$l, 1), 2)
            expect_equal(tail(curve[[2]], 1), unname(amounts[name]))
        }
        # Length bounds also select the same integrals on a weight axis.
        by_weight <- plotCDF(p, min_l = 1, max_l = 2,
                             resource = TRUE, biomass = FALSE,
                             normalise = FALSE, return_data = TRUE)
        expect_equal(by_weight[[2]], result[[2]][result$Species != "Total"])
        expect_equal(plotCDF(p, min_w = 999, max_w = 1000,
                             min_l = 1, max_l = 2, biomass = FALSE,
                             return_data = TRUE),
                     plotCDF(p, min_l = 1, max_l = 2, biomass = FALSE,
                             return_data = TRUE))
    }
})

test_that("CDF totals include hidden contributors and keep completed mass", {
    p <- NS_params_small
    p@initial_n[] <- 1
    p@initial_n_pp[] <- 1
    all <- plotCDF(p, biomass = FALSE, normalise = FALSE, total = TRUE,
                   resource = TRUE, size_axis = "l",
                   min_w = 0.2, max_w = 100, return_data = TRUE)
    hidden <- plotCDF(p, species = "Cod", biomass = FALSE, normalise = FALSE,
                      total = TRUE, resource = FALSE, background = FALSE,
                      size_axis = "l", min_w = 0.2, max_w = 100,
                      return_data = TRUE)
    total <- all[all$Species == "Total", ]
    expect_equal(hidden[hidden$Species == "Total", ], total)
    expect_true(all(diff(total[[2]]) >= -1e-10))
    expect_equal(tail(total[[2]], 1), 4 * (100 - 0.2))
})

test_that("CDF bounds validate, handle empty mass and reach every entry point", {
    p <- NS_params_small
    expect_error(plotCDF(p, min_w = NA), "min_w")
    expect_error(plotCDF(p, min_l = -1), "min_l")
    expect_error(plotCDF(p, max_w = c(1, 2)), "max_w")
    expect_error(plotCDF(p, min_w = 10, max_w = 1), "lower integration bound")
    expect_error(plotCDF(p, min_w = max(p@w) * 10), "do not overlap")
    p@initial_n["Cod", ] <- 0
    empty <- plotCDF(p, species = "Cod", return_data = TRUE)
    expect_true(all(is.na(empty[[2]])))
    empty <- plotCDF(p, species = "Cod", normalise = FALSE, return_data = TRUE)
    expect_true(all(empty[[2]] == 0))

    sim <- NS_sim_small
    expected <- plotCDF(getParams(sim, time_range = 1:2, geometric_mean = TRUE),
                        min_w = 1, max_w = 100, return_data = TRUE)
    expect_equal(plotCDF(sim, time_range = 1:2, geometric_mean = TRUE,
                         min_w = 1, max_w = 100, return_data = TRUE), expected)
    compared <- plotCDF2(sim, sim, time_range = 1:2, geometric_mean = TRUE,
                         min_w = 1, max_w = 100)
    expect_equal(compared$data[[2]][compared$data$Model == "First"], expected[[2]])
    expect_s3_class(plotlyCDF(sim, min_w = 1, max_w = 100), "plotly")
    expect_s3_class(plotlyCDF2(sim, sim, min_l = 1, max_l = 10), "plotly")
})

test_that("interactive wrappers preserve tooltip fields", {
    gp <- plotlyCDF(params, species = species, resource = FALSE)
    expect_s3_class(gp, "plotly")
    expect_equal(tooltip_fields(gp), c("Species", "w", "Cumulative proportion of biomass"))

    gp <- plotlyCDF2(params, params, species = species, resource = FALSE)
    expect_s3_class(gp, "plotly")
    expect_equal(tooltip_fields(gp), c("Legend", "w", "Cumulative proportion of biomass", "Model"))
})

test_that("plotCDF rejects conflicting biomass and power", {
    expect_error(plotCDF(params, power = 1, biomass = FALSE),
                 "use `biomass = FALSE` on its own")
})

test_that("CDFs and comparisons support length axes", {
    params_len <- params
    params_len@species_params$a <- 0.01
    params_len@species_params$b <- 3
    sim_len <- sim
    sim_len@params <- params_len
    sp_rows <- rep(1:2, times = length(params_len@w))
    lengths <- w2l(rep(params_len@w, each = 2),
                   params_len@species_params[sp_rows, ])
    llim <- stats::quantile(lengths, c(0.25, 0.75), names = FALSE)

    expect_true("l" %in% names(plotCDF(params_len, species = species,
                                       resource = FALSE, size_axis = "l",
                                       return_data = TRUE)))
    cdf_l_limited <- plotCDF(params_len, species = species, resource = FALSE,
                             size_axis = "l", llim = llim,
                             return_data = TRUE)
    expect_equal(cdf_l_limited,
                 plotCDF(params_len, species = species, resource = FALSE,
                         size_axis = "l", return_data = TRUE))
    y_var <- names(cdf_l_limited)[2]
    expect_equal(stats::aggregate(cdf_l_limited[[y_var]] ~ cdf_l_limited$Species,
                                  FUN = max)[[2]],
                 rep(1, length(unique(cdf_l_limited$Species))))

    expect_true("l" %in% names(plotCDF2(params_len, params_len,
                                        species = species,
                                        resource = FALSE,
                                        size_axis = "l")$data))
})

test_that("CDF length-axis totals are retained", {
    params <- NS_params_small
    sim <- NS_sim_small

    # plotCDF on length axis with total
    cdf_l <- plotCDF(params, size_axis = "l", total = TRUE, return_data = TRUE)
    expect_identical(names(cdf_l)[[1]], "l")
    expect_true("Total" %in% cdf_l$Species)
    expect_equal(max(cdf_l[["Cumulative proportion of biomass"]][cdf_l$Species == "Total"]), 1)

    # plotCDF2 with size_axis = "l" and total = TRUE
    p2 <- plotCDF2(params, params, size_axis = "l", total = TRUE)
    expect_s3_class(p2, "ggplot")
    expect_identical(p2$scales$get_scales("x")$name, "Length [cm]")
    expect_true("Total" %in% levels(p2$data$Legend))
})

test_that("size viewports preserve returned data and totals", {
    for (second in c(FALSE, TRUE)) {
        p <- NS_params_small
        second_order_w(p) <- c(bin_average = second)
        for (axis in c("w", "l")) {
            for (fun in list(plotCDF)) {
                full <- fun(p, species = "Cod", size_axis = axis,
                            return_data = TRUE)
                zoom <- fun(p, species = "Cod", size_axis = axis,
                            wlim = c(1, 100), llim = c(5, 20),
                            ylim = c(0.2, 0.8), return_data = TRUE)
                expect_equal(zoom, full)
            }
            for (fun in list(plotCDF2)) {
                full <- fun(p, p, total = TRUE, size_axis = axis)
                zoom <- fun(p, p, total = TRUE, size_axis = axis,
                            wlim = c(1, 100), llim = c(5, 20),
                            ylim = c(0.2, 0.8))
                expect_equal(zoom$data, full$data)
            }
        }
    }
})

test_that("CDF bounds need only be in order for the curves they affect", {
    # 10 cm is less than 100 g for each fish but not for the resource.
    cdf <- plotCDF(NS_params_small, min_l = 10, max_w = 100,
                   resource = FALSE, total = FALSE, return_data = TRUE)
    expect_setequal(unique(cdf$Species),
                    NS_params_small@species_params$species)
    expect_error(plotCDF(NS_params_small, min_l = 10, max_w = 100,
                         resource = TRUE, return_data = TRUE),
                 "lower integration bound must be less than")
    # A total sums the resource too, whichever curves are drawn.
    expect_error(plotCDF(NS_params_small, min_l = 10, max_w = 100,
                         resource = FALSE, total = TRUE, return_data = TRUE),
                 "lower integration bound must be less than")
})

test_that("the end of a CDF is the getBiomass() integral over its range", {
    for (second in c(FALSE, TRUE)) {
        p <- NS_params_small
        second_order_w(p) <- c(bin_average = second)
        w <- p@w
        lower <- 0.3 * w[5] + 0.7 * w[6]
        upper <- 0.5 * (w[15] + w[16])
        cdf <- plotCDF(p, return_data = TRUE, normalise = FALSE,
                       biomass = TRUE, background = FALSE,
                       min_w = lower, max_w = upper)
        ends <- c(tapply(cdf[[2]], cdf$Species, max)[p@species_params$species])
        expect_equal(unname(ends),
                     unname(getBiomass(p, min_w = lower, max_w = upper,
                                       exact_bounds = TRUE)))
    }
    # A model of edition 1 or 2 counts whole bins by default: a bound on a
    # bin edge ends the CDF there, while getBiomass(max_w = w[15]) still
    # includes the bin that starts at w[15], as documented in ?plotCDF.
    p <- NS_params_small
    w <- p@w
    cdf <- plotCDF(p, return_data = TRUE, normalise = FALSE, biomass = TRUE,
                   background = FALSE, min_w = w[5], max_w = w[15])
    ends <- c(tapply(cdf[[2]], cdf$Species, max)[p@species_params$species])
    expect_equal(unname(ends),
                 unname(getBiomass(p, min_w = w[5], max_w = w[14])))
    expect_true(all(getBiomass(p, min_w = w[5], max_w = w[15]) >=
                        ends * (1 - 1e-12)))
})
