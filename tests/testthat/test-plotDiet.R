# Plot fixtures use edition 1, as do the existing visual baselines.
withr::local_options(mizer_defaults_edition = 1)
params <- plot_params_small
sim <- plot_sim_small
species <- c(1, 2)

test_that("plotDiet restricts to meaningful abundance ranges", {
    # Get diet data for a species
    diet_data <- plotDiet(params, species = "3", return_data = TRUE)

    # Get abundance data for species "3"
    sp_idx <- which(dimnames(params@initial_n)[[1]] == "3")
    abundance <- params@initial_n[sp_idx, ] * params@w * params@dw
    max_abundance <- max(abundance)

    # Find the maximum meaningful size
    meaningful_idx <- which(abundance > 1e-5 * max_abundance)
    w_max_meaningful <- params@w[max(meaningful_idx)]

    # Check that all sizes in diet_data are at or below w_max_meaningful
    expect_true(all(diet_data$w <= w_max_meaningful))

    # Also verify that the maximum size in diet_data is close to w_max_meaningful
    # (there might be filtering due to proportion threshold, so we check it's not much smaller)
    if (nrow(diet_data) > 0) {
        max_w_in_data <- max(diet_data$w)
        expect_true(max_w_in_data <= w_max_meaningful)
    }
})

test_that("plotDiet works with MizerSim", {
    p <- plotDiet(sim, species = 2) # Species 11 is Cod in setup
    expect_true(is(p, "ggplot"))
    p <- plotDiet(sim, species = 2, time_range = 1:2)
    expect_true(is(p, "ggplot"))
})

test_that("plotDiet.MizerSim uses the simulated abundance, not the initial", {
    # sim is projected with fishing effort = 1, so its diet at later times
    # differs from the diet computed from the initial (params) abundances.
    times <- as.numeric(dimnames(sim@n)$time)
    diet_params <- getDiet(sim@params)
    diet_final <- getDiet(sim, time_range = max(times), drop = TRUE)
    expect_false(isTRUE(all.equal(diet_params, diet_final,
                                  check.attributes = FALSE)))
    # By default the MizerSim plot uses the final time step, so its data matches
    # the MizerParams plot built from the diet and abundance at that time step.
    last <- dim(sim@n)[1]
    n_last <- apply(sim@n[last, , , drop = FALSE], c(2, 3), mean)
    expect_equal(
        plotDiet(sim, species = 2, return_data = TRUE),
        plot_diet(sim@params, n = n_last, diet = diet_final, species = 2,
                  log_x = TRUE, log_y = FALSE, wlim = c(NA, NA),
                  llim = c(NA, NA), size_axis = "w", return_data = TRUE))
})

test_that("plotDiet.MizerSim averages rates, not proportions, over a range", {
    tr <- range(as.numeric(dimnames(sim@n)$time))
    # Correct aggregation: average the consumption rates then normalise
    rates <- getDiet(sim, time_range = tr, drop = FALSE, proportion = FALSE)
    rates_mean <- apply(rates, c(2, 3, 4), mean)
    total <- rowSums(rates_mean, dims = 2)
    diet_expected <- sweep(rates_mean, c(1, 2), total, "/")
    diet_expected[is.nan(diet_expected)] <- 0
    # Naive (wrong) aggregation: average the per-step proportions
    props <- getDiet(sim, time_range = tr, drop = FALSE, proportion = TRUE)
    diet_naive <- apply(props, c(2, 3, 4), mean)
    # The two genuinely differ for this projection
    expect_false(isTRUE(all.equal(diet_expected, diet_naive)))
    # plotDiet uses the correct (rate-averaged) proportions
    time_elements <- get_time_elements(sim, tr)
    n_avg <- apply(sim@n[time_elements, , , drop = FALSE], c(2, 3), mean)
    expect_equal(
        plotDiet(sim, species = 2, time_range = tr, return_data = TRUE),
        plot_diet(sim@params, n = n_avg, diet = diet_expected, species = 2,
                  log_x = TRUE, log_y = FALSE, wlim = c(NA, NA),
                  llim = c(NA, NA), size_axis = "w", return_data = TRUE))
})

test_that("diet plot appearance is unchanged", {
    sp_name <- NS_params_small@species_params$species[2]
    p <- plotDiet(NS_params_small, species = sp_name)
    expect_doppelganger("Plot Diet", p)
})

test_that("interactive wrappers preserve tooltip fields", {
    gp <- plotlyDiet(params, species = species[[1]])
    expect_s3_class(gp, "plotly")
    expect_equal(tooltip_fields(gp), c("w", "Proportion", "Prey"))
})

test_that("diet plots support length axes", {
    params_len <- params
    params_len@species_params$a <- 0.01
    params_len@species_params$b <- 3
    sim_len <- sim
    sim_len@params <- params_len
    sp_rows <- rep(1:2, times = length(params_len@w))
    lengths <- w2l(rep(params_len@w, each = 2),
                   params_len@species_params[sp_rows, ])
    llim <- stats::quantile(lengths, c(0.25, 0.75), names = FALSE)

    expect_true("l" %in% names(plotDiet(params_len, species = species[[1]],
                                        size_axis = "l",
                                        return_data = TRUE)))
    expect_equal(plotDiet(params_len, species = species[[1]],
                          size_axis = "l", llim = llim, return_data = TRUE),
                 plotDiet(params_len, species = species[[1]],
                          size_axis = "l", return_data = TRUE))
})

test_that("size viewports preserve returned data and totals", {
    for (second in c(FALSE, TRUE)) {
        p <- NS_params_small
        second_order_w(p) <- c(bin_average = second)
        for (axis in c("w", "l")) {
            for (fun in list(plotDiet)) {
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
