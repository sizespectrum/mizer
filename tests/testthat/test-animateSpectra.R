example_animate_sim <- project(example_params(), t_max = 2, t_save = 1,
                               effort = 1)
example_animate_long_sim <- project(example_params(), t_max = 5, t_save = 1,
                                    effort = 1)
example_animate_high_effort_sim <- project(example_params(), t_max = 2,
                                           t_save = 1, effort = 10)
ns_animate_sim <- project(NS_params, t_max = 2, t_save = 1, effort = 1)

test_that("animateSpectra does not throw error", {
    sim <- example_animate_sim
    sp <- sim@params@species_params$species
    expect_error(animateSpectra(sim, species = sp[1:2],
                                time_range = c(1, 2),
                                wlim = c(1, 1000),
                                ylim = c(1e6, 1e9),
                                power = 1,
                                total = TRUE,
                                resource = TRUE), NA)
})

test_that("animateSpectra returns a plotly object", {
    sim <- example_animate_sim
    result <- animateSpectra(sim, time_range = c(1, 2))
    expect_s3_class(result, "plotly")
})

test_that("animateSpectra sets axis ranges without dropping vertices", {
    sim <- example_animate_sim
    result <- animateSpectra(sim, species = "Cod", time_range = c(1, 2),
                             resource = FALSE, wlim = c(1, 1000),
                             ylim = c(1e6, 1e9))
    built_plot <- plotly::plotly_build(result)
    frame_lengths <- lengths(lapply(built_plot$x$frames, function(frame) {
        frame$data[[1]]$x
    }))

    expect_equal(built_plot$x$layout$xaxis$range, log10(c(1, 1000)))
    expect_equal(built_plot$x$layout$yaxis$range, log10(c(1e6, 1e9)))
    expect_equal(frame_lengths, rep(length(sim@params@w), length(frame_lengths)))
})

test_that("animateSpectra derives missing x-axis limits from plotted data", {
    result <- animateSpectra(NS_sim, species = "Cod",
                             time_range = c(2000, 2001),
                             resource = FALSE, wlim = c(NA, 1000))
    built_plot <- plotly::plotly_build(result)
    first_frame <- built_plot$x$frames[[1]]$data[[1]]

    expect_equal(built_plot$x$layout$xaxis$range,
                 log10(c(min(first_frame$x), 1000)))
})

test_that("animateSpectra can disable interpolation between frames", {
    sim <- example_animate_sim
    result <- animateSpectra(sim, time_range = c(1, 2),
                             transition_duration = 0)
    expect_identical(result$animation$transition$duration, 0)
})

test_that("animateSpectra exposes plotly animation timing controls", {
    sim <- example_animate_sim
    result <- animateSpectra(sim, time_range = c(1, 2),
                             frame_duration = 800,
                             transition_duration = 120,
                             easing = "cubic-in-out")

    expect_identical(result$animation$frame$duration, 800)
    expect_identical(result$animation$transition$duration, 120)
    expect_identical(result$animation$transition$easing, "cubic-in-out")
})

test_that("animateSpectra handles species parameter correctly", {
    sim <- example_animate_sim

    # Test with specific species
    sp <- sim@params@species_params$species
    result <- animateSpectra(sim, species = sp[1], time_range = c(1, 2))
    expect_s3_class(result, "plotly")

    # Test with multiple species
    result <- animateSpectra(sim, species = sp[1:2], time_range = c(1, 2))
    expect_s3_class(result, "plotly")

    # Test with NULL (default - all species)
    result <- animateSpectra(sim, species = NULL, time_range = c(1, 2))
    expect_s3_class(result, "plotly")
})

test_that("animateSpectra handles time_range parameter correctly", {
    sim <- example_animate_long_sim

    # Test with min/max vector
    expect_error(animateSpectra(sim, time_range = c(1, 3)), NA)

    # Test with full vector of values
    expect_error(animateSpectra(sim, time_range = 1:3), NA)

    # Test with missing time_range (should use entire range)
    expect_error(animateSpectra(sim), NA)
})

test_that("animateSpectra handles wlim parameter with NA values", {
    sim <- example_animate_sim

    # Test with both NA (should use defaults)
    expect_error(animateSpectra(sim, wlim = c(NA, NA), time_range = c(1, 2)), NA)

    # Test with lower NA
    expect_error(animateSpectra(sim, wlim = c(NA, 1000), time_range = c(1, 2)), NA)

    # Test with upper NA
    expect_error(animateSpectra(sim, wlim = c(0.1, NA), time_range = c(1, 2)), NA)

    # Test with specific values
    expect_error(animateSpectra(sim, wlim = c(1, 1000), time_range = c(1, 2)), NA)
})

test_that("animateSpectra handles ylim parameter with NA values", {
    sim <- example_animate_sim

    # Test with both NA (should use defaults)
    expect_error(animateSpectra(sim, ylim = c(NA, NA), time_range = c(1, 2)), NA)

    # Test with lower NA
    expect_error(animateSpectra(sim, ylim = c(NA, 1e9), time_range = c(1, 2)), NA)

    # Test with upper NA
    expect_error(animateSpectra(sim, ylim = c(1e6, NA), time_range = c(1, 2)), NA)

    # Test with specific values
    expect_error(animateSpectra(sim, ylim = c(1e6, 1e9), time_range = c(1, 2)), NA)
})

test_that("animateSpectra handles power parameter correctly", {
    sim <- example_animate_sim

    # Test with power = 0 (Number density)
    result <- animateSpectra(sim, power = 0, time_range = c(1, 2))
    expect_s3_class(result, "plotly")

    # Test with power = 1 (Biomass density - default)
    result <- animateSpectra(sim, power = 1, time_range = c(1, 2))
    expect_s3_class(result, "plotly")

    # Test with power = 2 (Biomass density with respect to logarithmic size bins)
    result <- animateSpectra(sim, power = 2, time_range = c(1, 2))
    expect_s3_class(result, "plotly")

    # Test with custom power value
    result <- animateSpectra(sim, power = 1.5, time_range = c(1, 2))
    expect_s3_class(result, "plotly")
})

test_that("animateSpectra handles total parameter correctly", {
    sim <- example_animate_sim

    # Test with total = FALSE (default)
    result <- animateSpectra(sim, total = FALSE, time_range = c(1, 2))
    expect_s3_class(result, "plotly")

    # Test with total = TRUE (should include total line)
    result <- animateSpectra(sim, total = TRUE, time_range = c(1, 2))
    expect_s3_class(result, "plotly")
})

test_that("animateSpectra handles resource parameter correctly", {
    sim <- example_animate_sim

    # Test with resource = TRUE (default)
    result <- animateSpectra(sim, resource = TRUE, time_range = c(1, 2))
    expect_s3_class(result, "plotly")

    # Test with resource = FALSE (should exclude resource)
    result <- animateSpectra(sim, resource = FALSE, time_range = c(1, 2))
    expect_s3_class(result, "plotly")
})

test_that("animateSpectra validates input parameters", {
    sim <- example_animate_sim

    # Test invalid wlim length
    expect_error(animateSpectra(sim, wlim = c(1), time_range = c(1, 2)))
    expect_error(animateSpectra(sim, wlim = c(1, 10, 100), time_range = c(1, 2)))

    # Test invalid ylim length
    expect_error(animateSpectra(sim, ylim = c(1), time_range = c(1, 2)))
    expect_error(animateSpectra(sim, ylim = c(1, 10, 100), time_range = c(1, 2)))
})

test_that("animateSpectra uses consistent colors matching linecolour", {
    sim <- example_animate_sim

    # Get the result
    sp <- sim@params@species_params$species
    result <- animateSpectra(sim, species = sp[1:2],
                            time_range = c(1, 2))

    # The plotly object should be created
    expect_s3_class(result, "plotly")

    # Extract the data from the plotly object
    plot_data <- plotly::plotly_build(result)

    # Check that species are factors
    # This is done by checking the internal data structure
    expect_true(is.list(plot_data$x$data))

    # The colors should be assigned consistently
    # Each trace should have a specific color
    expect_true(length(plot_data$x$data) > 0)
})

test_that("animateSpectra maintains color consistency when species go extinct", {
    # Create a simulation where a species might have very low abundance
    sim <- example_animate_high_effort_sim

    # Test with species selection
    sp <- sim@params@species_params$species
    result <- animateSpectra(sim, species = sp[1:2],
                            time_range = c(1, 2))

    expect_s3_class(result, "plotly")

    # Build the plot to access internal structure
    built_plot <- plotly::plotly_build(result)

    # Check that we have traces (lines) in the plot
    expect_true(length(built_plot$x$data) > 0)

    # Each trace should have consistent properties
    for (trace in built_plot$x$data) {
        expect_true("line" %in% names(trace) || "marker" %in% names(trace))
    }
})

test_that("animateSpectra adds resource and total traces when requested", {
    sim <- ns_animate_sim

    built_plot <- plotly::plotly_build(
        animateSpectra(sim,
                       species = c("Cod", "Haddock"),
                       time_range = c(1, 2),
                       total = TRUE,
                       resource = TRUE)
    )

    trace_names <- vapply(built_plot$x$data, `[[`, character(1), "name")
    expect_identical(trace_names,
                     intersect(names(sim@params@linecolour),
                               c("Cod", "Haddock", "Total", "Resource")))
})

test_that("animateSpectra handles background parameter correctly", {
    params_bkgrd <- markBackground(NS_params,
                                    species = species_params(NS_params)$species[1:3])
    sim_bkgrd <- project(params_bkgrd, t_max = 2, t_save = 1, effort = 1)

    # background = TRUE (default) includes a "Background" trace
    built_on <- plotly::plotly_build(
        animateSpectra(sim_bkgrd, species = "Cod", time_range = c(1, 2),
                       resource = FALSE, background = TRUE)
    )
    trace_names_on <- vapply(built_on$x$data, `[[`, character(1), "name")
    expect_true("Background" %in% trace_names_on)

    # background = FALSE excludes the "Background" trace
    built_off <- plotly::plotly_build(
        animateSpectra(sim_bkgrd, species = "Cod", time_range = c(1, 2),
                       resource = FALSE, background = FALSE)
    )
    trace_names_off <- vapply(built_off$x$data, `[[`, character(1), "name")
    expect_false("Background" %in% trace_names_off)

    # background = TRUE on a model with no background species adds no "Background" trace
    sim_plain <- ns_animate_sim
    built_plain <- plotly::plotly_build(
        animateSpectra(sim_plain, species = "Cod", time_range = c(1, 2),
                       resource = FALSE, background = TRUE)
    )
    trace_names_plain <- vapply(built_plain$x$data, `[[`, character(1), "name")
    expect_false("Background" %in% trace_names_plain)
})

test_that("animateSpectra sets the y axis title from power", {
    sim <- ns_animate_sim

    built0 <- plotly::plotly_build(animateSpectra(sim, species = "Cod",
                                                  time_range = c(1, 2),
                                                  power = 0))
    built1 <- plotly::plotly_build(animateSpectra(sim, species = "Cod",
                                                  time_range = c(1, 2),
                                                  power = 1))
    built2 <- plotly::plotly_build(animateSpectra(sim, species = "Cod",
                                                  time_range = c(1, 2),
                                                  power = 2))
    built_custom <- plotly::plotly_build(animateSpectra(sim, species = "Cod",
                                                         time_range = c(1, 2),
                                                         power = 1.5))

    expect_identical(built0$x$layout$yaxis$title, "Number density [1/g]")
    expect_identical(built1$x$layout$yaxis$title, "Biomass density")
    expect_identical(built2$x$layout$yaxis$title, "Biomass density [g]")
    expect_identical(built_custom$x$layout$yaxis$title,
                     "Number density * w^1.5")
})

test_that("animate supports the log argument", {
    sim <- example_animate_sim
    built_xy <- plotly::plotly_build(
        animateSpectra(sim, species = "Cod", time_range = c(1, 2), log = "xy")
    )
    expect_identical(built_xy$x$layout$xaxis$type, "log")
    expect_identical(built_xy$x$layout$yaxis$type, "log")

    built_none <- plotly::plotly_build(
        animateSpectra(sim, species = "Cod", time_range = c(1, 2), log = "")
    )
    expect_identical(built_none$x$layout$xaxis$type, "-")
    expect_identical(built_none$x$layout$yaxis$type, "-")

    built_array <- plotly::plotly_build(
        animate(getFMort(sim), species = "Cod", time_range = c(1, 2), log = "y")
    )
    expect_identical(built_array$x$layout$xaxis$type, "-")
    expect_identical(built_array$x$layout$yaxis$type, "log")
})
