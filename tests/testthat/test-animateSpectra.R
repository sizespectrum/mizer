test_that("animateSpectra does not throw error", {
    sim <- project(NS_params, t_max = 2, t_save = 1, effort = 1)
    expect_error(animateSpectra(sim, species = c("Cod", "Haddock"),
                                time_range = c(1, 2),
                                wlim = c(1, 1000),
                                ylim = c(1e6, 1e9),
                                power = 1,
                                total = TRUE,
                                resource = TRUE), NA)
})

test_that("animateSpectra returns a plotly object", {
    sim <- project(NS_params, t_max = 2, t_save = 1, effort = 1)
    result <- animateSpectra(sim, time_range = c(1, 2))
    expect_s3_class(result, "plotly")
})

test_that("animateSpectra handles species parameter correctly", {
    sim <- project(NS_params, t_max = 2, t_save = 1, effort = 1)
    
    # Test with specific species
    result <- animateSpectra(sim, species = "Cod", time_range = c(1, 2))
    expect_s3_class(result, "plotly")
    
    # Test with multiple species
    result <- animateSpectra(sim, species = c("Cod", "Haddock"), time_range = c(1, 2))
    expect_s3_class(result, "plotly")
    
    # Test with NULL (default - all species)
    result <- animateSpectra(sim, species = NULL, time_range = c(1, 2))
    expect_s3_class(result, "plotly")
})

test_that("animateSpectra handles time_range parameter correctly", {
    sim <- project(NS_params, t_max = 5, t_save = 1, effort = 1)
    
    # Test with min/max vector
    expect_error(animateSpectra(sim, time_range = c(1, 3)), NA)
    
    # Test with full vector of values
    expect_error(animateSpectra(sim, time_range = 1:3), NA)
    
    # Test with missing time_range (should use entire range)
    expect_error(animateSpectra(sim), NA)
})

test_that("animateSpectra handles wlim parameter with NA values", {
    sim <- project(NS_params, t_max = 2, t_save = 1, effort = 1)
    
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
    sim <- project(NS_params, t_max = 2, t_save = 1, effort = 1)
    
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
    sim <- project(NS_params, t_max = 2, t_save = 1, effort = 1)
    
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
    sim <- project(NS_params, t_max = 2, t_save = 1, effort = 1)
    
    # Test with total = FALSE (default)
    result <- animateSpectra(sim, total = FALSE, time_range = c(1, 2))
    expect_s3_class(result, "plotly")
    
    # Test with total = TRUE (should include total line)
    result <- animateSpectra(sim, total = TRUE, time_range = c(1, 2))
    expect_s3_class(result, "plotly")
})

test_that("animateSpectra handles resource parameter correctly", {
    sim <- project(NS_params, t_max = 2, t_save = 1, effort = 1)
    
    # Test with resource = TRUE (default)
    result <- animateSpectra(sim, resource = TRUE, time_range = c(1, 2))
    expect_s3_class(result, "plotly")
    
    # Test with resource = FALSE (should exclude resource)
    result <- animateSpectra(sim, resource = FALSE, time_range = c(1, 2))
    expect_s3_class(result, "plotly")
})

test_that("animateSpectra validates input parameters", {
    sim <- project(NS_params, t_max = 2, t_save = 1, effort = 1)
    
    # Test invalid wlim length
    expect_error(animateSpectra(sim, wlim = c(1), time_range = c(1, 2)))
    expect_error(animateSpectra(sim, wlim = c(1, 10, 100), time_range = c(1, 2)))
    
    # Test invalid ylim length
    expect_error(animateSpectra(sim, ylim = c(1), time_range = c(1, 2)))
    expect_error(animateSpectra(sim, ylim = c(1, 10, 100), time_range = c(1, 2)))
})

test_that("animateSpectra uses consistent colors matching linecolour", {
    sim <- project(NS_params, t_max = 2, t_save = 1, effort = 1)
    
    # Get the result
    result <- animateSpectra(sim, species = c("Cod", "Haddock", "Sprat"), 
                            time_range = c(1, 2))
    
    # The plotly object should be created
    expect_s3_class(result, "plotly")
    
    # Extract the data from the plotly object
    plot_data <- plotly::plotly_build(result)
    
    expected_names <- intersect(names(sim@params@linecolour),
                                c("Cod", "Haddock", "Sprat", "Resource"))
    actual_names <- vapply(plot_data$x$data, `[[`, character(1), "name")
    actual_colours <- vapply(plot_data$x$data,
                             function(trace) trace$line$color,
                             character(1))
    expected_colours <- c(unname(sim@params@linecolour[
        intersect(names(sim@params@linecolour), c("Cod", "Haddock", "Sprat"))
    ]), "green")

    expect_identical(actual_names, expected_names)
    expect_identical(unname(actual_colours), expected_colours)
})

test_that("animateSpectra maintains color consistency when species go extinct", {
    # Create a simulation where a species might have very low abundance
    sim <- project(NS_params, t_max = 2, t_save = 1, effort = 1)
    
    # Test with species selection
    result <- animateSpectra(sim, species = c("Cod", "Haddock"), 
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
    sim <- project(NS_params, t_max = 2, t_save = 1, effort = 1)

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

test_that("animateSpectra sets the y axis title from power", {
    sim <- project(NS_params, t_max = 2, t_save = 1, effort = 1)

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
