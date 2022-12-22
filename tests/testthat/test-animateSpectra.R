local_edition(3)
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
