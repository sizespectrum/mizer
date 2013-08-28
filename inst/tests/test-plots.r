context("Plotting methods")

test_that("plot",{
    data(NS_species_params_gears)
    data(inter)
    params <- MizerParams(NS_species_params_gears, inter)
    sim <- project(params, effort=1, t_max=40, dt = 1, t_save = 1)
    # Just check that it doesn't crash
    plot(sim)
    #expect_that(plot(sim),!throws_error())    
})

