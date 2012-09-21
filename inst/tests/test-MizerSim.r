context("MizerSim constructor dimension checks")

test_that("basic constructor sets dimensions properly",{
    # check dimension against t input arguments
    data(species_params_gears)
    data(inter)
    params <- MizerParams(species_params_gears, inter)
    t_max <- 50
    t_save <- 1
    sim <- MizerSim(params, t_max = t_max, t_save = t_save)
    expect_that(dim(sim@effort)[1], equals(t_max / t_save))
    t_max <- 20
    t_save <- 2
    sim <- MizerSim(params, t_max = t_max, t_save = t_save)
    expect_that(dim(sim@effort)[1], equals(t_max / t_save))
    expect_that(dim(sim@n)[1], equals(dim(sim@effort)[1] + 1))
    # Test for warning if t_max / t_save is not integer
    expect_error(MizerSim(params, t_max = 15, t_save = 2))
})



