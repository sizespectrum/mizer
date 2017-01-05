context("MizerSim constructor dimension checks")

test_that("basic constructor sets dimensions properly",{
    # check dimension against t input arguments
    data(NS_species_params_gears)
    data(inter)
    params <- MizerParams(NS_species_params_gears, inter)
    # Make MizerSims with t_max and t_save
    t_max <- 50
    t_save <- 1
    sim <- MizerSim(params, t_max = t_max, t_save = t_save)
    expect_that(dim(sim@effort)[1], equals(t_max / t_save))
    expect_that(dimnames(sim@effort)[[1]], is_identical_to(as.character(seq(from = t_save, to = t_max, by = t_save))))
    expect_that(all(dimnames(sim@effort)[[2]] %in% dimnames(params@selectivity)$gear), is_true())
    expect_that(dim(sim@n)[1], equals(1 + (t_max / t_save)))
    temp <- seq(from = t_save, to = t_max, by = t_save)
    expect_that(dimnames(sim@n)[[1]], is_identical_to(as.character(c(temp[1] - (temp[2]-temp[1]), temp))))

    t_max <- 20
    t_save <- 2
    sim <- MizerSim(params, t_max = t_max, t_save = t_save)
    expect_that(dim(sim@effort)[1], equals(t_max / t_save))
    expect_that(dimnames(sim@effort)[[1]], is_identical_to(as.character(seq(from = t_save, to = t_max, by = t_save))))
    expect_that(all(dimnames(sim@effort)[[2]] %in% dimnames(params@selectivity)$gear), is_true())
    expect_that(dim(sim@n)[1], equals(1 + (t_max / t_save)))
    temp <- seq(from = t_save, to = t_max, by = t_save)
    expect_that(dimnames(sim@n)[[1]], is_identical_to(as.character(c(temp[1] - (temp[2]-temp[1]), temp))))

    # Test for warning if t_max / t_save is not integer
    expect_error(MizerSim(params, t_max = 15, t_save = 2))

    # Make MizerSim using t_dimnames
    t_dimnames <- seq(from = 1990, to = 2000, by = 1)
    sim <- MizerSim(params, t_dimnames = t_dimnames)
    expect_that(dim(sim@effort)[1], equals(length(t_dimnames)))
    expect_that(dimnames(sim@effort)[[1]], is_identical_to(as.character(t_dimnames)))
    expect_that(all(dimnames(sim@effort)[[2]] %in% dimnames(params@selectivity)$gear), is_true())
    expect_that(dim(sim@n)[1], equals(1 + length(t_dimnames)))
    expect_that(dimnames(sim@n)[[1]], is_identical_to(as.character(c(t_dimnames[1] - (t_dimnames[2]-t_dimnames[1]), t_dimnames))))

    # Check error if t_dimnames is not numeric
    expect_that(MizerSim(params, t_dimnames = c("x","y","z")), throws_error())
    expect_that(MizerSim(params, t_dimnames = as.character(1:3)), throws_error())
})



