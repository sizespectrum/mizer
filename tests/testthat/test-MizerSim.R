params <- newMultispeciesParams(NS_species_params_gears, inter)

# basic constructor sets dimensions properly ----
test_that("basic constructor sets dimensions properly",{
    # check dimension against t input arguments
    # Make MizerSims with t_max and t_save
    t_max <- 5
    t_save <- 1
    sim <- MizerSim(params, t_max = t_max, t_save = t_save)
    expect_equal(dim(sim@effort)[1], 1 + t_max / t_save)
    expect_identical(dimnames(sim@effort)[[1]], 
                     as.character(seq(from = 0, to = t_max, by = t_save)))
    expect_setequal(dimnames(sim@effort)[[2]], 
                    dimnames(params@selectivity)$gear)
    expect_equal(dim(sim@n)[1], 1 + (t_max / t_save))
    expect_identical(dimnames(sim@n)[[1]], 
                     dimnames(sim@effort)[[1]])

    t_max <- 4
    t_save <- 2
    sim <- MizerSim(params, t_max = t_max, t_save = t_save)
    expect_equal(dim(sim@effort)[1], 1 + t_max / t_save)
    expect_identical(dimnames(sim@effort)[[1]], 
                     as.character(seq(from = 0, to = t_max, by = t_save)))
    expect_setequal(dimnames(sim@effort)[[2]], 
                    dimnames(params@selectivity)$gear)
    expect_equal(dim(sim@n)[1], 1 + (t_max / t_save))
    expect_identical(dimnames(sim@n)[[1]], 
                     dimnames(sim@effort)[[1]])

    # Make MizerSim using t_dimnames
    t_dimnames <- seq(from = 1990, to = 2000, by = 1)
    sim <- MizerSim(params, t_dimnames = t_dimnames)
    expect_equal(dim(sim@effort)[1], length(t_dimnames))
    expect_identical(dimnames(sim@effort)[[1]], as.character(t_dimnames))
    expect_setequal(dimnames(sim@effort)[[2]],
                    dimnames(params@selectivity)$gear)
    expect_equal(dim(sim@n)[1], length(t_dimnames))
    expect_identical(dimnames(sim@n)[[1]], as.character(t_dimnames))

    # Check error if t_dimnames is not numeric or not sorted
    expect_error(MizerSim(params, t_dimnames = c("x","y","z")),
                 "The t_dimnames argument must be numeric")
    expect_error(MizerSim(params, t_dimnames = as.character(1:3)),
                 "The t_dimnames argument must be numeric")
    expect_error(MizerSim(params, t_dimnames = 3:1),
                 "The t_dimnames argument should be increasing")
})

# upgradeSim leaves upgraded sim unchanged ----
test_that("upgradeSim leaves upgraded sim unchanged", {
    sim <- project(params, t_max = 0.1, t_save = 0.1)
    comment(sim) <- "sim"
    comment(sim@params) <- "params"
    comment(sim@params@intake_max) <- "intake_max"
    comment(sim@n) <- "n"
    expect_identical(sim, upgradeSim(sim))
})
