context("methods work for a single species data set")

# We choose the largest species for our single-species
params <- newMultispeciesParams(NS_species_params_gears[12, ])
n <- params@initial_n
npp <- params@initial_n_pp
effort <- array(abs(rnorm(10)), dim = c(10, 1))
sim1 <- project(params, effort = 1, t_max = 10)

test_that("project methods return arrays of correct dimension", {
    expect_length(dim(getEncounter(params, n, npp)), 2)
    expect_length(dim(getFeedingLevel(params, n, npp)), 2)
    expect_length(dim(getPredRate(params, n, npp)), 2)
    expect_length(dim(getPredMort(params, n, npp)), 2)
    expect_length(dim(getFMortGear(params, effort = 1)), 3)
    expect_length(dim(getFMortGear(params, effort = effort)), 4)
    expect_length(dim(getFMort(params, effort = 1)), 2)
    expect_length(dim(getFMort(params, effort = effort)), 3)
    expect_length(dim(getMort(params, n, npp, effort = 1)), 2)
    expect_length(dim(getEReproAndGrowth(params, n, npp)), 2)
    expect_length(dim(getERepro(params, n, npp)), 2)
    expect_length(dim(getEGrowth(params, n, npp)), 2)
})

test_that("summary methods return arrays of correct dimension", {
    expect_length(dim(get_size_range_array(params)), 2)
    expect_length(dim(getSSB(sim1)), 2)
    expect_length(dim(getBiomass(sim1)), 2)
    expect_length(dim(getN(sim1)), 2)
    expect_length(dim(getFMortGear(sim1)), 4)
    expect_length(dim(getYieldGear(sim1)), 3)
    expect_length(dim(getYield(sim1)), 2)
})

test_that("Can set up model with minimal information",{
    sp <- data.frame(species = "test")
    sp$w_inf <- 1000
    sp$k_vb <- 10
    params <- newMultispeciesParams(sp)
    sim <- project(params, t_max = 1)
})
