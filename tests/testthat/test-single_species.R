context("methods work for a single species data set")

data(inter)
data(NS_species_params_gears)
# We choose the largest species for our single-species
params <- MizerParams(NS_species_params_gears[12, ])
n <- params@initial_n
npp <- params@initial_n_pp
effort <- array(abs(rnorm(10)), dim = c(10, 1))
sim1 <- project(params, effort = 1, t_max = 10)

test_that("project methods return arrays of correct dimension", {
    expect_length(dim(getAvailEnergy(params, n, npp)), 2)
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
