# We choose the largest species for our single-species
params <- NS_params_cod_small
n <- params@initial_n
npp <- params@initial_n_pp
effort <- array(abs(rnorm(10)), dim = c(10, 1))
sim1 <- project(params, effort = 1, t_max = 2)

test_that("project methods return arrays of correct dimension", {
    expect_length(dim(getEncounter(params, n, npp)), 2)
    expect_length(dim(getFeedingLevel(params, n, npp)), 2)
    expect_length(dim(getFeedingLevel(sim1)), 3)
    expect_length(dim(getPredRate(params, n, npp)), 2)
    expect_length(dim(getPredMort(params, n, npp)), 2)
    expect_length(dim(getPredMort(sim1)), 2)
    expect_length(dim(getFMortGear(params, effort = 1)), 3)
    expect_length(dim(getFMortGear(sim1)), 4)
    expect_length(dim(getFMortGear(params, effort = effort)), 4)
    expect_length(dim(getFMort(params, effort = 1)), 2)
    expect_length(dim(getFMort(sim1, drop = FALSE)), 3)
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

test_that("Can set up model with minimal information", {
    sp <- data.frame(species = "test",
                     stringsAsFactors = FALSE)
    sp$w_max <- 1000
    sp$k_vb <- 10
    params <- newMultispeciesParams(sp, info_level = 0)
    expect_error(project(params, t_max = 1), NA)
})

# rate getters with a single species ----

# One species only ----
test_that("project function returns objects of correct dimension when community only has one species", {
    params <- newCommunityParams(z0 = 0.2, f0 = 0.7, alpha = 0.2)
    t_max <- 2
    sim <- project(params, t_max = t_max, effort = 0)
    n <- array(sim@n[t_max + 1, , ], dim = dim(sim@n)[2:3])
    dimnames(n) <- dimnames(sim@n)[2:3]
    n_pp <- sim@n_pp[1, ]
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)
    # MizerParams functions
    expect_equal(dim(getEncounter(params, n, n_pp)), c(1, no_w))
    expect_equal(dim(getFeedingLevel(params, n, n_pp)), c(1, no_w))
    expect_equal(dim(getPredRate(params, n, n_pp)), c(1, no_w_full))
    expect_equal(dim(getPredMort(params, n, n_pp)), c(1, no_w))
    expect_length(getResourceMort(params, n, n_pp), no_w_full)
    expect_equal(dim(getFMortGear(params, 0)), c(1, 1, no_w)) # 3D time x species x size
    expect_equal(dim(getFMortGear(params, matrix(c(0, 0), nrow = 2))),
                     c(2, 1, 1, no_w)) # 4D time x gear x species x size
    expect_equal(dim(getFMort(params, 0)), c(1, no_w)) # 2D species x size
    expect_equal(dim(getFMort(params, matrix(c(0, 0), nrow = 2))),
                     c(2, 1, no_w)) # 3D time x species x size
    expect_equal(dim(getMort(params, n, n_pp, effort = 0)), c(1, no_w))
    expect_equal(dim(getEReproAndGrowth(params, n, n_pp)), c(1, no_w))
    expect_equal(dim(getERepro(params, n, n_pp)), c(1, no_w))
    expect_equal(dim(getEGrowth(params, n, n_pp)), c(1, no_w))
    expect_length(getRDI(params, n, n_pp), 1)
    expect_length(getRDD(params, n, n_pp), 1)
    expect_equal(dim(getDiffusion(params, n, n_pp)), c(1, no_w))

    # MizerSim functions
    # time x species x size
    expect_equal(dim(getFeedingLevel(sim)), c(t_max + 1, 1, no_w))
    # time x species x size - default drop is TRUE, if called from
    # plots drop = FALSE
    expect_equal(dim(getPredMort(sim)), c(t_max + 1, no_w))
    # time x species x size
    expect_equal(dim(getPredMort(sim, drop = FALSE)), c(t_max + 1, 1, no_w))
    # time x gear x species x size
    expect_equal(dim(getFMortGear(sim)), c(t_max + 1, 1, 1, no_w))
    # time x species x size - note drop = TRUE
    expect_equal(dim(getFMort(sim)), c(t_max + 1, no_w))
    # time x species x size
    expect_equal(dim(getFMort(sim, drop = FALSE)), c(t_max + 1, 1, no_w))
})
