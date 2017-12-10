context("methods work for a single species data set")

test_that("project methods",{
    data(inter)
    data(NS_species_params_gears)
    # Multiple species params
    params <- MizerParams(NS_species_params_gears, inter=inter)
    # We choose the largest species for our single-species
    single_params <- MizerParams(NS_species_params_gears[12,])
    # Multiple species, single effort sim
    sim1 <- project(params, effort = 1, t_max = 10)
    n_mult <- sim1@n[11,,]
    n_single <- matrix(sim1@n[11,12,],nrow=1)
    dimnames(n_single) <- list(sp = "Sprat", w = dimnames(sim1@n)$w)
    npp <- sim1@n_pp[11,]
    # Multiple species, mult effort sim
    effort <- array(abs(rnorm(40)),dim=c(10,4))
    single_effort <- array(abs(rnorm(10)),dim=c(10,1))

    expect_that(length(dim(getPhiPrey(params,n_mult,npp))), is_identical_to(length(dim(getPhiPrey(single_params,n_single,npp)))))
    expect_that(length(dim(getFeedingLevel(params,n_mult,npp))), is_identical_to(length(dim(getFeedingLevel(single_params,n_single,npp)))))
    expect_that(length(dim(getPredRate(params,n_mult,npp))), is_identical_to(length(dim(getPredRate(single_params,n_single,npp)))))
    expect_that(length(dim(getM2(params,n_mult,npp))), is_identical_to(length(dim(getM2(single_params,n_single,npp)))))
    expect_that(length(getM2Background(params,n_mult,npp)), is_identical_to(length(getM2Background(single_params,n_single,npp))))
    expect_that(length(dim(getFMortGear(params,effort = 1))), is_identical_to(length(dim(getFMortGear(single_params,effort = 1)))))
    expect_that(length(dim(getFMortGear(params,effort = effort))), is_identical_to(length(dim(getFMortGear(single_params,effort = single_effort)))))
    expect_that(length(dim(getFMort(params,effort = 1))), is_identical_to(length(dim(getFMort(single_params,effort = 1)))))
    expect_that(length(dim(getFMort(params,effort = effort))), is_identical_to(length(dim(getFMort(single_params,effort = single_effort)))))
    expect_that(length(dim(getZ(params,n_mult,npp,effort=1))), is_identical_to(length(dim(getZ(single_params,n_single,npp, effort=1)))))
    expect_that(length(dim(getEReproAndGrowth(params,n_mult,npp))), is_identical_to(length(dim(getEReproAndGrowth(single_params,n_single,npp)))))
    expect_that(length(dim(getESpawning(params,n_mult,npp))), is_identical_to(length(dim(getESpawning(single_params,n_single,npp)))))
    expect_that(length(dim(getEGrowth(params,n_mult,npp))), is_identical_to(length(dim(getEGrowth(single_params,n_single,npp)))))
    expect_that(length(dim(getRDI(params,n_mult,npp))), is_identical_to(length(dim(getRDI(single_params,n_single,npp)))))
    expect_that(length(dim(getRDD(params,n_mult,npp))), is_identical_to(length(dim(getRDD(single_params,n_single,npp)))))
})

test_that("summary methods",{
    data(inter)
    data(NS_species_params_gears)
    # Multiple species params
    params <- MizerParams(NS_species_params_gears, inter=inter)
    single_params <- MizerParams(NS_species_params_gears[1,])
    # Multiple species, single effort sim
    sim1 <- project(params, effort = 1, t_max = 10)
    n_mult <- sim1@n[11,,]
    n_single <- matrix(sim1@n[11,1,],nrow=1)
    dimnames(n_single) <- list(sp = "Sprat", w = dimnames(sim1@n)$w)
    npp <- sim1@n_pp[11,]
    sim1s <- project(single_params, effort = 1, t_max = 10)

    expect_that(length(dim(get_size_range_array(params))), is_identical_to(length(dim(get_size_range_array(single_params)))))
    expect_that(length(dim(getSSB(sim1))), is_identical_to(length(dim(getSSB(sim1s)))))
    expect_that(length(dim(getBiomass(sim1))), is_identical_to(length(dim(getBiomass(sim1s)))))
    expect_that(length(dim(getN(sim1))), is_identical_to(length(dim(getN(sim1s)))))
    expect_that(length(dim(getFMortGear(sim1))), is_identical_to(length(dim(getFMortGear(sim1s)))))
    expect_that(length(dim(getYieldGear(sim1))), is_identical_to(length(dim(getYieldGear(sim1s)))))
    expect_that(length(dim(getYield(sim1))), is_identical_to(length(dim(getYield(sim1s)))))
    expect_that(length(getProportionOfLargeFish(sim1, threshold_w = 1)), is_identical_to(length(getProportionOfLargeFish(sim1s, threshold_w = 1))))
    expect_that(length(getMeanWeight(sim1, min_w = 1)), is_identical_to(length(getMeanWeight(sim1s, min_w = 1))))
    expect_that(length(getMeanMaxWeight(sim1, min_w = 1)), is_identical_to(length(getMeanMaxWeight(sim1s, min_w = 1))))
})
