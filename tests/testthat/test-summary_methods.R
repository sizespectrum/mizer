context("Summary methods")

## Initialisation ----
species_params <- NS_species_params_gears
species_params$pred_kernel_type <- "truncated_lognormal"
params <- newMultispeciesParams(species_params, inter, min_w_pp = 1e-12,
                                n = 2/3, p = 0.7, lambda = 2.8 - 2/3)
sim <- project(params, effort = 1, t_max = 10)
no_sp <- nrow(species_params)
no_w <- length(params@w)

# Random abundances
set.seed(0)
n <- abs(array(rnorm(no_w * no_sp), dim = c(no_sp, no_w)))
n_pp <- abs(rnorm(length(params@w_full)))

## get_size_range_array ----
test_that("get_size_range_array",{
    NS_species_params_gears[["a"]] <- 0.01
    NS_species_params_gears[["b"]] <- 3
    params <- newMultispeciesParams(NS_species_params_gears, inter)
    size_n <- get_size_range_array(params)
    expect_true(all(size_n))
    size_n <- get_size_range_array(params, min_w = 1)
    expect_true(!all(size_n[,which(params@w < 1)]))
    expect_true(all(size_n[,which(params@w >= 1)]))
    size_n <- get_size_range_array(params, max_w = 100)
    expect_true(all(size_n[,which(params@w <= 100)]))
    expect_true(!all(size_n[,which(params@w > 1)]))
    size_n <- get_size_range_array(params, min_w = 1, max_w = 100)
    expect_true(!all(size_n[,which(params@w > 100)]))
    expect_true(!all(size_n[,which(params@w < 1)]))
    expect_true(all(size_n[,which((params@w >= 1) & (params@w<=100))]))
    size_n <- get_size_range_array(params, min_l = 1)

    min_w <- params@species_params$a * 1 ^ params@species_params$b
    for (sp in 1:nrow(params@species_params)){ 
        expect_true(all(size_n[sp,which(params@w >= min_w[sp])]))
        expect_true(!all(size_n[sp,which(params@w < min_w[sp])]))
    }
    size_n <- get_size_range_array(params, max_l = 100)
    max_w <- params@species_params$a * 100 ^ params@species_params$b
    for (sp in 1:nrow(params@species_params)){ 
        expect_true(all(size_n[sp,which(params@w <= max_w[sp])]))
        expect_true(!all(size_n[sp,which(params@w > max_w[sp])]))
    }
    size_n <- get_size_range_array(params, min_l = 1, max_l = 100)
    min_w <- params@species_params$a * 1 ^ params@species_params$b
    max_w <- params@species_params$a * 100 ^ params@species_params$b
    for (sp in 1:nrow(params@species_params)){ 
        expect_true(all(size_n[sp,which((params@w <= max_w[sp]) & (params@w >= min_w[sp]))]))
        expect_true(!all(size_n[sp,which(params@w < min_w[sp])]))
        expect_true(!all(size_n[sp,which(params@w > max_w[sp])]))
    }
    size_n <- get_size_range_array(params, min_w = 1, max_l = 100)
    min_w <- rep(1,nrow(params@species_params))
    max_w <- params@species_params$a * 100 ^ params@species_params$b
    for (sp in 1:nrow(params@species_params)){ 
        expect_true(all(size_n[sp,which((params@w <= max_w[sp]) & (params@w >= min_w[sp]))]))
        expect_true(!all(size_n[sp,which(params@w < min_w[sp])]))
        expect_true(!all(size_n[sp,which(params@w > max_w[sp])]))
    }
    size_n <- get_size_range_array(params, min_l = 1, max_w = 100)
    min_w <- params@species_params$a * 1 ^ params@species_params$b
    max_w <- rep(100,nrow(params@species_params))
    for (sp in 1:nrow(params@species_params)){ 
        expect_true(all(size_n[sp,which((params@w <= max_w[sp]) & (params@w >= min_w[sp]))]))
        expect_true(!all(size_n[sp,which(params@w < min_w[sp])]))
        expect_true(!all(size_n[sp,which(params@w > max_w[sp])]))
    }
    expect_that(get_size_range_array(params, min_w = 1000, max_w = 1), throws_error())
    expect_that(get_size_range_array(params, min_l = 1000, max_l = 1), throws_error())
    expect_that(get_size_range_array(params, min_l = 1000, max_w = 1), throws_error())
    expect_that(get_size_range_array(params, min_w = 1000, max_l = 1), throws_error())
    # checking if fails if a and b not in species_params
    no_ab_params <- params
    no_ab_params@species_params <- params@species_params[,!(names(params@species_params) %in% c("a","b"))]
    expect_that(get_size_range_array(no_ab_params, min_l = 1, max_w = 100), throws_error())
})

# get_time_elements ----
test_that("get_time_elements",{
    sim <- project(params, effort=1, t_max=10, dt = 0.5, t_save = 0.5)
    expect_equal(length(get_time_elements(sim, as.character(3:4))),
                 dim(sim@n)[1])
    expect_equal(length(get_time_elements(sim, 3:4)),
                 dim(sim@n)[1])
    expect_that(sum(get_time_elements(sim,3:4)), equals(3))
    expect_that(sum(get_time_elements(sim,3:50)), throws_error())
    expect_that(which(get_time_elements(sim,seq(from=3,to=4,by = 0.1))), is_equivalent_to(c(7,8,9)))
    expect_that(length(get_time_elements(sim,seq(from=3,to=4,by = 0.1), slot_name="effort")), equals(dim(sim@effort)[1]))
    # What if real years are used
    effort <- array(1, dim = c(19,4), dimnames = list(year = seq(from = 1960, to = 1969, by = 0.5), gear = c("Industrial","Pelagic","Otter","Beam")))
    sim <- project(params, effort = effort, t_save = 0.5)
    expect_that(which(get_time_elements(sim,1965)), is_equivalent_to(11))
    expect_that(which(get_time_elements(sim,"1965")), is_equivalent_to(11))
    expect_that(which(get_time_elements(sim,1965:1969)), is_equivalent_to(11:19))
})


# getProportionOfLargeFish ----
test_that("getProportionOfLargeFish works",{
    sim <- project(params, effort = 1, t_max = 20, dt = 0.5, t_save = 0.5)
    # noddy test - using full range of sizes
    prop <- getProportionOfLargeFish(sim, threshold_w = 500)
    time_idx <- 40
    threshold_w <- sim@params@w > 500
    total_biomass <- sum(sweep(sim@n[time_idx,,],2, sim@params@w * sim@params@dw, "*"))
    larger_biomass <- sum(sweep(sim@n[time_idx,,],2, threshold_w * sim@params@w * sim@params@dw, "*"))
    expect_that(prop[time_idx] , is_equivalent_to(larger_biomass / total_biomass))
    # using a size range
    prop <- getProportionOfLargeFish(sim, min_w = 10, max_w = 5000, threshold_w = 500)
    range_w <- (sim@params@w >= 10) & (sim@params@w <= 5000)
    threshold_w <- sim@params@w > 500
    total_biomass <- sum(sweep(sim@n[time_idx,,],2, range_w * sim@params@w * sim@params@dw, "*"))
    larger_biomass <- sum(sweep(sim@n[time_idx,,],2, threshold_w * range_w * sim@params@w * sim@params@dw, "*"))
    expect_that(prop[time_idx] , is_equivalent_to(larger_biomass / total_biomass))
    # numeric test
    expect_known_value(prop, "values/getProportionOfLargeFish")
})


# check_species ----
test_that("check_species works",{
    sim <- project(params, effort=1, t_max=20, dt = 0.5, t_save = 0.5)
    expect_true(check_species(sim,c("Cod","Haddock")))
    expect_true(check_species(sim,c(10,11)))
    expect_that(check_species(sim,c("Arse","Balls")), throws_error())
    expect_that(check_species(sim,c(10,666)), throws_error())

})


# getMeanWeight ----
test_that("getMeanWeight works",{
    sim <- project(params, effort=1, t_max=20, dt = 0.5, t_save = 0.5)
    # all species, all size range
    total_biomass <- apply(sweep(sim@n, 3, sim@params@w * sim@params@dw, "*"),1,sum)
    total_n  <- apply(sweep(sim@n, 3, sim@params@dw, "*"),1,sum)
    mw1 <- total_biomass / total_n
    mw <- getMeanWeight(sim)
    expect_that(mw, equals(mw1))
    # select species
    species <- c("Cod","Haddock")
    total_biomass <- apply(sweep(sim@n[,species,], 3, sim@params@w * sim@params@dw, "*"),1,sum)
    total_n  <- apply(sweep(sim@n[,species,], 3, sim@params@dw, "*"),1,sum)
    mw2 <- total_biomass / total_n
    mw <- getMeanWeight(sim, species = species)
    expect_that(mw, equals(mw2))
    # select size range
    min_w <- 10
    max_w <- 10000
    size_n <- get_size_range_array(sim@params, min_w = min_w, max_w = max_w)
    total_biomass <- apply(sweep(sweep(sim@n, c(2,3), size_n, "*"), 3, sim@params@w * sim@params@dw, "*"),1,sum)
    total_n <- apply(sweep(sweep(sim@n, c(2,3), size_n, "*"), 3, sim@params@dw, "*"),1,sum)
    mw3 <- total_biomass / total_n
    mw <- getMeanWeight(sim, min_w = min_w, max_w=max_w)
    expect_that(mw, equals(mw3))
    # select size range and species
    total_biomass <- apply(sweep(sweep(sim@n, c(2,3), size_n, "*")[,species,], 3, sim@params@w * sim@params@dw, "*"),1,sum)
    total_n <- apply(sweep(sweep(sim@n, c(2,3), size_n, "*")[,species,], 3, sim@params@dw, "*"),1,sum)
    mw4 <- total_biomass / total_n
    mw <- getMeanWeight(sim, species=species, min_w = min_w, max_w=max_w)
    expect_that(mw, equals(mw4))
    # errors
    expect_that(getMeanWeight(sim,species=c("Dougal","Ted")), throws_error())
    # numeric test
    expect_known_value(mw, "values/getMeanWeight")
})

# getMeanMaxWeight ----
test_that("getMeanMaxWeight works", {
    expect_error(getMeanMaxWeight(sim, measure = NA),
                 "measure must be one of")
    expect_known_value(getMeanMaxWeight(sim, measure = "both"),
                       "values/getMeanMaxWeight")
})


# getYieldGear ----
test_that("getYieldGear works",{
    y <- getYieldGear(sim)
    # check dims
    expect_that(dim(y),equals(c(11,dim(params@catchability)[1],dim(params@catchability)[2])))
    # check a value and assume the others are right
    biomass <- sweep(sim@n,3,sim@params@w * sim@params@dw, "*")
    f_gear <- getFMortGear(sim)
    expect_that(sum((biomass*f_gear[,1,,])[1,1,]),equals(y[1,1,1]))
    # numeric test
    expect_known_value(y, "values/getYieldGear")
})


# getYield ----
test_that("getYield works",{
    y <- getYield(sim)
    # check dims
    expect_that(dim(y),equals(c(11,dim(params@catchability)[2])))
    # check a value and assume the others are right
    biomass <- sweep(sim@n,3,sim@params@w * sim@params@dw, "*")
    f <- getFMort(sim)
    expect_that(sum((f*biomass)[1,1,]),equals(y[1,1]))
    # numeric test
    expect_known_value(y, "values/getYield")
})


# getCommunitySlope ----
test_that("getCommunitySlope works",{
    slope_b <- getCommunitySlope(sim)
    # dims
    expect_that(dim(slope_b), equals(c(dim(sim@n)[1],3)))
    # sum biomasses
    biomass <- apply(sweep(sim@n,3,sim@params@w,"*"),c(1,3),sum)
    # r2, slope and intercept at last time step
    lm_res <- lm(log(biomass[dim(sim@n)[1],]) ~ log(sim@params@w))
    expect_that(slope_b[dim(sim@n)[1],"r2"],equals(summary(lm_res)$r.squared))
    expect_that(slope_b[dim(sim@n)[1],"slope"],equals(summary(lm_res)$coefficients[2,1]))
    expect_that(slope_b[dim(sim@n)[1],"intercept"],equals(summary(lm_res)$coefficients[1,1]))
    # Test just numbers not biomass
    slope_n <- getCommunitySlope(sim, biomass=FALSE)
    expect_that(dim(slope_n), equals(c(dim(sim@n)[1],3)))
    # sum numbers
    numbers <- apply(sim@n,c(1,3),sum)
    # r2, slope and intercept at last time step
    lm_res <- lm(log(numbers[dim(sim@n)[1],]) ~ log(sim@params@w))
    expect_that(slope_n[dim(sim@n)[1],"r2"],equals(summary(lm_res)$r.squared))
    expect_that(slope_n[dim(sim@n)[1],"slope"],equals(summary(lm_res)$coefficients[2,1]))
    expect_that(slope_n[dim(sim@n)[1],"intercept"],equals(summary(lm_res)$coefficients[1,1]))
    # Check the sizes
    slope_b2 <- slope_biomass <- getCommunitySlope(sim, min_w = 10, max_w = 10000)
    sizes <- (sim@params@w >= 10) & (sim@params@w <= 10000)
    biomass <- apply(sweep(sim@n,3,sim@params@w,"*"),c(1,3),sum)
    # r2, slope and intercept at last time step
    lm_res <- lm(log(biomass[dim(sim@n)[1],sizes]) ~ log(sim@params@w[sizes]))
    expect_that(slope_b2[dim(sim@n)[1],"r2"],equals(summary(lm_res)$r.squared))
    expect_that(slope_b2[dim(sim@n)[1],"slope"],equals(summary(lm_res)$coefficients[2,1]))
    expect_that(slope_b2[dim(sim@n)[1],"intercept"],equals(summary(lm_res)$coefficients[1,1]))
    # Check the species
    dem_species <- c("Dab","Whiting","Sole","Gurnard","Plaice","Haddock", "Cod","Saithe")
    slope_b3 <- getCommunitySlope(sim, species = dem_species)
    biomass <- apply(sweep(sim@n[,dem_species,],3,sim@params@w,"*"),c(1,3),sum)
    # r2, slope and intercept at last time step
    lm_res <- lm(log(biomass[dim(sim@n)[1],]) ~ log(sim@params@w))
    expect_that(slope_b3[dim(sim@n)[1],"r2"],equals(summary(lm_res)$r.squared))
    expect_that(slope_b3[dim(sim@n)[1],"slope"],equals(summary(lm_res)$coefficients[2,1]))
    expect_that(slope_b3[dim(sim@n)[1],"intercept"],equals(summary(lm_res)$coefficients[1,1]))
    # numeric test
    expect_known_value(slope_b3, "values/getCommunitySlope")
})


# getDiet ----
test_that("getDiet works with proportion = FALSE", {
    diet <- getDiet(params, n, n_pp, proportion = FALSE)
    expect_known_value(diet, "values/getDiet")
    # Check that summing over all species and resource gives 
    # total consumption
    consumption <- rowSums(diet, dims = 2)
    encounter <- getEncounter(params, n, n_pp)
    feeding_level <- getFeedingLevel(params, n, n_pp)
    expect_equivalent(consumption, encounter * (1 - feeding_level))
    # Check that using pred kernel instead of FFT gives the same result
    params <- setPredKernel(params, pred_kernel = getPredKernel(params))
    expect_equal(diet, getDiet(params, n, n_pp, proportion = FALSE))
})

test_that("getDiet works with proportion = TRUE", {
    diet <- getDiet(params)
    total <- rowSums(diet, dims = 2)
    ones <- total
    # Only check at sizes where there are actually fish
    ones[] <- as.numeric(params@initial_n > 0)
    expect_equal(total, ones)
})


# getSSB ----
test_that("getSSB works", {
    ssb <- getSSB(sim)
    expect_known_value(ssb, "values/getSSB")
})

# getBiomass ----
test_that("getBiomass works", {
    biomass <- getBiomass(sim)
    expect_known_value(biomass, "values/getBiomass")
})

# getN ----
test_that("getN works", {
    N <- getN(sim)
    expect_known_value(N, "values/getN")
})

# getGrowthCurves ----
test_that("getGrowthCurves works with MizerSim", {
    ps <- setInitialValues(params, sim)
    expect_identical(getGrowthCurves(sim),
                     getGrowthCurves(ps))
})

# summary ----
test_that("summary works", {
    expect_output(summary(params),
                  'An object of class "MizerParams"')
    expect_output(summary(sim),
                  'An object of class "MizerSim"')
})