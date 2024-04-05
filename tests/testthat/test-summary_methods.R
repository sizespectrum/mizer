context("Summary methods")

## Initialisation ----
species_params <- NS_species_params_gears
species_params$pred_kernel_type <- "truncated_lognormal"
params <- newMultispeciesParams(species_params, inter, min_w_pp = 1e-12,
                                n = 2/3, p = 0.7, lambda = 2.8 - 2/3,
                                initial_effort = 1)
sim <- project(params, t_max = 10)
no_sp <- nrow(species_params)
no_w <- length(params@w)

# Random abundances
set.seed(0)
n <- abs(array(rnorm(no_w * no_sp), dim = c(no_sp, no_w)))
n_pp <- abs(rnorm(length(params@w_full)))

## helper functions ----
test_that("getRetainProbGear works", {
    # Choose example size
    w50_idx <- 50
    w50 <- params@w[w50_idx]
    l50 <- w2l(w50, params)
    l25_default <- 1.5 * l50
    w25_default <- l2w(l25_default, params)
    w25_idx <- 60
    w25 <- params@w[w25_idx]
    l25 <- w2l(w25, params)
    
    # If retain_l50 is missing, the return value should be all zero
    # with same dimensions and names as the selectivity matrix
    rpones <- params@selectivity
    rpones[] <- 1
    expect_identical(getRetainProbGear(params), rpones)
    gear_params(params)$retain_l50 <- NA
    expect_identical(getRetainProbGear(params), rpones)
    
    # Set retain_l50 for Sprat, Industrial trawl
    sp_idx <- which(params@species_params$species == "Sprat")
    gear_params(params)["Sprat, Industrial", "retain_l50"] <- l50[sp_idx]
    retain_prob <- getRetainProbGear(params)
    # Expect 0.5 at l50
    expect_equal(retain_prob["Industrial", "Sprat", w50_idx], 0.5)
    # Expect 0.25 at l25
    # Check that it is above 0.25 in the previous size class and below in the next
    w25_default_idx <- sum(params@w <= w25_default[sp_idx])
    expect_gt(retain_prob["Industrial", "Sprat", w25_default_idx - 1], 0.25)
    expect_lt(retain_prob["Industrial", "Sprat", w25_default_idx + 1], 0.25)
    # Expect entries for other gear-species pairs to be one
    retain_prob["Industrial", "Sprat", ] <- 1
    expect_identical(retain_prob, rpones)
    
    # set also retain_l25 for Sprat, Industrial trawl
    gear_params(params)["Sprat, Industrial", "retain_l25"] <- l25[sp_idx]
    retain_prob <- getRetainProbGear(params)
    # Expect 0.5 at l50
    expect_equal(retain_prob["Industrial", "Sprat", w50_idx], 0.5)
    # Expect 0.25 at l25
    # Check that it is above 0.25 in the previous size class and below in the next
    w25_idx <- sum(params@w <= w25[sp_idx])
    expect_gt(retain_prob["Industrial", "Sprat", w25_idx - 1], 0.25)
    expect_lt(retain_prob["Industrial", "Sprat", w25_idx + 1], 0.25)
    
    # Expect errors for invalid parameters
    # l50 must be smaller than l25
    gear_params(params)["Sprat, Industrial", "retain_l25"] <- l50[sp_idx]
    expect_error(getRetainProbGear(params),
                 "The value for `retain_l25` must always be larger than")
    # l50 can not be negative
    gear_params(params)["Sprat, Industrial", "retain_l50"] <- -1
    expect_error(getRetainProbGear(params),
                 "The value for `retain_l50` must be non-negative")
})

## get_size_range_array ----
test_that("get_size_range_array works", {
    params@species_params[["a"]] <- 
        c(0.007, 0.001, 0.009, 0.002, 0.010, 0.006, 0.008, 0.004,
            0.007, 0.005, 0.005, 0.007)
    params@species_params[["b"]] <- 
        c(3.014, 3.320, 2.941, 3.429, 2.986, 3.080, 3.019, 3.198,
            3.101, 3.160, 3.173, 3.075)
    
    # no limits
    size_n <- get_size_range_array(params)
    expect_true(all(size_n))
    
    # specifying weights
    size_n <- get_size_range_array(params, min_w = 1)
    expect_true(!all(size_n[, which(params@w < 1)]))
    expect_true(all(size_n[, which(params@w >= 1)]))
    size_n <- get_size_range_array(params, max_w = 100)
    expect_true(all(size_n[, which(params@w <= 100)]))
    expect_true(!all(size_n[, which(params@w > 1)]))
    size_n <- get_size_range_array(params, min_w = 1, max_w = 100)
    expect_true(!all(size_n[, which(params@w > 100)]))
    expect_true(!all(size_n[, which(params@w < 1)]))
    expect_true(all(size_n[, which((params@w >= 1) & (params@w <= 100))]))
    
    # specifying lengths
    min_l <- 2
    size_n <- get_size_range_array(params, min_l = min_l)
    min_w <- params@species_params$a * min_l ^ params@species_params$b
    for (sp in seq_len(nrow(params@species_params))) { 
        expect_true(all(size_n[sp, which(params@w >= min_w[sp])]))
        expect_true(!all(size_n[sp, which(params@w < min_w[sp])]))
    }
    max_l <- 100
    size_n <- get_size_range_array(params, max_l = max_l)
    max_w <- params@species_params$a * max_l ^ params@species_params$b
    for (sp in seq_len(nrow(params@species_params))) { 
        expect_true(all(size_n[sp, which(params@w <= max_w[sp])]))
        expect_true(!all(size_n[sp, which(params@w > max_w[sp])]))
    }
    size_n <- get_size_range_array(params, min_l = min_l, max_l = max_l)
    min_w <- params@species_params$a * min_l ^ params@species_params$b
    max_w <- params@species_params$a * max_l ^ params@species_params$b
    for (sp in seq_len(nrow(params@species_params))) { 
        expect_true(all(size_n[sp, which((params@w <= max_w[sp]) & 
                                             (params@w >= min_w[sp]))]))
        expect_true(!all(size_n[sp, which(params@w < min_w[sp])]))
        expect_true(!all(size_n[sp, which(params@w > max_w[sp])]))
    }
    
    # mixed weights and lengths
    size_n <- get_size_range_array(params, min_w = 1, max_l = max_l)
    min_w <- rep(1, nrow(params@species_params))
    for (sp in seq_len(nrow(params@species_params))) { 
        expect_true(all(size_n[sp, which((params@w <= max_w[sp]) & 
                                             (params@w >= min_w[sp]))]))
        expect_true(!all(size_n[sp, which(params@w < min_w[sp])]))
        expect_true(!all(size_n[sp, which(params@w > max_w[sp])]))
    }
    size_n <- get_size_range_array(params, min_l = min_l, max_w = 100)
    max_w <- rep(100, nrow(params@species_params))
    for (sp in seq_len(nrow(params@species_params))) { 
        expect_true(all(size_n[sp, which((params@w <= max_w[sp]) & 
                                             (params@w >= min_w[sp]))]))
        expect_true(!all(size_n[sp, which(params@w < min_w[sp])]))
        expect_true(!all(size_n[sp, which(params@w > max_w[sp])]))
    }
    
    # Gives expected error messages
    expect_error(get_size_range_array(params, min_w = 1000, max_w = 1),
                 "min_w must be less than max_w")
    expect_error(get_size_range_array(params, min_l = 1000, max_l = 1),
                 "min_w must be less than max_w")
    expect_error(get_size_range_array(params, min_l = 1000, max_w = 1),
                 "min_w must be less than max_w")
    expect_error(get_size_range_array(params, min_w = 1000, max_l = 1),
                 "min_w must be less than max_w")
    expect_error(get_size_range_array(params, min_l = 1:4, max_w = 10),
                 "min_l must be a single number or a vector")
    expect_error(get_size_range_array(params, min_l = 1, max_l = 1:10),
                 "max_l must be a single number or a vector")
    expect_error(get_size_range_array(params, min_w = 1:4, max_w = 10),
                 "min_w and max_w must be a single number of a vector")
})


# getProportionOfLargeFish ----
test_that("getProportionOfLargeFish works", {
    sim <- project(params, effort = 1, t_max = 20, dt = 0.5, t_save = 0.5)
    # noddy test - using full range of sizes
    prop <- getProportionOfLargeFish(sim, threshold_w = 500)
    time_idx <- 40
    threshold_w <- sim@params@w > 500
    total_biomass <- sum(sweep(sim@n[time_idx, , ], 2,
                               sim@params@w * sim@params@dw, "*")
                         )
    larger_biomass <- sum(sweep(sim@n[time_idx, , ], 2,
                                threshold_w * sim@params@w * sim@params@dw, "*")
                          )
    expect_that(prop[time_idx], 
                is_equivalent_to(larger_biomass / total_biomass))
    # using a size range
    prop <- getProportionOfLargeFish(sim, min_w = 10, max_w = 5000,
                                     threshold_w = 500)
    range_w <- (sim@params@w >= 10) & (sim@params@w <= 5000)
    threshold_w <- sim@params@w > 500
    total_biomass <- sum(sweep(sim@n[time_idx,,],2, range_w * sim@params@w * sim@params@dw, "*"))
    larger_biomass <- sum(sweep(sim@n[time_idx,,],2, threshold_w * range_w * sim@params@w * sim@params@dw, "*"))
    expect_that(prop[time_idx] , is_equivalent_to(larger_biomass / total_biomass))
    # numeric test
    expect_known_value(prop, "values/getProportionOfLargeFish")
})

# getMeanWeight ----
test_that("getMeanWeight works",{
    sim <- project(params, t_max = 20, dt = 0.5, t_save = 0.5)
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
    expect_equal(getYieldGear(sim)[1, , ], 
                 getYieldGear(sim@params))
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
    expect_equal(getYield(sim)[1, ], getYield(sim@params))
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
test_that("getDiet works with additional components", {
    params <- NS_params
    e <- globalenv()
    e$test_dyn <- function(params, ...) {
        111
    }
    # switch off satiation for easier test of result
    species_params(params)$h <- Inf
    p <- setComponent(params, "test", 1, 
                      dynamics_fun = "test_dyn",
                      encounter_fun = "test_dyn")
    
    diet1 <- getDiet(params, proportion = FALSE)
    diet2 <- getDiet(p, proportion = FALSE)
    expect_identical(diet1[, , 1:13], diet2[, , 1:13])
    expect_identical(diet2[1, 1, 14], 111)
})


# getSSB ----
test_that("getSSB works", {
    ssb <- getSSB(sim)
    expect_known_value(ssb, "values/getSSB")
    expect_equal(getSSB(sim)[1, ], getSSB(sim@params))
})

# getBiomass ----
test_that("getBiomass works", {
    biomass <- getBiomass(sim)
    expect_known_value(biomass, "values/getBiomass")
    expect_equal(getBiomass(sim)[1, ], getBiomass(sim@params))
})

# getN ----
test_that("getN works", {
    N <- getN(sim)
    expect_known_value(N, "values/getN")
    expect_equal(getN(sim)[1, ], getN(sim@params))
})

# getGrowthCurves ----
test_that("getGrowthCurves works with MizerSim", {
    ps <- setInitialValues(params, sim)
    expect_identical(getGrowthCurves(sim),
                     getGrowthCurves(ps))
})

# summary ----
test_that("summary works", {
    # Check that it works also with nonstandard kernel
    params@species_params$ppmr_min <- 100
    params@species_params$ppmr_max <- 10000
    params@species_params$beta <- NULL
    params@species_params$sigma <- NULL
    species_params(params)$pred_kernel_type <- "box"
    expect_output(summary(params),
                  'An object of class "MizerParams"')
    sim <- project(params, t_max = 0.1)
    expect_output(summary(sim),
                  'An object of class "MizerSim"')
})
