## Initialisation ----
species_params <- NS_species_params_gears
species_params$pred_kernel_type <- "truncated_lognormal"
params <- newMultispeciesParams(species_params, inter, min_w_pp = 1e-12,
                                n = 2/3, p = 0.7, lambda = 2.8 - 2/3,
                                initial_effort = 1, info_level = 0)
sim <- project(params, t_max = 10)
no_sp <- nrow(species_params)
no_w <- length(params@w)

# Random abundances
set.seed(0)
n <- abs(array(rnorm(no_w * no_sp), dim = c(no_sp, no_w)))
n_pp <- abs(rnorm(length(params@w_full)))

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

    # Return all FALSE if no sizes in range
    expect_true(all(get_size_range_array(params, min_w = 1000, max_w = 1) == FALSE))
    expect_true(all(get_size_range_array(params, min_l = 1000, max_l = 1) == FALSE))
    expect_true(all(get_size_range_array(params, min_l = 1000, max_w = 1) == FALSE))
    expect_true(all(get_size_range_array(params, min_w = 1000, max_l = 1) == FALSE))

    # Expect errors
    expect_error(get_size_range_array(params, min_l = 1:4, max_w = 10),
                 "min_l must be a single number or a vector")
    expect_error(get_size_range_array(params, min_l = 1, max_l = 1:10),
                 "max_l must be a single number or a vector")
    expect_error(get_size_range_array(params, min_w = 1:4, max_w = 10),
                 "min_w and max_w must be a single number of a vector")
    # checking if fails if a and b not in species_params
    no_ab_params <- params
    no_ab_params@species_params$a[1] <- NA
    expect_error(get_size_range_array(no_ab_params, min_l = 1, max_w = 100),
                 "There must be no NAs in the species_params columns 'a' and 'b'")
    no_ab_params@species_params <-
        params@species_params[, !(names(params@species_params) %in% c("a", "b"))]
    expect_error(get_size_range_array(no_ab_params, min_l = 1, max_w = 100),
                 "pecies_params slot must have columns 'a' and 'b'")
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
    expect_equal(prop[time_idx], larger_biomass / total_biomass, ignore_attr = TRUE)
    # using a size range
    prop <- getProportionOfLargeFish(sim, min_w = 10, max_w = 5000,
                                     threshold_w = 500)
    range_w <- (sim@params@w >= 10) & (sim@params@w <= 5000)
    threshold_w <- sim@params@w > 500
    total_biomass <- sum(sweep(sim@n[time_idx,,],2, range_w * sim@params@w * sim@params@dw, "*"))
    larger_biomass <- sum(sweep(sim@n[time_idx,,],2, threshold_w * range_w * sim@params@w * sim@params@dw, "*"))
    expect_equal(prop[time_idx] , larger_biomass / total_biomass, ignore_attr = TRUE)
    # numeric test
    # expect_known_value(prop, "values/getProportionOfLargeFish")
    expect_snapshot(prop)
})

# getMeanWeight ----
test_that("getMeanWeight works",{
    sim <- project(params, t_max = 20, dt = 0.5, t_save = 0.5)
    # all species, all size range
    total_biomass <- apply(sweep(sim@n, 3, sim@params@w * sim@params@dw, "*"),1,sum)
    total_n  <- apply(sweep(sim@n, 3, sim@params@dw, "*"),1,sum)
    mw1 <- total_biomass / total_n
    mw <- getMeanWeight(sim)
    expect_equal(mw, mw1, ignore_attr = TRUE)
    # select species
    species <- c("Cod","Haddock")
    total_biomass <- apply(sweep(sim@n[,species,], 3, sim@params@w * sim@params@dw, "*"),1,sum)
    total_n  <- apply(sweep(sim@n[,species,], 3, sim@params@dw, "*"),1,sum)
    mw2 <- total_biomass / total_n
    mw <- getMeanWeight(sim, species = species)
    expect_equal(mw, mw2, ignore_attr = TRUE)
    # select size range
    min_w <- 10
    max_w <- 10000
    size_n <- get_size_range_array(sim@params, min_w = min_w, max_w = max_w)
    total_biomass <- apply(sweep(sweep(sim@n, c(2,3), size_n, "*"), 3, sim@params@w * sim@params@dw, "*"),1,sum)
    total_n <- apply(sweep(sweep(sim@n, c(2,3), size_n, "*"), 3, sim@params@dw, "*"),1,sum)
    mw3 <- total_biomass / total_n
    mw <- getMeanWeight(sim, min_w = min_w, max_w=max_w)
    expect_equal(mw, mw3, ignore_attr = TRUE)
    # select size range and species
    total_biomass <- apply(sweep(sweep(sim@n, c(2,3), size_n, "*")[,species,], 3, sim@params@w * sim@params@dw, "*"),1,sum)
    total_n <- apply(sweep(sweep(sim@n, c(2,3), size_n, "*")[,species,], 3, sim@params@dw, "*"),1,sum)
    mw4 <- total_biomass / total_n
    mw <- getMeanWeight(sim, species=species, min_w = min_w, max_w=max_w)
    expect_equal(mw, mw4, ignore_attr = TRUE)
    # numeric test
    # expect_known_value(mw, "values/getMeanWeight")
    expect_snapshot(mw)
})

# getMeanMaxWeight ----
test_that("getMeanMaxWeight works", {
    expect_error(getMeanMaxWeight(sim, measure = NA),
                 "measure must be one of")
    # expect_known_value(getMeanMaxWeight(sim, measure = "both"),
    #                    "values/getMeanMaxWeight")
    expect_snapshot(getMeanMaxWeight(sim, measure = "both"))

})


# getYieldGear ----
test_that("getYieldGear works",{
    y <- getYieldGear(sim)
    # check dims
    expect_equal(dim(y), c(11,dim(params@catchability)[1],dim(params@catchability)[2]), ignore_attr = TRUE)
    # check a value and assume the others are right
    biomass <- sweep(sim@n,3,sim@params@w * sim@params@dw, "*")
    f_gear <- getFMortGear(sim)
    expect_equal(sum((biomass*f_gear[,1,,])[1,1,]), y[1,1,1], ignore_attr = TRUE)
    # numeric test
    # expect_known_value(y, "values/getYieldGear")
    expect_snapshot(y)
    expect_equal(getYieldGear(sim)[1, , ],
                 getYieldGear(sim@params))
})


# getYield ----
test_that("getYield works",{
    y <- getYield(sim)
    # check dims
    expect_equal(dim(y), c(11,dim(params@catchability)[2]), ignore_attr = TRUE)
    # check a value and assume the others are right
    biomass <- sweep(sim@n,3,sim@params@w * sim@params@dw, "*")
    f <- getFMort(sim)
    expect_equal(sum((f*biomass)[1,1,]), y[1,1], ignore_attr = TRUE)
    # numeric test
    # expect_known_value(y, "values/getYield")
    expect_snapshot(y)
    expect_equal(getYield(sim)[1, ], getYield(sim@params))
})


# getCommunitySlope ----
test_that("getCommunitySlope works",{
    slope_b <- getCommunitySlope(sim)
    # dims
    expect_equal(dim(slope_b), c(dim(sim@n)[1],3), ignore_attr = TRUE)
    # sum biomasses
    biomass <- apply(sweep(sim@n,3,sim@params@w,"*"),c(1,3),sum)
    # r2, slope and intercept at last time step
    lm_res <- lm(log(biomass[dim(sim@n)[1],]) ~ log(sim@params@w))
    expect_equal(slope_b[dim(sim@n)[1],"r2"], summary(lm_res)$r.squared, ignore_attr = TRUE)
    expect_equal(slope_b[dim(sim@n)[1],"slope"], summary(lm_res)$coefficients[2,1], ignore_attr = TRUE)
    expect_equal(slope_b[dim(sim@n)[1],"intercept"], summary(lm_res)$coefficients[1,1], ignore_attr = TRUE)
    # Test just numbers not biomass
    slope_n <- getCommunitySlope(sim, biomass=FALSE)
    expect_equal(dim(slope_n),  c(dim(sim@n)[1],3), ignore_attr = TRUE)
    # sum numbers
    numbers <- apply(sim@n,c(1,3),sum)
    # r2, slope and intercept at last time step
    lm_res <- lm(log(numbers[dim(sim@n)[1],]) ~ log(sim@params@w))
    expect_equal(slope_n[dim(sim@n)[1],"r2"], summary(lm_res)$r.squared, ignore_attr = TRUE)
    expect_equal(slope_n[dim(sim@n)[1],"slope"], summary(lm_res)$coefficients[2,1], ignore_attr = TRUE)
    expect_equal(slope_n[dim(sim@n)[1],"intercept"], summary(lm_res)$coefficients[1,1], ignore_attr = TRUE)
    # Check the sizes
    slope_b2 <- slope_biomass <- getCommunitySlope(sim, min_w = 10, max_w = 10000)
    sizes <- (sim@params@w >= 10) & (sim@params@w <= 10000)
    biomass <- apply(sweep(sim@n,3,sim@params@w,"*"),c(1,3),sum)
    # r2, slope and intercept at last time step
    lm_res <- lm(log(biomass[dim(sim@n)[1],sizes]) ~ log(sim@params@w[sizes]))
    expect_equal(slope_b2[dim(sim@n)[1],"r2"], summary(lm_res)$r.squared, ignore_attr = TRUE)
    expect_equal(slope_b2[dim(sim@n)[1],"slope"], summary(lm_res)$coefficients[2,1], ignore_attr = TRUE)
    expect_equal(slope_b2[dim(sim@n)[1],"intercept"], summary(lm_res)$coefficients[1,1], ignore_attr = TRUE)
    # Check the species
    dem_species <- c("Dab","Whiting","Sole","Gurnard","Plaice","Haddock", "Cod","Saithe")
    slope_b3 <- getCommunitySlope(sim, species = dem_species)
    biomass <- apply(sweep(sim@n[,dem_species,],3,sim@params@w,"*"),c(1,3),sum)
    # r2, slope and intercept at last time step
    lm_res <- lm(log(biomass[dim(sim@n)[1],]) ~ log(sim@params@w))
    expect_equal(slope_b3[dim(sim@n)[1],"r2"], summary(lm_res)$r.squared, ignore_attr = TRUE)
    expect_equal(slope_b3[dim(sim@n)[1],"slope"], summary(lm_res)$coefficients[2,1], ignore_attr = TRUE)
    expect_equal(slope_b3[dim(sim@n)[1],"intercept"], summary(lm_res)$coefficients[1,1], ignore_attr = TRUE)
    # numeric test
    # expect_known_value(slope_b3, "values/getCommunitySlope")
    expect_snapshot(slope_b3)
})


# getDiet ----
test_that("getDiet works with proportion = FALSE", {
    diet <- getDiet(params, n, n_pp, proportion = FALSE)
    # expect_known_value(diet, "values/getDiet")
    expect_snapshot(diet)
    # Check that summing over all species and resource gives
    # total consumption
    consumption <- rowSums(diet, dims = 2)
    encounter <- getEncounter(params, n, n_pp)
    feeding_level <- getFeedingLevel(params, n, n_pp)
    expect_equal(consumption, encounter * (1 - feeding_level), ignore_attr = TRUE)
    # Check that using pred kernel instead of FFT gives the same result
    params <- setPredKernel(params, pred_kernel = getPredKernel(params))
    # Due to problem with fft on M1mac, skip this test on CRAN
    skip_on_cran()
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
    expect_identical(diet1[, , 1:14], diet2[, , 1:14])
    expect_identical(diet2[1, 1, 15], 111)
})


# getSSB ----
test_that("getSSB works", {
    ssb <- getSSB(sim)
    # expect_known_value(ssb, "values/getSSB")
    expect_snapshot(ssb)
    expect_equal(getSSB(sim)[1, ], getSSB(sim@params))
})

# getBiomass ----
test_that("getBiomass works", {
    biomass <- getBiomass(sim)
    # expect_known_value(biomass, "values/getBiomass")
    expect_snapshot(biomass)
    expect_equal(getBiomass(sim)[1, ], getBiomass(sim@params))
})

# getBiomass with biomass_cutoff ----
test_that("getBiomass works with biomass_cutoff", {
    # Add biomass_cutoff to species_params
    params_with_cutoff <- params
    params_with_cutoff@species_params$biomass_cutoff <- c(10, 20, 15, 5, 25, 8, 12, 18, 7, 9, 11, 14)

    # Create simulation with biomass_cutoff
    sim_with_cutoff <- project(params_with_cutoff, t_max = 10)

        # Test that biomass_cutoff is used when use_cutoff = TRUE
    biomass_with_cutoff <- getBiomass(sim_with_cutoff, use_cutoff = TRUE)

    # Test that use_cutoff = FALSE (default) ignores biomass_cutoff
    biomass_no_cutoff <- getBiomass(sim_with_cutoff, use_cutoff = FALSE)
    biomass_default <- getBiomass(sim_with_cutoff)
    expect_equal(biomass_no_cutoff, biomass_default)

            # Test that explicit size range arguments are ignored when use_cutoff = TRUE
    biomass_explicit <- getBiomass(sim_with_cutoff, use_cutoff = TRUE, min_w = 5, max_w = 1000)
    biomass_cutoff_used <- getBiomass(sim_with_cutoff, use_cutoff = TRUE)

    # These should be the same because explicit arguments are ignored when use_cutoff = TRUE
    expect_equal(biomass_explicit, biomass_cutoff_used)

    # Test that explicit size range arguments work when use_cutoff = FALSE
    biomass_explicit_no_cutoff <- getBiomass(sim_with_cutoff, use_cutoff = FALSE, min_w = 5, max_w = 1000)
    # This should be different from the biomass_cutoff result
    expect_false(all(biomass_explicit_no_cutoff == biomass_cutoff_used))

    # Test with some NA values in biomass_cutoff
    params_partial_cutoff <- params
    params_partial_cutoff@species_params$biomass_cutoff <- c(10, NA, 15, 5, NA, 8, 12, 18, 7, 9, 11, 14)
    sim_partial_cutoff <- project(params_partial_cutoff, t_max = 10)

    # Should work without error
    expect_no_error(getBiomass(sim_partial_cutoff))

    # Test with MizerParams object
    biomass_params <- getBiomass(params_with_cutoff)
    expect_equal(length(biomass_params), nrow(params_with_cutoff@species_params))

    # Test that use_cutoff = FALSE works with MizerParams
    biomass_params_no_cutoff <- getBiomass(params_with_cutoff, use_cutoff = FALSE)
    biomass_params_default <- getBiomass(params_with_cutoff)
    expect_equal(biomass_params_no_cutoff, biomass_params_default)
})

# getN ----
test_that("getN works", {
    N <- getN(sim)
    # expect_known_value(N, "values/getN")
    expect_snapshot(N)
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
