## Initialisation ----
species_params <- NS_species_params_gears
species_params$pred_kernel_type <- "truncated_lognormal"
params <- newMultispeciesParams(species_params, inter, min_w_pp = 1e-12,
                                n = 2/3, p = 0.7, lambda = 2.8 - 2/3,
                                initial_effort = 1, info_level = 0)
sim <- project(params, t_max = 2)
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

    expect_identical(names(dimnames(get_size_range_array(params))), c("sp", "w"))
})


# getYieldGear ----
test_that("getYieldGear works",{
    y <- getYieldGear(sim)
    # check dims
    expect_equal(dim(y), c(length(getTimes(sim)), dim(params@catchability)[1],
                           dim(params@catchability)[2]), ignore_attr = TRUE)
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

test_that("getYieldGear for params matches fishing mortality by gear times biomass", {
    biomass <- sweep(params@initial_n, 2, params@w * params@dw, "*")
    f_gear <- getFMortGear(params)
    expected <- apply(sweep(f_gear, c(2, 3), biomass, "*"), c(1, 2), sum)
    expect_equal(getYieldGear(params), expected)
})


# getYield ----
test_that("getYield works",{
    y <- getYield(sim)
    expect_true(is.ArrayTimeBySpecies(y))
    # check dims
    expect_equal(dim(y), c(length(getTimes(sim)), dim(params@catchability)[2]),
                 ignore_attr = TRUE)
    # check a value and assume the others are right
    biomass <- sweep(sim@n,3,sim@params@w * sim@params@dw, "*")
    f <- getFMort(sim)
    expect_equal(sum((f*biomass)[1,1,]), y[1,1], ignore_attr = TRUE)
    # numeric test
    # expect_known_value(y, "values/getYield")
    expect_snapshot(y)
    expect_equal(getYield(sim)[1, ], getYield(sim@params))
})

test_that("getYield for params matches fishing mortality times biomass", {
    biomass <- sweep(params@initial_n, 2, params@w * params@dw, "*")
    f <- getFMort(params, drop = FALSE)
    expected <- apply(f * biomass, 1, sum)
    expect_equal(getYield(params), expected)
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


# getTrophicLevel ----
test_that("getTrophicLevel returns matrix with correct structure", {
    tl <- getTrophicLevel(params, n, n_pp)
    expect_true(is.ArraySpeciesBySize(tl))
    expect_true(is.matrix(tl))
    expect_equal(dim(tl), c(no_sp, no_w))
    expect_equal(dimnames(tl), dimnames(params@initial_n))
    # All trophic levels >= 1 (since T = 1 + non-negative)
    expect_true(all(tl >= 1, na.rm = TRUE))
    # Trophic levels should be reasonable (< 10)
    expect_true(all(tl < 10, na.rm = TRUE))
})

test_that("getTrophicLevel gives same result with explicit pred_kernel", {
    tl1 <- getTrophicLevel(params, n, n_pp)
    # Force explicit pred_kernel storage
    params2 <- setPredKernel(params, pred_kernel = getPredKernel(params))
    tl2 <- getTrophicLevel(params2, n, n_pp)
    expect_equal(tl1, tl2, tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("getTrophicLevel increases along body size for apex predators", {
    tl <- getTrophicLevel(NS_params)
    # For Cod (apex predator), trophic level should increase with size
    cod_tl <- tl["Cod", ]
    cod_tl <- cod_tl[!is.na(cod_tl)]
    # Should be non-decreasing overall (allow small numerical fluctuations)
    expect_true(cod_tl[length(cod_tl)] >= cod_tl[1])
})

# getTrophicLevelBySpecies ----
test_that("getTrophicLevelBySpecies returns named vector", {
    tl_sp <- getTrophicLevelBySpecies(params, n, n_pp)
    expect_true(is.numeric(tl_sp))
    expect_equal(length(tl_sp), no_sp)
    expect_equal(names(tl_sp), params@species_params$species)
    expect_true(all(tl_sp >= 1, na.rm = TRUE))
})

test_that("getTrophicLevelBySpecies is consistent with getTrophicLevel", {
    tl <- getTrophicLevel(params, n, n_pp)
    tl_sp <- getTrophicLevelBySpecies(params, n, n_pp)
    # Species-level trophic level should be between min and max size-resolved tl
    for (i in seq_len(no_sp)) {
        tl_range <- range(tl[i, ], na.rm = TRUE)
        expect_gte(tl_sp[i], tl_range[1] - 1e-10)
        expect_lte(tl_sp[i], tl_range[2] + 1e-10)
    }
})


# getSSB ----
test_that("getSSB works", {
    ssb <- getSSB(sim)
    expect_true(is.ArrayTimeBySpecies(ssb))
    # expect_known_value(ssb, "values/getSSB")
    expect_snapshot(ssb)
    expect_equal(getSSB(sim)[1, ], getSSB(sim@params))
})

test_that("getSSB matches mature biomass formula for params and sim", {
    expected_params <- ((params@initial_n * params@maturity) %*%
                            (params@w * params@dw))[, , drop = TRUE]
    expected_sim <- apply(
        sweep(sweep(sim@n, c(2, 3), sim@params@maturity, "*"), 3,
              sim@params@w * sim@params@dw, "*"),
        c(1, 2), sum
    )

    expect_equal(getSSB(params), expected_params)
    expect_true(is.ArrayTimeBySpecies(getSSB(sim)))
    expect_equal(getSSB(sim), expected_sim, ignore_attr = TRUE)
})

# getBiomass ----
test_that("getBiomass works", {
    biomass <- getBiomass(sim)
    expect_true(is.ArrayTimeBySpecies(biomass))
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
    sim_with_cutoff <- project(params_with_cutoff, t_max = 2)

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
    sim_partial_cutoff <- project(params_partial_cutoff, t_max = 2)

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
    expect_true(is.ArrayTimeBySpecies(N))
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

test_that("getGrowthCurves percentage rescales by maximum weight", {
    curves <- getGrowthCurves(params, species = c("Cod", "Haddock"),
                              percentage = TRUE)
    raw <- getGrowthCurves(params, species = c("Cod", "Haddock"),
                           percentage = FALSE)
    w_max <- params@species_params$w_max[match(rownames(curves),
                                               params@species_params$species)]
    expected <- sweep(raw, 1, w_max, "/") * 100

    expect_equal(curves, expected)
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
