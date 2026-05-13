## Initialisation ----
species_params <- NS_species_params_gears
species_params$pred_kernel_type <- "truncated_lognormal"
params <- newMultispeciesParams(species_params, inter, min_w_pp = 1e-12,
                                n = 2/3, p = 0.7, lambda = 2.8 - 2/3,
                                initial_effort = 1, info_level = 0)
sim <- project(params, t_max = 2)

# getProportionOfLargeFish ----
test_that("getProportionOfLargeFish works", {
    sim <- project(params, effort = 1, t_max = 2, dt = 0.5, t_save = 0.5)
    # noddy test - using full range of sizes
    prop <- getProportionOfLargeFish(sim, threshold_w = 500)
    time_idx <- length(getTimes(sim))
    threshold_w <- sim@params@w > 500
    total_biomass <- sum(sweep(sim@n[time_idx, , ], 2,
                               sim@params@w * sim@params@dw, "*"))
    larger_biomass <- sum(sweep(sim@n[time_idx, , ], 2,
                                threshold_w * sim@params@w *
                                    sim@params@dw, "*"))
    expect_equal(prop[time_idx], larger_biomass / total_biomass,
                 ignore_attr = TRUE)
    # using a size range
    prop <- getProportionOfLargeFish(sim, min_w = 10, max_w = 5000,
                                     threshold_w = 500)
    range_w <- (sim@params@w >= 10) & (sim@params@w <= 5000)
    threshold_w <- sim@params@w > 500
    total_biomass <- sum(sweep(sim@n[time_idx, , ], 2,
                               range_w * sim@params@w * sim@params@dw, "*"))
    larger_biomass <- sum(sweep(sim@n[time_idx, , ], 2,
                                threshold_w * range_w * sim@params@w *
                                    sim@params@dw, "*"))
    expect_equal(prop[time_idx], larger_biomass / total_biomass,
                 ignore_attr = TRUE)
    expect_snapshot(prop)
})

test_that("getProportionOfLargeFish honours species, numbers, and threshold_l", {
    sim <- project(params, effort = 1, t_max = 2, dt = 0.5, t_save = 0.5)
    species <- c("Cod", "Haddock")
    sim@params@species_params$a <- 0.01
    sim@params@species_params$b <- 3
    threshold_l <- 10
    threshold_w <- 1

    by_length <- getProportionOfLargeFish(
        sim, species = species, threshold_w = threshold_w,
        threshold_l = threshold_l, biomass_proportion = FALSE
    )
    expected_large <- get_size_range_array(sim@params, max_l = threshold_l)[
        species, , drop = FALSE
    ]
    total_n <- apply(sweep(sim@n[, species, , drop = FALSE], 3,
                           sim@params@dw, "*"), 1, sum)
    upto_threshold_n <- apply(
        sweep(
            sweep(sim@n[, species, , drop = FALSE], c(2, 3), expected_large,
                  "*"),
            3, sim@params@dw, "*"
        ),
        1, sum
    )
    expect_equal(by_length, 1 - upto_threshold_n / total_n, ignore_attr = TRUE)
})

test_that("getProportionOfLargeFish works for MizerParams", {
    species <- c("Cod", "Haddock")
    expected_large <- get_size_range_array(params, min_w = 10, max_w = 500)[
        species, , drop = FALSE
    ]
    expected_total <- get_size_range_array(params, min_w = 10, max_w = 5000)[
        species, , drop = FALSE
    ]
    n <- params@initial_n[species, , drop = FALSE]
    total_biomass <- sum(n * expected_total * params@w * params@dw)
    upto_threshold_biomass <- sum(n * expected_large * params@w * params@dw)

    expect_equal(
        getProportionOfLargeFish(params, species = species, min_w = 10,
                                 max_w = 5000, threshold_w = 500),
        1 - upto_threshold_biomass / total_biomass
    )
})

# getMeanWeight ----
test_that("getMeanWeight works", {
    sim <- project(params, t_max = 2, dt = 0.5, t_save = 0.5)
    # all species, all size range
    total_biomass <- apply(sweep(sim@n, 3, sim@params@w * sim@params@dw, "*"),
                           1, sum)
    total_n <- apply(sweep(sim@n, 3, sim@params@dw, "*"), 1, sum)
    mw1 <- total_biomass / total_n
    mw <- getMeanWeight(sim)
    expect_equal(mw, mw1, ignore_attr = TRUE)
    # select species
    species <- sim@params@species_params$species[11:10]
    total_biomass <- apply(sweep(sim@n[, species, ], 3,
                                 sim@params@w * sim@params@dw, "*"), 1, sum)
    total_n <- apply(sweep(sim@n[, species, ], 3, sim@params@dw, "*"), 1, sum)
    mw2 <- total_biomass / total_n
    mw <- getMeanWeight(sim, species = species)
    expect_equal(mw, mw2, ignore_attr = TRUE)
    # select size range
    min_w <- 10
    max_w <- 10000
    size_n <- get_size_range_array(sim@params, min_w = min_w, max_w = max_w)
    total_biomass <- apply(
        sweep(sweep(sim@n, c(2, 3), size_n, "*"), 3,
              sim@params@w * sim@params@dw, "*"),
        1, sum
    )
    total_n <- apply(
        sweep(sweep(sim@n, c(2, 3), size_n, "*"), 3, sim@params@dw, "*"),
        1, sum
    )
    mw3 <- total_biomass / total_n
    mw <- getMeanWeight(sim, min_w = min_w, max_w = max_w)
    expect_equal(mw, mw3, ignore_attr = TRUE)
    # select size range and species
    total_biomass <- apply(
        sweep(sweep(sim@n, c(2, 3), size_n, "*")[, species, ], 3,
              sim@params@w * sim@params@dw, "*"),
        1, sum
    )
    total_n <- apply(
        sweep(sweep(sim@n, c(2, 3), size_n, "*")[, species, ], 3,
              sim@params@dw, "*"),
        1, sum
    )
    mw4 <- total_biomass / total_n
    mw <- getMeanWeight(sim, species = species, min_w = min_w, max_w = max_w)
    expect_equal(mw, mw4, ignore_attr = TRUE)
    expect_snapshot(mw)
})

test_that("getMeanWeight works for MizerParams", {
    species <- c("Cod", "Haddock")
    expected_n <- sum(getN(params, min_w = 10, max_w = 5000)[species])
    expected_biomass <- sum(getBiomass(params, min_w = 10, max_w = 5000)[species])

    expect_equal(
        getMeanWeight(params, species = species, min_w = 10, max_w = 5000),
        expected_biomass / expected_n
    )
})

# getMeanMaxWeight ----
test_that("getMeanMaxWeight works", {
    expect_error(getMeanMaxWeight(sim, measure = NA),
                 "measure must be one of")
    species <- c("Cod", "Haddock")
    n_species <- getN(sim)
    biomass_species <- getBiomass(sim)
    w_max <- sim@params@species_params$w_max
    mmw_numbers <- apply(sweep(n_species[, species, drop = FALSE], 2,
                               w_max[
                                   match(species,
                                         sim@params@species_params$species)
                               ], "*"), 1, sum) /
        apply(n_species[, species, drop = FALSE], 1, sum)
    mmw_biomass <- apply(sweep(biomass_species[, species, drop = FALSE], 2,
                               w_max[
                                   match(species,
                                         sim@params@species_params$species)
                               ], "*"), 1, sum) /
        apply(biomass_species[, species, drop = FALSE], 1, sum)
    expect_equal(getMeanMaxWeight(sim, species = species, measure = "numbers"),
                 mmw_numbers, ignore_attr = TRUE)
    expect_equal(getMeanMaxWeight(sim, species = species, measure = "biomass"),
                 mmw_biomass, ignore_attr = TRUE)
    expect_equal(getMeanMaxWeight(sim, species = species, measure = "both"),
                 cbind(mmw_numbers, mmw_biomass), ignore_attr = TRUE)
    expect_snapshot(getMeanMaxWeight(sim, measure = "both"))
})

test_that("getMeanMaxWeight works for MizerParams", {
    species <- c("Cod", "Haddock")
    n_species <- getN(params, min_w = 10, max_w = 5000)[species]
    biomass_species <- getBiomass(params, min_w = 10, max_w = 5000)[species]
    w_max <- params@species_params$w_max[
        match(species, params@species_params$species)
    ]
    expected_numbers <- sum(n_species * w_max) / sum(n_species)
    expected_biomass <- sum(biomass_species * w_max) / sum(biomass_species)

    expect_equal(
        getMeanMaxWeight(params, species = species, measure = "numbers",
                         min_w = 10, max_w = 5000),
        expected_numbers
    )
    expect_equal(
        getMeanMaxWeight(params, species = species, measure = "biomass",
                         min_w = 10, max_w = 5000),
        expected_biomass
    )
    expect_equal(
        getMeanMaxWeight(params, species = species, measure = "both",
                         min_w = 10, max_w = 5000),
        c(mmw_numbers = expected_numbers, mmw_biomass = expected_biomass)
    )
})

# getCommunitySlope ----
test_that("getCommunitySlope works", {
    slope_b <- getCommunitySlope(sim)
    # dims
    expect_equal(dim(slope_b), c(dim(sim@n)[1], 3), ignore_attr = TRUE)
    # sum biomasses
    biomass <- apply(sweep(sim@n, 3, sim@params@w, "*"), c(1, 3), sum)
    # r2, slope and intercept at last time step
    lm_res <- lm(log(biomass[dim(sim@n)[1], ]) ~ log(sim@params@w))
    expect_equal(slope_b[dim(sim@n)[1], "r2"], summary(lm_res)$r.squared,
                 ignore_attr = TRUE)
    expect_equal(slope_b[dim(sim@n)[1], "slope"],
                 summary(lm_res)$coefficients[2, 1], ignore_attr = TRUE)
    expect_equal(slope_b[dim(sim@n)[1], "intercept"],
                 summary(lm_res)$coefficients[1, 1], ignore_attr = TRUE)
    # Test just numbers not biomass
    slope_n <- getCommunitySlope(sim, biomass = FALSE)
    expect_equal(dim(slope_n), c(dim(sim@n)[1], 3), ignore_attr = TRUE)
    # sum numbers
    numbers <- apply(sim@n, c(1, 3), sum)
    # r2, slope and intercept at last time step
    lm_res <- lm(log(numbers[dim(sim@n)[1], ]) ~ log(sim@params@w))
    expect_equal(slope_n[dim(sim@n)[1], "r2"], summary(lm_res)$r.squared,
                 ignore_attr = TRUE)
    expect_equal(slope_n[dim(sim@n)[1], "slope"],
                 summary(lm_res)$coefficients[2, 1], ignore_attr = TRUE)
    expect_equal(slope_n[dim(sim@n)[1], "intercept"],
                 summary(lm_res)$coefficients[1, 1], ignore_attr = TRUE)
    # Check the sizes
    slope_b2 <- getCommunitySlope(sim, min_w = 10, max_w = 10000)
    sizes <- (sim@params@w >= 10) & (sim@params@w <= 10000)
    biomass <- apply(sweep(sim@n, 3, sim@params@w, "*"), c(1, 3), sum)
    # r2, slope and intercept at last time step
    lm_res <- lm(log(biomass[dim(sim@n)[1], sizes]) ~
                     log(sim@params@w[sizes]))
    expect_equal(slope_b2[dim(sim@n)[1], "r2"], summary(lm_res)$r.squared,
                 ignore_attr = TRUE)
    expect_equal(slope_b2[dim(sim@n)[1], "slope"],
                 summary(lm_res)$coefficients[2, 1], ignore_attr = TRUE)
    expect_equal(slope_b2[dim(sim@n)[1], "intercept"],
                 summary(lm_res)$coefficients[1, 1], ignore_attr = TRUE)
    # Check the species
    dem_species <- sim@params@species_params$species[5:12]
    slope_b3 <- getCommunitySlope(sim, species = dem_species)
    biomass <- apply(sweep(sim@n[, dem_species, ], 3, sim@params@w, "*"),
                     c(1, 3), sum)
    # r2, slope and intercept at last time step
    lm_res <- lm(log(biomass[dim(sim@n)[1], ]) ~ log(sim@params@w))
    expect_equal(slope_b3[dim(sim@n)[1], "r2"], summary(lm_res)$r.squared,
                 ignore_attr = TRUE)
    expect_equal(slope_b3[dim(sim@n)[1], "slope"],
                 summary(lm_res)$coefficients[2, 1], ignore_attr = TRUE)
    expect_equal(slope_b3[dim(sim@n)[1], "intercept"],
                 summary(lm_res)$coefficients[1, 1], ignore_attr = TRUE)
    expect_snapshot(slope_b3)
})

test_that("getCommunitySlope works for MizerParams", {
    expected_total_n <- colSums(
        params@initial_n[c("Cod", "Haddock"), , drop = FALSE] *
            get_size_range_array(params, min_w = 10, max_w = 5000)[
                c("Cod", "Haddock"), , drop = FALSE
            ]
    )
    expected_total_n <- expected_total_n * params@w
    expected_total_n[expected_total_n <= 0] <- NA
    expected_fit <- summary(lm(log(expected_total_n) ~ log(params@w)))

    expect_equal(
        getCommunitySlope(params, species = c("Cod", "Haddock"),
                          min_w = 10, max_w = 5000),
        data.frame(
            slope = expected_fit$coefficients[2, 1],
            intercept = expected_fit$coefficients[1, 1],
            r2 = expected_fit$r.squared
        )
    )
})
