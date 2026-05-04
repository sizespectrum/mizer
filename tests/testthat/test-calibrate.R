test_that("calibrateBiomass works", {
    params <- example_params()
    # Does nothing when no observed biomass
    expect_identical(calibrateBiomass(params), params)
    species_params(params)$biomass_observed <- NA
    expect_identical(calibrateBiomass(params), params)
    # Does nothing if observed already equals model
    species_params(params)$biomass_cutoff <- 1e-4
    species_params(params)$biomass_observed <-
        rowSums(sweep(params@initial_n, 2, params@w * params@dw, "*"))
    expect_unchanged(calibrateBiomass(params), params)
    # Even if only partially observed
    species_params(params)$biomass_observed[1:2] <- NA
    expect_unchanged(calibrateBiomass(params), params)
    # If we double the observations, we get twice the abundance
    species_params(params)$biomass_observed <-
        species_params(params)$biomass_observed * 2
    params2 <- calibrateBiomass(params)
    expect_equal(params2@initial_n, params@initial_n * 2)
    # We don't need to check other slots because this function uses
    # `scaleModel()` which is unit-tested separately.
})


test_that("calibrateNumber works", {
    params <- example_params()
    # Does nothing when no observed Number
    expect_identical(calibrateNumber(params), params)
    species_params(params)$number_observed <- NA
    expect_identical(calibrateNumber(params), params)
    # Does nothing if observed already equals model
    species_params(params)$number_cutoff <- 1e-4
    species_params(params)$number_observed <-
        rowSums(sweep(params@initial_n, 2, params@dw, "*"))
    expect_unchanged(calibrateNumber(params), params)
    # Even if only partially observed
    species_params(params)$number_observed[1:2] <- NA
    expect_unchanged(calibrateNumber(params), params)
    # If we double the observations, we get twice the abundance
    species_params(params)$number_observed <-
        species_params(params)$number_observed * 2
    params2 <- calibrateNumber(params)
    expect_equal(params2@initial_n, params@initial_n * 2)
    # We don't need to check other slots because this function uses
    # `scaleModel()` which is unit-tested separately.
})

test_that("calibrateBiomass and calibrateNumber honour cutoffs", {
    params <- NS_params
    cutoff <- params@w[10]

    species_params(params)$biomass_cutoff <- cutoff
    species_params(params)$biomass_observed <-
        c(rep(NA, 3),
          rowSums((params@initial_n * params@w * params@dw)[, params@w >= cutoff])[4:12] * 2)
    observed_total <- sum(species_params(params)$biomass_observed, na.rm = TRUE)
    model_total <- sum(vapply(which(!is.na(species_params(params)$biomass_observed)),
                              function(i) {
                                  sum((params@initial_n[i, ] * params@w * params@dw)[params@w >= cutoff])
                              }, numeric(1)))
    biomass_scaled <- calibrateBiomass(params)
    expect_equal(biomass_scaled@initial_n,
                 params@initial_n * observed_total / model_total)

    params <- NS_params
    species_params(params)$number_cutoff <- cutoff
    species_params(params)$number_observed <-
        c(rep(NA, 2),
          rowSums((params@initial_n * params@dw)[, params@w >= cutoff])[3:12] * 1.5)
    observed_total <- sum(species_params(params)$number_observed, na.rm = TRUE)
    model_total <- sum(vapply(which(!is.na(species_params(params)$number_observed)),
                              function(i) {
                                  sum((params@initial_n[i, ] * params@dw)[params@w >= cutoff])
                              }, numeric(1)))
    number_scaled <- calibrateNumber(params)
    expect_equal(number_scaled@initial_n,
                 params@initial_n * observed_total / model_total)
})

test_that("calibrateYield works and warns about deprecation", {
    params <- NS_params
    expect_warning(expect_identical(calibrateYield(params), params),
                   "deprecated")

    species_params(params)$yield_observed <- NA
    expect_warning(expect_identical(calibrateYield(params), params),
                   "deprecated")

    species_params(params)$yield_observed <-
        rowSums(sweep(params@initial_n, 2, params@w * params@dw, "*") *
                    getFMort(params))
    expect_warning(expect_unchanged(calibrateYield(params), params),
                   "deprecated")

    species_params(params)$yield_observed <- species_params(params)$yield_observed * 2
    params2 <- suppressWarnings(calibrateYield(params))
    expect_equal(params2@initial_n, params@initial_n * 2)
})

test_that("scaleModel does not change dynamics.", {
    factor <- 10
    sim <- project(NS_params, t_max = 1)
    params2 <- scaleModel(NS_params, factor)
    sim2 <- project(params2, t_max = 1)
    expect_equal(sim2@n[1, , ], sim@n[1, , ] * factor)
    expect_equal(sim2@n[2, , ], sim@n[2, , ] * factor)
})

test_that("scaleModel rescales documented state and parameter slots", {
    factor <- 4
    params <- NS_params
    params@species_params$R_max <- 1:nrow(params@species_params)

    scaled <- scaleModel(params, factor)

    expect_equal(scaled@cc_pp, params@cc_pp * factor)
    expect_equal(scaled@search_vol, params@search_vol / factor)
    expect_equal(scaled@species_params$gamma, params@species_params$gamma / factor)
    expect_equal(scaled@species_params$R_max, params@species_params$R_max * factor)
    expect_equal(scaled@initial_n, params@initial_n * factor)
    expect_equal(scaled@initial_n_pp, params@initial_n_pp * factor)
    expect_equal(scaled@resource_params$kappa, params@resource_params$kappa * factor)
    expect_equal(scaled@sc, params@sc * factor)
})

test_that("scaleModel renames deprecated r_max column", {
    params <- NS_params
    params@species_params$r_max <- 1:nrow(params@species_params)
    params@species_params$R_max <- NULL

    scaled <- scaleModel(params, 2)
    expect_false("r_max" %in% names(scaled@species_params))
    expect_equal(scaled@species_params$R_max, 2 * (1:nrow(params@species_params)))
})

test_that("scaleModel updates `time_modified`", {
    params <- scaleModel(NS_params, factor = 2)
    expect_false(identical(params@time_modified, NS_params@time_modified))
})
