test_that("matchBiomasses works", {
    params <- setBevertonHolt(NS_params_small)

    # Does nothing when no observed biomass
    expect_identical(matchBiomasses(params), params)
    species_params(params)$biomass_observed <- NA
    expect_identical(matchBiomasses(params), params)

    # Does nothing if observed already equals model
    species_params(params)$biomass_cutoff <- 1e-4
    biomass_actual <- getBiomass(params, use_cutoff = TRUE)
    species_params(params)$biomass_observed <- biomass_actual
    expect_unchanged(matchBiomasses(params), params)

    # Even if only partially observed
    species_params(params)$biomass_observed[1] <- NA
    expect_unchanged(matchBiomasses(params), params)

    # If we double the observations, we get twice the abundance
    species <- 1:2
    species_params(params)$biomass_observed <-
        species_params(params)$biomass_observed * 2
    params2 <- matchBiomasses(params, species)
    expect_equal(params2@initial_n[2, ], params@initial_n[2, ] * 2)
    # but unobserved species don't change
    expect_equal(params2@initial_n[1, ], params@initial_n[1, ])
    # and unselected species don't change
    expect_equal(params2@initial_n[3, ], params@initial_n[3, ])
})

test_that("matchNumbers works", {
    params <- setBevertonHolt(NS_params_small)

    # Does nothing when no observed numbers
    expect_identical(matchNumbers(params), params)
    species_params(params)$number_observed <- NA
    expect_unchanged(matchNumbers(params), params)
    # Does nothing if observed already equals model
    species_params(params)$number_cutoff <- 1e-4
    number_actual <-
        rowSums(sweep(params@initial_n, 2, params@dw, "*"))
    species_params(params)$number_observed <- number_actual
    expect_unchanged(matchNumbers(params), params)
    # Even if only partially observed
    species_params(params)$number_observed[1] <- NA
    expect_unchanged(matchNumbers(params), params)

    # If we double the observations, we get twice the abundance
    species <- 1:2
    species_params(params)$number_observed <-
        species_params(params)$number_observed * 2
    params2 <- matchNumbers(params, species)
    expect_equal(params2@initial_n[2, ], params@initial_n[2, ] * 2)
    # but unobserved species don't change
    expect_equal(params2@initial_n[1, ], params@initial_n[1, ])
    # and unselected species don't change
    expect_equal(params2@initial_n[3, ], params@initial_n[3, ])
})

test_that("matchBiomasses and matchNumbers fail for unreachable cutoff", {
    params <- setBevertonHolt(NS_params_small)
    species_params(params)$biomass_observed <- 0
    species_params(params)$biomass_observed[1] <- 1
    species_params(params)$biomass_cutoff <- 0
    species_params(params)$biomass_cutoff[1] <- 1e20
    expect_error(matchBiomasses(params, 1),
                 "Sprat does not grow up to the biomass_cutoff size of 1e\\+20 grams")

    species_params(params)$number_observed <- 0
    species_params(params)$number_observed[1] <- 1
    species_params(params)$number_cutoff <- 0
    species_params(params)$number_cutoff[1] <- 1e20
    expect_error(matchNumbers(params, 1),
                 "Sprat does not grow up to the number_cutoff size of 1e\\+20 grams")
})
