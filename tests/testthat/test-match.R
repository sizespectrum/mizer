test_that("matchBiomasses works", {
    params <- setBevertonHolt(NS_params)

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
    species_params(params)$biomass_observed[1:5] <- NA
    expect_unchanged(matchBiomasses(params), params)

    # If we double the observations, we get twice the abundance
    species <- 1:9
    species_params(params)$biomass_observed <-
        species_params(params)$biomass_observed * 2
    params2 <- matchBiomasses(params, species)
    expect_equal(params2@initial_n[6:9, ], params@initial_n[6:9, ] * 2)
    # but unobserved species don't change
    expect_equal(params2@initial_n[1:5, ], params@initial_n[1:5, ])
    # and unselected species don't change
    expect_equal(params2@initial_n[10:12, ], params@initial_n[10:12, ])

    # Throws an error if biomass_cutoff > w_max
    params@species_params$biomass_cutoff[6] <- 1e16
    expect_error(matchBiomasses(params),
                 "Whiting does not grow up to the biomass_cutoff")
})

test_that("matchNumbers works", {
    params <- setBevertonHolt(NS_params)

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
    species_params(params)$number_observed[1:5] <- NA
    expect_unchanged(matchNumbers(params), params)

    # If we double the observations, we get twice the abundance
    species <- 1:9
    species_params(params)$number_observed <-
        species_params(params)$number_observed * 2
    params2 <- matchNumbers(params, species)
    expect_equal(params2@initial_n[6:9, ], params@initial_n[6:9, ] * 2)
    # but unobserved species don't change
    expect_equal(params2@initial_n[1:5, ], params@initial_n[1:5, ])
    # and unselected species don't change
    expect_equal(params2@initial_n[10:12, ], params@initial_n[10:12, ])
})
