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

test_that("matchYields works", {
    params <- setBevertonHolt(NS_params_small)
    initial_effort(params) <- setNames(rep(1, length(initial_effort(params))),
                                       names(initial_effort(params)))
    yield_actual <- rowSums(sweep(params@initial_n, 2, params@w * params@dw, "*") *
                                getFMort(params))

    species_params(params)$yield_observed <- NA_real_
    expect_message(
        params_na <- suppressWarnings(matchYields(params, info_level = 3)),
        "The following species have no yield observations"
    )
    expect_equal(params_na@initial_n, params@initial_n)

    species_params(params)$yield_observed <- 0
    species_params(params)$yield_observed[3] <- 2 * yield_actual[3]
    expect_warning(params2 <- matchYields(params, 3, info_level = 0), "deprecated")
    expect_equal(params2@initial_n[3, ], 2 * params@initial_n[3, ])
    expect_equal(params2@initial_n[-3, ], params@initial_n[-3, ])
})

test_that("matchYields updates `time_modified`", {
    p <- NS_params_small
    species_params(p)$yield_observed <- 1
    p2 <- suppressWarnings(matchYields(p))
    expect_false(identical(p2@time_modified, p@time_modified))
})

# The steady-state notice ----

test_that("the match functions say that they moved the model off steady state", {
    p <- NS_params_small
    species_params(p)$biomass_observed <- as.numeric(getBiomass(p)) * 2
    expect_message(matchBiomasses(p), "moved it off its steady state")

    species_params(p)$number_observed <- as.numeric(getN(p)) * 2
    expect_message(matchNumbers(p), "moved it off its steady state")
})

test_that("the steady-state notice is silenced by info_level", {
    p <- NS_params_small
    species_params(p)$biomass_observed <- as.numeric(getBiomass(p)) * 2
    expect_no_message(matchBiomasses(p, info_level = 0))
    withr::local_options(mizer_info_level = 0)
    expect_no_message(matchBiomasses(p))
})

test_that("the calibrate functions do not claim to break the steady state", {
    # They apply one overall scaling factor, which is an exact symmetry of the
    # model, so the steady state survives them untouched. If this ever stops
    # being true they need the notice that the match functions carry.
    p <- suppressMessages(
        steady(NS_params_small, tol = 1e-6, t_max = 500,
               progress_bar = FALSE, info_level = 0))
    before <- steady_biomass_drift(p)
    species_params(p)$biomass_observed <- as.numeric(getBiomass(p)) * 2
    expect_equal(steady_biomass_drift(suppressMessages(calibrateBiomass(p))),
                 before, tolerance = 1e-6)
    expect_equal(steady_biomass_drift(suppressMessages(scaleModel(p, 2))),
                 before, tolerance = 1e-6)
})
