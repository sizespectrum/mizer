# get_initial_n ----
edition1_initial_n_params <- withr::with_options(
    list(mizer_defaults_edition = 1),
    newMultispeciesParams(NS_species_params_gears_small, inter_small, info_level = 0)
)

test_that("get_initial_n is working properly in edition 1", {
    old <- getOption("mizer_defaults_edition")
    on.exit(options(mizer_defaults_edition = old), add = TRUE)
    options(mizer_defaults_edition = 1)
    params <- edition1_initial_n_params
    n <- get_initial_n(params)
    no_sp <- nrow(params@species_params)
    for (i in 1:no_sp) {
        expect_true(all(n[i, params@w > params@species_params$w_max[i]] == 0))
        expect_true(all(n[i, params@w < params@species_params$w_min[i]] == 0))
    }
    # Check slope of all species is the same
    slopes <- rep(NA, no_sp)
    for (i in 1:no_sp) {
        n_idx <- which(n[i, ] != 0)
        slopes[i] <- (log(n[i, min(n_idx)]) - log(n[i, max(n_idx)])) / 
            (log(params@w[min(n_idx)]) - log(params@w[max(n_idx)]))
    }
    expect_equal(slopes, rep(slopes[1], no_sp))
    # Check that slopes = slope0
})

test_that("get_initial_n validates params and honours n0_mult in edition 1", {
    expect_error(get_initial_n(1), "params argument must of type MizerParams")

    old <- getOption("mizer_defaults_edition")
    on.exit(options(mizer_defaults_edition = old), add = TRUE)
    options(mizer_defaults_edition = 1)

    params <- edition1_initial_n_params
    n1 <- get_initial_n(params, n0_mult = 1)
    n2 <- get_initial_n(params, n0_mult = 2)
    expect_equal(n2, 2 * n1, ignore_attr = TRUE)
})

test_that("get_initial_n keeps the egg size class populated in edition 1", {
    # Issue #610: zeroing every size class below `w_min` also emptied the class
    # that contains `w_min`, so the species had no abundance at `w_min_idx`.
    old <- getOption("mizer_defaults_edition")
    on.exit(options(mizer_defaults_edition = old), add = TRUE)
    options(mizer_defaults_edition = 1)

    sp <- NS_species_params_small
    # The egg sizes of the last two species do not lie on grid points
    sp$w_min <- c(1e-3, 1e-2, 1e-1)
    params <- newMultispeciesParams(sp, no_w = 20, info_level = 0)
    expect_true(all(params@w[params@w_min_idx[2:3]] < sp$w_min[2:3]))

    n <- get_initial_n(params)
    for (i in seq_len(nrow(params@species_params))) {
        idx <- params@w_min_idx[[i]]
        # The size class holding the egg size is populated ...
        expect_gt(n[i, idx], 0)
        # ... and the classes below it are empty
        expect_true(all(n[i, seq_len(idx - 1)] == 0))
        expect_true(all(n[i, params@w > params@species_params$w_max[i]] == 0))
    }
})
