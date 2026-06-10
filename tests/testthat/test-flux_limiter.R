# Tests for the `flux_limiter` option of project() and the steady-state
# machinery. The limiter folds a frozen flux-limited (TVD) high-order advective
# flux into the implicit transport operator. See the numerical-details vignette.

test_that("flux_limiter = 'none' reproduces the upwind projection", {
    sim_default <- project(NS_params_small, t_max = 1, dt = 0.1,
                           progress_bar = FALSE)
    sim_none <- project(NS_params_small, t_max = 1, dt = 0.1,
                        flux_limiter = "none", progress_bar = FALSE)
    expect_identical(sim_none@n, sim_default@n)

    # The limiter must actually change the projection.
    sim_vl <- project(NS_params_small, t_max = 1, dt = 0.1,
                      flux_limiter = "van_leer", progress_bar = FALSE)
    expect_false(isTRUE(all.equal(sim_vl@n, sim_none@n)))
})

test_that("project validates the flux_limiter argument", {
    expect_error(project(NS_params_small, t_max = 1, flux_limiter = "bogus",
                         progress_bar = FALSE))
})

test_that("flux_limiter is stored and restored on append", {
    sim <- project(NS_params_small, t_max = 1, flux_limiter = "van_leer",
                   progress_bar = FALSE)
    expect_identical(sim@sim_params$flux_limiter, "van_leer")
    # Appending without specifying flux_limiter reuses the stored value.
    sim2 <- project(sim, t_max = 1, progress_bar = FALSE)
    expect_identical(sim2@sim_params$flux_limiter, "van_leer")
})

test_that("flux limiter keeps abundances non-negative", {
    for (m in c("euler", "predictor_corrector", "tr_bdf2")) {
        sim <- project(NS_params_small, t_max = 2, dt = 0.1, method = m,
                       flux_limiter = "van_leer", progress_bar = FALSE)
        expect_true(all(sim@n >= 0))
    }
})

test_that("flux limiter steady state is preserved by all methods", {
    # Build the limiter's own steady state and hold reproduction constant so the
    # only thing being tested is whether the projection scheme keeps that state
    # steady. This mirrors the predictor-corrector test in
    # test-getRequiredRDD.R, but with the flux limiter switched on everywhere.
    params <- single_sp_params
    species <- params@species_params$species[1]
    n <- params@species_params[species, "n"]
    ext_diffusion(params)[species, ] <- 0.1 * params@w^(n + 1)

    params <- steadySingleSpecies(params, species = species,
                                  flux_limiter = "van_leer")
    params@species_params$constant_reproduction <-
        getRequiredRDD(params, flux_limiter = "van_leer")
    params <- setRateFunction(params, "RDD", "constantRDD")

    idx <- params@w_min_idx[species]
    n0 <- initialN(params)[species, ]
    for (m in c("euler", "predictor_corrector", "tr_bdf2")) {
        sim <- project(params, t_max = 1, t_save = 1, method = m,
                       flux_limiter = "van_leer", progress_bar = FALSE)
        nf <- finalN(sim)[species, ]
        # Egg density and the whole spectrum stay put.
        expect_equal(nf[idx], n0[idx], tolerance = 1e-8)
        expect_lt(max(abs(nf - n0)) / max(n0), 1e-6)
    }
})

test_that("steadySingleSpecies flux_limiter default matches 'none'", {
    a <- steadySingleSpecies(single_sp_params)
    b <- steadySingleSpecies(single_sp_params, flux_limiter = "none")
    # Compare the abundance slots directly (initialN() carries a params
    # attribute whose time_modified timestamp differs between calls).
    expect_identical(a@initial_n, b@initial_n)
})

test_that("getRequiredRDD default matches flux_limiter = 'none'", {
    expect_identical(getRequiredRDD(single_sp_params),
                     getRequiredRDD(single_sp_params, flux_limiter = "none"))
})
