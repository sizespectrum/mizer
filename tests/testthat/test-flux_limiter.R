# Tests for the flux-limited advective flux, controlled by the `flux_limiter`
# entry of the `second_order_w` slot. The limiter folds a frozen flux-limited
# (TVD) high-order advective flux into the implicit transport operator. See the
# numerical-details vignette.

test_that("flux_limiter_scheme reads the second_order_w slot", {
    expect_identical(flux_limiter_scheme(NS_params_small), "none")
    p <- NS_params_small
    second_order_w(p) <- c(flux_limiter = TRUE)
    expect_identical(flux_limiter_scheme(p), "van_leer")
})

test_that("the flux limiter changes the projection when switched on", {
    sim_none <- project(NS_params_small, t_max = 1, dt = 0.1,
                        progress_bar = FALSE)
    p_vl <- NS_params_small
    second_order_w(p_vl) <- c(flux_limiter = TRUE)
    sim_vl <- project(p_vl, t_max = 1, dt = 0.1, progress_bar = FALSE)
    # With the limiter off the result is the plain upwind projection; switching
    # it on must change the spectrum.
    expect_false(isTRUE(all.equal(sim_vl@n, sim_none@n)))
})

test_that("the flux limiter scheme is recorded and carried over on append", {
    p_vl <- NS_params_small
    second_order_w(p_vl) <- c(flux_limiter = TRUE)
    sim <- project(p_vl, t_max = 1, progress_bar = FALSE)
    expect_identical(sim@sim_params$flux_limiter, "van_leer")
    # The scheme is a property of the params, so appending carries it over.
    sim2 <- project(sim, t_max = 1, progress_bar = FALSE)
    expect_identical(sim2@sim_params$flux_limiter, "van_leer")
})

test_that("flux limiter keeps abundances non-negative", {
    p_vl <- NS_params_small
    second_order_w(p_vl) <- c(flux_limiter = TRUE)
    for (m in c("euler", "predictor_corrector", "tr_bdf2")) {
        sim <- project(p_vl, t_max = 2, dt = 0.1, method = m,
                       progress_bar = FALSE)
        expect_true(all(sim@n >= 0))
    }
})

test_that("flux limiter steady state is preserved by all methods", {
    # Build the limiter's own steady state and hold reproduction constant so the
    # only thing being tested is whether the projection scheme keeps that state
    # steady. This mirrors the predictor-corrector test in
    # test-getRequiredRDD.R, but with the flux limiter switched on everywhere.
    params <- single_sp_params
    second_order_w(params) <- c(flux_limiter = TRUE)
    species <- params@species_params$species[1]
    n <- params@species_params[species, "n"]
    ext_diffusion(params)[species, ] <- 0.1 * params@w^(n + 1)

    params <- steadySingleSpecies(params, species = species)
    params@species_params$constant_reproduction <- getRequiredRDD(params)
    params <- setRateFunction(params, "RDD", "constantRDD")

    idx <- params@w_min_idx[species]
    n0 <- initialN(params)[species, ]
    for (m in c("euler", "predictor_corrector", "tr_bdf2")) {
        sim <- project(params, t_max = 1, t_save = 1, method = m,
                       progress_bar = FALSE)
        nf <- finalN(sim)[species, ]
        # Egg density and the whole spectrum stay put.
        expect_equal(nf[idx], n0[idx], tolerance = 1e-8)
        expect_lt(max(abs(nf - n0)) / max(n0), 1e-6)
    }
})

test_that("steadySingleSpecies honours the slot", {
    # Switching the limiter on via the slot changes the interior steady state.
    p_vl <- single_sp_params
    second_order_w(p_vl) <- c(flux_limiter = TRUE)

    a <- steadySingleSpecies(single_sp_params)
    b <- steadySingleSpecies(p_vl)
    expect_false(isTRUE(all.equal(a@initial_n, b@initial_n)))
    # getRequiredRDD reads the same slot, but its value is scheme-independent
    # here: the van Leer limiter reverts to pure upwind at the non-smooth
    # recruitment boundary, so the egg-cell coefficients (and hence the required
    # reproduction) are identical with the limiter on or off. The consistency of
    # getRequiredRDD with the limited projection is covered by the steady-state
    # preservation test above.
    expect_identical(getRequiredRDD(single_sp_params), getRequiredRDD(p_vl))
})
