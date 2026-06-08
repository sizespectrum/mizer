test_that("scaleRates with factor = 1 leaves params unchanged", {
    # validParams is called internally, so use it as the base to avoid
    # spurious differences from columns it adds (e.g. D_ext)
    params <- suppressMessages(validParams(NS_params_small))
    expect_unchanged(scaleRates(params, 1), params)
})

test_that("scaleRates validates factor", {
    params <- NS_params_small
    expect_error(scaleRates(params, 0))
    expect_error(scaleRates(params, -1))
    expect_error(scaleRates(params, "a"))
})

test_that("scaleRates scales all rate slots by factor", {
    params <- NS_params_small
    f <- 3
    scaled <- scaleRates(params, f)

    expect_equal(scaled@search_vol,    params@search_vol    * f)
    expect_equal(scaled@intake_max,    params@intake_max    * f)
    expect_equal(scaled@metab,         params@metab         * f)
    expect_equal(scaled@mu_b,          params@mu_b          * f)
    expect_equal(scaled@ext_encounter, params@ext_encounter * f)
    expect_equal(scaled@ext_diffusion, params@ext_diffusion * f)
    expect_equal(scaled@catchability,  params@catchability  * f)
    expect_equal(scaled@rr_pp,         params@rr_pp         * f)
})

test_that("scaleRates scales species_params columns by factor", {
    params <- NS_params_small
    f <- 3
    scaled <- scaleRates(params, f)
    sp  <- params@species_params
    sp2 <- scaled@species_params

    for (col in c("gamma", "h", "ks", "k", "z0", "E_ext")) {
        if (col %in% names(sp)) {
            expect_equal(sp2[[col]], sp[[col]] * f,
                         label = paste0("species_params$", col))
        }
    }
    # R_max may be Inf by default; Inf * f == Inf
    if ("R_max" %in% names(sp)) {
        expect_equal(sp2$R_max, sp$R_max * f)
    }
})

test_that("scaleRates scales given_species_params columns by factor", {
    params <- NS_params_small
    f <- 3
    scaled <- scaleRates(params, f)
    gsp  <- params@given_species_params
    gsp2 <- scaled@given_species_params

    for (col in c("gamma", "h", "ks", "k", "z0", "R_max")) {
        if (col %in% names(gsp)) {
            expect_equal(gsp2[[col]], gsp[[col]] * f,
                         label = paste0("given_species_params$", col))
        }
    }
})

test_that("scaleRates scales gear_params catchability by factor", {
    params <- NS_params_small
    f <- 3
    scaled <- scaleRates(params, f)
    expect_equal(scaled@gear_params$catchability,
                 params@gear_params$catchability * f)
})

test_that("scaleRates does not add absent species_params columns", {
    params <- NS_params_small
    # z0pre is absent from NS_params_small and should not be invented by scaleRates
    expect_false("z0pre" %in% names(params@species_params))
    scaled <- scaleRates(params, 2)
    expect_false("z0pre" %in% names(scaled@species_params))
})

test_that("scaleRates scales ext_diffusion and D_ext when set", {
    params <- example_params()
    f <- 4
    scaled <- scaleRates(params, f)
    expect_equal(scaled@ext_diffusion, params@ext_diffusion * f)
    if ("D_ext" %in% names(params@species_params)) {
        expect_equal(scaled@species_params$D_ext,
                     params@species_params$D_ext * f)
    }
})

test_that("scaleRates produces time-rescaled dynamics", {
    # Scaling all rates by f is equivalent to running the original model for
    # f times as long with f times larger time steps.
    # project(scaleRates(p, f), t_max=T, dt=dt) should equal
    # project(p, t_max=f*T, dt=f*dt) since each step covers the same
    # biological time.
    f <- 2
    params <- NS_params_small
    params_scaled <- scaleRates(params, f)

    sim_orig   <- project(params,        t_max = 1,   dt = 0.1,  t_save = 1)
    sim_scaled <- project(params_scaled, t_max = 0.5, dt = 0.05, t_save = 0.5)

    # Compare raw arrays to avoid the params attribute carried by ArraySpeciesBySize
    expect_equal(sim_scaled@n[dim(sim_scaled@n)[1], , ],
                 sim_orig@n[dim(sim_orig@n)[1], , ],
                 tolerance = 1e-10)
    expect_equal(sim_scaled@n_pp[dim(sim_scaled@n_pp)[1], ],
                 sim_orig@n_pp[dim(sim_orig@n_pp)[1], ],
                 tolerance = 1e-10)
})
