test_that("project_simple matches project for one saved step", {
    params <- NS_params
    dt <- 0.1
    t_save <- 1
    steps <- as.integer(t_save / dt)
    effort <- getInitialEffort(params)
    res_fn <- get(params@resource_dynamics)
    other_fns <- lapply(params@other_dynamics, get)
    rates_fns <- lapply(params@rates_funcs, get)
    start <- list(
        n = params@initial_n,
        n_pp = params@initial_n_pp,
        n_other = params@initial_n_other
    )
    out <- project_simple(
        params,
        n = start$n,
        n_pp = start$n_pp,
        n_other = start$n_other,
        t = 0,
        dt = dt,
        steps = steps,
        effort = effort,
        resource_dynamics_fn = res_fn,
        other_dynamics_fns = other_fns,
        rates_fns = rates_fns
    )
    sim <- project(params, t_max = t_save, dt = dt, t_save = t_save, progress_bar = FALSE)
    expected_n <- params@initial_n
    expected_n[] <- N(sim)[2, , ]
    expect_equal(out$n, expected_n, tolerance = 1e-12)
    expect_equal(out$n_pp, NResource(sim)[2, ], tolerance = 1e-12)
})

test_that("project_simple returns rates from the final update step", {
    params <- NS_params
    dt <- 0.1
    steps <- 2
    effort <- getInitialEffort(params)
    res_fn <- get(params@resource_dynamics)
    other_fns <- lapply(params@other_dynamics, get)
    rates_fns <- lapply(params@rates_funcs, get)
    out <- project_simple(
        params,
        n = params@initial_n,
        n_pp = params@initial_n_pp,
        n_other = params@initial_n_other,
        t = 0,
        dt = dt,
        steps = steps,
        effort = effort,
        resource_dynamics_fn = res_fn,
        other_dynamics_fns = other_fns,
        rates_fns = rates_fns
    )

    step1 <- project_simple(
        params,
        n = params@initial_n,
        n_pp = params@initial_n_pp,
        n_other = params@initial_n_other,
        t = 0,
        dt = dt,
        steps = 1,
        effort = effort,
        resource_dynamics_fn = res_fn,
        other_dynamics_fns = other_fns,
        rates_fns = rates_fns
    )
    expected_rates <- getRates(
        params,
        n = step1$n,
        n_pp = step1$n_pp,
        n_other = step1$n_other,
        effort = effort,
        t = dt
    )

    expect_identical(names(out$rates), names(expected_rates))
    expect_equal(out$rates$encounter, expected_rates$encounter)
    expect_equal(out$rates$rdd, expected_rates$rdd)
})

test_that("project_simple accepts predictor-corrector method", {
    params <- NS_params
    effort <- getInitialEffort(params)
    res_fn <- get(params@resource_dynamics)
    other_fns <- lapply(params@other_dynamics, get)
    rates_fns <- lapply(params@rates_funcs, get)

    out <- project_simple(
        params,
        n = params@initial_n,
        n_pp = params@initial_n_pp,
        n_other = params@initial_n_other,
        t = 0,
        dt = 0.1,
        steps = 1,
        effort = effort,
        resource_dynamics_fn = res_fn,
        other_dynamics_fns = other_fns,
        rates_fns = rates_fns,
        method = "predictor_corrector"
    )

    expect_true(all(is.finite(out$n)))
    expect_true(all(out$n >= 0))
})
