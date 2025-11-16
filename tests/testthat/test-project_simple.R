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
    expect_equal(out$n, N(sim)[2, , ], tolerance = 1e-12)
    expect_equal(out$n_pp, NResource(sim)[2, ], tolerance = 1e-12)
})


