test_that("steady works", {
    expect_message(params <- newTraitParams(no_sp = 4, no_w = 30, R_factor = Inf,
                                            n = 2/3, lambda = 2 + 3/4 - 2/3,
                                            max_w_inf = 1e3, min_w = 1e-4,
                                            w_pp_cutoff = 10, ks = 4),
                   "Increased no_w to 36")
    params@species_params$gamma[2] <- 2000
    params <- setSearchVolume(params)
    p <- steady(params, t_per = 2)
    expect_known_value(getRDI(p), "values/steady")
    # and works the same when returning sim
    sim <- steady(params, t_per = 2, return_sim = TRUE)
    expect_is(sim, "MizerSim")
    expect_known_value(getRDI(sim@params), "values/steady")
})
