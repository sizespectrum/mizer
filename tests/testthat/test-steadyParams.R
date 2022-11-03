local_edition(3)
test_that("steadyParams works", {
    params <- newTraitParams(no_sp = 4, no_w = 30,
                             n = 2/3, lambda = 2 + 3/4 - 2/3,
                             max_w_inf = 1e3, min_w = 1e-4,
                             w_pp_cutoff = 10, ks = 4,
                             reproduction_level = 0) |>
        suppressMessages()
    species_params(params)$gamma[2] <- 2000
    p <- steadyParams(params, t_per = 2) |>
        suppressMessages()
    expect_snapshot_value(getRDD(p), style = "deparse")
    # and works the same when returning sim
    sim <- steadyParams(params, t_per = 2, return_sim = TRUE) |>
        suppressMessages()
    expect_s4_class(sim, "MizerSim")
    expect_snapshot_value(getRDD(sim@params), style = "deparse")
})

test_that("steadyParams() preserves reproduction function", {
    params <- NS_params
    params@rates_funcs$RDD <- "noRDD"
    p2 <- steadyParams(params, t_per = 1, t_max = 1) |> suppressMessages()
    expect_equal(params@rates_funcs$RDD, "noRDD")
})
test_that("steadyParams() preserves erepro", {
    params <- NS_params
    species_params(params)$R_max <- 1.01 * species_params(params)$R_max
    p2 <- steadyParams(params, t_per = 1, preserve = "erepro") |>
        suppressMessages()
    expect_equal(p2@species_params$erepro, params@species_params$erepro)
})
test_that("steady() preserves reproduction level", {
    params <- NS_params
    species_params(params)$R_max <- 1.01 * species_params(params)$R_max
    p2 <- steadyParams(params, t_per = 1, preserve = "reproduction_level") |>
        suppressMessages()
    expect_equal(getReproductionLevel(p2), getReproductionLevel(params))
})
test_that("steady() preserves R_max", {
    params <- NS_params
    species_params(params)$erepro <- 1.01 * species_params(params)$erepro
    expect_warning(p2 <- steadyParams(params, t_per = 1, preserve = "R_max") |>
                       suppressMessages(),
                   "The following species require an unrealistic reproductive")
    expect_equal(p2@species_params$R_max, params@species_params$R_max)
})
