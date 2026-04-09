test_that("idxFinalT matches final time index and results", {
    idx <- idxFinalT(NS_sim)
    expect_equal(idx, length(getTimes(NS_sim)))
    expect_equal(N(NS_sim)[idx, , ], finalN(NS_sim), ignore_attr = TRUE)
    expect_identical(NResource(NS_sim)[idx, ], finalNResource(NS_sim))
})

test_that("getEffort and getParams expose stored simulation metadata", {
    sim <- project(NS_params, t_max = 1, effort = 2, progress_bar = FALSE)

    expect_identical(getParams(sim), sim@params)
    expect_identical(getEffort(sim), sim@effort)
    expect_true(all(getEffort(sim) == 2))
})

test_that("N, NResource and getTimes expose stored arrays and times", {
    sim <- project(NS_params, t_max = 1, t_save = 0.5, progress_bar = FALSE)

    expect_identical(N(sim), sim@n)
    expect_identical(NResource(sim), sim@n_pp)
    expect_identical(getTimes(sim), c(0, 0.5, 1))
})
