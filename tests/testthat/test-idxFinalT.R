test_that("idxFinalT matches final time index and results", {
    idx <- idxFinalT(NS_sim)
    expect_equal(idx, length(getTimes(NS_sim)))
    expect_identical(N(NS_sim)[idx, , ], finalN(NS_sim))
    expect_identical(NResource(NS_sim)[idx, ], finalNResource(NS_sim))
})

test_that("getEffort and getParams expose stored simulation metadata", {
    sim <- project(NS_params, t_max = 1, effort = 2, progress_bar = FALSE)

    expect_identical(getParams(sim), sim@params)
    expect_identical(getEffort(sim), sim@effort)
    expect_true(all(getEffort(sim) == 2))
})

