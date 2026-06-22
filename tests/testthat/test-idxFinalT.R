test_that("idxFinalT matches final time index and results", {
    idx <- idxFinalT(NS_sim_small)
    expect_equal(idx, length(getTimes(NS_sim_small)))
    expect_equal(N(NS_sim_small)[idx, , ], finalN(NS_sim_small), ignore_attr = TRUE)
    expect_identical(NResource(NS_sim_small)[idx, ], finalNResource(NS_sim_small))
})

test_that("getEffort works", {
    sim <- project(NS_params_small, t_max = 1, effort = 2, progress_bar = FALSE)
    expect_identical(getEffort(sim), sim@effort)
    expect_true(all(getEffort(sim) == 2))
})

test_that("N, NResource and getTimes expose stored arrays and times", {
    sim <- project(NS_params_small, t_max = 1, t_save = 0.5, progress_bar = FALSE)

    expect_equal(unclass(N(sim)), sim@n, ignore_attr = TRUE)
    expect_equal(NResource(sim), sim@n_pp, ignore_attr = TRUE)
    expect_identical(getTimes(sim), c(0, 0.5, 1))
})

test_that("MizerSim accessors are registered as S3 methods", {
    expect_identical(utils::getS3method("validSim", "MizerSim"), validSim.MizerSim)
    expect_identical(utils::getS3method("N", "MizerSim"), N.MizerSim)
    expect_identical(utils::getS3method("NResource", "MizerSim"), NResource.MizerSim)
    expect_identical(utils::getS3method("finalN", "MizerSim"), finalN.MizerSim)
    expect_identical(utils::getS3method("finalNResource", "MizerSim"), finalNResource.MizerSim)
    expect_identical(utils::getS3method("idxFinalT", "MizerSim"), idxFinalT.MizerSim)
    expect_identical(utils::getS3method("getTimes", "MizerSim"), getTimes.MizerSim)
    expect_identical(utils::getS3method("getEffort", "MizerSim"), getEffort.MizerSim)
    expect_identical(utils::getS3method("getParams", "MizerSim"), getParams.MizerSim)
})
