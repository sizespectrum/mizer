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

    expect_equal(unclass(N(sim)), sim@n, ignore_attr = TRUE)
    expect_identical(NResource(sim), sim@n_pp)
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
