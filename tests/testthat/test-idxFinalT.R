test_that("idxFinalT matches final time index and results", {
    idx <- idxFinalT(NS_sim)
    expect_equal(idx, length(getTimes(NS_sim)))
    expect_identical(N(NS_sim)[idx, , ], finalN(NS_sim))
    expect_identical(NResource(NS_sim)[idx, ], finalNResource(NS_sim))
})


