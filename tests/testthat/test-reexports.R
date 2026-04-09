test_that("melt is re-exported from reshape2", {
    expect_identical(melt, reshape2::melt)
})
