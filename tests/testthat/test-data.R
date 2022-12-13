test_that("in NS_params linecolours are all distinct", {
    expect_equal(length(unique(NS_params@linecolour)),
                 length(NS_params@linecolour))
})
