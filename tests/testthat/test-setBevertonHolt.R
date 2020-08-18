# setBevertonHolt works ----
test_that("setBevertonHolt works", {
    params <- NS_params
    rdd <- getRDD(params)
    expect_error(setBevertonHolt(params, c(1, rep(2, 10), 1)),
                 "Elements 1, 12 of R_factor > 1 are not true")
    R_factor <- c(Inf, rep(5, 11))
    params <- setBevertonHolt(params, R_factor)
    expect_equivalent(params@species_params$R_max, rdd * R_factor)
    expect_equal(getRDD(params), rdd)
    # An R_factor of Inf means no density dependence
    expect_equal(getRDD(params)[[1]], getRDI(params)[[1]])
})

