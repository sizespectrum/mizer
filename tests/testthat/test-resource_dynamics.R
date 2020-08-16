# resource_semichemostat ----
test_that("resource_semichemostat works", {
    params <- NS_params
    x <- resource_semichemostat(params, 
                                n = params@initial_n,
                                n_pp = params@initial_n_pp,
                                n_other = params@initial_n_other,
                                rates = getRates(params),
                                t = 0, dt = 0.1)
    # result should not be above carrying capacity
    expect_true(all(x <= params@cc_pp))
})

# resource_constant ----
test_that("resource_constant works", {
    params <- NS_params
    x <- resource_constant(params,
                           n = params@initial_n,
                           n_pp = params@initial_n_pp,
                           n_other = params@initial_n_other,
                           rates = getRates(params),
                           t = 0, dt = 0.1)
    # result should not be above carrying capacity
    expect_equal(x, params@initial_n_pp)
})