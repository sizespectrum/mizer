

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

test_that("We can set and get resource parameters", {
    params <- NS_params
    # get
    rp <- resource_params(params)
    expect_identical(rp, params@resource_params)
    # set
    resource_params(params)$test <- "hi"
    expect_identical(params@resource_params$test, "hi")
})
