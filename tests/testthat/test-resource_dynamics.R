

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

    rp2 <- params@resource_params
    rp2$lambda <- rp2$lambda + 0.1
    resource_params(params) <- rp2
    expect_identical(resource_params(params), rp2)
})

test_that("Deprecated resource getters warn and delegate to accessors", {
    params <- NS_params
    expect_warning(expect_identical(getResourceDynamics(params),
                                    resource_dynamics(params)),
                   "deprecated")
    expect_warning(expect_identical(getResourceLevel(params),
                                    resource_level(params)),
                   "deprecated")
    expect_warning(expect_identical(getResourceRate(params),
                                    resource_rate(params)),
                   "deprecated")
    expect_warning(expect_identical(getResourceCapacity(params),
                                    resource_capacity(params)),
                   "deprecated")
})

test_that("resource_params<- updates `time_modified`", {
    params <- NS_params
    resource_params(params) <- resource_params(NS_params)
    expect_false(identical(params@time_modified, NS_params@time_modified))
})
