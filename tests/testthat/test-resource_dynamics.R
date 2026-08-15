

# resource_constant ----
test_that("resource_constant works", {
    params <- NS_params_small
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
    params <- NS_params_small
    # get
    rp <- resource_params(params)
    expect_identical(rp, params@resource_params)
    # set
    resource_params(params)$test <- "hi"
    expect_identical(params@resource_params$test, "hi")

    rp2 <- params@resource_params
    rp2$kappa <- 1e13
    rp2$lambda <- rp2$lambda + 0.1
    resource_params(params) <- rp2
    expect_identical(resource_params(params), rp2)
})

test_that("Deprecated resource getters warn and delegate to accessors", {
    params <- NS_params_small
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
    params <- NS_params_small
    resource_params(params) <- resource_params(NS_params_small)
    expect_false(identical(params@time_modified, NS_params_small@time_modified))
})

test_that("resource_params<- refreshes calculated gamma and q (#497)", {
    params <- newMultispeciesParams(NS_species_params_small, no_w = 20,
                                    info_level = 0)
    gamma0 <- species_params(params)$gamma
    q0 <- species_params(params)$q

    # Protect the first species explicitly. The other entries remain
    # calculated even though the given columns now exist.
    params@given_species_params$gamma <- c(gamma0[[1]], NA, NA)
    params@given_species_params$q <- c(q0[[1]], NA, NA)

    resource_params(params)$lambda <- resource_params(params)$lambda + 0.15

    expect_equal(species_params(params)$gamma[[1]], gamma0[[1]])
    gamma_ratio <- as.numeric(species_params(params)$gamma[-1] / gamma0[-1])
    expect_false(isTRUE(all.equal(gamma_ratio, rep(1, 2))))
    expect_equal(species_params(params)$q[[1]], q0[[1]])
    expect_equal(species_params(params)$q[-1], q0[-1] + 0.15)

    params2 <- newMultispeciesParams(NS_species_params_small, no_w = 20,
                                     info_level = 0)
    gamma0 <- species_params(params2)$gamma
    q0 <- species_params(params2)$q
    resource_params(params2)$kappa <- 2 * resource_params(params2)$kappa

    expect_equal(species_params(params2)$gamma, gamma0 / 2,
                 ignore_attr = TRUE)
    expect_equal(species_params(params2)$q, q0, ignore_attr = TRUE)
})
