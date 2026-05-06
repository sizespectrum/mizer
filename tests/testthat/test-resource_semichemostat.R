test_that("resource_semichemostat preserves steady state", {
    # Set resource parameters so that we are at steady state
    params <- newTraitParams()
    params <- setResource(params,
                          resource_capacity = 2 * initialNResource(params),
                          resource_dynamics = "resource_semichemostat")
    x <- resource_semichemostat(params, 
                                n = params@initial_n,
                                n_pp = params@initial_n_pp,
                                n_other = params@initial_n_other,
                                rates = getRates(params),
                                t = 0, dt = 0.1,
                                resource_capacity = params@cc_pp,
                                resource_rate = params@rr_pp)
    expect_equal(x, params@initial_n_pp,
                 tolerance = 1e-15,
                 ignore_attr = TRUE)
})


test_that("resource_semichemostat evolves towards steady state", {
    # This might go wrong if numerical errors are too large. This in turn
    # might happen when the rates are small. To achieve that we
    # reduce plankton mortality to a small value by scaling gamma down by
    # a large factor s and simultaneously scaling up the plankton abundance
    # so as to keep the income of the fish constant.
    s <- 1e10
    params <- newTraitParams()
    species_params(params)$gamma <- species_params(params)$gamma / s
    initialNResource(params) <- initialNResource(params) * s
    # Set resource parameters so that we are at steady state
    params <- setResource(params,
                          resource_capacity = 2 * initialNResource(params),
                          resource_rate = getResourceMort(params),
                          balance = FALSE)
    # Now perturb a bit away from steady state
    n_pp_pert <- params@initial_n_pp * (1 + 1e-10)
    # run the dynamics for one time step and check that the perturbation has
    # not grown
    x <- resource_semichemostat(params, 
                                n = params@initial_n,
                                n_pp = n_pp_pert,
                                n_other = params@initial_n_other,
                                rates = getRates(params),
                                t = 0, dt = 0.1,
                                resource_capacity = params@cc_pp,
                                resource_rate = params@rr_pp)
    expect_lte(max(abs(x - n_pp_pert)), 0)
})


test_that("balance_resource_semichemostat works", {
    # semichemostat rate
    sc <- function(params) {
        r <- params@rr_pp
        c <- params@cc_pp
        N <- params@initial_n_pp
        mu <- getResourceMort(params)
        N_steady <- r * c / (mu + r)
        sel <- (mu + r) > 0
        all.equal(N_steady[sel], N[sel])
    }
    
    params <- NS_params
    
    # setting rate
    rate <- getResourceMort(params)
    params <- setResource(params, resource_rate = rate,
                          resource_dynamics = "resource_semichemostat")
    expect_identical(params@rr_pp, rate)
    expect_true(sc(params))
    
    # setting capacity
    capacity <- params@cc_pp * 2
    p1 <- setResource(params, resource_capacity = capacity)
    expect_identical(p1@cc_pp, capacity)
    expect_true(sc(p1))
})

test_that("balance_resource_semichemostat validates balancing inputs", {
    params <- NS_params

    rate <- getResourceMort(params)
    rate[1] <- 0
    expect_error(balance_resource_semichemostat(params,
                                                resource_rate = rate,
                                                resource_capacity = NULL),
                 "resource rate is zero while the resource mortality is not")

    expect_error(balance_resource_semichemostat(
        params,
        resource_rate = NULL,
        resource_capacity = initialNResource(params) * 0.9
    ), "capacity is less than the current abundance")
})

test_that("balance_resource_semichemostat nudges capacity to avoid division by zero", {
    params <- NS_params
    capacity <- initialNResource(params)
    death <- getResourceMort(params) * capacity != 0

    expect_warning(
        balanced <- balance_resource_semichemostat(
            params,
            resource_rate = NULL,
            resource_capacity = capacity
        ),
        "division by zero"
    )

    expect_true(all(is.finite(balanced$resource_rate)))
    expect_true(all(balanced$resource_capacity[death] > capacity[death]))
})

test_that("balance_resource_semichemostat keeps current rate when unidentifiable", {
    params <- newTraitParams()
    initialN(params)[] <- 0
    keep <- params@rr_pp[1]
    capacity <- initialNResource(params)

    balanced <- balance_resource_semichemostat(params,
                                               resource_rate = NULL,
                                               resource_capacity = capacity)

    expect_identical(balanced$resource_rate[1], keep)
    expect_identical(balanced$resource_capacity[1], capacity[1])
})
