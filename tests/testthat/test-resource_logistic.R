# resource_logistic ----
test_that("resource_logistic preserves steady state", {
    # Set resource parameters so that we are at steady state
    params <- newTraitParams()
    params <- setResourceLogistic(
        params,
        resource_capacity = 2 * initialNResource(params)
        )
    x <- resource_logistic(params, 
                           n = params@initial_n,
                           n_pp = params@initial_n_pp,
                           n_other = params@initial_n_other,
                           rates = getRates(params),
                           t = 0, dt = 0.1,
                           resource_capacity = params@cc_pp,
                           resource_rate = params@rr_pp)
    expect_equal(x, params@initial_n_pp,
                 tolerance = 1e-15,
                 check.attributes = FALSE)
})


test_that("resource_logistic evolves towards steady state", {
    # This might go wrong if numerical errors are too large. This in turn
    # might happen when the rates are small. To achieve that we
    # reduce plankton mortality to a small value by scaling gamma down by
    # a large factor s and simultaneously scaling up the plankton abundance
    # so as to keep the income of the fish constant.
    s <- 1e10
    params <- newTraitParams()
    species_params(params)$gamma = species_params(params)$gamma / s
    initialNResource(params) <- initialNResource(params) * s
    # Set resource parameters so that we are at steady state
    params <- setResourceLogistic(params,
                          resource_capacity = 2 * initialNResource(params))
    # Now perturb a bit away from steady state
    n_pp_pert <- params@initial_n_pp * (1 + 1e-10)
    # run the dynamics for one time step and check that the perturbation has
    # not grown
    x <- resource_logistic(params, 
                           n = params@initial_n,
                           n_pp = n_pp_pert,
                           n_other = params@initial_n_other,
                           rates = getRates(params),
                           t = 0, dt = 0.1,
                           resource_capacity = params@cc_pp,
                           resource_rate = params@rr_pp)
    expect_lte(max(abs(x - n_pp_pert)), 0)
})