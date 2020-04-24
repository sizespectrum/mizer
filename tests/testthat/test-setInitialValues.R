params <- NS_params
sim <- project(params, t_max = 0.1, t_save = 0.1)

test_that("We can set and get initial values", {
    expect_identical(initialN(sim), initialN(params))
    expect_identical(initialNResource(sim), initialNResource(params))
    initialN(params) <- params@metab
    expect_identical(initialN(params), params@metab)
    initialNResource(params) <- params@cc_pp
    expect_identical(initialNResource(params), params@cc_pp)
    params <- setInitialValues(params, sim)
    expect_identical(finalN(sim), initialN(params))
    expect_identical(finalNResource(sim), initialNResource(params))
})

test_that("Can set initial values in a model with a single species", {
    species_params <- NS_species_params[1, ]
    params <- newMultispeciesParams(species_params)
    sim <- project(params, t_max=0.1)
    p <- setInitialValues(params, sim)
    expect_identical(finalN(sim), initialN(params))
})
