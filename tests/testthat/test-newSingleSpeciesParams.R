# newSingleSpeciesParams ----
test_that("newSingleSpeciesParams works", {
    params <- newSingleSpeciesParams()
    no_w <- length(params@w)
    expect_equal(dim(params@initial_n), c(1, no_w))
    expect_equal(params@w[1], params@species_params$w_min)
    expect_equal(params@w[no_w], params@species_params$w_max)
    sim <- project(params, t_max = 1)
    expect_equal(sim@n[1, 1, ], sim@n[2, 1, ])
})
