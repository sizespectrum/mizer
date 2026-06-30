# newSingleSpeciesParams ----
test_that("newSingleSpeciesParams works", {
    params <- newSingleSpeciesParams()
    no_w <- length(params@w)
    expect_equal(dim(params@initial_n), c(1, no_w))
    expect_equal(params@w[1], params@species_params$w_min, ignore_attr = TRUE)
    expect_equal(params@w[no_w], params@species_params$w_max, ignore_attr = TRUE)
    sim <- project(params, t_max = 1)
    expect_equal(sim@n[1, 1, ], sim@n[2, 1, ])
})

test_that("newSingleSpeciesParams documents and applies grid and deprecation behaviour", {
    expect_message(
        params <- newSingleSpeciesParams(w_max = 100, w_min = 0.001, no_w = 2),
        "Increased no_w to 26"
    )
    expect_equal(length(w(params)), 26)

    expect_snapshot_warning(params2 <- newSingleSpeciesParams(R_factor = Inf))
    expect_equal(unname(getReproductionLevel(params2)), rep(0, nrow(species_params(params2))))
})

test_that("Sets given_species_params", {
    # Calling `given_species_params<-()` should not make a change
    params <- newSingleSpeciesParams()
    p2 <- params
    given_species_params(p2) <- given_species_params(p2)
    expect_unchanged(p2, params)
})
