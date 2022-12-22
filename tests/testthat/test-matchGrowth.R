local_edition(3)

test_that("matchGrowth only affects selected species", {
    params <- matchGrowth(NS_params, species = "Cod")
    # Haddock unaffected
    expect_identical(params@initial_n["Haddock", ], 
                     NS_params@initial_n["Haddock", ])
    # but Cod changed
    expect_gt(params@initial_n["Cod", 100], 
              NS_params@initial_n["Cod", 100])
    # and changes again when called again
    params2 <- matchGrowth(params, species = "Cod")
    expect_lt(params2@initial_n["Cod", 100], 
                     params@initial_n["Cod", 100])
})

test_that("matchGrowth is idempotent on single species", {
    ss <- newSingleSpeciesParams()
    ss2 <- matchGrowth(ss)
    expect_equal(ss@initial_n, ss2@initial_n)
})

test_that("matchGrowth `keep` argument works", {
    params <- matchGrowth(NS_params, species = 1:2)
    expect_equal(params@initial_n[1, 1], NS_params@initial_n[1, 1])
    params <- matchGrowth(NS_params, species = 3, keep = "biomass")
    expect_equal(getBiomass(params)[3], getBiomass(NS_params)[3])
    params <- matchGrowth(NS_params, species = 3, keep = "number")
    expect_equal(getN(params)[3], getN(NS_params)[3])
    expect_gt(getBiomass(params)[3], getBiomass(NS_params)[3])
})

test_that("matchGrowth does nothing when no info is given", {
    params <- NS_params
    params@species_params$k_vb <- NULL
    params2 <- matchGrowth(params, species = "Cod")
    expect_identical(params2@initial_n["Cod", ], 
                     params@initial_n["Cod", ])
})
