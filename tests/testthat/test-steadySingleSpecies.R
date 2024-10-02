test_that("steadySingleSpecies only affects selected species", {
    params <- steadySingleSpecies(NS_params, species = "Cod")
    # Haddock unaffected
    expect_identical(params@initial_n["Haddock", ], 
                     NS_params@initial_n["Haddock", ])
    # but Cod changed
    expect_gt(params@initial_n["Cod", 100], 
              NS_params@initial_n["Cod", 100])
    # Test that steadySingleSpecies updates time_modified
    expect_false(identical(NS_params@time_modified, params@time_modified))
})

test_that("steadySingleSpecies is idempotent on single-species model", {
    ss <- newSingleSpeciesParams()
    ss2 <- steadySingleSpecies(ss)
    expect_unchanged(ss, ss2)
})

test_that("steadySingleSpecies `keep` argument works", {
    params <- steadySingleSpecies(NS_params, species = 1:2)
    expect_equal(params@initial_n[1, 1], NS_params@initial_n[1, 1])
    params <- steadySingleSpecies(NS_params, species = 3, keep = "biomass")
    expect_equal(getBiomass(params)[3], getBiomass(NS_params)[3])
    params <- steadySingleSpecies(NS_params, species = 3, keep = "number")
    expect_equal(getN(params)[3], getN(NS_params)[3])
    expect_gt(getBiomass(params)[3], getBiomass(NS_params)[3])
})
