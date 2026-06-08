test_that("matchGrowth only affects selected species", {
    sp <- NS_params_small@species_params$species
    species1 <- sp[3]
    species2 <- sp[2]
    params <- matchGrowth(NS_params_small, species = species1)
    # Herring unaffected
    expect_identical(params@initial_n[species2, ],
                     NS_params_small@initial_n[species2, ])
    # but Cod changed
    expect_gt(params@initial_n[species1, 10],
              NS_params_small@initial_n[species1, 10])
    # and changes again when called again
    params2 <- matchGrowth(params, species = species1)
    expect_lt(params2@initial_n[species1, 10],
                     params@initial_n[species1, 10])
})

test_that("matchGrowth is idempotent on single species", {
    ss <- newSingleSpeciesParams()
    ss2 <- matchGrowth(ss)
    expect_equal(ss@initial_n, ss2@initial_n)
})

test_that("matchGrowth `keep` argument works", {
    params <- matchGrowth(NS_params_small, species = 1:2)
    expect_equal(params@initial_n[1, 1], NS_params_small@initial_n[1, 1])
    params <- matchGrowth(NS_params_small, species = 3, keep = "biomass")
    expect_equal(getBiomass(params)[3], getBiomass(NS_params_small)[3])
    params <- matchGrowth(NS_params_small, species = 3, keep = "number")
    expect_equal(getN(params)[3], getN(NS_params_small)[3])
    expect_false(isTRUE(all.equal(getBiomass(params)[3], getBiomass(NS_params_small)[3])))
})

test_that("matchGrowth does nothing when no info is given", {
    params <- NS_params_small
    params@species_params$k_vb <- NULL
    sp_name <- params@species_params$species[3]
    params2 <- matchGrowth(params, species = sp_name)
    expect_identical(params2@initial_n[sp_name, ], 
                     params@initial_n[sp_name, ])
})

test_that("matchGrowth rescales rates and species parameters by age ratio", {
    params <- NS_params_small
    i <- which(species_params(params)$species == "Cod")
    sp <- set_species_param_default(species_params(params), "age_mat", NA)
    if (all(c("k_vb", "w_inf") %in% names(sp))) {
        sp <- set_species_param_default(sp, "age_mat", age_mat_vB(params))
    }
    factor <- age_mat(params)[[i]] / sp$age_mat[[i]]

    params2 <- matchGrowth(params, species = i)
    expect_equal(params2@search_vol[i, ], params@search_vol[i, ] * factor,
                 ignore_attr = TRUE)
    expect_equal(params2@intake_max[i, ], params@intake_max[i, ] * factor,
                 ignore_attr = TRUE)
    expect_equal(params2@metab[i, ], params@metab[i, ] * factor,
                 ignore_attr = TRUE)
    expect_equal(params2@ext_encounter[[i]], params@ext_encounter[[i]] * factor)
    expect_equal(species_params(params2)$gamma[[i]],
                 species_params(params)$gamma[[i]] * factor)
    expect_equal(species_params(params2)$h[[i]],
                 species_params(params)$h[[i]] * factor)
})
