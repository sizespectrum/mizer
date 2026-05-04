test_that("steadySingleSpecies only affects abundance of selected species", {
    params1 <- NS_params
    sp_names <- params1@species_params$species
    species1 <- sp_names[11]
    species2 <- sp_names[10]
    # make sure it is not in steady state
    params1@initial_n[,50:80] <- params1@initial_n[,50:80] * 2

    params2 <- steadySingleSpecies(params1, species = species1) |>
        suppressWarnings()
    # Haddock unaffected
    expect_identical(params2@initial_n[species2, ],
                     params1@initial_n[species2, ])
    # but Cod changed
    expect_lt(params2@initial_n[species1, 100],
              params1@initial_n[species1, 100])
    # Test that steadySingleSpecies updates time_modified
    expect_false(identical(params1@time_modified, params2@time_modified))
    # Nothing else changed
    params2@initial_n <- params1@initial_n
    params2@time_modified <- params1@time_modified
    expect_identical(params1, params2)
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

test_that("steadySingleSpecies produces steady state with diffusion", {
    # Use a single-species model so that changing the species abundance does
    # not affect its own growth and mortality rates via self-predation.
    params <- newSingleSpeciesParams()
    species <- params@species_params$species[1]
    n <- params@species_params[species, "n"]
    d <- 0.1 * params@w^(n + 1)
    ext_diffusion(params)[species, ] <- d

    # Increase minimum size to test boundary condition
    params@w_min_idx[species] <- 10
    params@species_params[species, "w_min"] <- params@w[10]

    params <- steadySingleSpecies(params, species = species)

    # Now that we have the steady state, we can use setBevertonHolt() to
    # set the reproduction parameters to values that are consistent with it.
    suppressWarnings(params <- setBevertonHolt(params, reproduction_level = 0.5))

    # And then the steady state should be preserved by project()
    sim <- project(params, t_max = 5)
    initial_n <- params@initial_n[species, ]
    final_n <- finalN(sim)[species, ]
    rel_error <- abs(initial_n - final_n) / initial_n
    # Ignore indices where initial_n is very small/zero to avoid division by zero or numerical noise
    valid_idx <- initial_n > 1e-20
    max_rel_error <- max(rel_error[valid_idx], na.rm = TRUE)
    expect_lt(max_rel_error, 1e-10)
})

test_that("steadySingleSpecies errors when growth stops before maturity", {
    # Create a simple params object
    params <- newSingleSpeciesParams()

    # Artificially set growth rate to zero before maturity
    # by setting very high metabolic rate
    params@metab[1, ] <- params@metab[1, ] * 1000

    expect_error(steadySingleSpecies(params, species = 1),
                 "cannot grow to maturity")
})
