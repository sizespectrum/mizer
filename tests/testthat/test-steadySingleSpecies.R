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

test_that("steadySingleSpecies produces steady state with diffusion", {
    params <- NS_params
    
    # Enable diffusion for Cod
    species <- "Cod"
    p <- params@species_params[species, "p"]
    d <- 0.1 * params@w^p
    diffusion(params)[species, ] <- d
    
    # Keep original params to calculate rates that steadySingleSpecies used
    params_orig <- params
    
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

test_that("steadySingleSpecies warns when growth stops after maturity", {
    # Create a simple params object
    params <- newSingleSpeciesParams()

    # Get the species and find indices for maturity and max size
    w_mat <- params@species_params$w_mat[1]
    w_max <- params@species_params$w_max[1]
    w_mat_idx <- sum(params@w <= w_mat)
    w_max_idx <- sum(params@w <= w_max)

    # Increase metabolic rate significantly after maturity
    if (w_mat_idx < length(params@w)) {
        params@metab[1, (w_mat_idx + 1):length(params@w)] <-
            params@metab[1, (w_mat_idx + 1):length(params@w)] * 1000
    }

    expect_warning(steadySingleSpecies(params, species = 1),
                   "has zero growth rate after maturity size")
})
