test_that("steadySingleSpecies only affects abundance of selected species", {
    params1 <- NS_params_small
    sp_names <- params1@species_params$species
    species1 <- sp_names[3]
    species2 <- sp_names[2]
    # make sure it is not in steady state
    params1@initial_n[, 10:15] <- params1@initial_n[, 10:15] * 2

    params2 <- steadySingleSpecies(params1, species = species1) |>
        suppressWarnings()
    # Herring unaffected
    expect_identical(params2@initial_n[species2, ],
                     params1@initial_n[species2, ])
    # but Cod changed
    expect_lt(params2@initial_n[species1, 20],
              params1@initial_n[species1, 20])
    # Test that steadySingleSpecies updates time_modified
    expect_false(identical(params1@time_modified, params2@time_modified))
    # Nothing else changed
    params2@initial_n <- params1@initial_n
    params2@time_modified <- params1@time_modified
    expect_identical(params1, params2)
})

test_that("steadySingleSpecies is idempotent on single-species model", {
    ss <- single_sp_params
    ss2 <- steadySingleSpecies(ss)
    expect_unchanged(ss, ss2)
})

test_that("steadySingleSpecies `keep` argument works", {
    params <- steadySingleSpecies(NS_params_small, species = 1:2)
    expect_equal(params@initial_n[1, 1], NS_params_small@initial_n[1, 1])
    params <- steadySingleSpecies(NS_params_small, species = 3, keep = "biomass")
    expect_equal(getBiomass(params)[3], getBiomass(NS_params_small)[3])
    params <- steadySingleSpecies(NS_params_small, species = 3, keep = "number")
    expect_equal(getN(params)[3], getN(NS_params_small)[3])
    expect_false(isTRUE(all.equal(getBiomass(params)[3], getBiomass(NS_params_small)[3])))
})

test_that("steadySingleSpecies produces steady state with diffusion", {
    # Use a single-species model so that changing the species abundance does
    # not affect its own growth and mortality rates via self-predation.
    params <- single_sp_params
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

    sim_pc <- project(params, t_max = 5, method = "predictor-corrector")
    final_n_pc <- finalN(sim_pc)[species, ]
    rel_error_pc <- abs(initial_n - final_n_pc) / initial_n
    max_rel_error_pc <- max(rel_error_pc[valid_idx], na.rm = TRUE)
    expect_lt(max_rel_error_pc, 1e-10)
})

test_that("steadySingleSpecies holds abundance at zero above w_max", {
    # example_params() has diffusion on species 1, which would otherwise carry
    # density above w_max.
    params <- example_params()
    params <- steadySingleSpecies(params)

    w_top <- mizer:::support_top_idx(params)
    no_w <- length(params@w)
    checked_any <- FALSE
    for (sp in seq_len(nrow(params@species_params))) {
        if (w_top[sp] < no_w) {
            checked_any <- TRUE
            expect_true(all(params@initial_n[sp, (w_top[sp] + 1):no_w] == 0))
        }
    }
    # Make sure the assertion above was actually exercised.
    expect_true(checked_any)
})

test_that("steadySingleSpecies errors when growth stops before maturity", {
    # Create a simple params object
    params <- single_sp_params

    # Artificially set growth rate to zero before maturity
    # by setting very high metabolic rate
    params@metab[1, ] <- params@metab[1, ] * 1000

    expect_error(steadySingleSpecies(params, species = 1),
                 "cannot grow to maturity")
})
