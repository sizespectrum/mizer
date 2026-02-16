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
    
    # Since N is fixed at the boundary (and potentially inconsistent with R_dd),
    # projecting with `project()` (which uses R_dd) will immediately change the boundary 
    # and propagate effects.
    # Instead, we verify that the calculated N satisfies the steady state transport equation
    # with the fixed boundary condition.
    
    # We need to access internal functions
    get_transport_coefs <- mizer:::get_transport_coefs
    
    # Calculate coefficients
    # steadySingleSpecies uses rates from the *original* state
    # and dt=1
    growth <- getEGrowth(params_orig)
    mort <- getMort(params_orig)
    rates <- list(e_growth = growth, mort = mort)
    dt <- 1
    
    coefs <- get_transport_coefs(params_orig, params_orig@initial_n, params_orig@initial_n_pp, 
                                 params_orig@initial_n_other, rates, dt)
                                 
    # Check residual for Cod
    sp <- species
    n <- params@initial_n[sp, ]
    
    a <- coefs$a[sp, ]
    b <- coefs$b[sp, ] - 1 # Adjust for steady state
    c <- coefs$c[sp, ]
    
    # Boundary correction used in steadySingleSpecies
    w_min_idx <- params@w_min_idx[sp]
    # In steadySingleSpecies we set b[w_min_idx] = 1 and c[w_min_idx] = 0.
    # We should replicate that here to verify the interior.
    b[w_min_idx] <- 1
    c[w_min_idx] <- 0
    a[w_min_idx] <- 0
    
    # Calculate residual A*N_{i-1} + B*N_i + C*N_{i+1}
    # For interior points within the solved range
    w_max_idx <- sum(params@w <= params@species_params[sp, "w_max"])
    
    residuals <- numeric(w_max_idx)
    # Start from w_min_idx + 1 (the first interior node)
    # End at w_max_idx (the last solved node)
    for (i in (w_min_idx + 1):w_max_idx) {
        # Note: n[i+1] will be 0 if i = w_max_idx, which is consistent with the solver 
        # (assuming 0 flux from above or just absorbing boundary)
        val_next <- if (i < length(n)) n[i+1] else 0
        residuals[i] <- a[i] * n[i-1] + b[i] * n[i] + c[i] * val_next
    }
    
    # Check max residual relative to N
    # We exclude the boundary point because we fixed it explicitly.
    
    valid_range <- (w_min_idx + 1):w_max_idx
    max_rel_resid <- max(abs(residuals[valid_range]) / n[valid_range])
    
    expect_lt(max_rel_resid, 1e-10)
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
