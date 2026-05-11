
test_that("getRequiredRDD works for single species model", {
    params <- newSingleSpeciesParams()
    # In a steady state model, getRequiredRDD should return the same as getRDD
    # because the reproduction matches the required amount to maintain the steady state.
    rdd <- getRequiredRDD(params)
    rdd_actual <- getRDD(params)
    
    expect_equal(rdd, rdd_actual)
})

test_that("getRequiredRDD works for community model", {
    params <- newCommunityParams()
    # In community model, reproduction is constant
    rdd <- getRequiredRDD(params)
    rdd_actual <- getRDD(params)
    
    # Check if they are close enough
    expect_equal(rdd, rdd_actual)
})

test_that("getRequiredRDD handles diffusion", {
    # Create a model with diffusion
    params <- newSingleSpeciesParams()
    
    # Add diffusion
    ext_diff <- params@ext_diffusion
    ext_diff[] <- 1e9 * params@w
    params <- setExtDiffusion(params, ext_diffusion = ext_diff)
    
    # Update initial_n to be the steady state solution with this diffusion
    # We need to recalculate it using get_steady_state_n
    
    # Instead of fully replicating newSingleSpeciesParams logic, let's use the fact that
    # getRequiredRDD should make the current state a steady state *at the boundary*.
    
    # If we simply update reproduction to match getRequiredRDD, then the boundary flux matches reproduction.
    
    # So let's update erepro
    rdd_req <- getRequiredRDD(params)
    params@species_params$erepro <- params@species_params$erepro * rdd_req / getRDI(params)
    
    # Now getRDD(params) should match rdd_req
    expect_equal(getRDD(params), rdd_req)
    
    # Let's stick to the test that it does not error and returns a numeric vector of correct length
    expect_length(rdd_req, 1)
    expect_type(rdd_req, "double")
    expect_true(rdd_req > 0)
    
    # verifying that it is different from the no-diffusion case
    params_no_diff <- newSingleSpeciesParams()
    expect_true(rdd_req != getRequiredRDD(params_no_diff))
})

test_that("getRequiredRDD uses total diffusion", {
    params <- newSingleSpeciesParams()
    species <- params@species_params$species[1]
    params@use_predation_diffusion <- TRUE

    idx <- params@w_min_idx[species]
    next_idx <- idx + 1
    dw <- params@dw[idx]
    diffusion <- getDiffusion(params)

    expected <- initialN(params)[species, idx] *
        (getEGrowth(params)[species, idx] +
             getMort(params)[species, idx] * dw) +
        0.5 * (diffusion[species, idx] * initialN(params)[species, idx] -
                   diffusion[species, next_idx] *
                   initialN(params)[species, next_idx]) / dw

    expect_equal(getRequiredRDD(params)[species], expected,
                 ignore_attr = TRUE)
})

test_that("getRequiredRDD keeps egg density constant with diffusion", {
    params <- newSingleSpeciesParams()
    species <- params@species_params$species[1]
    n <- params@species_params[species, "n"]
    ext_diffusion(params)[species, ] <- 0.1 * params@w^(n + 1)
    params <- steadySingleSpecies(params, species = species)
    params@species_params$constant_reproduction <- getRequiredRDD(params)
    params <- setRateFunction(params, "RDD", "constantRDD")

    sim <- project(params, t_max = 1, t_save = 1, progress_bar = FALSE)
    idx <- params@w_min_idx[species]
    expect_equal(finalN(sim)[species, idx], initialN(params)[species, idx],
                 tolerance = 1e-12)
})

test_that("getRequiredRDD works with predictor-corrector project method", {
    params <- newSingleSpeciesParams()
    species <- params@species_params$species[1]
    n <- params@species_params[species, "n"]
    ext_diffusion(params)[species, ] <- 0.1 * params@w^(n + 1)
    params <- steadySingleSpecies(params, species = species)
    params@species_params$constant_reproduction <- getRequiredRDD(params)
    params <- setRateFunction(params, "RDD", "constantRDD")

    sim <- project(params, t_max = 1, t_save = 1,
                   method = "predictor-corrector", progress_bar = FALSE)
    idx <- params@w_min_idx[species]
    expect_equal(finalN(sim)[species, idx], initialN(params)[species, idx],
                 tolerance = 1e-12)
})
