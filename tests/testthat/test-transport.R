test_that("get_transport_coefs works correctly", {
    params <- newTraitParams(no_sp = 2)
    # Force different w_min to test zeroing logic
    params@species_params$w_min[2] <- 0.01 
    params@w_min_idx[2] <- which.min(abs(params@w - 0.01))
    
    n <- params@initial_n
    n[] <- 1 # Set n to 1 to make S checking easier
    
    dt <- 0.1
    recruitment_flux <- c(10, 20)
    
    # We need to access the internal function
    get_transport_coefs <- mizer:::get_transport_coefs
    getEGrowth <- mizer:::getEGrowth
    getMort <- mizer:::getMort
    
    coefs <- get_transport_coefs(params, n, getEGrowth(params), getMort(params), dt, recruitment_flux)
    
    # Check dimensions
    expect_equal(dim(coefs$a), dim(n))
    expect_equal(dim(coefs$b), dim(n))
    expect_equal(dim(coefs$c), dim(n))
    expect_equal(dim(coefs$S), dim(n))
    
    # Check zeroing out below w_min_idx
    w_min_idx_2 <- params@w_min_idx[2]
    expect_true(w_min_idx_2 > 1) # Ensure we are actually testing something interesting
    
    expect_true(all(coefs$a[2, 1:(w_min_idx_2 - 1)] == 0))
    expect_true(all(coefs$b[2, 1:(w_min_idx_2 - 1)] == 0))
    expect_true(all(coefs$c[2, 1:(w_min_idx_2 - 1)] == 0))
    expect_true(all(coefs$S[2, 1:(w_min_idx_2 - 1)] == 0))
    
    # Check boundary condition at w_min_idx
    # S should include recruitment flux
    # S[i, j_start] = n[i, j_start] + R[i] * dt / dw[j_start]
    
    # Species 1
    j_start_1 <- params@w_min_idx[1]
    expected_S_1 <- n[1, j_start_1] + recruitment_flux[1] * dt / params@dw[j_start_1]
    expect_equivalent(coefs$S[1, j_start_1], expected_S_1)
    
    # Species 2
    j_start_2 <- params@w_min_idx[2]
    expected_S_2 <- n[2, j_start_2] + recruitment_flux[2] * dt / params@dw[j_start_2]
    expect_equivalent(coefs$S[2, j_start_2], expected_S_2)
    
    # Check that 'a' at boundary is 0
    expect_equivalent(coefs$a[1, j_start_1], 0)
    expect_equivalent(coefs$a[2, j_start_2], 0)
    
    # Check 'b' at boundary
    # b_j_start = 1 + dt*mu + dt/dw * (g + D/(2*dw))
    # Note: we need to calculate expected values using the same inputs
    g <- getEGrowth(params)
    mu <- getMort(params)
    dw <- params@dw
    
    # Species 2
    expected_b_2 <- 1 + dt * mu[2, j_start_2] + 
        (dt / dw[j_start_2]) * (g[2, j_start_2] + 0.5 * params@diffusion[2, j_start_2] / dw[j_start_2])
        
    expect_equivalent(coefs$b[2, j_start_2], expected_b_2)
})
