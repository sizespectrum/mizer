test_that("getFlux works correctly", {
    params <- newTraitParams(no_sp = 2)
    # Force different w_min to test zeroing logic
    params@species_params$w_min[2] <- 0.01 
    params@w_min_idx[2] <- which.min(abs(params@w - 0.01))
    
    n <- params@initial_n
    n[] <- 1 # Set n to 1 to make checking easier
    
    t <- 0
    g <- getEGrowth(params, n = n, t = t)
    d <- params@diffusion
    dw <- params@dw
    rdd <- getRDD(params, n = n, t = t)
    
    flux <- getFlux(params, n = n, t = t)
    
    # Check dimensions
    expect_equal(dim(flux), dim(n))
    
    # Check zeroing out below w_min_idx
    w_min_idx_2 <- params@w_min_idx[2]
    expect_true(w_min_idx_2 > 1) 
    
    expect_true(all(flux[2, 1:(w_min_idx_2 - 1)] == 0))
    
    # Check boundary condition at w_min_idx
    # flux[i, j_start] = Rdd[i]
    
    # Species 1
    j_start_1 <- params@w_min_idx[1]
    expect_equal(flux[1, j_start_1], rdd[1], ignore_attr = TRUE)
    
    # Species 2
    j_start_2 <- params@w_min_idx[2]
    expect_equal(flux[2, j_start_2], rdd[2], ignore_attr = TRUE)
    
    # Check general calculation for some j > j_start
    # J_{i,j} = g_{i, j-1} N_{i, j-1} - 1/2 * (d_{i, j} N_{i, j} - d_{i, j-1} N_{i, j-1}) / dw_{j-1}
    j <- j_start_2 + 5
    expected_flux_2_j <- g[2, j - 1] * n[2, j - 1] - 0.5 * (d[2, j] * n[2, j] - d[2, j - 1] * n[2, j - 1]) / dw[j - 1]
    
    expect_equal(flux[2, j], expected_flux_2_j, ignore_attr = TRUE)
})
