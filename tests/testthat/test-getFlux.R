test_that("getFlux works correctly", {
    params <- trait_params_2sp
    # Force different w_min to test zeroing logic
    params@species_params$w_min[2] <- 0.01 
    params@w_min_idx[2] <- which.min(abs(params@w - 0.01))
    
    n <- params@initial_n
    n[] <- 1 # Set n to 1 to make checking easier
    
    t <- 0
    g <- getEGrowth(params, n = n, t = t)
    d <- getDiffusion(params, n = n, t = t)
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

test_that("getFlux uses total diffusion", {
    params <- single_sp_params
    params@use_predation_diffusion <- TRUE
    species <- params@species_params$species[1]

    n <- params@initial_n
    t <- 0
    g <- getEGrowth(params, n = n, t = t)
    d <- getDiffusion(params, n = n, t = t)
    flux <- getFlux(params, n = n, t = t)

    j <- params@w_min_idx[species] + 5
    expected <- g[species, j - 1] * n[species, j - 1] -
        0.5 * (d[species, j] * n[species, j] -
                   d[species, j - 1] * n[species, j - 1]) /
        params@dw[j - 1]

    expect_equal(flux[species, j], expected, ignore_attr = TRUE)
})

test_that("getFlux power argument scales the flux by a power of weight", {
    params <- trait_params_2sp
    n <- params@initial_n

    flux0 <- getFlux(params, n = n)
    flux1 <- getFlux(params, n = n, power = 1)
    flux2 <- getFlux(params, n = n, power = 2)

    # power = 0 is the default
    expect_identical(flux0, getFlux(params, n = n, power = 0))

    expect_equal(flux1[], sweep(flux0[], 2, params@w, "*"),
                 ignore_attr = TRUE)
    expect_equal(flux2[], sweep(flux0[], 2, params@w^2, "*"),
                 ignore_attr = TRUE)

    # units are updated
    expect_equal(attr(flux0, "units"), "1/year")
    expect_equal(attr(flux1, "units"), "g/year")
    expect_equal(attr(flux2, "units"), "g^2/year")

    expect_error(getFlux(params, n = n, power = "a"))
})

test_that("getFlux power argument works for MizerSim", {
    params <- trait_params_2sp
    sim <- project(params, t_max = 1, progress_bar = FALSE)

    flux0 <- getFlux(sim)
    flux1 <- getFlux(sim, power = 1)

    expect_equal(flux1[], sweep(flux0[], 3, params@w, "*"),
                 ignore_attr = TRUE)
    expect_equal(attr(flux1, "units"), "g/year")
})
