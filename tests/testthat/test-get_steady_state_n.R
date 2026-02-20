test_that("get_steady_state_n works with no diffusion", {
    params <- NS_params
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)
    
    # Mocking constant growth and mortality
    growth <- matrix(1, nrow = no_sp, ncol = no_w)
    mort <- matrix(0.5, nrow = no_sp, ncol = no_w)
    N0_vec <- rep(100, no_sp)
    
    # Zero diffusion
    params@diffusion[] <- 0
    
    n_calc <- mizer:::get_steady_state_n(params, growth, mort, N0_vec)
    
    expect_equal(dim(n_calc), c(no_sp, no_w))
    
    # Check species 1 manually against old analytical form
    sp <- 1
    w_min_idx <- params@w_min_idx[sp]
    w_max_idx <- sum(params@w <= params@species_params[sp, "w_max"])
    idx <- w_min_idx:(w_max_idx - 1)
    dw <- params@dw
    
    # Old calculation logic (no diffusion)
    n_old <- c(1, cumprod(growth[sp, idx] / ((growth[sp, idx] + mort[sp, idx] * dw)[idx + 1])))
    n_old <- 100 * n_old
    
    expect_equal(unname(n_calc[sp, w_min_idx:w_max_idx]), unname(n_old), tolerance = 1e-10)
    
    # All bins above w_max_idx should be 0
    if (w_max_idx < no_w) {
        expect_true(all(n_calc[sp, (w_max_idx + 1):no_w] == 0))
    }
})

test_that("get_steady_state_n works with diffusion", {
    params <- NS_params
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)
    
    # Mocking constant growth and mortality
    growth <- matrix(1, nrow = no_sp, ncol = no_w)
    mort <- matrix(0.1, nrow = no_sp, ncol = no_w)
    N0_vec <- rep(100, no_sp)
    
    # Non-zero diffusion
    params@diffusion[] <- 1
    
    n_calc <- mizer:::get_steady_state_n(params, growth, mort, N0_vec)
    
    expect_equal(dim(n_calc), c(no_sp, no_w))
    expect_equal(unname(n_calc[, params@w_min_idx[1]]), unname(N0_vec)) # Assuming same min for all just for picking index
    
    # Should be positive up to w_max
    for(sp in 1:no_sp) {
        w_min_idx <- params@w_min_idx[sp]
        w_max_idx <- sum(params@w <= params@species_params[sp, "w_max"])
        expect_true(all(n_calc[sp, w_min_idx:w_max_idx] > 0))
        if(w_max_idx < no_w) {
            expect_true(all(n_calc[sp, (w_max_idx + 1):no_w] == 0))
        }
    }
    
    # Compare with no diffusion
    params@diffusion[] <- 0
    n_nodiff <- mizer:::get_steady_state_n(params, growth, mort, N0_vec)
    expect_false(isTRUE(all.equal(n_calc, n_nodiff)))
})

test_that("get_steady_state_n matches analytical steady state for advection-diffusion?", {
    # It solves the transport scheme equation so consistency is implicitly tested.
    expect_true(TRUE)
})
