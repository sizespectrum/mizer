test_that("get_steady_state_n works with no diffusion", {
    # Simple case: constant growth, constant mortality
    # dN/dw = - (mu/g) N
    # N(w) = N0 * exp(- mu/g * (w - w0))
    # But we work with bins.
    
    # 3 bins
    growth <- c(1, 1, 1)
    mort <- c(0.5, 0.5, 0.5)
    dw <- c(1, 1, 1)
    # diffusion 0
    idx <- 1:2
    N0 <- 100
    
    # Analytical/recursive result from old implementation
    # n_exact <- c(1, cumprod(growth[idx] / ((growth + mort * dw)[idx + 1]))) * N0
    # idx=1: g[1] / (g[2] + mort[2]*dw[2]) = 1 / (1 + 0.5*1) = 1/1.5 = 2/3
    # idx=2: g[2] / (g[3] + mort[3]*dw[3]) = 1 / 1.5 = 2/3
    # n[1] = 100
    # n[2] = 100 * 2/3
    # n[3] = 100 * 2/3 * 2/3 = 100 * 4/9
    
    n_expected <- c(100, 100 * 2/3, 100 * 4/9)
    
    # Using new implementation with default diffusion (0)
    n_calc <- mizer:::get_steady_state_n(growth, mort, dw, idx = idx, N0 = N0)
    
    expect_equal(n_calc, n_expected)
})

test_that("get_steady_state_n works with diffusion", {
    # 3 bins
    # Diffusion dominates or mixes
    # If we have strong diffusion, distribution should flatten or change
    
    growth <- c(1, 1, 1)
    mort <- c(0.1, 0.1, 0.1)
    dw <- c(1, 1, 1)
    diffusion <- c(1, 1, 1) # Non-zero
    idx <- 1:2
    N0 <- 100
    
    n_calc <- mizer:::get_steady_state_n(growth, mort, dw, diffusion, idx, N0)
    
    # Check it runs and returns vector of correct length
    expect_length(n_calc, 3)
    expect_equal(n_calc[1], N0) # Boundary condition
    expect_true(all(n_calc > 0)) # Should be positive
    
    # Compare with no diffusion
    n_nodiff <- mizer:::get_steady_state_n(growth, mort, dw, diffusion = rep(0, 3), idx, N0)
    expect_false(isTRUE(all.equal(n_calc, n_nodiff)))
})

test_that("get_steady_state_n matches analytical steady state for advection-diffusion?", {
    # Hard to test exact analytical without solving ODE, but we can check consistency
    # or rely on the fact it solves the system we defined.
})
