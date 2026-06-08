test_that("get_steady_state_n works with no diffusion", {
    params <- NS_params_small
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)

    # Mocking constant growth and mortality
    growth <- matrix(1, nrow = no_sp, ncol = no_w)
    mort <- matrix(0.5, nrow = no_sp, ncol = no_w)
    N0_vec <- rep(100, no_sp)

    # Zero diffusion
    params@ext_diffusion[] <- 0

    n_calc <- mizer:::get_steady_state_n(params, growth, mort,
                                         params@ext_diffusion, N0_vec)

    expect_equal(dim(n_calc), c(no_sp, no_w))

    # Check species 1 manually against old analytical form
    sp <- 1
    w_min_idx <- params@w_min_idx[sp]
    w_max_idx <- sum(params@w <= params@species_params[sp, "w_max"])
    idx <- w_min_idx:(w_max_idx - 1)
    dw <- params@dw

    # Old calculation logic (no diffusion)
    n_old <- c(1, cumprod(growth[sp, idx] / (growth[sp, idx + 1] + mort[sp, idx + 1] * dw[idx + 1])))
    n_old <- 100 * n_old

    expect_equal(unname(n_calc[sp, w_min_idx:w_max_idx]), unname(n_old), tolerance = 1e-10)


})

test_that("get_steady_state_n works with diffusion", {
    params <- example_params()
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)

    # Mocking constant growth and mortality
    growth <- matrix(1, nrow = no_sp, ncol = no_w)
    mort <- matrix(0.1, nrow = no_sp, ncol = no_w)
    N0_vec <- rep(100, no_sp)

    D <- getDiffusion(params)
    n_calc <- mizer:::get_steady_state_n(params, growth, mort, D, N0_vec)

    expect_equal(dim(n_calc), c(no_sp, no_w))

    # Should be positive up to w_max
    for(sp in 1:no_sp) {
        w_min_idx <- params@w_min_idx[sp]
        w_max_idx <- sum(params@w <= params@species_params[sp, "w_max"])
        expect_equal(unname(n_calc[sp, w_min_idx]), unname(N0_vec[sp]))
    }

    # Compare with no diffusion
    zero_D <- D
    zero_D[] <- 0
    n_nodiff <- mizer:::get_steady_state_n(params, growth, mort, zero_D,
                                           N0_vec)
    expect_false(isTRUE(all.equal(n_calc, n_nodiff)))
})

test_that("get_steady_state_n uses supplied diffusion array", {
    params <- example_params()
    params@ext_diffusion[] <- 0
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)

    growth <- matrix(1, nrow = no_sp, ncol = no_w)
    mort <- matrix(0.1, nrow = no_sp, ncol = no_w)
    N0_vec <- rep(100, no_sp)
    D <- params@ext_diffusion
    D[1, ] <- 0.1 * params@w^(params@species_params$n[1] + 1)

    n_with_D <- mizer:::get_steady_state_n(params, growth, mort, D, N0_vec)
    D[] <- 0
    n_without_D <- mizer:::get_steady_state_n(params, growth, mort, D, N0_vec)

    expect_false(isTRUE(all.equal(n_with_D, n_without_D)))
})
