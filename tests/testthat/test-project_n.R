test_that("project_n follows the documented one-step update", {
    params <- newMultispeciesParams(NS_species_params_gears[3, , drop = FALSE],
                                    info_level = 0)
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)
    idx <- 2:no_w
    w_min_idx_array_ref <- (params@w_min_idx - 1) * no_sp + seq_len(no_sp)

    n <- matrix(seq_len(no_w), nrow = no_sp,
                dimnames = dimnames(params@initial_n))
    e_growth <- matrix(seq(0.2, 0.2 * no_w, by = 0.2), nrow = no_sp,
                       dimnames = dimnames(params@initial_n))
    mort <- matrix(seq(0.05, 0.05 * no_w, by = 0.05), nrow = no_sp,
                   dimnames = dimnames(params@initial_n))
    diffusion <- matrix(0, nrow = no_sp, ncol = no_w,
                        dimnames = dimnames(params@initial_n))
    r <- list(e_growth = e_growth, mort = mort, rdd = 3, diffusion = diffusion)
    dt <- 0.1
    a <- matrix(0, nrow = no_sp, ncol = no_w)
    b <- matrix(0, nrow = no_sp, ncol = no_w)
    c <- matrix(0, nrow = no_sp, ncol = no_w)
    S <- matrix(0, nrow = no_sp, ncol = no_w)

    result <- project_n(params, r, n, dt, a, b, c, S, idx, w_min_idx_array_ref,
                        no_sp, no_w)

    a_expected <- a
    b_expected <- b
    S_expected <- S
    a_expected[, idx] <- sweep(-e_growth[, idx - 1, drop = FALSE] * dt, 2,
                               params@dw[idx], "/")
    b_expected[] <- 1 + sweep(e_growth * dt, 2, params@dw, "/") + mort * dt
    S_expected[, idx] <- n[, idx, drop = FALSE]
    expected <- n
    expected[w_min_idx_array_ref] <-
        (n[w_min_idx_array_ref] + r$rdd * dt /
             params@dw[params@w_min_idx]) /
        b_expected[w_min_idx_array_ref]
    for (j in (params@w_min_idx[1] + 1):no_w) {
        expected[1, j] <-
            (S_expected[1, j] - a_expected[1, j] * expected[1, j - 1]) /
            b_expected[1, j]
    }

    expect_equal(result, expected)
})

test_that("project_n_2 follows Crank-Nicolson update for fixed rates", {
    params <- newMultispeciesParams(NS_species_params_gears[3, , drop = FALSE],
                                    info_level = 0)
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)
    idx <- 2:no_w
    w_min_idx_array_ref <- (params@w_min_idx - 1) * no_sp + seq_len(no_sp)

    n <- matrix(seq_len(no_w), nrow = no_sp,
                dimnames = dimnames(params@initial_n))
    e_growth <- matrix(seq(0.2, 0.2 * no_w, by = 0.2), nrow = no_sp,
                       dimnames = dimnames(params@initial_n))
    mort <- matrix(seq(0.05, 0.05 * no_w, by = 0.05), nrow = no_sp,
                   dimnames = dimnames(params@initial_n))
    diffusion <- matrix(0, nrow = no_sp, ncol = no_w,
                        dimnames = dimnames(params@initial_n))
    r <- list(e_growth = e_growth, mort = mort, rdd = 3, diffusion = diffusion)
    dt <- 0.1
    a <- matrix(0, nrow = no_sp, ncol = no_w)
    b <- matrix(0, nrow = no_sp, ncol = no_w)
    c <- matrix(0, nrow = no_sp, ncol = no_w)
    S <- matrix(0, nrow = no_sp, ncol = no_w)

    result <- project_n_2(params, r, n, dt, a, b, c, S, idx,
                          w_min_idx_array_ref, no_sp, no_w)

    coefs <- mizer:::get_transport_coefs(params, n, e_growth, mort, dt / 2,
                                         recruitment_flux = r$rdd,
                                         d = diffusion)
    rhs <- 2 * n
    matrix_n <- coefs$b * n
    matrix_n[, idx] <- matrix_n[, idx] +
        coefs$a[, idx] * n[, idx - 1, drop = FALSE]
    matrix_n[, -no_w] <- matrix_n[, -no_w] +
        coefs$c[, -no_w] * n[, -1, drop = FALSE]
    rhs <- rhs - matrix_n
    rhs[w_min_idx_array_ref] <- rhs[w_min_idx_array_ref] +
        dt * r$rdd / params@dw[params@w_min_idx]

    expected <- mizer:::project_n_loop(n, coefs$a, coefs$b, coefs$c, rhs,
                                       params@w_min_idx)

    expect_equal(result, expected)
})

test_that("project_n_2 uses midpoint rates from provisional rates", {
    params <- newMultispeciesParams(NS_species_params_gears[3, , drop = FALSE],
                                    info_level = 0)
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)
    idx <- 2:no_w
    w_min_idx_array_ref <- (params@w_min_idx - 1) * no_sp + seq_len(no_sp)

    n <- matrix(seq_len(no_w), nrow = no_sp,
                dimnames = dimnames(params@initial_n))
    e_growth <- matrix(seq(0.2, 0.2 * no_w, by = 0.2), nrow = no_sp,
                       dimnames = dimnames(params@initial_n))
    mort <- matrix(seq(0.05, 0.05 * no_w, by = 0.05), nrow = no_sp,
                   dimnames = dimnames(params@initial_n))
    diffusion <- matrix(0.01 * seq_len(no_w), nrow = no_sp, ncol = no_w,
                        dimnames = dimnames(params@initial_n))
    r <- list(e_growth = e_growth, mort = mort, rdd = 3, diffusion = diffusion)
    r_hat <- list(e_growth = 2 * e_growth, mort = 2 * mort, rdd = 5,
                  diffusion = 2 * diffusion)
    r_mid <- list(e_growth = 0.5 * (r$e_growth + r_hat$e_growth),
                  mort = 0.5 * (r$mort + r_hat$mort),
                  rdd = 0.5 * (r$rdd + r_hat$rdd),
                  diffusion = 0.5 * (r$diffusion + r_hat$diffusion))
    dt <- 0.1
    a <- b <- c <- S <- matrix(0, nrow = no_sp, ncol = no_w)

    result_hat <- project_n_2(params, r, n, dt, a, b, c, S, idx,
                              w_min_idx_array_ref, no_sp, no_w,
                              r_hat = r_hat)
    result_mid <- project_n_2(params, r, n, dt, a, b, c, S, idx,
                              w_min_idx_array_ref, no_sp, no_w,
                              r_mid = r_mid)
    predicted_n <- project_n(params, r, n, dt, a, b, c, S, idx,
                             w_min_idx_array_ref, no_sp, no_w)
    seen_n <- NULL
    rates_fns <- list(Rates = function(params, n, n_pp, n_other, t, effort,
                                       rates_fns, ...) {
        seen_n <<- n
        r_hat
    })
    result_rates_fns <- project_n_2(params, r, n, dt, a, b, c, S, idx,
                                    w_min_idx_array_ref, no_sp, no_w,
                                    rates_fns = rates_fns,
                                    n_pp = params@initial_n_pp,
                                    n_other = list(), t = 2, effort = 1)

    expect_equal(result_hat, result_mid)
    expect_equal(result_rates_fns, result_mid)
    expect_equal(seen_n, predicted_n)
})

test_that("project_n with nonzero diffusion produces a different result", {
    params <- newMultispeciesParams(NS_species_params_gears[3, , drop = FALSE],
                                    info_level = 0)
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)
    idx <- 2:no_w
    w_min_idx_array_ref <- (params@w_min_idx - 1) * no_sp + seq_len(no_sp)

    n <- matrix(seq_len(no_w), nrow = no_sp,
                dimnames = dimnames(params@initial_n))
    e_growth <- matrix(seq(0.2, 0.2 * no_w, by = 0.2), nrow = no_sp,
                       dimnames = dimnames(params@initial_n))
    mort <- matrix(seq(0.05, 0.05 * no_w, by = 0.05), nrow = no_sp,
                   dimnames = dimnames(params@initial_n))
    dt <- 0.1
    a <- b <- c <- S <- matrix(0, nrow = no_sp, ncol = no_w)

    r_zero <- list(e_growth = e_growth, mort = mort, rdd = 3,
                   diffusion = matrix(0, nrow = no_sp, ncol = no_w,
                                      dimnames = dimnames(params@initial_n)))
    r_diff <- list(e_growth = e_growth, mort = mort, rdd = 3,
                   diffusion = matrix(0.01 * seq_len(no_w), nrow = no_sp, ncol = no_w,
                                      dimnames = dimnames(params@initial_n)))

    result_zero <- project_n(params, r_zero, n, dt, a, b, c, S, idx,
                             w_min_idx_array_ref, no_sp, no_w)
    result_diff <- project_n(params, r_diff, n, dt, a, b, c, S, idx,
                             w_min_idx_array_ref, no_sp, no_w)

    expect_false(identical(result_zero, result_diff))
    expect_snapshot_value(result_diff, style = 'json2', tolerance = 1e-5)
})

test_that("project_n with predation diffusion stays non-negative and finite", {
    params_d <- NS_params
    params_d@use_predation_diffusion <- TRUE
    r <- getRates(params_d)
    # confirm predation diffusion is active
    expect_true(any(r$diffusion > 0))

    no_sp <- nrow(params_d@species_params)
    no_w <- length(params_d@w)
    idx <- 2:no_w
    w_min_idx_array_ref <- (params_d@w_min_idx - 1) * no_sp + seq_len(no_sp)
    a <- b <- c <- S <- matrix(0, nrow = no_sp, ncol = no_w)

    n_new <- project_n(params_d, r, params_d@initial_n, dt = 0.1,
                       a, b, c, S, idx, w_min_idx_array_ref, no_sp, no_w)
    expect_true(all(is.finite(n_new)))
    expect_true(all(n_new >= 0))
})
