test_that("project_n follows the documented one-step update", {
    params <- newMultispeciesParams(NS_species_params_gears[12, , drop = FALSE],
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

test_that("project_n with nonzero diffusion produces a different result", {
    params <- newMultispeciesParams(NS_species_params_gears[12, , drop = FALSE],
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
