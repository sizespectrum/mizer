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
