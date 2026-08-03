test_that("project_n follows the documented one-step update", {
    params <- NS_params_cod_small
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
    params <- NS_params_cod_small
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
    params <- NS_params_cod_small
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
    params <- NS_params_cod_small
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
    params_d <- NS_params_small
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

# project_n_tr_bdf2 ----

# Tests for the L-stable TR-BDF2 time stepper (method = "tr_bdf2")

test_that("tr_bdf2 runs and keeps densities finite and non-negative", {
    sim <- project(NS_params_small, t_max = 2, dt = 0.1,
                   progress_bar = FALSE, method = "tr_bdf2")
    expect_s4_class(sim, "MizerSim")
    expect_true(all(is.finite(sim@n)))
    expect_true(all(sim@n >= 0))
    expect_equal(getSimParams(sim)$method, "tr_bdf2")
})

test_that("tr_bdf2 preserves a steady state", {
    ps <- suppressMessages(steady(NS_params_small, progress_bar = FALSE))
    sim_t <- project(ps, t_max = 5, dt = 0.1, method = "tr_bdf2",
                     progress_bar = FALSE)
    sim_e <- project(ps, t_max = 5, dt = 0.1, method = "euler",
                     progress_bar = FALSE)
    # From the same fixed point both methods stay put and agree closely, for
    # both the consumer spectrum and the resource.
    rel_diff <- max(abs(finalN(sim_t) - finalN(sim_e)) /
                        (finalN(sim_e) + 1e-30))
    expect_lt(rel_diff, 1e-3)
    rel_diff_pp <- max(abs(finalNResource(sim_t) - finalNResource(sim_e)) /
                           (finalNResource(sim_e) + 1e-30))
    expect_lt(rel_diff_pp, 1e-3)
})

test_that("second-order methods advance the resource with midpoint mortality", {
    # The resource update should use the midpoint resource mortality
    # (r + r_hat)/2 rather than the start-of-step value. We replicate one step
    # by hand and check the resource matches the midpoint update exactly.
    p <- NS_params_small
    skip_if_not(length(p@other_dynamics) == 0)
    # Make the resource genuinely dynamic over the step so that the midpoint and
    # start-of-step updates differ: slow the resource and perturb the state.
    p@rr_pp[] <- p@rr_pp * 0.3
    initialN(p)[] <- initialN(p) * 3
    n0 <- initialN(p)
    npp0 <- initialNResource(p)
    nother0 <- initialNOther(p)
    dt <- 0.4

    rates_fns <- projectRateFunctions(p)
    res_fn <- get(p@resource_dynamics)
    no_sp <- nrow(p@species_params)
    no_w <- length(p@w)
    idx <- 2:no_w
    wref <- (p@w_min_idx - 1) * no_sp + (1:no_sp)
    zero <- matrix(0, no_sp, no_w)

    r <- rates_fns$Rates(p, n = n0, n_pp = npp0, n_other = nother0,
                         t = 0, effort = p@initial_effort, rates_fns = rates_fns)
    npp_hat <- res_fn(p, n = n0, n_pp = npp0, n_other = nother0, rates = r,
                      t = 0, dt = dt, resource_rate = p@rr_pp,
                      resource_capacity = p@cc_pp)
    n_hat <- project_n(p, r, n0, dt, zero, zero, zero, zero,
                       idx, wref, no_sp, no_w)
    r_hat <- rates_fns$Rates(p, n = n_hat, n_pp = npp_hat, n_other = nother0,
                             t = dt, effort = p@initial_effort,
                             rates_fns = rates_fns)
    r_mid <- average_rates(r, r_hat)

    expected_npp <- res_fn(p, n = n0, n_pp = npp0, n_other = nother0,
                           rates = r_mid, t = 0, dt = dt,
                           resource_rate = p@rr_pp,
                           resource_capacity = p@cc_pp)
    startrate_npp <- res_fn(p, n = n0, n_pp = npp0, n_other = nother0,
                            rates = r, t = 0, dt = dt,
                            resource_rate = p@rr_pp,
                            resource_capacity = p@cc_pp)
    # The midpoint correction genuinely changes the resource here (per-bin
    # relative difference, over the bins where the resource is present).
    present <- npp0 > 0
    rel_change <- max(abs(expected_npp[present] - startrate_npp[present]) /
                          abs(startrate_npp[present]))
    expect_gt(rel_change, 1e-4)

    for (m in c("tr_bdf2", "predictor_corrector")) {
        sim <- project(p, t_max = dt, dt = dt, t_save = dt, method = m,
                       progress_bar = FALSE)
        expect_equal(as.numeric(finalNResource(sim)),
                     as.numeric(expected_npp))
    }
})

test_that("tr_bdf2 is second order in time and beats euler", {
    relerr <- function(x, ref) {
        xf <- finalN(x)
        rf <- finalN(ref)
        sqrt(sum((xf - rf)^2)) / sqrt(sum(rf^2))
    }
    p <- NS_params_small
    # Perturb the initial state so the run is dominated by time-stepping error.
    initialN(p)[] <- initialN(p) * 1.5
    t_max <- 1
    reference <- project(p, dt = 0.2 / 2^6, t_max = t_max,
                         method = "tr_bdf2", progress_bar = FALSE)

    dt_values <- c(0.2, 0.1, 0.05)
    err_tr <- numeric(length(dt_values))
    err_eu <- numeric(length(dt_values))
    for (i in seq_along(dt_values)) {
        dt <- dt_values[i]
        tr <- project(p, dt = dt, t_max = t_max, t_save = t_max,
                      method = "tr_bdf2", progress_bar = FALSE)
        eu <- project(p, dt = dt, t_max = t_max, t_save = t_max,
                      method = "euler", progress_bar = FALSE)
        err_tr[i] <- relerr(tr, reference)
        err_eu[i] <- relerr(eu, reference)
    }

    # TR-BDF2 is more accurate than Euler at every step size.
    expect_true(all(err_tr < err_eu))
    # Halving dt reduces the TR-BDF2 error by clearly more than the factor ~2
    # expected from a first-order method (second order gives ~4).
    ratios <- err_tr[-length(err_tr)] / err_tr[-1]
    expect_true(all(ratios > 3))
})

test_that("tr_bdf2 damps stiff modes that make Crank-Nicolson oscillate", {
    # Frozen-rate pure-decay problem: no growth, no diffusion, no recruitment,
    # large mortality and a large time step. Crank-Nicolson rings (amplification
    # factor near -1); the L-stable TR-BDF2 decays quickly.
    p <- NS_params_small
    no_sp <- nrow(p@species_params)
    no_w <- length(p@w)
    zero <- matrix(0, no_sp, no_w)
    r <- list(e_growth = zero, diffusion = zero,
              mort = matrix(50, no_sp, no_w), rdd = rep(0, no_sp))

    idx <- 2:no_w
    wref <- (p@w_min_idx - 1) * no_sp + (1:no_sp)
    dt <- 0.5

    trajectory <- function(fn) {
        nn <- p@initial_n
        out <- numeric(6)
        out[1] <- nn[3, 10]
        for (k in 2:6) {
            nn <- fn(p, r, nn, dt, zero, zero, zero, zero,
                     idx, wref, no_sp, no_w)
            out[k] <- nn[3, 10]
        }
        out
    }
    tr <- trajectory(project_n_tr_bdf2)
    pc <- trajectory(project_n_2)

    # TR-BDF2 amplitude decays monotonically and is small after a few steps.
    expect_true(all(diff(abs(tr)) < 0))
    expect_lt(abs(tr[6]) / abs(tr[1]), 0.01)
    # Crank-Nicolson barely decays: amplitude stays large (rings).
    expect_gt(min(abs(pc)) / abs(pc[1]), 0.3)
})
