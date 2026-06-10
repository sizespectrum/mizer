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
    # From the same fixed point both methods stay put and agree closely.
    rel_diff <- max(abs(finalN(sim_t) - finalN(sim_e)) /
                        (finalN(sim_e) + 1e-30))
    expect_lt(rel_diff, 1e-3)
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
