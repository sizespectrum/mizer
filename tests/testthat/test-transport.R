test_that("get_transport_coefs works correctly", {
    params <- trait_params_2sp
    # Force different w_min to test zeroing logic
    params@species_params$w_min[2] <- 0.01 
    params@w_min_idx[2] <- which.min(abs(params@w - 0.01))
    
    n <- params@initial_n
    n[] <- 1 # Set n to 1 to make S checking easier
    
    dt <- 0.1
    recruitment_flux <- c(10, 20)
    
    # We need to access the internal function
    get_transport_coefs <- mizer:::get_transport_coefs
    getEGrowth <- mizer:::getEGrowth
    getMort <- mizer:::getMort
    
    d <- params@ext_diffusion
    coefs <- get_transport_coefs(params, n, getEGrowth(params), getMort(params), dt, recruitment_flux, d)
    
    # Check dimensions
    expect_equal(dim(coefs$a), dim(n))
    expect_equal(dim(coefs$b), dim(n))
    expect_equal(dim(coefs$c), dim(n))
    expect_equal(dim(coefs$S), dim(n))
    
    # Check zeroing out below w_min_idx
    w_min_idx_2 <- params@w_min_idx[2]
    expect_true(w_min_idx_2 > 1) # Ensure we are actually testing something interesting
    
    expect_true(all(coefs$a[2, 1:(w_min_idx_2 - 1)] == 0))
    expect_true(all(coefs$b[2, 1:(w_min_idx_2 - 1)] == 0))
    expect_true(all(coefs$c[2, 1:(w_min_idx_2 - 1)] == 0))
    expect_true(all(coefs$S[2, 1:(w_min_idx_2 - 1)] == 0))
    
    # Check boundary condition at w_min_idx
    # S should include recruitment flux
    # S[i, j_start] = n[i, j_start] + R[i] * dt / dw[j_start]
    
    # Species 1
    j_start_1 <- params@w_min_idx[1]
    expected_S_1 <- n[1, j_start_1] + recruitment_flux[1] * dt / params@dw[j_start_1]
    expect_equal(coefs$S[1, j_start_1], expected_S_1, ignore_attr = TRUE)
    
    # Species 2
    j_start_2 <- params@w_min_idx[2]
    expected_S_2 <- n[2, j_start_2] + recruitment_flux[2] * dt / params@dw[j_start_2]
    expect_equal(coefs$S[2, j_start_2], expected_S_2, ignore_attr = TRUE)
    
    # Check that 'a' at boundary is 0
    expect_equal(coefs$a[1, j_start_1], 0, ignore_attr = TRUE)
    expect_equal(coefs$a[2, j_start_2], 0, ignore_attr = TRUE)
    
    # Check 'b' at boundary
    # b_j_start = 1 + dt*mu + dt/dw * (g + D/(2*dw))
    # Note: we need to calculate expected values using the same inputs
    g <- getEGrowth(params)
    mu <- getMort(params)
    dw <- params@dw
    
    # Species 2
    expected_b_2 <- 1 + dt * mu[2, j_start_2] +
        (dt / dw[j_start_2]) * (g[2, j_start_2] + 0.5 * d[2, j_start_2] / dw[j_start_2])
        
    expect_equal(coefs$b[2, j_start_2], expected_b_2, ignore_attr = TRUE)
})

# flux_limiter_scheme ----

# Tests for the flux-limited advective flux, controlled by the `flux_limiter`
# entry of the `second_order_w` slot. The limiter folds a frozen flux-limited
# (TVD) high-order advective flux into the implicit transport operator. See the
# numerical-details vignette.

test_that("flux_limiter_scheme reads the second_order_w slot", {
    expect_identical(flux_limiter_scheme(NS_params_small), "none")
    p <- NS_params_small
    second_order_w(p) <- c(flux = TRUE)
    expect_identical(flux_limiter_scheme(p), "van_leer")
})

test_that("the flux limiter changes the projection when switched on", {
    sim_none <- project(NS_params_small, t_max = 1, dt = 0.1,
                        progress_bar = FALSE)
    p_vl <- NS_params_small
    second_order_w(p_vl) <- c(flux = TRUE)
    sim_vl <- project(p_vl, t_max = 1, dt = 0.1, progress_bar = FALSE)
    # With the limiter off the result is the plain upwind projection; switching
    # it on must change the spectrum.
    expect_false(isTRUE(all.equal(sim_vl@n, sim_none@n)))
})

test_that("the flux limiter scheme is recorded and carried over on append", {
    p_vl <- NS_params_small
    second_order_w(p_vl) <- c(flux = TRUE)
    sim <- project(p_vl, t_max = 1, progress_bar = FALSE)
    expect_identical(sim@sim_params$flux_limiter, "van_leer")
    # The scheme is a property of the params, so appending carries it over.
    sim2 <- project(sim, t_max = 1, progress_bar = FALSE)
    expect_identical(sim2@sim_params$flux_limiter, "van_leer")
})

test_that("flux limiter keeps abundances non-negative", {
    p_vl <- NS_params_small
    second_order_w(p_vl) <- c(flux = TRUE)
    for (m in c("euler", "predictor_corrector", "tr_bdf2")) {
        sim <- project(p_vl, t_max = 2, dt = 0.1, method = m,
                       progress_bar = FALSE)
        expect_true(all(sim@n >= 0))
    }
})

test_that("flux limiter steady state is preserved by all methods", {
    # Build the limiter's own steady state and hold reproduction constant so the
    # only thing being tested is whether the projection scheme keeps that state
    # steady. This mirrors the predictor-corrector test in
    # test-getRequiredRDD.R, but with the flux limiter switched on everywhere.
    params <- single_sp_params
    second_order_w(params) <- c(flux = TRUE)
    species <- params@species_params$species[1]
    n <- params@species_params[species, "n"]
    ext_diffusion(params)[species, ] <- 0.1 * params@w^(n + 1)

    params <- steadySingleSpecies(params, species = species)
    params@species_params$constant_reproduction <- getRequiredRDD(params)
    params <- setRateFunction(params, "RDD", "constantRDD")

    idx <- params@w_min_idx[species]
    n0 <- initialN(params)[species, ]
    for (m in c("euler", "predictor_corrector", "tr_bdf2")) {
        sim <- project(params, t_max = 1, t_save = 1, method = m,
                       progress_bar = FALSE)
        nf <- finalN(sim)[species, ]
        # Egg density and the whole spectrum stay put.
        expect_equal(nf[idx], n0[idx], tolerance = 1e-8)
        expect_lt(max(abs(nf - n0)) / max(n0), 1e-6)
    }
})

test_that("the second-order scheme is second order in the size step", {
    skip_on_cran()
    # Power-law transport problem with the exact Bessel-kernel solution from the
    # analytic-test vignette. The default upwind scheme is first order in the
    # grid spacing; the second-order scheme (advective flux + bin-averaged
    # sinks) is second order. The finite-volume cells are the bins, so the
    # second-order field N_j is the cell average, compared at the cell centre
    # sqrt(w_j w_{j+1}); the first-order field is compared at the node.
    p <- 0.7; A <- 1; B <- 0.5; K <- 0.1
    assign("so_growth", function(params, ...) matrix(A * params@w^p, nrow = 1),
           envir = globalenv())
    on.exit(rm(list = "so_growth", envir = globalenv()), add = TRUE)
    N_analytic <- function(w, t, w0) {
        U <- A - 0.5 * K; V <- 0.5 * K * (1 - p); b <- B / (1 - p)
        nu <- sqrt((U / V)^2 + 4 * b / V)
        x <- w^(1 - p) / (1 - p); x0 <- w0^(1 - p) / (1 - p)
        z <- 2 * sqrt(x * x0) / (V * t)
        exp(-log(V * t) + (U / (2 * V)) * log(x / x0) - (x + x0) / (V * t) +
                z + log(besselI(z, nu, expon.scaled = TRUE))) * w^(-p)
    }
    w0 <- 1e-2; t_start <- 0.1; t_end <- 1.5

    err <- function(no_w, second_order) {
        pr <- newMultispeciesParams(
            data.frame(species = "Test", w_max = 1000, w_mat = 100,
                       n = p, z0 = 0, z_ext = B, d = p - 1, D_ext = K),
            no_w = no_w, min_w = 1e-3, info_level = 0)
        second_order_w(pr) <- second_order
        pr <- setRateFunction(pr, "EGrowth", "so_growth")
        pr@interaction[] <- 0
        pr <- setExtMort(pr); pr <- setExtDiffusion(pr)
        pr@species_params$constant_reproduction <- 0
        pr <- setRateFunction(pr, "RDD", "constantRDD")
        # Cell average (= cell-centre value) for the second-order scheme, node
        # value for the first-order scheme.
        wv <- w(pr)
        refw <- if (second_order) wv * sqrt(wv[2] / wv[1]) else wv
        initialN(pr) <- matrix(N_analytic(refw, t_start, w0), nrow = 1)
        sim <- project(pr, t_max = t_end - t_start, dt = 0.005,
                       t_save = t_end - t_start, method = "tr_bdf2",
                       progress_bar = FALSE)
        num <- finalN(sim)[1, ]; ana <- N_analytic(refw, t_end, w0)
        mask <- ana > max(ana) * 1e-3
        sqrt(sum((num[mask] - ana[mask])^2 * pr@dw[mask]) /
             sum(ana[mask]^2 * pr@dw[mask]))
    }

    no_ws <- c(200, 400, 800)
    dx <- log(1000 / 1e-3) / no_ws
    slope <- function(e) coef(lm(log(e) ~ log(dx)))[2]
    ord_upwind <- slope(sapply(no_ws, err, second_order = FALSE))
    ord_second <- slope(sapply(no_ws, err, second_order = TRUE))
    expect_lt(ord_upwind, 1.2)            # first order
    expect_gt(ord_second, 1.8)            # second order
})

test_that("steadySingleSpecies honours the slot", {
    # Switching the limiter on via the slot changes the interior steady state.
    p_vl <- single_sp_params
    second_order_w(p_vl) <- c(flux = TRUE)

    a <- steadySingleSpecies(single_sp_params)
    b <- steadySingleSpecies(p_vl)
    expect_false(isTRUE(all.equal(a@initial_n, b@initial_n)))
    # getRequiredRDD reads the same slot. The second-order scheme discretises the
    # transport step (including the egg-cell balance) on the log-size grid, so the
    # required reproduction differs from the first-order upwind value. Its
    # consistency with the limited projection is covered by the steady-state
    # preservation test above.
    expect_false(isTRUE(all.equal(getRequiredRDD(single_sp_params),
                                  getRequiredRDD(p_vl))))
})
