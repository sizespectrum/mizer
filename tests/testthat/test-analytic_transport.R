library(mizer)

# Parameters
p <- 0.7
A <- 1
B <- 0.5
K <- 0.1

# Helper functions defined at top level for visibility
# We assign to global environment so that mizer can find them with get()
assign("start_growth", function(params, ...) {
    matrix(A * params@w^p, nrow = 1, byrow = TRUE)
}, envir = globalenv())

assign("start_mort", function(params, ...) {
    matrix(B * params@w^(p - 1), nrow = 1, byrow = TRUE)
}, envir = globalenv())

assign("N_analytic", function(w, t, w0, t0, K_eff = K) {
    U <- A - 0.5 * K_eff
    V <- 0.5 * K_eff * (1 - p)
    b <- B / (1 - p)
    nu <- sqrt((U/V)^2 + 4 * b / V)
    dt <- t - t0
    x <- w^(1 - p) / (1 - p)
    x0 <- w0^(1 - p) / (1 - p)
    z <- 2 * sqrt(x * x0) / (V * dt)
    bessel_scaled <- besselI(z, nu, expon.scaled = TRUE)
    log_N_tilde <- -log(V * dt) + (U / (2 * V)) * log(x / x0) -
        (x + x0) / (V * dt) + z + log(bessel_scaled)
    N_tilde <- exp(log_N_tilde)
    N <- N_tilde * w^(-p)
    return(N)
}, envir = globalenv())

p_euler <- 1
A_euler <- 1
B_euler <- 0.5
K_euler <- 0.1
dt_euler <- 0.01

assign("start_growth_euler", function(params, ...) {
    matrix(A_euler * params@w^p_euler, nrow = 1, byrow = TRUE)
}, envir = globalenv())

assign("start_mort_euler", function(params, ...) {
    matrix(B_euler * params@w^(p_euler - 1), nrow = 1, byrow = TRUE)
}, envir = globalenv())

assign("N_lognormal", function(w, t, w0, t0, K_eff) {
    dt <- t - t0
    x <- log(w)
    mean_x <- log(w0) + (A_euler - 0.5 * K_eff) * dt
    
    exp(-B_euler * dt) *
        dnorm(x, mean = mean_x, sd = sqrt(K_eff * dt)) / w
}, envir = globalenv())

# 1. Steady State Test
test_that("Exact steady state is maintained", {
    # Calculate lambda for initial condition
    a_quad <- K
    b_quad <- -(2 * A - K)
    c_quad <- -2 * B
    det <- b_quad^2 - 4 * a_quad * c_quad
    x <- (-b_quad - sqrt(det)) / (2 * a_quad)
    lambda <- p - x
    
    # Setup Params
    params <- newMultispeciesParams(data.frame(species = "Test",
                                               w_max = 1000,
                                               w_mat = 100),
                                    no_w = 1000, min_w = 1e-3,
                                    info_level = 0)
    initialN(params) <- matrix(params@w^(-lambda), nrow = 1, byrow = TRUE)
    
    params <- setRateFunction(params, "EGrowth", "start_growth")
    params <- setRateFunction(params, "Mort", "start_mort")
    params@ext_diffusion[1, ] <- K * params@w^(p + 1)
    # RDD (Constant Flux)
    params@species_params$constant_reproduction <- getRequiredRDD(params)
    params <- setRateFunction(params, "RDD", "constantRDD")
    
    
    # Run
    sim <- project(params, t_max = 1, dt = 0.1, method = "predictor-corrector")
    
    # Compare
    n0 <- initialN(params)[1, ]
    n1 <- finalN(sim)[1, ]
    # Exclude boundaries
    idx <- 1:980
    rel_err <- abs(n1[idx] - n0[idx]) / n0[idx]
    
    expect_lt(max(rel_err), 0.01)
})

# 2. Time Dependent Test
test_that("Exact time-dependent solution is followed", {
    # Setup Params
    params <- newMultispeciesParams(data.frame(species = "Test",
                                               w_max = 1000,
                                               w_mat = 100),
                                    no_w = 1000, min_w = 1e-3,
                                    info_level = 0)
    
    params <- setRateFunction(params, "EGrowth", "start_growth")
    params <- setRateFunction(params, "Mort", "start_mort")
    params@ext_diffusion[1, ] <- K * params@w^(p + 1)
    # Set RDD to 0
    params@species_params$constant_reproduction <- 0
    params <- setRateFunction(params, "RDD", "constantRDD")
    
    w0 <- 1e-2
    t_start <- 0.1
    t_end <- 2
    
    initial_n <- N_analytic(params@w, t_start, w0, 0)
    initialN(params) <- matrix(initial_n, nrow = 1, byrow = TRUE)
    
    sim <- project(params, t_max = t_end - t_start, dt = 0.05,
                   t_save = t_end - t_start,
                   method = "predictor-corrector")
    
    final_n_num <- finalN(sim)[1, ]
    final_n_ana <- N_analytic(params@w, t_end, w0, 0)
    
    # Metrics
    total_n_num <- sum(final_n_num * params@dw)
    total_n_ana <- sum(final_n_ana * params@dw)
    rel_err_total <- abs(total_n_num - total_n_ana) / total_n_ana
    
    peak_idx_num <- which.max(final_n_num)
    peak_idx_ana <- which.max(final_n_ana)
    peak_w_num <- params@w[peak_idx_num]
    peak_w_ana <- params@w[peak_idx_ana]
    rel_err_peak_loc <- abs(peak_w_num - peak_w_ana) / peak_w_ana
    
    peak_val_num <- max(final_n_num)
    peak_val_ana <- max(final_n_ana)
    rel_err_peak_val <- abs(peak_val_num - peak_val_ana) / peak_val_ana
    
    expect_lt(rel_err_total, 0.005)
    expect_lt(rel_err_peak_loc, 0.05)
    expect_lt(rel_err_peak_val, 0.04)
})

test_that("Predictor-corrector agreement improves with spatial numerical diffusion", {
    params <- newMultispeciesParams(data.frame(species = "Test",
                                               w_max = 1000,
                                               w_mat = 100),
                                    no_w = 1000, min_w = 1e-3,
                                    info_level = 0)
    
    params <- setRateFunction(params, "EGrowth", "start_growth")
    params <- setRateFunction(params, "Mort", "start_mort")
    params@ext_diffusion[1, ] <- K * params@w^(p + 1)
    params@species_params$constant_reproduction <- 0
    params <- setRateFunction(params, "RDD", "constantRDD")
    
    w0 <- 1e-2
    t_start <- 0.1
    t_end <- 2
    
    initial_n <- N_analytic(params@w, t_start, w0, 0)
    initialN(params) <- matrix(initial_n, nrow = 1, byrow = TRUE)
    
    sim <- project(params, t_max = t_end - t_start, dt = 0.05,
                   t_save = t_end - t_start,
                   method = "predictor-corrector")
    
    final_n_num <- finalN(sim)[1, ]
    final_n_ana <- N_analytic(params@w, t_end, w0, 0)
    
    beta_grid <- params@w[2] / params@w[1]
    K_eff_pc <- K + A * (beta_grid - 1)
    final_n_ana_pc <- N_analytic(params@w, t_end, w0, 0, K_eff_pc)
    
    rel_err_total <- abs(sum(final_n_num * params@dw) -
                             sum(final_n_ana * params@dw)) /
        sum(final_n_ana * params@dw)
    rel_err_total_pc <- abs(sum(final_n_num * params@dw) -
                                sum(final_n_ana_pc * params@dw)) /
        sum(final_n_ana_pc * params@dw)
    
    peak_w_num <- params@w[which.max(final_n_num)]
    peak_w_ana <- params@w[which.max(final_n_ana)]
    peak_w_ana_pc <- params@w[which.max(final_n_ana_pc)]
    rel_err_peak_loc <- abs(peak_w_num - peak_w_ana) / peak_w_ana
    rel_err_peak_loc_pc <- abs(peak_w_num - peak_w_ana_pc) / peak_w_ana_pc
    
    rel_err_peak_val <- abs(max(final_n_num) - max(final_n_ana)) /
        max(final_n_ana)
    rel_err_peak_val_pc <- abs(max(final_n_num) - max(final_n_ana_pc)) /
        max(final_n_ana_pc)
    
    expect_lt(rel_err_total_pc, rel_err_total)
    expect_lte(rel_err_peak_loc_pc, rel_err_peak_loc)
    expect_lt(rel_err_peak_val_pc, rel_err_peak_val)
    
    expect_lt(rel_err_total_pc, 0.001)
    expect_lt(rel_err_peak_loc_pc, 1e-6)
    expect_lt(rel_err_peak_val_pc, 0.03)
})

test_that("Euler method follows solution with spatial and temporal numerical diffusion", {
    params <- newMultispeciesParams(data.frame(species = "Test",
                                               w_max = 1000,
                                               w_mat = 100),
                                    no_w = 1000, min_w = 1e-3,
                                    info_level = 0)
    
    params <- setRateFunction(params, "EGrowth", "start_growth_euler")
    params <- setRateFunction(params, "Mort", "start_mort_euler")
    params@ext_diffusion[1, ] <- K_euler * params@w^2
    params@species_params$constant_reproduction <- 0
    params <- setRateFunction(params, "RDD", "constantRDD")
    
    beta_grid <- params@w[2] / params@w[1]
    K_num <- A_euler * (beta_grid - 1) + A_euler^2 * dt_euler
    K_eff <- K_euler + K_num
    
    w0 <- 1e-2
    t_start <- 0.1
    t_end <- 1
    
    initial_n <- N_lognormal(params@w, t_start, w0, 0, K_eff)
    initialN(params) <- matrix(initial_n, nrow = 1, byrow = TRUE)
    
    sim <- project(params, t_max = t_end - t_start, dt = dt_euler,
                   t_save = t_end - t_start, t_start = t_start,
                   method = "euler", progress_bar = FALSE)
    
    final_n_num <- finalN(sim)[1, ]
    final_n_ana <- N_lognormal(params@w, t_end, w0, 0, K_eff)
    
    rel_err_total <- abs(sum(final_n_num * params@dw) -
                             sum(final_n_ana * params@dw)) /
        sum(final_n_ana * params@dw)
    
    peak_w_num <- params@w[which.max(final_n_num)]
    peak_w_ana <- params@w[which.max(final_n_ana)]
    rel_err_peak_loc <- abs(peak_w_num - peak_w_ana) / peak_w_ana
    
    rel_err_peak_val <- abs(max(final_n_num) - max(final_n_ana)) /
        max(final_n_ana)
    
    expect_lt(rel_err_total, 0.005)
    expect_lt(rel_err_peak_loc, 0.005)
    expect_lt(rel_err_peak_val, 0.05)
})

# Cleanup
rm(list = c("start_growth", "start_mort", "N_analytic",
            "start_growth_euler", "start_mort_euler", "N_lognormal"),
   envir = globalenv())
