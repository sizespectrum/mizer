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

assign("constant_rdd", function(rdi, species_params, params, ...) {
    w_min <- 1e-3
    # Recalculate lambda from A, B, K, p
    a_quad <- K
    b_quad <- -(2 * A - K)
    c_quad <- -2 * B
    det <- b_quad^2 - 4 * a_quad * c_quad
    x <- (-b_quad - sqrt(det)) / (2 * a_quad)
    lambda <- p - x
    
    J_min <- w_min^(p - lambda) * (A - 0.5 * K * (p + 1 - lambda))
    structure(rep(J_min, length(rdi)), names = names(rdi))
}, envir = globalenv())

assign("N_analytic", function(w, t, w0, t0, params) {
    U <- A - 0.5 * K 
    V <- 0.5 * K * (1 - p)
    b <- B / (1 - p)
    nu <- sqrt((U/V)^2 + 4 * b / V)
    dt <- t - t0
    x <- w^(1 - p) / (1 - p)
    x0 <- w0^(1 - p) / (1 - p)
    z <- 2 * sqrt(x * x0) / (V * dt)
    bessel_scaled <- besselI(z, nu, expon.scaled = TRUE)
    log_N_tilde <- -log(V * dt) + (U / (2 * V)) * log(x / x0) - (x + x0) / (V * dt) + z + log(bessel_scaled)
    N_tilde <- exp(log_N_tilde)
    N <- N_tilde * w^(-p)
    return(N)
}, envir = globalenv())

assign("time_dep_rdd", function(rdi, species_params, params, t, ...) {
    t0 <- 0
    w0 <- 10
    U <- A - 0.5 * K 
    V <- 0.5 * K * (1 - p)
    b <- B / (1 - p)
    nu <- sqrt((U/V)^2 + 4 * b / V)
    dt <- t - t0
    if (dt <= 0) return(structure(rep(0, length(rdi)), names = names(rdi)))
    
    w_min <- min(params@w)
    x <- w_min^(1 - p) / (1 - p)
    x0 <- w0^(1 - p) / (1 - p)
    z <- 2 * sqrt(x * x0) / (V * dt)
    
    I_nu <- besselI(z, nu, expon.scaled = TRUE)
    I_nu_plus_1 <- besselI(z, nu + 1, expon.scaled = TRUE)
    
    ratio <- if (I_nu == 0) 0 else I_nu_plus_1 / I_nu
    log_G <- -log(V * dt) + (U / (2 * V)) * log(x / x0) - (x + x0) / (V * dt) + z + log(I_nu)
    G <- exp(log_G)
    
    term <- U/2 + x/dt - (V * z / 2) * ratio - (V * nu / 2)
    J <- G * term
    
    structure(rep(J, length(rdi)), names = names(rdi))
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
                                               w_inf = 1000, 
                                               w_max = 1000,
                                               w_mat = 100, 
                                               beta = 100, 
                                               sigma = 1, 
                                               k_vb = 0.1), 
                                    no_w = 1000, min_w = 1e-3,
                                    info_level = 0)
    
    params <- setRateFunction(params, "EGrowth", "start_growth")
    params <- setRateFunction(params, "Mort", "start_mort")
    params <- setRateFunction(params, "RDD", "constant_rdd")
    params@diffusion[1, ] <- K * params@w^(p + 1)
    params <- setResource(params, resource_dynamics = "resource_constant")
    initialNResource(params) <- 0
    
    initialN(params) <- matrix(params@w^(-lambda), nrow = 1, byrow = TRUE)
    
    # Run
    sim <- project(params, t_max = 1, dt = 0.001)
    
    # Compare
    n0 <- initialN(params)[1, ]
    n1 <- finalN(sim)[1, ]
    # Exclude boundaries
    idx <- 10:(length(params@w) - 10)
    rel_err <- abs(n1[idx] - n0[idx]) / n0[idx]
    
    expect_lt(max(rel_err), 0.05)
})

# 2. Time Dependent Test
test_that("Exact time-dependent solution is followed", {
    # Setup Params
    params <- newMultispeciesParams(data.frame(species = "Test", 
                                               w_inf = 1000, 
                                               w_max = 1000,
                                               w_mat = 100, 
                                               beta = 100, 
                                               sigma = 1, 
                                               k_vb = 0.1), 
                                    no_w = 1000, min_w = 1e-3,
                                    info_level = 0)
    
    params <- setRateFunction(params, "EGrowth", "start_growth")
    params <- setRateFunction(params, "Mort", "start_mort")
    params <- setRateFunction(params, "RDD", "time_dep_rdd")
    params@diffusion[1, ] <- K * params@w^(p + 1)
    params <- setResource(params, resource_dynamics = "resource_constant")
    initialNResource(params) <- 0
    
    w0 <- 10
    t_start <- 1
    t_end <- 2
    
    initial_n <- N_analytic(params@w, t_start, w0, 0, params)
    initialN(params) <- matrix(initial_n, nrow = 1, byrow = TRUE)
    
    sim <- project(params, t_max = t_end - t_start, dt = 0.001)
    
    final_n_num <- finalN(sim)[1, ]
    final_n_ana <- N_analytic(params@w, t_end, w0, 0, params)
    
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
    
    expect_lt(rel_err_total, 0.05)
    expect_lt(rel_err_peak_loc, 0.05)
    expect_lt(rel_err_peak_val, 0.4)
})

# Cleanup
rm(list = c("start_growth", "start_mort", "constant_rdd", "N_analytic", "time_dep_rdd"), envir = globalenv())
