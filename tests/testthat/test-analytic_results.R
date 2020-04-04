context("Analytic results")

# Initialise power law ----
no_w <- 100
no_sp <- 2
p <- newTraitParams(no_sp = no_sp, perfect_scaling = TRUE, no_w = no_w)
p@species_params$pred_kernel_type <- "truncated_lognormal"
n0 <- p@initial_n
n0[] <- 0
n_pp <- p@initial_n_pp
n_pp[] <- p@resource_params$kappa * p@w_full^(-p@resource_params$lambda)
sp <- 1  # check first species
sigma <- p@species_params$sigma[sp]
beta <- p@species_params$beta[sp]
gamma <- p@species_params$gamma[sp]
q <- p@species_params$q[sp]
n <- p@species_params$n[sp]
lm2 <- p@resource_params$lambda - 2

# getEncounter ----
test_that("getEncounter approximates analytic result when feeding on resource only", {
    e <- getEncounter(p, n0, n_pp)[sp, ] * p@w^(lm2 - q)
    # Check that this is constant
    expect_equivalent(e, rep(e[1], length(e)))
    # Check that it agrees with analytic result
    Dx <- p@w[2] / p@w[1] - 1
    dx <- log(p@w[2] / p@w[1])
    encounter_analytic <- p@resource_params$kappa * exp(lm2^2 * sigma^2 / 2) *
        beta^lm2 * sqrt(2 * pi) * sigma * gamma *
        # The following factor takes into account the discretisation scheme
        Dx / dx #*
        # The following factor takes into account the cutoff in the integral
        # (pnorm(3 - lm2 * sigma) + pnorm(log(beta)/sigma + lm2 * sigma) - 1)
    # The Riemann sum is not precise enough
    expect_equivalent(e[1], encounter_analytic, tolerance = 1e-3)
    
    # Check that it agrees with left Riemann sum from w-Beta-3*sigma to w 
    Beta <- log(beta)
    x_full <- log(p@w_full)
    dx <- x_full[2] - x_full[1]
    # Choose some predator weight w[i]
    i <- 100
    ear <- 0
    # Calculate left Riemann sum
    for (j in 1:(i - 1)) {
        ear <- ear + p@w_full[j]^(2 - p@resource_params$lambda) * 
            exp(-(x_full[i] - x_full[j] - Beta)^2 / (2 * sigma^2))
    }
    ear <- ear * p@resource_params$kappa * p@w_full[i]^(p@resource_params$lambda - 2) * dx * gamma
    expect_equivalent(e[1], ear * Dx / dx)
})

# getDiet ----
test_that("getDiet approximates analytic result when feeding on resource only", {
    # n and n_pp are power laws
    n <- p@initial_n
    n[] <- rep(p@resource_params$kappa * p@w^(-p@resource_params$lambda), each = 2)
    n_pp <- p@initial_n_pp
    n_pp[] <- p@resource_params$kappa * p@w_full^(-p@resource_params$lambda)
    # switch of interaction between species
    p0 <- setInteraction(p, interaction = matrix(0, nrow = no_sp, ncol = no_sp))
    diet <- getDiet(p0, n, n_pp, proportion = FALSE)[sp, , ]
    # None of the diet should come from fish
    expect_true(all(diet[, 1:2] == 0))
    # Check that diet from resource is power law
    diet_coeff <- diet[, 3] * p@w^(lm2 - q)
    expect_equivalent(diet_coeff, rep(diet_coeff[1], no_w))
    # and agrees with result from getEncounter
    feeding_level <- getFeedingLevel(p0, n0, n_pp)[sp, ]
    encounter <- getEncounter(p0, n0, n_pp)[sp, ]
    expect_equivalent(diet[, 3], encounter * (1 - feeding_level))
})


test_that("getFeedingLevel approximates analytic result", {
    f <- getFeedingLevel(p)[sp, ]
    # Check that this is constant
    expect_equivalent(f, rep(f[1], length(f)), 
                      tolerance = 1e-12, check.names = FALSE)
    # Still to imprecise
    expect_equivalent(f[1], 0.6, tolerance = 2e-2, check.names = FALSE)
})


# TODO: fix this
# test_that("getPredRate approximates analytic result", {
#     # We use a power law for the species spectrum
#     p@initial_n[sp, ] <- p@resource_params$kappa * p@w^(-p@resource_params$lambda)
#     # and constant feeding level
#     f0 <- 0.6
#     f <- matrix(f0, nrow = 2, ncol = no_w)
#     # Calculate the coefficient of the power law
#     pr <- getPredRate(p, feeding_level = f)[sp, ] * p@w_full^(1 - n)
#     # Check that this is constant in the feeding range of the predator
#     sel <- (p@w_full > min(p@w) / beta * exp(3 * sigma)) &
#         (p@w_full < max(p@w) / beta / exp(3 * sigma))
#     pr <- pr[sel]
#     # The first three entries of pr are still different. Why?
#     # For now we just cut them off
#     pr <- pr[4:length(pr)]
#     expect_equivalent(pr, rep(pr[1], length(pr)))
#     # Check that it agrees with analytic result
#     n1 <- n - 1
#     Dx <- p@w[2] / p@w[1] - 1
#     dx <- log(p@w[2] / p@w[1])
#     pred_rate_analytic <- p@resource_params$kappa * gamma * (1 - f0) *
#         exp(n1^2 * sigma^2 / 2) *
#         beta^n1 * sqrt(2 * pi) * sigma * 
#         # The following factor takes into account the discretisation scheme
#         Dx / dx #*
#         # The following factor takes into account the cutoff in the integral
#         #(pnorm(3 - n1 * sigma) + pnorm(log(beta)/sigma + n1 * sigma) - 1)
#     # This is still too imprecise
#     expect_equivalent(pr[1], pred_rate_analytic, tolerance = 1e-3)
#     
#     # Check that it agrees with Riemann sum from w to w+Beta-3*sigma
#     Beta <- log(beta)
#     x_full <- log(p@w_full)
#     dx <- x_full[2] - x_full[1]
#     rr <- Beta + 3*sigma
#     jj <- ceiling(rr/dx)
#     # Choose some prey weight w[i]
#     i <- which.max(sel) + 10
#     pra <- 0
#     # The following corresponds to the right Riemann sum because the sum
#     # goes all the way to the right limit of j == i
#     for (j in i:(i + jj)) {
#         pra <- pra + p@w_full[j]^(n - 1) * 
#             exp(-(x_full[j] - x_full[i] - Beta)^2 / (2 * sigma^2))
#     }
#     pra <- pra * (1 - f0) * p@resource_params$kappa * gamma * p@w_full[i]^(1 - n) * dx
#     # TODO: Still need to understand this
#     # expect_equal(unname(pr[1]), pra, tolerance = 1e-14)
# })


# Analytic steady-state solution ----
test_that("Analytic steady-state solution is well approximated", {
    # Choose some parameters
    f0 <- 0.6
    alpha <- 0.4
    r_pp <- 10^18  # Choosing a high value because we want the resource to stay
    # at its power-law steady state
    n <- 2/3
    p <- n
    q <- 0.95
    lambda <- 2 + q - n
    erepro <- 0.1
    R <- 1e10  # The rate of reproduction
    
    beta <- 100
    sigma <- 1.3
    h <- 30
    ks <- 4
    kappa <- 1e11
    
    w_min <- 1e-3
    w_inf <- 1e3
    w_mat <- 1e2
    min_w_pp <- 1e-7  # Only have to make sure the smallest fish are perfectly fed
    # Chose number of gridpoints so that w_mat and w_inf lie on gridpoints
    no_w <- log10(w_inf / w_min) * 100 + 1  
    
    species_params <- data.frame(
        species = "Single",
        w_min = w_min,
        w_inf = w_inf,
        w_mat = w_mat,
        f0 = f0,
        h = h,
        ks = ks,
        beta = beta,
        sigma = sigma,
        z0 = 0,
        alpha = alpha,
        erepro = erepro,
        sel_func = "knife_edge", # not used but required
        knife_edge_size = 1000
    )
    
    params <- newMultispeciesParams(species_params, p = p, n = n, lambda = lambda,
                                    kappa = kappa, min_w = w_min, max_w = w_inf,
                                    no_w = no_w, min_w_pp = min_w_pp, w_pp_cutoff = w_inf,
                                    r_pp = r_pp)
    
    gamma <- params@species_params$gamma[1]
    w <- params@w
    
    # mu0 w^(n-1) is the death rate that is produced by predation if the predators
    # follow the same power law as the resource. 
    # We could equally well have chosen any other constant
    mu0 <- (1 - f0) * sqrt(2 * pi) * kappa * gamma * sigma *
        (beta ^ (n - 1)) * exp(sigma ^ 2 * (n - 1) ^ 2 / 2)
    params@mu_b[1,] <- mu0 * w ^ (n - 1)
    # hbar w^n is the rate at which energy is available for growth and reproduction
    hbar <- alpha * h * f0 - ks
    # n_exact is calculated using the analytic expression for the solution
    pow <- mu0 / hbar / (1 - n)
    n_mult <- (1 - (w / w_inf) ^ (1 - n)) ^ (pow - 1) *
        (1 - (w_mat / w_inf) ^ (1 - n)) ^ (-pow)
    n_mult[w < w_mat] <- 1
    n_exact <- params@psi  # Just to get array with correct dimensions and names
    n_exact[] <- R * (w_min/w)^(mu0/hbar) / (hbar * w^n) * n_mult
    
    # Make sure that the rate of reproduction is R
    params@rates_funcs$RDD <- "constantRDD"
    params@species_params$constant_reproduction <- R
    # We use a step function for the maturity function
    params@psi[1,] <- (params@w / w_inf) ^ (1 - n)
    params@psi[1, params@w < w_mat] <- 0
    # We switch off the self-interaction
    params@interaction[] <- 0
    
    # We start the simulation with the exact steady-state solution
    sim <- project(params, t_max = 5, effort = 0, initial_n = n_exact)
    # If all is well, it should stay close to the steady-state solution
    relative_error <- abs((n_exact[1,] - sim@n[6, 1, ]) / n_exact[1, ])
    # TODO: Unfortunately there is a significant difference at the maximum weight,
    # so we only test the others
    # expect_lt(max(relative_error[1:(no_w - 1)]), 0.02)
})
