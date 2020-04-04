context("Wrapper functions for trait and community models")

# Multiple gears work correctly in trait-based model ----
test_that("Multiple gears work correctly in trait-based model", {
    # Check multiple gears are working properly
    min_w_inf <- 10
    max_w_inf <- 1e5
    no_sp <- 10
    w_inf <- 10^seq(from = log10(min_w_inf), 
                    to = log10(max_w_inf), 
                    length = no_sp)
    knife_edges <- w_inf * 0.05
    params <- newTraitParams(no_sp = no_sp, 
                             min_w_inf = min_w_inf, 
                             max_w_inf = max_w_inf, 
                             knife_edge_size = knife_edges)
    expect_identical(params@gear_params$knife_edge_size, 
                     knife_edges)
    # All gears fire
    sim1 <- project(params, t_max = 10, effort = 1)
    fmg <- getFMortGear(sim1)
    for (i in 1:no_sp) {
        expect_true(all(fmg[10,1,i,params@w < knife_edges[i]] == 0))
        expect_true(all(fmg[10,1,i,params@w >= knife_edges[i]] == 1))
    }
    # Only the 4th gear fires
    params <- newTraitParams(no_sp = no_sp, 
                             min_w_inf = min_w_inf, 
                             max_w_inf = max_w_inf, 
                             knife_edge_size = knife_edges, 
                             gear_names = 1:no_sp)
    effort <- c(0,0,0,1,0,0,0,0,0,0)
    names(effort) <- 1:no_sp
    sim2 <- project(params, t_max = 10, effort = effort)
    fmg <- getFMortGear(sim2)
    expect_true(all(fmg[10, c(1:3,5:10),c(1:3,5:10),] == 0))
    expect_true(all(fmg[10, 4, 4, params@w < knife_edges[4]] == 0))
    expect_true(all(fmg[10, 4, 4, params@w >= knife_edges[4]] == 1))
    
})

# Scaling model is set up correctly ----
test_that("Scaling model is set up correctly", {
    p <- newTraitParams(perfect_scaling = TRUE, sigma = 1,
                        n = 2/3, lambda = 2 + 3/4 - 2/3)
    sim <- project(p, t_max = 5)
    
    # Check some dimensions
    no_sp <- length(p@species_params$species)
    expect_equal(no_sp, 11)
    
    # Check against analytic results
    sp <- 6  # check middle species
    gamma <- p@species_params$gamma[[sp]]
    sigma <- p@species_params$sigma[[sp]]
    beta <- p@species_params$beta[[sp]]
    alpha <- p@species_params$alpha[[sp]]
    h <- p@species_params$h[[sp]]
    ks <- p@species_params$ks[[sp]]
    f0 <- 0.6
    n <- p@species_params$n[[sp]]
    mu0 <- (1 - f0) * sqrt(2 * pi) * 
        p@resource_params$kappa * gamma * sigma *
        (beta ^ (n - 1)) * exp(sigma ^ 2 * (n - 1) ^ 2 / 2)
    hbar <- alpha * h * f0 - ks
    # Check encounter rate
    lm2 <- p@resource_params$lambda - 2
    q <- p@species_params$q[[sp]]
    e <- getEncounter(p, p@initial_n, p@initial_n_pp)[sp, ] * p@w^(lm2 - q)
    ae <- gamma * p@resource_params$kappa * exp(lm2^2 * sigma^2 / 2) *
        beta^lm2 * sqrt(2 * pi) * sigma * 
        # The following factor takes into account the cutoff in the integral
        (pnorm(3 - lm2 * sigma) + pnorm(log(beta)/sigma + lm2 * sigma) - 1)
    # TODO: not precise enough yet
    expect_equivalent(e, rep(ae, length(e)), tolerance = 1e-1)
    # Check feeding level
    f <- getFeedingLevel(p, p@initial_n, p@initial_n_pp)[sp, ]
    names(f) <- NULL
    expect_equal(f, rep(f[1], length(f)), tolerance = 1e-14)
    # Death rate
    mu <- getMort(p, p@initial_n, p@initial_n_pp, effort = 0)[sp, ]
    mumu <- mu  # To set the right names
    mumu[] <- mu0 * p@w^(p@species_params$n[[sp]] - 1)
    expect_equal(mu, mumu, tolerance = 0.2)
    # Growth rate
    g <- getEGrowth(p, p@initial_n, p@initial_n_pp)[sp, ]
    gg <- g  # To set the right names
    gg[] <- hbar * p@w^p@species_params$n[[sp]] * (1 - p@psi[sp, ])
    # TODO: not precise enough yet
    expect_equal(g, gg, tolerance = 1e-4)
    
    # Check that community is perfect power law
    expect_identical(p@sc, colSums(p@initial_n))
    total <- p@initial_n_pp
    fish_idx <- (length(p@w_full) - length(p@w) + 1):length(p@w_full)
    total[fish_idx] <- total[fish_idx] + p@sc
    total <- total * p@w_full^p@resource_params$lambda
    expected <- rep(p@resource_params$kappa, length(p@w_full))
    expect_equivalent(total, expected, tolerance = 1e-15, check.names = FALSE)
    
    # All erepros should be equal
    expect_equal(p@species_params$erepro, rep(p@species_params$erepro[1], no_sp))
    
    # Check that total biomass changes little (relatively)
    bm <- getBiomass(sim)
    expect_lt(max(abs(bm[1, ] - bm[6, ])), 1.3e-4)
})


# setRmax works ----
test_that("setRmax works", {
    params <- NS_params
    params@rates_funcs$RDD <- "noRDD"
    rdd <- getRDD(params)
    R_factor <- 5
    params <- setRmax(params, R_factor)
    expect_equivalent(params@species_params$R_max, rdd * R_factor)
    expect_equal(getRDD(params), rdd)
})
