# Trait based ----
test_that("Providing gamma overrules f0 in newTraitParams()", {
    gamma <- 2000
    params <- newTraitParams(f0 = 0.4, gamma = gamma)
    expect_identical(params@species_params$gamma,
                     rep(gamma, nrow(params@species_params)))
})
# * check a few messages ----
test_that("newTraitParams produces errors and messages", {
    expect_error(newTraitParams(ext_mort_prop = 2),
                 "ext_mort_prop must be a number between 0 and 1")
    expect_error(newTraitParams(reproduction_level = 1),
                   "The reproduction level must be smaller than 1 and non-negative.")
    expect_error(newTraitParams(min_w = -1),
                 "The smallest egg size min_w must be greater than zero.")
    expect_error(newTraitParams(min_w_max = 10^4),
                 "The maximum size of the smallest species min_w_max must be smaller than")
})
# * Multiple gears work correctly in trait-based model ----
test_that("Multiple gears work correctly in trait-based model", {
    # Check multiple gears are working properly
    min_w_max <- 10
    max_w_max <- 1e5
    no_sp <- 10
    w_max <- 10^seq(from = log10(min_w_max), 
                    to = log10(max_w_max), 
                    length = no_sp)
    knife_edges <- w_max * 0.05
    params <- newTraitParams(no_sp = no_sp, 
                             min_w_max = min_w_max, 
                             max_w_max = max_w_max, 
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
                             min_w_max = min_w_max, 
                             max_w_max = max_w_max, 
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

# * Scaling model is set up correctly ----
test_that("Scaling model is set up correctly", {
    (p <- newTraitParams(perfect_scaling = TRUE, sigma = 1,
                         n = 2/3, lambda = 2 + 3/4 - 2/3)) |>
        expect_message("Note: Negative resource abundances")
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
    expect_equal(e, rep(ae, length(e)), tolerance = 1e-1, ignore_attr = TRUE)
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
    expect_equal(total, expected, tolerance = 1e-15, ignore_attr = TRUE)
    
    # All erepros should be equal
    expect_equal(p@species_params$erepro, rep(p@species_params$erepro[1], no_sp))
    
    # Check that total biomass changes little (relatively)
    bm <- getBiomass(sim)
    expect_lt(max(abs(bm[1, ] - bm[6, ])), 1.3e-4)
})

# Community ----
test_that("newCommunityParams works", {
    newCommunityParams(z0 = 0.05, f0 = 0.5) |>
        expect_warning(NA)
})
