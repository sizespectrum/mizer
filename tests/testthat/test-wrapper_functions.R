context("Wrapper functions for trait and community models")

# Scaling model is set up correctly ----
test_that("Scaling model is set up correctly", {
    skip("Does not work yet")
    p <- set_scaling_model(perfect = TRUE, sigma=1)
    sim <- project(p, t_max = 5)
    
    # Check some dimensions
    no_sp <- length(p@species_params$species)
    expect_equal(no_sp, 11)
    
    # Check against analytic results
    sp <- 6  # check middle species
    gamma <- p@species_params$gamma[sp]
    sigma <- p@species_params$sigma[sp]
    beta <- p@species_params$beta[sp]
    alpha <- p@species_params$alpha[sp]
    h <- p@species_params$h[sp]
    ks <- p@species_params$ks[sp]
    mu0 <- (1 - p@f0) * sqrt(2 * pi) * p@kappa * gamma * sigma *
        (beta ^ (p@n - 1)) * exp(sigma ^ 2 * (p@n - 1) ^ 2 / 2)
    hbar <- alpha * h * p@f0 - ks
    # Check encounter rate
    lm2 <- p@lambda - 2
    e <- getEncounter(p, p@initial_n, p@initial_n_pp)[sp, ] * p@w^(lm2 - p@q)
    ae <- gamma * p@kappa * exp(lm2^2 * sigma^2 / 2) *
        beta^lm2 * sqrt(2 * pi) * sigma * 
        # The following factor takes into account the cutoff in the integral
        (pnorm(3 - lm2 * sigma) + pnorm(log(beta)/sigma + lm2 * sigma) - 1)
    expect_equal(ea, rep(ae, length(ea)), tolerance = 1e-15)
    # Check feeding level
    f <- getFeedingLevel(p, p@initial_n, p@initial_n_pp)[sp, ]
    names(f) <- NULL
    expect_equal(f, rep(f[1], length(f)), tolerance = 1e-14)
    # Death rate
    mu <- getMort(p, p@initial_n, p@initial_n_pp, effort = 0)[sp, ]
    mumu <- mu  # To set the right names
    mumu[] <- mu0 * p@w^(p@n - 1)
    expect_equal(mu, mumu, tolerance = 1e-15)
    # Growth rate
    g <- getEGrowth(p, p@initial_n, p@initial_n_pp)[sp, ]
    gg <- g  # To set the right names
    gg[] <- hbar * p@w^p@n * (1 - p@psi[sp, ])
    expect_equal(g, gg)
    
    # Check that community is perfect power law
    expect_identical(p@sc, colSums(p@initial_n))
    total <- p@initial_n_pp
    fish_idx <- (length(p@w_full) - length(p@w) + 1):length(p@w_full)
    total[fish_idx] <- total[fish_idx] + p@sc
    total <- total * p@w_full^p@lambda
    expected <- rep(p@kappa, length(p@w_full))
    expect_equivalent(total, expected, tolerance = 1e-15, check.names = FALSE)
    
    # All erepros should be equal
    expect_equal(p@species_params$erepro, rep(p@species_params$erepro[1], no_sp))
    
    # Check that total biomass changes little (relatively)
    bm <- getBiomass(sim)
    expect_lt(max(abs(bm[1, ] - bm[6, ])), 4*10^(-5))
})

# retuneAbundance() reproduces scaling model ----
test_that("retuneAbundance() reproduces scaling model", {
    # This numeric test failed on Solaris and without long doubles. So for now
    # skipping it on CRAN
    skip_on_cran()
    p <- set_scaling_model()
    initial_n <- p@initial_n
    p@initial_n[5, ] <- 5 * p@initial_n[5, ]
    retune <- rep(TRUE, length(p@A))
    pr <- retuneAbundance(p, retune)
    expect_lt(max(abs(initial_n - pr@initial_n)), 2e-11)
})

# addSpecies works when adding a second identical species ----
test_that("addSpecies works when adding a second identical species", {
    p <- set_scaling_model()
    no_sp <- length(p@A)
    p <- setBackground(p)
    species_params <- p@species_params[5,]
    species_params$species = "new"
    SSB <- sum(p@initial_n[5, ] * p@w * p@dw * p@psi[5, ])
    # Adding species 5 again at half its background biomass should lead two
    # two copies of the species each with half the biomass
    pa <- addSpecies(p, species_params, SSB = SSB/2, rfac=Inf)
    pa@initial_n[5, ] <- pa@initial_n[5, ] + pa@initial_n[no_sp+1, ]
    expect_lt(max(abs(p@initial_n - pa@initial_n[1:no_sp, ])), 1)
    expect_lt(max(abs(p@initial_n[5,] - pa@initial_n[5, ])), 0.7)
})

# Multiple gears work correctly in trait-based model ----
test_that("Multiple gears work correctly in trait-based model", {
    # Check multiple gears are working properly
    min_w_inf <- 10
    max_w_inf <- 1e5
    no_sp <- 10
    w_inf <- 10^seq(from=log10(min_w_inf), to = log10(max_w_inf), length=no_sp)
    knife_edges <- w_inf * 0.05
    params <- set_trait_model(no_sp = no_sp, min_w_inf = min_w_inf, max_w_inf = max_w_inf, knife_edge_size = knife_edges)
    expect_that(params@species_params$knife_edge_size, is_identical_to(knife_edges))
    # All gears fire
    sim1 <- project(params, t_max = 10, effort = 1)
    fmg <- getFMortGear(sim1)
    for (i in 1:no_sp){
        expect_that(all(fmg[10,1,i,params@w < knife_edges[i]] == 0), is_true())
        expect_that(all(fmg[10,1,i,params@w >= knife_edges[i]] == 1), is_true())
    }
    # Only the 4th gear fires
    params <- set_trait_model(no_sp = no_sp, min_w_inf = min_w_inf, max_w_inf = max_w_inf, knife_edge_size = knife_edges, gear_names = 1:no_sp)
    effort <- c(0,0,0,1,0,0,0,0,0,0)
    names(effort) = 1:no_sp
    sim2 <- project(params, t_max = 10, effort = effort)
    fmg <- getFMortGear(sim2)
    expect_that(all(fmg[10,c(1:3,5:10),c(1:3,5:10),] == 0), is_true())
    expect_that(all(fmg[10,4,4,params@w < knife_edges[4]] == 0), is_true())
    expect_that(all(fmg[10,4,4,params@w >= knife_edges[4]] == 1), is_true())
    
})
