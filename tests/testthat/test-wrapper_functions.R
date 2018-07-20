context("Wrapper functions for trait and community models")

test_that("scaling model is set up correctly", {
  p <- set_scaling_model(perfect = Inf)
  sim <- project(p, t_max=5, effort = 0)
  
  # Check some dimensions
  no_sp <- length(p@species_params$species)
  expect_equal(no_sp, 11)
  
  # Check growth and death rates
  sp <- 6  # check middle species
  gamma <- p@species_params$gamma[sp]
  sigma <- p@species_params$sigma[sp]
  beta <- p@species_params$beta[sp]
  alpha <- p@species_params$alpha[sp]
  h <- p@species_params$h[sp]
  ks <- p@species_params$ks[sp]
  mu0 <- (1 - p@f0) * sqrt(2 * pi) * p@kappa * gamma * sigma *
      (beta ^ (n - 1)) * exp(sigma ^ 2 * (n - 1) ^ 2 / 2)
  hbar <- alpha * h * f0 - ks
  # Death rate
  mu <- getZ(p, p@initial_n, p@initial_n_pp, effort = 0)[sp, ]
  mumu <- mu  # To set the right names
  mumu[] <- mu0 * w^(n-1)
  expect_equal(mu, mumu)
  # Growth rate
  g <- getEGrowth(p, p@initial_n, p@initial_n_pp)[1, ]
  gg <- g  # To set the right names
  gg[] <- hbar * w^n * (1-p@psi[1, ])
  expect_equal(g, gg)
  
  # Check that community is perfect power law
  expect_equal(p@sc, colSums(p@initial_n))
  total <- p@initial_n_pp
  fish_idx <- (length(p@w_full)-length(p@w)+1):length(p@w_full)
  total[fish_idx] <- total[fish_idx] + p@sc
  expected <- total  # To set the names
  expected[] <- p@kappa * p@w_full ^ (-p@lambda)
  expect_equal(total, expected)
  
  # All erepros should be equal
  expect_equal(p@species_params$erepro, rep(p@species_params$erepro[1], no_sp))
  
  # Check that total biomass changes little (relatively)
  bm <- getBiomass(sim)
  expect_lt(max(abs(bm[1, ]-bm[6, ])), 4*10^(-5))
})

test_that("retune_abundance reproduces scaling model", {
    # This numeric test failed on Solaris and without long doubles. So for now
    # skipping it on CRAN
    skip_on_cran()
    p <- set_scaling_model()
    initial_n <- p@initial_n
    p@initial_n[5, ] <- 5 * p@initial_n[5, ]
    retune <- rep(TRUE, length(p@A))
    pr <- retune_abundance(p, retune)
    expect_lt(max(abs(initial_n - pr@initial_n)), 2e-11)
})

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

test_that("trait-based model multiple gears",{
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
