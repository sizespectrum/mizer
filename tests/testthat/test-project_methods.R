context("Methods used in project")


# Initialise --------------------------------------------------------------

# North sea
data(NS_species_params_gears)
data(inter)
params <- MizerParams(NS_species_params_gears, inter)
no_gear <- dim(params@catchability)[1]
no_sp <- dim(params@catchability)[2]
no_w <- length(params@w)
no_w_full <- length(params@w_full)
sim <- project(params, effort=1, t_max=20, dt = 0.5, t_save = 0.5)

# Random abundances
n <- abs(array(rnorm(no_w * no_sp), dim = c(no_sp, no_w)))
n_full <- abs(rnorm(no_w_full))

# Power law
p <- set_scaling_model(no_sp = 2, lambda = 1.5, perfect = TRUE, no_w = 100)
p@initial_n[] <- 0
p@initial_n_pp[] <- p@kappa * p@w_full^(-p@lambda)
sp <- 1  # check first species
sigma <- p@species_params$sigma[sp]
beta <- p@species_params$beta[sp]
gamma <- p@species_params$gamma[sp]
lm2 <- p@lambda - 2


# getAvailEnergy --------------------------------------------------------------

test_that("getAvailEnergy approximates analytic result", {
    ea <- getAvailEnergy(p, p@initial_n, p@initial_n_pp)[sp, ] * p@w^(lm2)
    # Check that this is constant
    expect_equivalent(ea, rep(ea[1], length(ea)), tolerance = 1e-14)
    # Check that it agrees with analytic result
    Dx <- p@w[2] / p@w[1] - 1
    dx <- log(p@w[2] / p@w[1])
    avail_energy_analytic <- p@kappa * exp(lm2^2 * sigma^2 / 2) *
        beta^lm2 * sqrt(2 * pi) * sigma * 
        # The following factor takes into account the cutoff in the integral
        (pnorm(3 - lm2 * sigma) + pnorm(log(beta)/sigma + lm2 * sigma) - 1) *
        Dx / dx
    expect_equivalent(ea[1], avail_energy_analytic, tolerance = 1e-6)
    
    # Check that it agrees with left Riemann sum from w-Beta-3*sigma to w 
    Beta <- log(beta)
    x_full <- log(p@w_full)
    dx <- x_full[2] - x_full[1]
    rr <- Beta + 3*sigma
    jj <- ceiling(rr/dx)
    # Choose some predator weight w[i]
    i <- jj + 100
    ear <- 0
    # Calculate left Riemann sum
    for (j in (i - jj + 1):(i - 1)) {
        ear <- ear + p@w_full[j]^(2 - p@lambda) * 
            exp(-(x_full[i] - x_full[j] - Beta)^2 / (2 * sigma^2))
    }
    ear <- ear * p@kappa * p@w_full[i]^(p@lambda - 2) * dx
    expect_equivalent(ea[1], ear * Dx / dx, tolerance = 1e-14)
})


# getFeedingLevel -----------------------------------------

test_that("getFeedingLevel for MizerParams", {
    fl <- getFeedingLevel(params, n, n_full)
    # test dim
    expect_identical(dim(fl), c(no_sp, no_w))
    # A crap test - just returns what's already in the method
    avail_energy <- getAvailEnergy(params, n = n, n_pp = n_full)
    encount <- params@search_vol * avail_energy
    f <- encount / (encount + params@intake_max)
    expect_identical(fl, f)
    # passing in avail_energy gives the same as not
    fl1 <- getFeedingLevel(params, n, n_full)
    avail_energy <- getAvailEnergy(params, n, n_full)
    fl2 <- getFeedingLevel(params, n, n_full, avail_energy = avail_energy)
    expect_identical(fl1, fl2)
    # calling with avail_energy of wrong dimension gives error
    avail_energy = matrix(rnorm(10 * (no_sp - 1)), ncol = 10, nrow = no_sp - 1)
    expect_error(getFeedingLevel(params, n, n_full, avail_energy = avail_energy),
                 'avail_energy argument must have dimensions: no\\. species \\(12\\) x no. size bins \\(100\\)')
})

test_that("getFeedingLevel for MizerSim", {
    time_range <- 15:20
    expect_length(dim(getFeedingLevel(sim, time_range = time_range)), 3)
    time_range <- 20
    expect_length(dim(getFeedingLevel(sim, time_range = time_range)), 3)
    expect_identical(
        getFeedingLevel(sim, time_range = time_range)[1, , ],
        getFeedingLevel(sim@params, sim@n[as.character(time_range), , ], 
                        sim@n_pp[as.character(time_range), ])
    )
})

test_that("getFeedingLevel approximates analytic result", {
    skip("This is still too imprecise.")
    f <- getFeedingLevel(p, p@initial_n, p@initial_n_pp)[sp, ]
    # Check that this is constant
    expect_equivalent(f, rep(f[1], length(f)), 
                      tolerance = 1e-12, check.names = FALSE)
    expect_equivalent(f[1], 0.6, tolerance = 1.2e-5, check.names = FALSE)
})


# getPredRate -------------------------------------------------------------

test_that("getPredRate approximates analytic result", {
    skip("Still need to understand this.")
    # We use a power law for the species spectrum
    p@initial_n[sp, ] <- p@kappa * p@w^(-p@lambda)
    # and constant feeding level
    f0 <- 0.6
    f <- matrix(f0, nrow = 2, ncol = no_w)
    # Calculate the coefficient of the power law
    pr <- getPredRate(p, p@initial_n, p@initial_n_pp)[sp, ] * p@w_full^(1 - p@n)
    # Check that this is constant
    expect_equal(pr, rep(pr[1], length(pr)), tolerance = 1e-20)
    # Check that it agrees with analytic result
    pred_rate_analytic <- p@kappa * gamma * (1 - f0)
        exp(lm2^2 * sigma^2 / 2) *
        beta^lm2 * sqrt(2 * pi) * sigma * 
        # The following factor takes into account the cutoff in the integral
        (pnorm(3 - lm2 * sigma) + pnorm(log(beta)/sigma + lm2 * sigma) - 1)
    expect_equal(pr[1], pred_rate_analytic, tolerance = 1e-6)
    
    # Check that it agrees with Riemann sum from w to w+Beta-3*sigma
    Beta <- log(beta)
    x_full <- log(p@w_full)
    dx <- x_full[2] - x_full[1]
    rr <- Beta + 3*sigma
    jj <- ceiling(rr/dx)
    # Choose some prey weight w[i]
    i <- 100
    ear <- 0
    # The following corresponds to the right Riemann sum because the sum
    # goes all the way to the right limit of j == i
    for (j in i:(i + jj)) {
        ear <- ear + p@w_full[j]^(p@n - 1) * 
            exp(-(x_full[j] - x_full[i] - Beta)^2 / (2 * sigma^2))
    }
    ear <- ear * (1 - f0) * p@kappa * gamma * p@w_full[i]^(1 - p@n) * dx
    expect_equal(unname(pr[i]), ear, tolerance = 1e-14)
})


# getPredMort -------------------------------------------------------------------

test_that("getPredMort for MizerParams", {
    # Randomize selectivity and catchability for proper test
    params@catchability[] <-
        runif(prod(dim(params@catchability)), min = 0, max = 1)
    params@selectivity[] <-
        runif(prod(dim(params@selectivity)), min = 0, max = 1)
    # Two methods:
    # Params + pred_rate
    # Params + n + n_pp
    
    pred_rate <- getPredRate(params, n, n_full)
    m21 <- getPredMort(params, pred_rate = pred_rate)
    m22 <- getPredMort(params, n, n_full)
    # Test dims
    expect_identical(dim(m21), c(no_sp, no_w))
    expect_identical(dim(m21), c(no_sp, no_w))
    expect_equal(m22[1, ], m21[1, ])
    
    # Look at numbers in a single prey
    w_offset <- no_w_full - no_w
    ##@@ With the new fft based definition of pred_rate, we can just set pred_total equal to pred_rate
    pred_total <- pred_rate
    m2temp <- rep(NA, no_w)
    sp <- runif(1, min = 1, max = no_sp)
    for (i in 1:no_w) {
        m2temp[i] <- sum(params@interaction[, sp] * pred_total[, w_offset + i])
    }
    expect_equal(m2temp, m21[sp, ], check.names = FALSE)
})

test_that("getPredMort for MizerSim", {
    time_range <- 15:20
    expect_length(dim(getPredMort(sim, time_range = time_range)), 3)
    time_range <- 20
    expect_length(dim(getPredMort(sim, time_range = time_range)), 2)
    ##expect_that(getPredMort(sim, time_range=time_range), equals(getPredMort(sim@params, sim@n[as.character(time_range),,], sim@n_pp[as.character(time_range),])))
    aq1 <- getPredMort(sim, time_range = time_range)
    aq2 <- getPredMort(sim@params, sim@n[as.character(time_range), , ],
                 sim@n_pp[as.character(time_range), ])
    
    ttot <- 0
    for (i in (1:dim(aq1)[1])) {
        ttot <- ttot + sum(aq1[i, ] != aq2[i, ])
    }
    
    expect_equal(ttot, 0)
})

test_that("interaction is right way round in getPredMort method", {
    inter[, "Dab"] <- 0  # Dab not eaten by anything
    params <- MizerParams(NS_species_params_gears, inter)
    m2 <- getPredMort(params, get_initial_n(params), params@cc_pp)
    expect_true(all(m2["Dab", ] == 0))
})


# getPlanktonMort ---------------------------------------------------------

test_that("getPlanktonMort", {
    m2 <- getPlanktonMort(params, n, n_full)
    # test dim
    expect_length(m2, no_w_full)
    # Check number in final prey size group
    m22 <- colSums(getPredRate(params, n, n_full))
    expect_identical(m22, m2)
    # Passing in pred_rate gives the same
    pr <- getPredRate(params, n, n_full)
    m2b1 <- getPlanktonMort(params, n, n_full)
    m2b2 <- getPlanktonMort(params, n, n_full, pred_rate = pr)
    expect_identical(m2b1, m2b2)
})


# getFmortGear ------------------------------------------------------------

test_that("getFmortGear", {
    # Two methods:
    # MizerParams + numeric
    # MizerParams + matrix
    # Randomize selectivity and catchability for proper test
    params@catchability[] <-
        runif(prod(dim(params@catchability)), min = 0, max = 1)
    params@selectivity[] <-
        runif(prod(dim(params@selectivity)), min = 0, max = 1)
    no_gear <- dim(params@catchability)[1]
    no_sp <- dim(params@catchability)[2]
    no_w <- length(params@w)
    # Single numeric
    effort_num1 <- runif(1, min = 0.1, max = 1)
    # Numeric vector
    effort_num2 <- runif(no_gear, min = 0.1, max = 1)
    # Matrix (or 2D  array) - here with 7 timesteps
    effort_mat <-
        array(runif(no_gear * 7, min = 0.1, max = 1), dim = c(7, no_gear))
    # Call both methods with different effort inputs
    f1 <- getFMortGear(params, effort_num1)
    f2 <- getFMortGear(params, effort_num2)
    f3 <- getFMortGear(params, effort_mat)
    # Check dimnames are right
    expect_named(dimnames(f1), c("gear", "sp", "w"))
    expect_named(dimnames(f2), c("gear", "sp", "w"))
    expect_named(dimnames(f3)[2:4], c("gear", "sp", "w"))
    expect_identical(dim(f3), c(dim(effort_mat)[1], no_gear, no_sp, no_w))
    # check fails if effort is not right size
    bad_effort <- rep(effort_num1, no_gear - 1)
    expect_error(getFMortGear(params, bad_effort))
    # Check contents of output
    widx <- round(runif(1, min = 1, max = no_w))
    sp <- round(runif(1, min = 1, max = no_sp))
    gear <- round(runif(1, min = 1, max = no_gear))
    expect_identical(f1[gear, sp, widx],
                     effort_num1 * params@catchability[gear, sp] * 
                         params@selectivity[gear, sp, widx])
    expect_identical(f2[gear, sp, widx],
                     effort_num2[gear] * params@catchability[gear, sp] * 
                         params@selectivity[gear, sp, widx])
    expect_identical(f3[, gear, sp, widx],
                     effort_mat[, gear] * params@catchability[gear, sp] * 
                         params@selectivity[gear, sp, widx])
})


# getFMort ----------------------------------------------------------------

test_that("getFMort", {
    effort1 <- 0.5
    effort2 <- rep(effort1, no_gear)
    effort3 <- array(effort1, dim = c(7, no_gear))
    f1 <- getFMort(params, effort1)
    f2 <- getFMort(params, effort2)
    f3 <- getFMort(params, effort3)
    # check that length of dims is right
    expect_identical(dim(f1), c(no_sp, no_w))
    expect_identical(dim(f2), c(no_sp, no_w))
    expect_identical(dim(f3), c(dim(effort3)[1], no_sp, no_w))
    # Check dimnames are right
    expect_named(dimnames(f3)[2:3], c("sp", "w"))
    # check fails if effort is not right size
    expect_error(getFMort(params, c(1, 2)))
    # check contents of output
    fmg1 <- getFMortGear(params, effort1)
    fmg2 <- getFMortGear(params, effort2)
    fmg3 <- getFMortGear(params, effort3)
    fmg11 <- array(0, dim = c(no_sp, no_w))
    fmg22 <- array(0, dim = c(no_sp, no_w))
    fmg33 <- array(0, dim = c(dim(effort3)[1], no_sp, no_w))
    for (i in 1:no_gear) {
        fmg11 <- fmg11 + fmg1[i, , ]
        fmg22 <- fmg22 + fmg2[i, , ]
        fmg33 <- fmg33 + fmg3[, i, , ]
    }
    expect_equal(f1, fmg11)
    expect_equal(f2, fmg22)
    expect_equal(f3, fmg33)
})


# getMort --------------------------------------------------------------------

test_that("getMort", {
    no_gear <- dim(params@catchability)[1]
    effort1 <- 0.5
    effort2 <- rep(effort1, no_gear)
    z <- getMort(params, n, n_full, effort2)
    # test dim
    expect_identical(dim(z), c(no_sp, no_w))
    # Look at numbers in species 1
    f <- getFMort(params, effort2)
    m2 <- getPredMort(params, n, n_full)
    z1 <- f[1, ] + m2[1, ] + params@species_params$z0[1]
    expect_equal(z1, z[1, ], check.names = FALSE)
    # Passing in M2 gives the same
    m2 <- getPredMort(params, n, n_full)
    z1 <- getMort(params, n, n_full, effort = effort2)
    z2 <- getMort(params, n, n_full, effort = effort2, m2 = m2)
    expect_identical(z1, z2)
})


# getEReproAndGrowth ------------------------------------------------------

test_that("getEReproAndGrowth", {
    erg <- getEReproAndGrowth(params, n, n_full)
    expect_true(all(erg >= 0))
    # test dim
    expect_identical(dim(erg), c(no_sp, no_w))
    # Check number in final prey size group
    f <- getFeedingLevel(params, n = n, n_pp = n_full)
    e <-  (f[1, ] * params@intake_max[1, ]) * params@species_params$alpha[1]
    e <- e - params@metab[1, ]
    e[e < 0] <- 0 # Do not allow negative growth
    expect_identical(e, erg[1, ])
    # Adding feeding level gives the same result
    f <- getFeedingLevel(params, n = n, n_pp = n_full)
    erg1 <- getEReproAndGrowth(params, n, n_full)
    erg2 <- getEReproAndGrowth(params, n, n_full, feeding_level = f)
    expect_identical(erg1, erg2)
})


# getERepro ------------------------------------------------------------

test_that("getERepro", {
    n <- 1e6 * abs(array(rnorm(no_w * no_sp), dim = c(no_sp, no_w)))
    es <- getERepro(params, n, n_full)
    # test dim
    expect_identical(dim(es), c(no_sp, no_w))
    e <- getEReproAndGrowth(params, n = n, n_pp = n_full)
    e_repro <- params@psi * e
    expect_identical(es, e_repro)
    e_growth <- getEGrowth(params, n, n_full)
    expect_identical(e_growth, e - es)
    # Including ESpawningAndGrowth gives the same
    e <- getEReproAndGrowth(params, n = n, n_pp = n_full)
    es1 <- getERepro(params, n, n_full)
    es2 <- getERepro(params, n, n_full, e = e)
    expect_identical(es1, es2)
})


# getRDI ------------------------------------------------------------------

test_that("getRDI", {
    n <- 1e6 * abs(array(rnorm(no_w * no_sp), dim = c(no_sp, no_w)))
    sex_ratio <- 0.5
    rdi <- getRDI(params, n, n_full, sex_ratio = sex_ratio)
    # test dim
    expect_length(rdi, no_sp)
    # test values
    e_repro <- getERepro(params, n = n, n_pp = n_full)
    e_repro_pop <- apply(sweep(e_repro * n, 2, params@dw, "*"), 1, sum)
    rdix <- sex_ratio * (e_repro_pop * params@species_params$erepro) / 
        params@w[params@w_min_idx]
    expect_equal(rdix, rdi, tolerance = 1e-15, check.names = FALSE)
    # Including ESpawning is the same
    e_repro <- getERepro(params, n = n, n_pp = n_full)
    rdi1 <- getRDI(params, n, n_full, sex_ratio = sex_ratio)
    rdi2 <- getRDI(params, n, n_full, sex_ratio = sex_ratio, 
                   e_repro = e_repro)
    expect_identical(rdi1, rdi2)
})


# getRDD ------------------------------------------------------------------

test_that("getRDD", {
    rdd <- getRDD(params, n, n_full)
    expect_length(rdd, no_sp)
    rdi <- getRDI(params, n, n_full)
    rdd2 <- getRDD(params, n, n_full, rdi = rdi)
    expect_identical(rdd, rdd2)
    rdd2 <- params@srr(rdi = rdi, species_params = params@species_params)
    expect_identical(rdd, rdd2)
})


# getEGrowth --------------------------------------------------------------

test_that("getEGrowth is working", {
    n <- 1e6 * abs(array(rnorm(no_w * no_sp), dim = c(no_sp, no_w)))
    e_repro <- getERepro(params, n = n, n_pp = n_full)
    e <- getEReproAndGrowth(params, n = n, n_pp = n_full)
    eg1 <- getEGrowth(params, n = n, n_pp = n_full)
    eg2 <- getEGrowth(params, n = n, n_pp = n_full, e = e, 
                      e_repro = e_repro)
    expect_identical(eg1, eg2)
    expect_identical(e - e_repro, eg1)
})


# general -----------------------------------------------------------------

test_that("Test that fft based integrator gives similar result as old code", {
    # Make it harder by working with kernels that need a lot of cutoff
    species_params <- NS_species_params_gears
    species_params$sigma[3] <- 3
    species_params$beta <- species_params$beta / 100
    # and use different egg sizes
    species_params$w_min <- seq(0.001, 1, length.out = no_sp)
    params <- MizerParams(species_params, inter, no_w = 30)
    # create a second params object that does not use fft
    params2 <- params
    params2@ft_pred_kernel_e <- array()
    params2@ft_pred_kernel_p <- array()
    # Test available energy integral
    afft <- getAvailEnergy(params, params@initial_n, params@initial_n_pp)
    a <- getAvailEnergy(params2, params@initial_n, params@initial_n_pp)
    # Only check values at fish sizes
    fish <- outer(1:no_sp, 1:no_w, function(i, a) a >= params@w_min_idx[i])
    expect_equivalent(afft[fish], a[fish], tolerance = 2e-15)
    # Test available energy integral
    pfft <- getPredRate(params, params@initial_n, params@initial_n_pp)
    p <- getPredRate(params2, params@initial_n, params@initial_n_pp)
    expect_equivalent(pfft, p, tolerance = 1e-15)
})

test_that("project methods return objects of correct dimension when community only has one species",{
    params <- set_community_model(z0 = 0.2, f0 = 0.7, alpha = 0.2, recruitment = 4e7)
    t_max <- 50
    sim <- project(params, t_max=t_max, effort = 0)
    n <- array(sim@n[t_max+1,,],dim=dim(sim@n)[2:3])
    dimnames(n) <- dimnames(sim@n)[2:3]
	n_pp <- sim@n_pp[1,]
    nw <- length(params@w)
    # MizerParams methods
    expect_that(dim(getAvailEnergy(params,n,n_pp)), equals(c(1,nw)))
    expect_that(dim(getFeedingLevel(params,n,n_pp)), equals(c(1,nw)))
    expect_that(dim(getPredRate(params,n,n_pp)), equals(c(1,length(params@w_full))))
    expect_that(dim(getPredMort(params,n,n_pp)), equals(c(1,nw)))
    expect_that(length(getPlanktonMort(params,n,n_pp)), equals(length(params@w_full)))
    expect_that(dim(getFMortGear(params,0)), equals(c(1,1,nw))) # 3D time x species x size
    expect_that(dim(getFMortGear(params,matrix(c(0,0),nrow=2))), equals(c(2,1,1,nw))) # 4D time x gear x species x size
    expect_that(dim(getFMort(params,0)), equals(c(1,nw))) # 2D species x size
    expect_that(dim(getFMort(params,matrix(c(0,0),nrow=2))), equals(c(2,1,nw))) # 3D time x species x size
    expect_that(dim(getMort(params,n,n_pp,0)), equals(c(1,nw)))
    expect_that(dim(getEReproAndGrowth(params,n,n_pp)), equals(c(1,nw)))
    expect_that(dim(getERepro(params,n,n_pp)), equals(c(1,nw)))
    expect_that(dim(getEGrowth(params,n,n_pp)), equals(c(1,nw)))
    expect_that(length(getRDI(params,n,n_pp)), equals(1))
    expect_that(length(getRDD(params,n,n_pp)), equals(1))

    # MizerSim methods
    expect_that(dim(getFeedingLevel(sim)), equals(c(t_max+1,1,nw))) # time x species x size
    expect_that(dim(getPredMort(sim)), equals(c(t_max+1,nw))) # time x species x size - default drop is TRUE, if called from plots drop = FALSE
    expect_that(dim(getPredMort(sim, drop=FALSE)), equals(c(t_max+1,1,nw))) # time x species x size 
    expect_that(dim(getFMortGear(sim)), equals(c(t_max+1,1,1,nw))) # time x gear x species x size
    expect_that(dim(getFMort(sim)), equals(c(t_max+1,nw))) # time x species x size - note drop = TRUE
    expect_that(dim(getFMort(sim, drop=FALSE)), equals(c(t_max+1,1,nw))) # time x species x size 

})
