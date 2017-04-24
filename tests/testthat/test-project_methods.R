context("methods used in project")

test_that("getFmortGear",{
    # Two methods:
    # MizerParams + numeric
    # MizerParams + matrix
    data(NS_species_params_gears)
    data(inter)
    params <- MizerParams(NS_species_params_gears, inter)
    # Randomize selectivity and catchability for proper test
    params@catchability[] <- runif(prod(dim(params@catchability)), min=0, max=1)
    params@selectivity[] <- runif(prod(dim(params@selectivity)), min=0, max=1)
    no_gear <- dim(params@catchability)[1]
    no_sp <- dim(params@catchability)[2]
    no_w <- length(params@w)
    # Single numeric
    effort_num1 <- runif(1, min=0.1, max=1)
    # Numeric vector
    effort_num2 <- runif(no_gear, min=0.1, max=1)
    # Matrix (or 2D  array) - here with 7 timesteps
    effort_mat <- array(runif(no_gear*7, min=0.1, max=1),dim=c(7,no_gear))
    # Call both methods with different effort inputs
    f1 <- getFMortGear(params,effort_num1)
    f2 <- getFMortGear(params,effort_num2)
    f3 <- getFMortGear(params,effort_mat)
    # check that number of dims are right
    expect_that(length(dim(f1)), equals(3))
    expect_that(length(dim(f2)), equals(3))
    expect_that(length(dim(f3)), equals(4))
    # check that length of dims is right
    expect_that(dim(f1), equals(c(no_gear,no_sp,no_w)))
    expect_that(dim(f2), equals(c(no_gear,no_sp,no_w)))
    expect_that(dim(f3), equals(c(dim(effort_mat)[1],no_gear,no_sp,no_w)))
    # Check dimnames are right
    expect_equal(names(dimnames(f1)),c("gear","sp","w"))
    expect_equal(names(dimnames(f2)),c("gear","sp","w"))
    expect_equal(names(dimnames(f3))[2:4],c("gear","sp","w"))
    # check fails if effort is not right size
    bad_effort <- rep(effort_num1, no_gear-1)
    expect_error(getFMortGear(params, bad_effort))
    # Check contents of output
    widx <- round(runif(1, min=1, max=no_w))
    sp <- round(runif(1, min=1, max=no_sp))
    gear <- round(runif(1, min=1, max=no_gear))
    expect_equal(effort_num1 * params@catchability[gear,sp] * params@selectivity[gear,sp,widx], f1[gear,sp,widx])
    expect_equal(effort_num2[gear] * params@catchability[gear,sp] * params@selectivity[gear,sp,widx], f2[gear,sp,widx])
    expect_equal(effort_mat[,gear] * params@catchability[gear,sp] * params@selectivity[gear,sp,widx], unname(f3[,gear,sp,widx]))
})

test_that("getFMort",{
    data(NS_species_params_gears)
    data(inter)
    params <- MizerParams(NS_species_params_gears, inter)
    no_gear <- dim(params@catchability)[1]
    no_sp <- dim(params@catchability)[2]
    no_w <- length(params@w)
    effort1 <- 0.5
    effort2 <- rep(effort1,no_gear)
    effort3 <- array(effort1,dim=c(7,no_gear))
    f1 <- getFMort(params,effort1)
    f2 <- getFMort(params,effort2)
    f3 <- getFMort(params,effort3)
    # check that number of dims are right
    expect_that(length(dim(f1)), equals(2))
    expect_that(length(dim(f2)), equals(2))
    expect_that(length(dim(f3)), equals(3))
    # check that length of dims is right
    expect_that(dim(f1), equals(c(no_sp,no_w)))
    expect_that(dim(f2), equals(c(no_sp,no_w)))
    expect_that(dim(f3), equals(c(dim(effort3)[1],no_sp,no_w)))
    # Check dimnames are right
    expect_that(names(dimnames(f3))[2], equals("sp"))
    expect_that(names(dimnames(f3))[3], equals("w"))
    #expect_that(dimnames(f3)$gear
    # check fails if effort is not right size
    expect_error(getFMort(params,c(1,2)))
    # check contents of output
    fmg1 <- getFMortGear(params,effort1)
    fmg2 <- getFMortGear(params,effort2)
    fmg3 <- getFMortGear(params,effort3)
    fmg11 <- array(0,dim=c(no_sp,no_w))
    fmg22 <- array(0,dim=c(no_sp,no_w))
    fmg33 <- array(0,dim=c(dim(effort3)[1],no_sp,no_w))
    for (i in 1:no_gear){
	fmg11 <- fmg11 + fmg1[i,,]
	fmg22 <- fmg22 + fmg2[i,,]
	fmg33 <- fmg33 + fmg3[,i,,]
    }
    expect_that(fmg11, is_equivalent_to(f1))
    expect_that(fmg22, is_equivalent_to(f2))
    expect_that(fmg33, is_equivalent_to(f3))
})


test_that("getPhiPrey",{
    data(NS_species_params_gears)
    data(inter)
    params <- MizerParams(NS_species_params_gears, inter)
    no_sp <- dim(params@catchability)[2]
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)
    n <- abs(array(rnorm(no_w * no_sp), dim = c(no_sp, no_w)))
    n_full <- abs(rnorm(no_w_full))
    pp <- getPhiPrey(params,n,n_full)
    # test dim
    expect_that(dim(pp), equals(c(no_sp,no_w)))
    
    # Test fft based integrator works ok
    
    n <- abs(matrix(rnorm(no_sp*no_w),nrow = no_sp,ncol = no_w))
    n_pp <- abs(rnorm(no_w_full))
    
    ppp <- getPhiPrey(params,n,n_pp)
    
    object <- params
    w0 <- object@w[1]
    Beta <- log(object@species_params$beta)
    sigma <- object@species_params$sigma
    wFull <- object@w_full
    xFull <- log(wFull)
    xFull <- xFull - xFull[1]
    dx <- xFull[2]-xFull[1]
    idx_sp <- (length(object@w_full) - length(object@w) + 1):length(object@w_full)
    fullEnergyMat <- matrix(0, nrow = dim(n)[1], ncol=length(object@w))
    fishEaten <- rep(0, length.out = length(wFull))
    for(i in 1:(dim(n)[1])) {
        fishEaten[idx_sp] <- (object@interaction %*% n)[i, ]
        f2 <- wFull*wFull*(n_pp + fishEaten)
        fullEnergy <- dx*Re(fft(object@ft_pred_kernel_e[i, ]*fft(f2), inverse=TRUE)/length(object@ft_pred_kernel_e[i,]))
        fullEnergyMat[i,] <- fullEnergy[idx_sp]
    }
    expect_that(fullEnergyMat[1,], is_equivalent_to(ppp[1,]))
})

test_that("getFeedingLevel for MizerParams",{
    data(NS_species_params_gears)
    data(inter)
    params <- MizerParams(NS_species_params_gears, inter)
    no_sp <- dim(params@catchability)[2]
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)
    n <- abs(array(rnorm(no_w * no_sp), dim = c(no_sp, no_w)))
    n_full <- abs(rnorm(no_w_full))
    fl <- getFeedingLevel(params,n,n_full)
    # test dim
    expect_that(dim(fl), equals(c(no_sp,no_w)))
    # A crap test - just returns what's already in the method
    phi_prey <- getPhiPrey(params, n=n, n_pp=n_full)
    encount <- params@search_vol * phi_prey
    f <- encount/(encount + params@intake_max)
    expect_that(fl, is_equivalent_to(f))
    # passing in phi_prey gives the same as not
    fl1 <- getFeedingLevel(params,n,n_full)
    phiprey <- getPhiPrey(params,n,n_full)
    fl2 <- getFeedingLevel(params,n,n_full,phi_prey=phi_prey)
    expect_that(fl1, is_identical_to(fl2))
    expect_that(getFeedingLevel(params,n,n_full,phi_prey=matrix(rnorm(10*(no_sp-1)),ncol=10,nrow=no_sp-1)), throws_error())
})

test_that("getFeedingLevel for MizerSim",{
    data(NS_species_params_gears)
    data(inter)
    params <- MizerParams(NS_species_params_gears, inter)
    sim <- project(params, effort=1, t_max=20, dt = 0.5, t_save = 0.5)
    time_range <- 15:20
    expect_that(length(dim(getFeedingLevel(sim, time_range=time_range))), equals(3))
    time_range <- 20
    expect_that(length(dim(getFeedingLevel(sim, time_range=time_range))), equals(3))
    expect_that(getFeedingLevel(sim, time_range=time_range)[1,,], equals(getFeedingLevel(sim@params, sim@n[as.character(time_range),,], sim@n_pp[as.character(time_range),])))
})

# This tests pred rate is calculated as is was in project_methods
 test_that("getPredRate",{
     # First we calculate pr using getPredRate
     data(NS_species_params_gears)
     data(inter)
     params <- MizerParams(NS_species_params_gears, inter)
     no_sp <- dim(params@catchability)[2]
     no_w <- length(params@w)
     no_w_full <- length(params@w_full)
     n <- abs(array(rnorm(no_w * no_sp), dim = c(no_sp, no_w)))
     n_full <- abs(rnorm(no_w_full))
     pr <- getPredRate(params,n,n_full)
     expect_that(dim(pr), equals(c(no_sp,no_w_full)))
     # Next we calculate output, which uses the same procedure as getPredRate, directly
     noSpecies <- dim(params@interaction)[1]
     wfull <- params@w_full
     xfull <- log(wfull)
     xfull <- xfull - xfull[1]
     dx <- xfull[2]-xfull[1]
     idx_sp <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)
     no_P <- length(params@ft_pred_kernel_p[1,])
     Q <- matrix(0, nrow = noSpecies, ncol = no_P )
     feeding_level <- getFeedingLevel(params, n=n, n_pp=n_full)
     Q[, idx_sp] <- sweep((1-feeding_level)*params@search_vol*n, 2, params@w, "*")
     mortLonger <- dx*Re(t(mvfft(t(params@ft_pred_kernel_p) * mvfft(t(Q)), inverse=TRUE)))/no_P
     output <- mortLonger[, 1:length(params@w_full), drop = FALSE]
     # Test that output matches pr
     expect_that(all(output==pr), equals(TRUE))
 })

test_that("getM2 for MizerParams",{
    data(NS_species_params_gears)
    data(inter)
    params <- MizerParams(NS_species_params_gears, inter)
    # Randomize selectivity and catchability for proper test
    params@catchability[] <- runif(prod(dim(params@catchability)), min=0, max=1)
    params@selectivity[] <- runif(prod(dim(params@selectivity)), min=0, max=1)
    no_sp <- dim(params@catchability)[2]
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)
    n <- abs(array(rnorm(no_w * no_sp), dim = c(no_sp, no_w)))
    n_full <- abs(rnorm(no_w_full))
    # Two methods:
    # Params + pred_rate
    # Params + n + n_pp
    
    pred_rate <- getPredRate(params,n,n_full)
    m21 <- getM2(params,pred_rate=pred_rate)
    m22 <- getM2(params,n,n_full)
    # Test dims
    expect_equal(dim(m21), c(no_sp,no_w))
    expect_equal(dim(m21), c(no_sp,no_w))
    # We have gone back to using the old test below since it is compatible with fftclean2
    expect_equal(sum(m22[1,]!=m21[1,]), 0)
    ########################################
    # Look at numbers in a single prey
    w_offset <- no_w_full - no_w
    ##@@ With the new fft based definition of pred_rate, we can just set pred_total equal to pred_rate
    pred_total <- pred_rate
    m2temp <- rep(NA,no_w)
    sp <- runif(1, min=1, max=no_sp)
    for (i in 1:no_w){
        m2temp[i] <- sum(params@interaction[,sp] * pred_total[,w_offset+i])
    }
    expect_equal(m2temp, unname(m21[sp,]))
    expect_equal(m2temp, unname(m21[sp,]))
})

test_that("getM2 for MizerSim",{
    data(NS_species_params_gears)
    data(inter)
    params <- MizerParams(NS_species_params_gears, inter)
    sim <- project(params, effort=1, t_max=20, dt = 0.5, t_save = 0.5)
    time_range <- 15:20
    expect_that(length(dim(getM2(sim, time_range=time_range))), equals(3))
    time_range <- 20
    expect_that(length(dim(getM2(sim, time_range=time_range))), equals(2))
    ##expect_that(getM2(sim, time_range=time_range), equals(getM2(sim@params, sim@n[as.character(time_range),,], sim@n_pp[as.character(time_range),])))
    aq1 <- getM2(sim, time_range=time_range)
    aq2 <- getM2(sim@params, sim@n[as.character(time_range),,], sim@n_pp[as.character(time_range),])
    
    ttot <- 0
    for (i in (1:dim(aq1)[1])){
        ttot <- ttot + sum(aq1[i,]!=aq2[i,])    
    }
    
    expect_that(ttot, equals(0))
})


test_that("getM2Background",{
    data(NS_species_params_gears)
    data(inter)
    params <- MizerParams(NS_species_params_gears, inter)
    no_sp <- dim(params@catchability)[2]
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)
    n <- abs(array(rnorm(no_w * no_sp), dim = c(no_sp, no_w)))
    n_full <- abs(rnorm(no_w_full))
    m2 <- getM2Background(params,n,n_full)
    # test dim
    expect_that(length(m2), equals(no_w_full))
    # Check number in final prey size group
    m22 <- colSums(getPredRate(params,n,n_full))
    expect_that(m22, is_equivalent_to(m2))
    # Passing in pred_rate gives the same
    pr <- getPredRate(params,n,n_full)
    m2b1 <- getM2Background(params,n,n_full)
    m2b2 <- getM2Background(params,n,n_full, pred_rate=pr)
    expect_that(m2b1, is_identical_to(m2b2))
})


test_that("getZ",{
    data(NS_species_params_gears)
    data(inter)
    params <- MizerParams(NS_species_params_gears, inter)
    no_sp <- dim(params@catchability)[2]
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)
    no_gear <- dim(params@catchability)[1]
    n <- abs(array(rnorm(no_w * no_sp), dim = c(no_sp, no_w)))
    n_full <- abs(rnorm(no_w_full))
    effort1 <- 0.5
    effort2 <- rep(effort1,no_gear)
    z <- getZ(params,n,n_full,effort2)
    # test dim
    expect_that(dim(z), equals(c(no_sp,no_w)))
    # Look at numbers in species 1
    f <- getFMort(params,effort2)
    m2 <- getM2(params,n,n_full)
    z1 <- f[1,] + m2[1,] + params@species_params$z0[1]
    expect_that(z1, is_equivalent_to(z[1,]))
    # Passing in M2 gives the same
    m2 <- getM2(params,n,n_full)
    z1 <- getZ(params,n,n_full,effort=effort2)
    z2 <- getZ(params,n,n_full,effort=effort2, m2=m2)
    expect_that(z1, is_identical_to(z2))
})

test_that("getEReproAndGrowth",{
    data(NS_species_params_gears)
    data(inter)
    params <- MizerParams(NS_species_params_gears, inter)
    no_sp <- dim(params@catchability)[2]
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)
    n <- 1e6*abs(array(rnorm(no_w * no_sp), dim = c(no_sp, no_w)))
    n_full <- abs(rnorm(no_w_full))
    erg <- getEReproAndGrowth(params,n,n_full)
    expect_that(all(erg>=0), is_true())
    # test dim
    expect_that(dim(erg), equals(c(no_sp,no_w)))
    # Check number in final prey size group
    f <- getFeedingLevel(params, n=n, n_pp=n_full)
    e <- (f[1,] * params@intake_max[1,]) * params@species_params$alpha[1]
    e <- e - params@std_metab[1,] - params@activity[1,]
    e[e<0] <- 0 # Do not allow negative growth
    expect_that(e, is_equivalent_to(erg[1,]))
    # Adding feeding level gives the same result
    f <- getFeedingLevel(params, n=n, n_pp=n_full)
    erg1 <- getEReproAndGrowth(params,n,n_full)
    erg2 <- getEReproAndGrowth(params,n,n_full, feeding_level=f)
    expect_that(erg1, is_identical_to(erg2))
})


test_that("getESpawning",{
    data(NS_species_params_gears)
    data(inter)
    params <- MizerParams(NS_species_params_gears, inter)
    no_sp <- dim(params@catchability)[2]
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)
    n <- 1e6*abs(array(rnorm(no_w * no_sp), dim = c(no_sp, no_w)))
    n_full <- abs(rnorm(no_w_full))
    es <- getESpawning(params,n,n_full)
    # test dim
    expect_that(dim(es), equals(c(no_sp,no_w)))
    e <- getEReproAndGrowth(params,n=n,n_pp=n_full)
    e_spawning <- params@psi * e 
    expect_that(es, is_equivalent_to(e_spawning))
    e_growth <- getEGrowth(params,n,n_full)
    expect_that(e_growth, is_equivalent_to(e - es))
    # Including ESpawningAndGrowth gives the same
    e <- getEReproAndGrowth(params,n=n,n_pp=n_full)
    es1 <- getESpawning(params,n,n_full)
    es2 <- getESpawning(params,n,n_full, e=e)
    expect_that(es1, is_identical_to(es2))
})


test_that("getRDI",{
    data(NS_species_params_gears)
    data(inter)
    params <- MizerParams(NS_species_params_gears, inter)
    no_sp <- dim(params@catchability)[2]
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)
    n <- 1e6*abs(array(rnorm(no_w * no_sp), dim = c(no_sp, no_w)))
    n_full <- abs(rnorm(no_w_full))
    sex_ratio <- 0.5
    rdi <- getRDI(params,n,n_full, sex_ratio=sex_ratio)
    # test dim
    expect_that(length(rdi), equals(no_sp))
    # test values
    e_spawning <- getESpawning(params, n=n, n_pp=n_full)
    e_spawning_pop <- apply(sweep(e_spawning*n,2,params@dw,"*"),1,sum)
    rdix <- sex_ratio*(e_spawning_pop*params@species_params$erepro)/params@w[params@species_params$w_min_idx] 
    expect_that(c(rdix), is_equivalent_to(c(rdi)))
    rdd <- getRDD(params,n,n_full, sex_ratio=sex_ratio)
    expect_that(length(rdd), equals(no_sp))
    rdd2 <- params@srr(rdi = rdi, species_params = params@species_params)
    expect_that(rdd, is_equivalent_to(rdd))
    # Including ESpawning is the same
    e_spawning <- getESpawning(params, n=n, n_pp=n_full)
    rdi1 <- getRDI(params,n,n_full, sex_ratio=sex_ratio)
    rdi2 <- getRDI(params,n,n_full, sex_ratio=sex_ratio, e_spawning=e_spawning)
    expect_that(rdi1, is_identical_to(rdi2))
})

test_that("interaction is right way round in getM2 method",{
    data(NS_species_params_gears)
    data(inter)
    inter[,"Dab"] <- 0 # Dab not eaten by anything
    params <- MizerParams(NS_species_params_gears, inter)
    m2 <- getM2(params,get_initial_n(params),params@cc_pp)
    expect_that(all(m2["Dab",] == 0), is_true())
})

test_that("getEGrowth is working",{
    data(NS_species_params_gears)
    data(inter)
    params <- MizerParams(NS_species_params_gears, inter)
    no_sp <- dim(params@catchability)[2]
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)
    n <- 1e6*abs(array(rnorm(no_w * no_sp), dim = c(no_sp, no_w)))
    n_full <- abs(rnorm(no_w_full))
    e_spawning <- getESpawning(params, n=n, n_pp=n_full)
    e <- getEReproAndGrowth(params, n=n, n_pp=n_full)
    eg1 <- getEGrowth(params, n=n, n_pp=n_full)
    eg2 <- getEGrowth(params, n=n, n_pp=n_full, e=e, e_spawning=e_spawning)
    expect_that(eg1, is_identical_to(eg2))
    expect_that(e-e_spawning, is_identical_to(eg1))
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
    expect_that(dim(getPhiPrey(params,n,n_pp)), equals(c(1,nw)))
    expect_that(dim(getFeedingLevel(params,n,n_pp)), equals(c(1,nw)))
    expect_that(dim(getPredRate(params,n,n_pp)), equals(c(1,length(params@w_full))))
    expect_that(dim(getM2(params,n,n_pp)), equals(c(1,nw)))
    expect_that(length(getM2Background(params,n,n_pp)), equals(length(params@w_full)))
    expect_that(dim(getFMortGear(params,0)), equals(c(1,1,nw))) # 3D time x species x size
    expect_that(dim(getFMortGear(params,matrix(c(0,0),nrow=2))), equals(c(2,1,1,nw))) # 4D time x gear x species x size
    expect_that(dim(getFMort(params,0)), equals(c(1,nw))) # 2D species x size
    expect_that(dim(getFMort(params,matrix(c(0,0),nrow=2))), equals(c(2,1,nw))) # 3D time x species x size
    expect_that(dim(getZ(params,n,n_pp,0)), equals(c(1,nw)))
    expect_that(dim(getEReproAndGrowth(params,n,n_pp)), equals(c(1,nw)))
    expect_that(dim(getESpawning(params,n,n_pp)), equals(c(1,nw)))
    expect_that(dim(getEGrowth(params,n,n_pp)), equals(c(1,nw)))
    expect_that(length(getRDI(params,n,n_pp)), equals(1))
    expect_that(length(getRDD(params,n,n_pp)), equals(1))

    # MizerSim methods
    expect_that(dim(getFeedingLevel(sim)), equals(c(t_max+1,1,nw))) # time x species x size
    expect_that(dim(getM2(sim)), equals(c(t_max+1,nw))) # time x species x size - default drop is TRUE, if called from plots drop = FALSE
    expect_that(dim(getM2(sim, drop=FALSE)), equals(c(t_max+1,1,nw))) # time x species x size 
    expect_that(dim(getFMortGear(sim)), equals(c(t_max,1,1,nw))) # time x gear x species x size
    expect_that(dim(getFMort(sim)), equals(c(t_max,nw))) # time x species x size - note drop = TRUE
    expect_that(dim(getFMort(sim, drop=FALSE)), equals(c(t_max,1,nw))) # time x species x size 


})
