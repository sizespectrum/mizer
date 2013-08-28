context("methods used in project")

test_that("getFmortGear",{
    data(NS_species_params_gears)
    data(inter)
    params <- MizerParams(NS_species_params_gears, inter)
    no_gear <- dim(params@catchability)[1]
    no_sp <- dim(params@catchability)[2]
    no_w <- length(params@w)
    effort1 <- 0.5
    effort2 <- rep(effort1,no_gear)
    effort3 <- array(effort1,dim=c(7,no_gear))
    f1 <- getFMortGear(params,effort1)
    f2 <- getFMortGear(params,effort2)
    f3 <- getFMortGear(params,effort3)
    # check that number of dims are right
    expect_that(length(dim(f1)), equals(3))
    expect_that(length(dim(f2)), equals(3))
    expect_that(length(dim(f3)), equals(4))
    # check that length of dims is right
    expect_that(dim(f1), equals(c(no_gear,no_sp,no_w)))
    expect_that(dim(f2), equals(c(no_gear,no_sp,no_w)))
    expect_that(dim(f3), equals(c(dim(effort3)[1],no_gear,no_sp,no_w)))
    # Check dimnames are right
    expect_that(names(dimnames(f3))[2], equals("gear"))
    expect_that(names(dimnames(f3))[3], equals("sp"))
    expect_that(names(dimnames(f3))[4], equals("w"))
    #expect_that(dimnames(f3)$gear
    # check fails if effort is not right size
    expect_error(getFMortGear(params,c(1,2)))
    # check contents of output
    expect_that(f1[1,1,], equals(params@catchability[1,1] * params@selectivity[1,1,] * effort1[1])) 
    expect_that(f1[no_gear,no_sp,], equals(params@catchability[no_gear,no_sp] * params@selectivity[no_gear,no_sp,] * effort1[1])) 
    expect_that(f2[1,1,1], equals(params@catchability[1,1] * params@selectivity[1,1,1] * effort2[1])) 
    expect_that(f2[no_gear,no_sp,no_w], equals(params@catchability[no_gear,no_sp] * params@selectivity[no_gear,no_sp,no_w] * effort2[no_gear])) 
    expect_that(f3[1,1,1,], equals(params@catchability[1,1] * params@selectivity[1,1,] * effort3[1,1])) 
    expect_that(f3[dim(effort3)[1],no_gear,no_sp,], equals(params@catchability[no_gear,no_sp] * params@selectivity[no_gear,no_sp,] * effort3[dim(effort3)[1],no_gear])) 
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
    # Test numbers are right
    # Hideous - doing it by hand for first predator - should be the same
    n <- abs(matrix(rnorm(no_sp*no_w),nrow = no_sp,ncol = no_w))
    n_pp <- abs(rnorm(no_w_full))
    neff1 <- rep(0,length(params@w))
    for (i in 1:no_w)
        neff1[i] <- neff1[i] + sum(params@interaction[1,] * n[,i] * params@w[i] * params@dw[i])
    w_offset <- no_w_full - no_w
    pks <- rep(NA,length(params@w))
    pkpp <- rep(NA,length(params@w))
    for (i in 1:length(params@w)){
        pks[i] <- sum(params@pred_kernel[1,i,(w_offset+1):no_w_full] * neff1)
        pkpp[i] <- sum(n_pp * params@w_full * params@dw_full * params@pred_kernel[1,i,])
    }
    pp1 <- pks + pkpp
    pp <- getPhiPrey(params,n,n_pp)
    expect_that(pp1, is_equivalent_to(pp[1,]))
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

test_that("getPredRate",{
    data(NS_species_params_gears)
    data(inter)
    params <- MizerParams(NS_species_params_gears, inter)
    no_sp <- dim(params@catchability)[2]
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)
    n <- abs(array(rnorm(no_w * no_sp), dim = c(no_sp, no_w)))
    n_full <- abs(rnorm(no_w_full))
    pr <- getPredRate(params,n,n_full)
    # test dim
    expect_that(dim(pr), equals(c(no_sp,no_w,no_w_full)))
    # Look at numbers in predator 1
    n_total <- n[1,] * params@dw
    fl <- getFeedingLevel(params, n=n, n_pp=n_full)
    prr <- (1-fl[1,])*params@search_vol[1,]*n_total
    prr1 <- array(NA,dim=c(no_w,no_w_full))
    for(i in 1:no_w_full)
        prr1[,i] <- prr * params@pred_kernel[1,,i]
    expect_that(pr[1,,], is_equivalent_to(prr1))
    # Passing in feeding_level should yield the same
    pr1 <- getPredRate(params,n=n,n_pp=n_full)
    fl <- getFeedingLevel(params, n=n, n_pp=n_full)
    pr2 <- getPredRate(params,n=n,n_pp=n_full, feeding_level=fl)
    expect_that(pr1, is_identical_to(pr2))
    expect_that(getPredRate(params,n,n_full,feeding_level=matrix(rnorm(10*(no_sp-1)),ncol=10,nrow=no_sp-1)), throws_error())
})

test_that("getM2 for MizerParams",{
    data(NS_species_params_gears)
    data(inter)
    params <- MizerParams(NS_species_params_gears, inter)
    no_sp <- dim(params@catchability)[2]
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)
    n <- abs(array(rnorm(no_w * no_sp), dim = c(no_sp, no_w)))
    n_full <- abs(rnorm(no_w_full))
    m2 <- getM2(params,n,n_full)
    # test dim
    expect_that(dim(m2), equals(c(no_sp,no_w)))
    # Look at numbers in prey 1
    pr <- getPredRate(params,n,n_full)
    w_offset <- no_w_full - no_w
    m2_temp <- array(0,dim=c(no_sp,no_w_full))
    # sum over predator sizes to give total predation rate of each predator on each prey size
    for (i in 1:no_w)
        m2_temp <- m2_temp + pr[,i,]
    m22 <- rep(NA,no_w)
    for (i in 1:no_w)
        m22[i] <- sum(params@interaction[,1] * m2_temp[,w_offset+i])
    expect_that(m22, is_equivalent_to(m2[1,]))
    # Passing in pred_rate yields the same
    m21 <- getM2(params,n,n_full)
    pr <- getPredRate(params,n,n_full)
    m22 <- getM2(params,n,n_full, pred_rate=pr)
    expect_that(m21, is_identical_to(m22))
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
    expect_that(getM2(sim, time_range=time_range), equals(getM2(sim@params, sim@n[as.character(time_range),,], sim@n_pp[as.character(time_range),])))
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
    pr <- getPredRate(params,n,n_full)
    m22 <- rep(NA,no_w_full)
    for (i in 1:no_w_full)
        m22[i] <- sum(pr[,,i])
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
    expect_that(dim(getPredRate(params,n,n_pp)), equals(c(1,nw,length(params@w_full))))
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
