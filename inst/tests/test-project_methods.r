context("methods used in project")

test_that("getFmortGear",{
    data(species_params_gears)
    data(inter)
    params <- MizerParams(species_params_gears, inter)
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
    data(species_params_gears)
    data(inter)
    params <- MizerParams(species_params_gears, inter)
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
    data(species_params_gears)
    data(inter)
    params <- MizerParams(species_params_gears, inter)
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
    data(species_params_gears)
    data(inter)
    params <- MizerParams(species_params_gears, inter)
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
})

test_that("getFeedingLevel for MizerSim",{
    data(species_params_gears)
    data(inter)
    params <- MizerParams(species_params_gears, inter)
    sim <- project(params, effort=1, t_max=20, dt = 0.5, t_save = 0.5)
    time_range <- 15:20
    expect_that(length(dim(getFeedingLevel(sim, time_range=time_range))), equals(3))
    time_range <- 20
    expect_that(length(dim(getFeedingLevel(sim, time_range=time_range))), equals(2))
    expect_that(getFeedingLevel(sim, time_range=time_range), equals(getFeedingLevel(sim@params, sim@n[as.character(time_range),,], sim@n_pp[as.character(time_range),])))
})

test_that("getPredRate",{
    data(species_params_gears)
    data(inter)
    params <- MizerParams(species_params_gears, inter)
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
})

test_that("getM2",{
    data(species_params_gears)
    data(inter)
    params <- MizerParams(species_params_gears, inter)
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
})
test_that("getM2 for MizerSim",{
    data(species_params_gears)
    data(inter)
    params <- MizerParams(species_params_gears, inter)
    sim <- project(params, effort=1, t_max=20, dt = 0.5, t_save = 0.5)
    time_range <- 15:20
    expect_that(length(dim(getM2(sim, time_range=time_range))), equals(3))
    time_range <- 20
    expect_that(length(dim(getM2(sim, time_range=time_range))), equals(2))
    expect_that(getM2(sim, time_range=time_range), equals(getM2(sim@params, sim@n[as.character(time_range),,], sim@n_pp[as.character(time_range),])))
})


test_that("getM2Background",{
    data(species_params_gears)
    data(inter)
    params <- MizerParams(species_params_gears, inter)
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
})


test_that("getZ",{
    data(species_params_gears)
    data(inter)
    params <- MizerParams(species_params_gears, inter)
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
})

test_that("getEForReproAndGrowth",{
    data(species_params_gears)
    data(inter)
    params <- MizerParams(species_params_gears, inter)
    no_sp <- dim(params@catchability)[2]
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)
    n <- abs(array(rnorm(no_w * no_sp), dim = c(no_sp, no_w)))
    n_full <- abs(rnorm(no_w_full))
    erg <- getEForReproAndGrowth(params,n,n_full)
    expect_that(all(erg>=0), is_true())
    # test dim
    expect_that(dim(erg), equals(c(no_sp,no_w)))
    # Check number in final prey size group
    f <- getFeedingLevel(params, n=n, n_pp=n_full)
    e <- (f[1,] * params@intake_max[1,]) * params@species_params$alpha[1]
    e <- e - params@std_metab[1,] - params@activity[1,]
    e[e<0] <- 0 # Do not allow negative growth
    expect_that(e, is_equivalent_to(erg[1,]))
})


test_that("getESpawning",{
    data(species_params_gears)
    data(inter)
    params <- MizerParams(species_params_gears, inter)
    no_sp <- dim(params@catchability)[2]
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)
    n <- abs(array(rnorm(no_w * no_sp), dim = c(no_sp, no_w)))
    n_full <- abs(rnorm(no_w_full))
    es <- getESpawning(params,n,n_full)
    # test dim
    expect_that(dim(es), equals(c(no_sp,no_w)))
    e <- getEForReproAndGrowth(params,n=n,n_pp=n_full)
    e_spawning <- params@psi * e 
    expect_that(es, is_equivalent_to(e_spawning))
    e_growth <- getEGrowth(params,n,n_full)
    expect_that(e_growth, is_equivalent_to(e - es))
})


test_that("getRDI",{
    data(species_params_gears)
    data(inter)
    params <- MizerParams(species_params_gears, inter)
    no_sp <- dim(params@catchability)[2]
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)
    n <- abs(array(rnorm(no_w * no_sp), dim = c(no_sp, no_w)))
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
})

test_that("interaction is right way round in getM2 method",{
    data(species_params_gears)
    data(inter)
    inter[,"Dab"] <- 0 # Dab not eaten by anything
    params <- MizerParams(species_params_gears, inter)
    m2 <- getM2(params,getInitialN(params),params@cc_pp)
    expect_that(all(m2["Dab",] == 0), is_true())
})



