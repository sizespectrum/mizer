context("project method")

test_that("time dimension is dealt with properly",{
    data(species_params_gears)
    data(inter)
    params <- MizerParams(species_params_gears, inter)
    no_gear <- dim(params@catchability)[1]
    no_sp <- dim(params@catchability)[2]
    max_t_effort <- 10
    effort <- array(abs(rnorm(max_t_effort*no_gear)),dim=c(max_t_effort,no_gear))
    expect_error(project(params,effort=1,t_save=3,dt=2))
    expect_error(project(params,effort=1,t_max=7,dt=2))
    sim <- project(params,t_max=10,t_save=2)
    expect_that(dim(sim@n)[1], equals(6))
    expect_that(dim(sim@effort)[1], equals(5))
    sim <- project(params,t_max=10,t_save=2, dt=0.5)
    expect_that(dim(sim@n)[1], equals(6))
    expect_that(dim(sim@effort)[1], equals(5))
    sim <- project(params,t_max=5,t_save=0.5, dt=0.5)
    expect_that(dim(sim@n)[1], equals(11))
    expect_that(dim(sim@effort)[1], equals(10))
})

test_that("Can pass in initial species",{
    data(species_params_gears)
    data(inter)
    params <- MizerParams(species_params_gears, inter)
    no_gear <- dim(params@catchability)[1]
    no_sp <- dim(params@catchability)[2]
    max_t_effort <- 10
    effort <- array(abs(rnorm(max_t_effort*no_gear)),dim=c(max_t_effort,no_gear))


})

test_that("get_initial_n is working properly",{
    data(species_params_gears)
    data(inter)
    params <- MizerParams(species_params_gears, inter)
    n <- get_initial_n(params)
    no_sp <- nrow(params@species_params)
    for(i in 1:no_sp){
        expect_that(all(n[i,params@w > params@species_params$w_inf[i]] == 0), is_true())
        expect_that(all(n[i,params@w < params@species_params$w_min[i]] == 0), is_true())
    }
    # Check slope of all species is the same
    slopes <- rep(NA, no_sp)
    for(i in 1:no_sp){
        n_idx <- which(n[i,] != 0)
        slopes[i] <- (log(n[i,min(n_idx)]) - log(n[i,max(n_idx)])) / (log(params@w[min(n_idx)]) - log(params@w[max(n_idx)]))
    }
    expect_that(slopes, equals(rep(slopes[1],no_sp)))
    # Check that slopes = slope0
})

test_that("w_min array reference is working OK",{
    data(species_params_gears)
    data(inter)
    species_params_gears$w_min <- 0.001
    species_params_gears$w_min[1] <- 1
    params <- MizerParams(species_params_gears, inter)
    sim <- project(params, effort=1, t_max=5)
    expect_that(all(sim@n[6,1,1:(sim@params@species_params$w_min_idx[1]-1)]==0),is_true())
})


