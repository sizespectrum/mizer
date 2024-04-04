params <- NS_params
sim <- project(params, t_max = 0.1, t_save = 0.1, effort = 1)

test_that("We can set and get initial values from sim object", {
    no_t <- dim(sim@effort)[[1]]
    gear_names <- names(params@initial_effort)
    expect_identical(initialN(sim), initialN(params))
    expect_identical(initialNResource(sim), initialNResource(params))
    initialN(params) <- params@metab
    expect_identical(initialN(params), params@metab)
    initialNResource(params) <- params@cc_pp
    expect_identical(initialNResource(params), params@cc_pp)
    params <- setInitialValues(params, sim)
    expect_identical(finalN(sim), initialN(params))
    expect_identical(finalNResource(sim), initialNResource(params))
    expect_identical(sim@effort[no_t, ], params@initial_effort)
    expect_named(params@initial_effort, gear_names)
    names(params@initial_effort) <- NULL
    expect_error(setInitialValues(params, sim),
                 "The gears in the simulation in `sim` have different names")
    params@initial_effort <- 1
    expect_error(setInitialValues(params, sim),
                 "The number of gears in the simulation in `sim` is different")
})

test_that("setInitialValues gives correct errors", {
    params1 <- newMultispeciesParams(NS_species_params, no_w = 20, info_level = 0)
    sim <- project(params1, t_max = 2, dt = 1)
    params2 <- newMultispeciesParams(NS_species_params, no_w = 30, info_level = 0)
    expect_error(setInitialValues(params2, sim),
                 "The consumer size spectrum of the simulation in `sim` has a different size from that in `params`")
    
    params3 <- newMultispeciesParams(NS_species_params, no_w = 20,
                                     min_w_pp = 1e-4, info_level = 0)
    expect_error(setInitialValues(params3, sim),
                 "The resource size spectrum of the simulation in `sim` has a different size from that in `params`.")
    params4 <- setComponent(params1, "test",
                            initial_value = 0, dynamics_fun = "sum")
    expect_error(setInitialValues(params4, sim),
                 "The number of other components in the simulation in `sim` is different from that in `params`.")
    gear_params(params1)$gear[1] <- "test"
    expect_error(setInitialValues(params1, sim),
                 "The number of gears in the simulation in `sim` is different from that in `params`.")
    gear_params(params1)$gear <- "test"
    expect_error(setInitialValues(params1, sim),
                 "The gears in the simulation in `sim` have different names from those in `params`.")
    
})

test_that("Can set initial values in a model with a single species", {
    species_params <- NS_species_params[1, ]
    params <- newMultispeciesParams(species_params, info_level = 0)
    sim <- project(params, t_max = 0.1, t_save = 0.1)
    p <- setInitialValues(params, sim)
    expect_identical(finalN(sim), initialN(p))
})

test_that("Can set initial values in a model with a single other component", {
    e <- globalenv()
    e$test_dyn <- function(params, ...) {
        111
    }
    params <- setComponent(params, 
                           component = "test",
                           initial_value = 1,
                           dynamics_fun = "test_dyn")
    sim <- project(params, t_max = 0.1, t_save = 0.1)
    p <- setInitialValues(params, sim)
    expect_identical(initialNOther(p), list(test = 111))
})

test_that("setInitialValues averages correctly over time range", {
    time_sel <- c(24:33)
    time_range <- getTimes(NS_sim)[time_sel]
    # arithmetic mean
    params <- setInitialValues(NS_sim@params, NS_sim, time_range = time_range)
    expected_n <- mean(NS_sim@n[time_sel, 1, 50])
    expect_equal(params@initial_n[1, 50], expected_n)
    # geometric mean
    params <- setInitialValues(NS_sim@params, NS_sim, time_range = time_range,
                               geometric_mean = TRUE)
    expected_n <- exp(mean(log(NS_sim@n[time_sel, 1, 50])))
    expect_equal(params@initial_n[1, 50], expected_n)
})
