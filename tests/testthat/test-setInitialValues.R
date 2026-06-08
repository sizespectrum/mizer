params <- NS_params
sim <- project(params, t_max = 0.1, t_save = 0.1, effort = 1)

test_that("setInitialValues is deprecated", {
    expect_warning(setInitialValues(params, sim),
                   "deprecated")
})

test_that("We can set and get initial values from sim object", {
    no_t <- dim(sim@effort)[[1]]
    gear_names <- names(params@initial_effort)
    expect_identical(initialN(sim), initialN(params))
    expect_identical(initialNResource(sim), initialNResource(params))
    initialN(params) <- params@metab
    expect_equal(initialN(params), params@metab, ignore_attr = TRUE)
    initialNResource(params) <- params@cc_pp
    expect_identical(initialNResource(params), params@cc_pp)
    params <- suppressWarnings(setInitialValues(params, sim))
    expect_equal(finalN(sim), initialN(params), ignore_attr = TRUE)
    expect_identical(finalNResource(sim), initialNResource(params))
    expect_identical(sim@effort[no_t, ], params@initial_effort)
    expect_named(params@initial_effort, gear_names)
    names(params@initial_effort) <- NULL
    expect_error(suppressWarnings(setInitialValues(params, sim)),
                 "The gears in the simulation in `sim` have different names")
    params@initial_effort <- 1
    expect_error(suppressWarnings(setInitialValues(params, sim)),
                 "The number of gears in the simulation in `sim` is different")
})

test_that("initialN and initialNResource setters validate dimnames and values", {
    params <- NS_params
    new_n <- initialN(params)
    dimnames(new_n)[[1]][1] <- "wrong"
    expect_warning(initialN(params) <- new_n,
                   "The dimnames do not match. I will ignore them.")
    expect_equal(initialN(params), NS_params@initial_n, ignore_attr = TRUE)

    new_n_pp <- initialNResource(params)
    names(new_n_pp)[1] <- "wrong"
    initialNResource(params) <- new_n_pp
    expect_identical(initialNResource(params), NS_params@initial_n_pp)

    expect_error(initialN(params) <- (-1) * initialN(params))
    expect_error(initialNResource(params) <- (-1) * initialNResource(params))
})

test_that("setInitialValues gives correct errors", {
    params1 <- newMultispeciesParams(NS_species_params, no_w = 20, info_level = 0)
    sim <- project(params1, t_max = 2, dt = 1)
    params2 <- newMultispeciesParams(NS_species_params, no_w = 30, info_level = 0)
    expect_error(suppressWarnings(setInitialValues(params2, sim)),
                 "The consumer size spectrum of the simulation in `sim` has a different size from that in `params`")

    params3 <- newMultispeciesParams(NS_species_params, no_w = 20,
                                     min_w_pp = 1e-4, info_level = 0)
    expect_error(suppressWarnings(setInitialValues(params3, sim)),
                 "The resource size spectrum of the simulation in `sim` has a different size from that in `params`.")
    params4 <- setComponent(params1, "test",
                            initial_value = 0, dynamics_fun = "sum")
    expect_error(suppressWarnings(setInitialValues(params4, sim)),
                 "The number of other components in the simulation in `sim` is different from that in `params`.")
    gear_params(params1)$gear[1] <- "test"
    expect_error(suppressWarnings(setInitialValues(params1, sim)),
                 "The number of gears in the simulation in `sim` is different from that in `params`.")
    gear_params(params1)$gear <- "test"
    expect_error(suppressWarnings(setInitialValues(params1, sim)),
                 "The gears in the simulation in `sim` have different names from those in `params`.")

})

test_that("Can set initial values in a model with a single species", {
    species_params <- NS_species_params[1, ]
    params <- newMultispeciesParams(species_params, info_level = 0)
    sim <- project(params, t_max = 0.1, t_save = 0.1)
    p <- suppressWarnings(setInitialValues(params, sim))
    expect_equal(finalN(sim), initialN(p), ignore_attr = TRUE)
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
    p <- suppressWarnings(setInitialValues(params, sim))
    expect_identical(initialNOther(p), list(test = 111))
})

test_that("setInitialValues averages correctly over time range", {
    time_sel <- c(2:4)
    time_range <- getTimes(NS_sim)[time_sel]
    # arithmetic mean
    params <- suppressWarnings(
        setInitialValues(NS_sim@params, NS_sim, time_range = time_range))
    expected_n <- mean(NS_sim@n[time_sel, 1, 10])
    expect_equal(params@initial_n[1, 10], expected_n)
    # geometric mean
    params <- suppressWarnings(
        setInitialValues(NS_sim@params, NS_sim, time_range = time_range,
                         geometric_mean = TRUE))
    expected_n <- exp(mean(log(NS_sim@n[time_sel, 1, 10])))
    expect_equal(params@initial_n[1, 10], expected_n)
})

test_that("setInitialValues updates `time_modified`", {
    params2 <- suppressWarnings(setInitialValues(params, sim))
    expect_false(identical(params2@time_modified, params@time_modified))
})

test_that("initialN<- updates `time_modified`", {
    p <- params
    initialN(p) <- params@initial_n
    expect_false(identical(p@time_modified, params@time_modified))
})

test_that("initialNResource<- updates `time_modified`", {
    p <- params
    initialNResource(p) <- params@initial_n_pp
    expect_false(identical(p@time_modified, params@time_modified))
})
