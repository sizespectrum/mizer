params <- NS_params
test_dyn <- function(params, ...) {
    111
}

# setRateFunction works ----
test_that("setRateFunction works", {
    expect_error(setRateFunction(params, rate = "wrong", fun = "sum"),
               "The `rate` argument must be one of")
    xxx <- "xx"
    expect_error(setRateFunction(params, rate = "Mort", fun = "xxx"),
                 "`fun` should be a function")
    p <- setRateFunction(params, "Mort", "mizerMort")
    expect_identical(params@rates_funcs, p@rates_funcs)
    p <- setRateFunction(params, rate = "EGrowth", fun = "test_dyn")
    expect_identical(p@rates_funcs[["EGrowth"]], "test_dyn")
    r <- mizerRates(p, n = initialN(p), 
                    n_pp = initialNResource(p), 
                    n_other = initialNOther(p),
                    effort = 0,
                    rates_fns = lapply(p@rates_funcs, get))
    expect_identical(r$e_growth, test_dyn(params))
})

# getRateFunction works ----
test_that("getRateFunction works", {
    expect_error(getRateFunction(params, "test"),
        "The `rate` argument must be one of")
    all <- getRateFunction(params)
    expect_type(all , "list")
    expect_named(all)
    expect_identical(all$Rates, "mizerRates")
    expect_identical(getRateFunction(params, rate = "Mort"), "mizerMort")
})

# other_params ----
test_that("We can set and get other params", {
    expect_length(other_params(params), 0)
    expect_error(other_params(params) <- 5,
                 "other_params should be a named list")
    expect_error(other_params(params) <- list(5),
                 "other_params should be a named list")
    other_params(params)$test <- 5
    expect_identical(other_params(params)$test, 5)
    other_params(params)$test <- NULL
    expect_length(other_params(params), 0)
    expect_null(other_params(params)$test)
})

# components ----
test_that("We can set, get and remove components", {
    expect_error(setComponent(params, "test", 1),
                 '"dynamics_fun" is missing')
    expect_error(setComponent(params, "test", 1, 
                      dynamics_fun = "test_dyn",
                      encounter_fun = "test_dyn",
                      mort_fun = "test_dyn",
                      component_params = list(2, 3)),
                 "`component_params` needs to be NULL or a named list.")
    p <- setComponent(params, "test", 1, 
                      dynamics_fun = "test_dyn",
                      encounter_fun = "test_dyn",
                      mort_fun = "test_dyn",
                      component_params = list(a = 2))
    comp <- getComponent(p, "test")
    expect_mapequal(comp, list(initial_value = 1,
                               encounter_fun = "test_dyn",
                               dynamics_fun = "test_dyn",
                               mort_fun = "test_dyn",
                               component_params = list(a = 2)))
    all <- getComponent(p)
    expect_type(all, "list")
    expect_identical(all$test, comp)
    p2 <- setComponent(p, "test2", 2, 
                      dynamics_fun = "test_dyn",
                      mort_fun = "test_dyn")
    all2 <- getComponent(p2)
    expect_length(all2, 2)
    expect_length(all2$test2, 5)
    expect_null(all2$test2$encounter_fun)
    p <- setComponent(p2, "test2", 1, "test_dyn", mort_fun = NULL)
    expect_null(getComponent(p, "test2")$mort_fun)
    expect_error(removeComponent(p2, "test3"),
                 "There is no component named test3")
    expect_error(getComponent(p2, "test3"),
                 "There is no component named test3")
    p1 <- removeComponent(p2, "test")
    d <- getComponent(p1, "test2")
    expect_length(p1@other_dynamics, 1)
    expect_length(p1@other_encounter, 0)
})

# initial values ----
test_that("We can set and get initial values", {
    p <- setComponent(params, "test", 1, 
                      dynamics_fun = "test_dyn")
    expect_identical(initialNOther(p), list(test = 1))
    p <- setComponent(p, "test", list(a = 1, b = 2), 
                      dynamics_fun = "test_dyn")
    expect_identical(initialNOther(p), list(test = list(a = 1, b = 2)))
    initialNOther(p)$test <- 3
    expect_identical(initialNOther(p), list(test = 3))
    expect_error(initialNOther(p)$test2 <- 3,
                 "The following components do not exist: test2")
    p <- setComponent(p, "test2", 2, 
                      dynamics_fun = "test_dyn")
    expect_identical(initialNOther(p), list(test = 3, test2 = 2))
    expect_error(initialNOther(p) <- list(test = 4),
                 "Missing values for components test2")
    initialNOther(p)$test <- 4
    expect_identical(initialNOther(p)$test, 4)
})

# encounter and mortality functions are called ----
test_that("encounter and mortality functions are called", {
    e <- getEncounter(params)
    m <- getMort(params)
    p <- setComponent(params, "test", 1, 
                      dynamics_fun = "test_dyn",
                      encounter_fun = "test_dyn",
                      mort_fun = "test_dyn")
    expect_identical(getEncounter(p), e + test_dyn(params))
    p <- setComponent(params, "test", 1, 
                      dynamics_fun = "test_dyn",
                      mort_fun = "test_dyn")
    expect_identical(getMort(p), m + test_dyn(params))
})

test_that("We can access simulation results", {
    p <- setComponent(params, "test", 1, 
                      dynamics_fun = "test_dyn",
                      encounter_fun = "test_dyn",
                      mort_fun = "test_dyn")
    sim <- project(p, t_max = 0.2, t_save = 0.1)
    expect_identical(finalNOther(sim), list("test" = test_dyn(p)))
    expect_identical(NOther(sim)[2, ], list(test_dyn(p)))
    expect_identical(NOther(sim)[1, ], list(1))
})
