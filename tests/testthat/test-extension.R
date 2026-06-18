params <- NS_params_small
e <- globalenv()
e$test_dyn <- function(params, ...) {
    111
}
e$test_egrowth <- function(params, ...) {
    array(111, dim = c(nrow(params@species_params), length(params@w)))
}
e$nt <- function(params, t, ...) {
    params@initial_n * t
}
e$resource_encounter <- function(params, n, n_pp, n_other, ...) {
    mizerEncounter(params, n = n, n_pp = n_other$resource, ...)
}
e$semichemostat <- function(params, n_other, rates, dt, component, ...) {
    c <- params@other_params[[component]]
    interaction <- params@species_params$interaction_resource
    mort <- as.vector(interaction  %*% rates$pred_rate)
    tmp <- c$rate * c$capacity / (c$rate + mort)
    return(tmp - (tmp - n_other[[component]]) * exp(-(c$rate + mort) * dt))
}

# setRateFunction works ----
test_that("setRateFunction works", {
    expect_error(setRateFunction(params, rate = "wrong", fun = "sum"),
               "The `rate` argument must be one of")
    xxx <- "xx"
    expect_error(setRateFunction(params, rate = "Mort", fun = "xxx"),
                 "There is no function ")
    p <- setRateFunction(params, "Mort", "mizerMort")
    expect_identical(params@rates_funcs, p@rates_funcs)
    p <- setRateFunction(params, rate = "EGrowth", fun = "test_egrowth")
    expect_identical(p@rates_funcs[["EGrowth"]], "test_egrowth")
    r <- mizerRates(p, n = initialN(p),
                    n_pp = initialNResource(p),
                    n_other = initialNOther(p),
                    effort = 0,
                    rates_fns = lapply(p@rates_funcs, get))
    expect_true(all(r$e_growth == 111))
})

test_that("setRateFunction validates RDI functions with diffusion", {
    assign("rdi_needs_diffusion",
           function(params, diffusion, ...) {
               rowSums(diffusion)
           }, envir = .GlobalEnv)
    withr::defer(rm("rdi_needs_diffusion", envir = .GlobalEnv))

    p <- setRateFunction(params, "RDI", "rdi_needs_diffusion")

    expect_identical(p@rates_funcs[["RDI"]], "rdi_needs_diffusion")
    expect_equal(getRDI(p), rowSums(getDiffusion(p)), ignore_attr = TRUE)
})

test_that("Time is passed correctly to rate functions", {
    params@rates_funcs$Encounter <- "nt"
    expect_equal(getEncounter(params, t = 2), nt(params, 2),
                 ignore_attr = TRUE)
    params@rates_funcs$FeedingLevel <- "nt"
    expect_equal(getFeedingLevel(params, time_range = 2), nt(params, 2),
                 ignore_attr = TRUE)
    
    gears <- unique(gear_params(params)$gear)
    effort <- array(0, dim = c(3, length(gears)),
                    dimnames = list(time = 2020:2022,
                                    gear = gears))
    sim <- project(params, effort = effort, dt = 1)
    expect_equal(getFeedingLevel(sim, time_range = 2021:2022)[1, , ],
                 nt(params, 2021), ignore_attr = TRUE)
    #TODO: extend this
})

# getRateFunction works ----
test_that("getRateFunction works", {
    expect_error(getRateFunction(params, "test"),
        "The `rate` argument must be one of")
    all <- getRateFunction(params)
    expect_type(all, "list")
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
    other_params(params) <- list(test = 5, other = 6)
    expect_identical(other_params(params), list(test = 5, other = 6))
    params2 <- setComponent(NS_params_small, "component", 1,
                            dynamics_fun = "test_dyn",
                            component_params = list(a = 2))
    expect_null(other_params(params2))
    expect_identical(params2@other_params$component, list(a = 2))
    other_params(params) <- structure(list(), names = character())
    expect_length(other_params(params), 0)
    expect_null(other_params(params)$test)
})

# components ----
test_that("We can set, get and remove components", {
    expect_error(setComponent(params, "test", 1),
                 '"dynamics_fun" is missing')
    p <- setComponent(params, "test", 1, 
                      dynamics_fun = "test_dyn",
                      encounter_fun = "test_dyn",
                      mort_fun = "test_dyn",
                      component_params = list(a = 2))
    expect_identical(p@other_dynamics, list(test = "test_dyn"))
    expect_identical(p@other_encounter, list(test = "test_dyn"))
    expect_identical(p@other_mort, list(test = "test_dyn"))
    expect_identical(p@other_params, list(test = list(a = 2)))
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
    p3 <- setComponent(p, "test2", 3, "test_dyn")
    expect_identical(getComponent(p3, "test2")$initial_value, 3)
    expect_null(getComponent(p3, "test2")$mort_fun)
    expect_null(getComponent(p3, "test2")$encounter_fun)
    expect_error(removeComponent(p2, "test3"),
                 "There is no component named test3")
    expect_null(getComponent(p2, "test3"))
    p1 <- removeComponent(p2, "test")
    d <- getComponent(p1, "test2")
    expect_length(p1@other_dynamics, 1)
    expect_length(p1@other_encounter, 0)
    expect_null(p1@initial_n_other[["test"]])
    expect_null(p1@other_params[["test"]])
})

test_that("We can set component colours and linetypes", {
    p <- setComponent(params, "test", 1,
                      dynamics_fun = "test_dyn",
                      colour = "orange",
                      linetype = "dashed")
    expect_identical(getColours(p)[["test"]], "orange")
    expect_identical(getLinetypes(p)[["test"]], "dashed")

    p2 <- suppressWarnings(
        setComponent(p, "test2", 1,
                     dynamics_fun = "test_dyn",
                     colour = "not-a-colour",
                     linetype = "not-a-linetype")
    )
    expect_false("test2" %in% names(getColours(p2)))
    expect_false("test2" %in% names(getLinetypes(p2)))

    p3 <- setComponent(params, "test3", 1,
                       dynamics_fun = "test_dyn")
    expect_identical(getColours(p3)[["test3"]], "grey")
    expect_identical(getLinetypes(p3)[["test3"]], "solid")
})

# initial values ----
test_that("We can set and get initial values for MizerParams", {
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
    # test that we can get initial values from MizerSim object
    sim <- project(p, t_max = 0.2, t_save = 0.1)
    expect_identical(initialNOther(sim)$test, 4)
})

# encounter and mortality functions are called ----
test_that("encounter and mortality functions are called", {
    e <- getEncounter(params)
    m <- getMort(params)
    p <- setComponent(params, "test", 1, 
                      dynamics_fun = "test_dyn",
                      encounter_fun = "test_dyn",
                      mort_fun = "test_dyn")
    expect_equal(getEncounter(p), e + 111, ignore_attr = TRUE)
    p <- setComponent(params, "test", 1, 
                      dynamics_fun = "test_dyn",
                      mort_fun = "test_dyn")
    expect_equal(getMort(p), m + 111, ignore_attr = TRUE)
})

test_that("We can access simulation results", {
    p <- setComponent(params, "test", 1, 
                      dynamics_fun = "test_dyn",
                      encounter_fun = "test_dyn",
                      mort_fun = "test_dyn")
    sim <- project(p, t_max = 0.2, t_save = 0.1)
    expect_identical(finalNOther(sim), list("test" = 111))
    expect_identical(NOther(sim)[2, ], list(111))
    expect_identical(NOther(sim)[1, ], list(1))
    expect_identical(dimnames(NOther(sim))$component, "test")
})

test_that("component can mimic resource", {
    params <- NS_params_small
    initialNResource(params) <- initialNResource(params) / 10
    component_params <- list(capacity = params@cc_pp,
                             rate = params@rr_pp)
    params2 <- params %>% 
        setComponent("resource",
                     initial_value = initialNResource(params),
                     dynamics_fun = "semichemostat",
                     component_params = component_params) %>% 
        setRateFunction("Encounter", "resource_encounter")
    sim <- project(params, t_max = 0.1, t_save = 0.1, effort = 0)
    sim2 <- project(params2, t_max = 0.1, t_save = 0.1, effort = 0)
    expect_equal(finalNResource(sim2), finalNOther(sim2)$resource,
                 ignore_attr = TRUE)
    expect_equal(finalNResource(sim), finalNResource(sim2),
                 ignore_attr = TRUE)
    expect_equal(finalN(sim), finalN(sim2), ignore_attr = "params")
})

# time_modified ----
test_that("setRateFunction updates `time_modified`", {
    p <- setRateFunction(params, "Encounter", "mizerEncounter")
    expect_false(identical(p@time_modified, params@time_modified))
})

test_that("other_params<- updates `time_modified`", {
    p <- params
    other_params(p) <- list(test = 1)
    expect_false(identical(p@time_modified, params@time_modified))
})

test_that("setComponent updates `time_modified`", {
    p <- setComponent(params, "test", 1, dynamics_fun = "test_dyn")
    expect_false(identical(p@time_modified, params@time_modified))
})

test_that("removeComponent updates `time_modified`", {
    p <- setComponent(params, "test", 1, dynamics_fun = "test_dyn")
    p2 <- removeComponent(p, "test")
    expect_false(identical(p2@time_modified, p@time_modified))
})

test_that("initialNOther<- updates `time_modified`", {
    p <- setComponent(params, "test", 1, dynamics_fun = "test_dyn")
    p2 <- p
    initialNOther(p2)$test <- 2
    expect_false(identical(p2@time_modified, p@time_modified))
})

# .MizerSim rate getters preserve n_other component names ----
e$test_encounter_n_other <- function(params, n, n_pp, n_other, t = 0, ...) {
    if (is.null(n_other[["my_comp"]])) stop("n_other component names were lost!")
    mizerEncounter(params, n, n_pp, n_other, t, ...)
}
e$test_predrate_n_other <- function(params, n, n_pp, n_other, t = 0, ...) {
    if (is.null(n_other[["my_comp"]])) stop("n_other component names were lost!")
    mizerPredRate(params, n, n_pp, n_other, t, ...)
}
e$test_fmort_n_other <- function(params, n, n_pp, n_other, effort, t, ...) {
    if (is.null(n_other[["my_comp"]])) stop("n_other component names were lost!")
    mizerFMort(params, n, n_pp, n_other, effort, t, ...)
}

test_that("getFeedingLevel.MizerSim preserves n_other component names with time_range", {
    p <- NS_params_small |>
        setComponent("my_comp", 1, dynamics_fun = "test_dyn") |>
        setRateFunction("Encounter", "test_encounter_n_other")
    sim <- project(p, t_max = 3, dt = 1, t_save = 1)
    expect_no_error(getFeedingLevel(sim, time_range = 2:3))
})

test_that("getPredMort.MizerSim preserves n_other component names with time_range", {
    p <- NS_params_small |>
        setComponent("my_comp", 1, dynamics_fun = "test_dyn") |>
        setRateFunction("PredRate", "test_predrate_n_other")
    sim <- project(p, t_max = 3, dt = 1, t_save = 1)
    expect_no_error(getPredMort(sim, time_range = 2:3))
})

test_that("getFMort.MizerSim preserves n_other component names with time_range", {
    p <- NS_params_small |>
        setComponent("my_comp", 1, dynamics_fun = "test_dyn") |>
        setRateFunction("FMort", "test_fmort_n_other")
    sim <- project(p, t_max = 3, dt = 1, t_save = 1)
    expect_no_error(getFMort(sim, time_range = 2:3))
})

test_that("Other components are advanced once per time step, not per save step", {
    # A component whose dynamics add `dt` at every time step. After projecting
    # for one year it must equal 1 regardless of how often results are saved.
    e$ticker_dyn <- function(params, n_other, component, dt, ...) {
        n_other[[component]] + dt
    }
    p <- setComponent(NS_params_small, "ticker", initial_value = 0,
                      dynamics_fun = "ticker_dyn")

    sim <- project(p, t_max = 2, t_save = 1, dt = 0.1, progress_bar = FALSE)
    expect_equal(sim@n_other[[2, "ticker"]], 1)
    expect_equal(sim@n_other[[3, "ticker"]], 2)

    # The result at a given time must not depend on the save frequency.
    sim_fine <- project(p, t_max = 1, t_save = 0.1, dt = 0.1,
                        progress_bar = FALSE)
    expect_equal(sim_fine@n_other[[11, "ticker"]], 1)
})

test_that("Other components get a corrector step under second-order methods", {
    # A component that integrates a time-varying, rate-derived signal. Because
    # its dynamics read `rates`, the corrector (midpoint rates) differs from
    # the predictor (start-of-step rates), so the second-order methods must be
    # more accurate than Euler. This guards against the other components being
    # left at first-order accuracy while the resource and consumers are
    # advanced to second order.
    e$accum_dyn <- function(params, n_other, component, rates, dt, ...) {
        n_other[[component]] + dt * sum(rates$mort)
    }
    make_p <- function() {
        p <- NS_params_small
        p@initial_n[] <- p@initial_n * 1.5   # perturb to create a transient
        setComponent(p, "accum", initial_value = 0, dynamics_fun = "accum_dyn")
    }
    accum <- function(method, dt) {
        project(make_p(), t_max = 2, t_save = 2, dt = dt, method = method,
                progress_bar = FALSE)@n_other[[2, "accum"]]
    }
    ref <- accum("euler", 0.0005)            # fine-step reference

    euler_err <- abs(accum("euler", 0.1) - ref)
    for (method in c("predictor_corrector", "tr_bdf2")) {
        # Generous margin: in practice the second-order error is several times
        # smaller, but we only need to confirm the corrector is taking effect.
        expect_lt(abs(accum(method, 0.1) - ref), euler_err / 2)
    }
})
