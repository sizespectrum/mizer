test_that("get_time_elements", {
    params <- NS_params_small
    sim <- project(params, effort = 1, t_max = 4, dt = 0.5, t_save = 0.5)
    expect_error(get_time_elements(sim, time_range = 0.1),
                 "The time range does not contain any simulation results.")
    expect_identical(get_time_elements(sim, as.character(3:4)),
                     get_time_elements(sim, 3:4))
    expect_identical(length(get_time_elements(sim, 3:4)),
                     dim(sim@n)[1])
    expect_identical(sum(get_time_elements(sim, 3:4)), 3L)
    expect_error(get_time_elements(sim, 3:50), 
                 "Time range is outside the time range of the model")
    expect_equal(which(get_time_elements(sim, seq(3, 4, by = 0.1))), 
                      c(7, 8, 9), ignore_attr = TRUE)
    # What if real years are used
    effort <- array(1, dim = c(19, 3),
                    dimnames = list(year = seq(1960, 1969, by = 0.5), 
                                    gear = c("Industrial", "Pelagic",
                                             "Otter")
                                    )
                    )
    sim <- project(params, effort = effort, t_save = 0.5)
    expect_equal(which(get_time_elements(sim, 1965)), 11, ignore_attr = TRUE)
    expect_equal(which(get_time_elements(sim, "1965")), 11, ignore_attr = TRUE)
    expect_equal(which(get_time_elements(sim, 1965:1969)), 11:19, ignore_attr = TRUE)
})

# get_sim_rate_time_elements ----

test_that("get_sim_rate_time_elements selects all saved times when time_range is missing", {
    sim <- NS_sim_small
    all_times <- get_sim_rate_time_elements(sim)
    expect_true(all(all_times))
    expect_length(all_times, dim(sim@n)[1])
})

test_that("get_sim_rate_time_elements delegates to get_time_elements", {
    sim <- NS_sim_small
    expect_identical(get_sim_rate_time_elements(sim, c(1, 2)),
                     get_time_elements(sim, c(1, 2)))
})

# get_sim_rate_slice ----

test_that("get_sim_rate_slice rebuilds the state at one saved time step", {
    sim <- NS_sim_small
    slice <- get_sim_rate_slice(sim, 2)
    expect_named(slice, c("n", "n_pp", "n_other", "effort", "t"))
    expect_equal(slice$n, sim@n[2, , ])
    expect_equal(slice$n_pp, sim@n_pp[2, ])
    expect_equal(slice$effort, sim@effort[2, ])
    expect_identical(slice$t, as.numeric(dimnames(sim@n)$time[[2]]))
})

test_that("get_sim_rate_slice keeps n a matrix for a one-species simulation", {
    sim <- project(NS_params_cod_small, t_max = 2, t_save = 1)
    slice <- get_sim_rate_slice(sim, 2)
    # the point of rebuilding `n` explicitly: it must not collapse to a vector
    expect_length(dim(slice$n), 2)
    expect_identical(dim(slice$n), c(1L, length(sim@params@w)))
    expect_identical(dimnames(slice$n), dimnames(sim@n)[2:3])
})

test_that("get_sim_rate_slice restores the names of n_other", {
    sim <- NS_sim_small
    slice <- get_sim_rate_slice(sim, 1)
    expect_identical(names(slice$n_other), dimnames(sim@n_other)$component)
})

# needed_rates ----

test_that("needed_rates returns a rate with no dependencies unchanged", {
    expect_identical(needed_rates("Encounter"), "Encounter")
})

test_that("needed_rates takes the transitive closure of the dependencies", {
    expect_setequal(needed_rates("EGrowth"),
                    c("Encounter", "FeedingLevel", "EReproAndGrowth",
                      "ERepro", "EGrowth"))
    expect_setequal(needed_rates("RDD"),
                    c("Encounter", "FeedingLevel", "EReproAndGrowth", "ERepro",
                      "EGrowth", "Diffusion", "PredRate", "PredMort", "FMort",
                      "Mort", "RDI", "RDD"))
})

test_that("needed_rates orders each rate after everything it depends on", {
    for (target in names(.rate_dependencies)) {
        need <- needed_rates(target)
        for (i in seq_along(need)) {
            deps <- .rate_dependencies[[need[[i]]]]
            expect_true(all(match(deps, need) < i),
                        info = paste(need[[i]], "in the closure of", target))
        }
    }
})

test_that("needed_rates lists each rate only once for repeated targets", {
    need <- needed_rates(c("EGrowth", "EGrowth", "ERepro"))
    expect_identical(anyDuplicated(need), 0L)
    expect_setequal(need, needed_rates("EGrowth"))
})

# mizer_rates_subset ----

test_that("mizer_rates_subset agrees with mizerRates on the requested targets", {
    params <- NS_params_small
    rates_fns <- projectRateFunctions(params)
    slice <- get_sim_rate_slice(NS_sim_small, 2)

    full <- mizerRates(params, n = slice$n, n_pp = slice$n_pp,
                       n_other = slice$n_other, t = slice$t,
                       effort = slice$effort, rates_fns = rates_fns)
    subset <- mizer_rates_subset(params, n = slice$n, n_pp = slice$n_pp,
                                 n_other = slice$n_other, t = slice$t,
                                 effort = slice$effort, rates_fns = rates_fns,
                                 targets = c("EGrowth", "Mort"))
    for (nm in names(subset)) {
        expect_equal(subset[[nm]], full[[nm]], ignore_attr = TRUE, info = nm)
    }
})

test_that("mizer_rates_subset skips the rates that are not needed", {
    params <- NS_params_small
    rates_fns <- projectRateFunctions(params)
    slice <- get_sim_rate_slice(NS_sim_small, 2)
    subset <- mizer_rates_subset(params, n = slice$n, n_pp = slice$n_pp,
                                 n_other = slice$n_other, t = slice$t,
                                 effort = slice$effort, rates_fns = rates_fns,
                                 targets = "FeedingLevel")
    # FeedingLevel needs only Encounter, so nothing downstream is calculated
    expect_setequal(names(subset), c("encounter", "feeding_level"))
})
