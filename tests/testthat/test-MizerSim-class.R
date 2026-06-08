params <- NS_params_default_small
short_ns_sim <- project(NS_params_small, t_max = 0.2, t_save = 0.1)

# basic constructor sets dimensions properly ----
test_that("basic constructor sets dimensions properly", {
    # check dimension against t input arguments
    # Make MizerSims with t_max and t_save
    t_max <- 5
    t_save <- 1
    sim <- MizerSim(params, t_max = t_max, t_save = t_save)
    expect_equal(dim(sim@effort)[1], 1 + t_max / t_save)
    expect_identical(dimnames(sim@effort)[[1]], 
                     as.character(seq(from = 0, to = t_max, by = t_save)))
    expect_setequal(dimnames(sim@effort)[[2]], 
                    dimnames(params@selectivity)$gear)
    expect_equal(dim(sim@n)[1], 1 + (t_max / t_save))
    expect_identical(dimnames(sim@n)[[1]], 
                     dimnames(sim@effort)[[1]])

    t_max <- 4
    t_save <- 2
    sim <- MizerSim(params, t_max = t_max, t_save = t_save)
    expect_equal(dim(sim@effort)[1], 1 + t_max / t_save)
    expect_identical(dimnames(sim@effort)[[1]], 
                     as.character(seq(from = 0, to = t_max, by = t_save)))
    expect_setequal(dimnames(sim@effort)[[2]], 
                    dimnames(params@selectivity)$gear)
    expect_equal(dim(sim@n)[1], 1 + (t_max / t_save))
    expect_identical(dimnames(sim@n)[[1]], 
                     dimnames(sim@effort)[[1]])

    # Make MizerSim using t_dimnames
    t_dimnames <- seq(from = 1990, to = 2000, by = 1)
    sim <- MizerSim(params, t_dimnames = t_dimnames)
    expect_equal(dim(sim@effort)[1], length(t_dimnames))
    expect_identical(dimnames(sim@effort)[[1]], as.character(t_dimnames))
    expect_setequal(dimnames(sim@effort)[[2]],
                    dimnames(params@selectivity)$gear)
    expect_equal(dim(sim@n)[1], length(t_dimnames))
    expect_identical(dimnames(sim@n)[[1]], as.character(t_dimnames))

    # Check error if t_dimnames is not numeric or not sorted
    expect_error(MizerSim(params, t_dimnames = c("x", "y", "z")),
                 "The t_dimnames argument must be numeric")
    expect_error(MizerSim(params, t_dimnames = as.character(1:3)),
                 "The t_dimnames argument must be numeric")
    expect_error(MizerSim(params, t_dimnames = 3:1),
                 "The t_dimnames argument should be increasing")
})

test_that("validSim works", {
    sim <- short_ns_sim
    sim@n[3, 1, 1] <- Inf
    expect_warning(simt <- validSim(sim),
                   "The simulation failed to work beyond time = 0.1")
    expect_equal(dim(simt@n), c(2, nrow(NS_params_small@species_params), length(NS_params_small@w)))
    expect_equal(dim(simt@n_pp), c(2, length(NS_params_small@w_full)))
    expect_equal(dim(simt@effort), c(2, length(unique(gear_params(NS_params_small)$gear))))
    sim@n[2, 2, 2] <- NaN
    expect_warning(simt <- validSim(sim),
                   "The simulation failed to work beyond time = 0")
    expect_equal(dim(simt@n), c(1, nrow(NS_params_small@species_params), length(NS_params_small@w)))
})

test_that("validSim also validates embedded params", {
    sim <- short_ns_sim
    sim@params@species_params$w_min[1] <- 1e-10
    expect_warning(sim2 <- validSim(sim), "smaller than the minimum")
    expect_equal(sim2@params@w_min_idx[[1]], 1)
})

# getParams, finalParams, initialParams ----

test_that("getParams() returns final time step by default", {
    sim <- short_ns_sim
    p <- getParams(sim)
    expect_s4_class(p, "MizerParams")
    no_t <- dim(sim@n)[1]
    expect_equal(p@initial_n, sim@n[no_t, , ], ignore_attr = TRUE)
    expect_equal(p@initial_n_pp, sim@n_pp[no_t, ], ignore_attr = TRUE)
    expect_equal(p@initial_effort, sim@effort[no_t, ], ignore_attr = TRUE)
})

test_that("getParams() with a single time_range selects that time step", {
    t <- getTimes(NS_sim_small)[3]
    p <- getParams(NS_sim_small, time_range = t)
    idx <- which(as.numeric(dimnames(NS_sim_small@n)$time) == t)
    expect_equal(p@initial_n, NS_sim_small@n[idx, , ], ignore_attr = TRUE)
})

test_that("getParams() averages arithmetic mean over time range", {
    time_sel <- c(2:4)
    time_range <- getTimes(NS_sim_small)[time_sel]
    p <- getParams(NS_sim_small, time_range = time_range)
    expect_equal(p@initial_n[1, 10], mean(NS_sim_small@n[time_sel, 1, 10]))
})

test_that("getParams() averages geometric mean over time range", {
    time_sel <- c(2:4)
    time_range <- getTimes(NS_sim_small)[time_sel]
    p <- getParams(NS_sim_small, time_range = time_range, geometric_mean = TRUE)
    expect_equal(p@initial_n[1, 10],
                 exp(mean(log(NS_sim_small@n[time_sel, 1, 10]))))
})

test_that("getParams() updates time_modified", {
    p <- getParams(short_ns_sim)
    expect_false(identical(p@time_modified,
                           short_ns_sim@params@time_modified))
})

test_that("finalParams() matches getParams() with no time_range", {
    p_get  <- getParams(NS_sim_small)
    p_final <- finalParams(NS_sim_small)
    expect_equal(p_final@initial_n, p_get@initial_n)
    expect_equal(p_final@initial_n_pp, p_get@initial_n_pp)
    expect_equal(p_final@initial_effort, p_get@initial_effort)
})

test_that("initialParams() returns values from the first time step", {
    sim <- short_ns_sim
    p <- initialParams(sim)
    expect_equal(p@initial_n, sim@n[1, , ], ignore_attr = TRUE)
    expect_equal(p@initial_n_pp, sim@n_pp[1, ], ignore_attr = TRUE)
    expect_equal(p@initial_effort, sim@effort[1, ], ignore_attr = TRUE)
})
