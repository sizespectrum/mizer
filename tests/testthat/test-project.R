params <- NS_params_default_small

# time dimension ----
test_that("time dimension is dealt with properly", {

    # Effort is a single numeric
    t_max <- 2
    t_save <- 1
    dt <- 0.1
    sim <- project(params, t_max = t_max, t_save = t_save, dt = dt, effort = 1)
    expect_identical(names(dimnames(sim@effort)), c("time", "gear"))
    expect_equal(dim(sim@effort)[1], 
                 length(seq(from = 0, to = t_max, by = t_save)))
    expect_equal(dim(sim@n)[1], length(seq(from = 0, to = t_max, by = t_save)))
    expect_identical(dimnames(sim@effort)[[1]], 
                     as.character(seq(from = 0, to = t_max, by = t_save)))
    expect_identical(dimnames(sim@n)[[1]], 
                     as.character(seq(from = 0, to = t_max, by = t_save)))
    dt <- 0.5
    t_save <- 2
    sim <- project(params, t_max = t_max, t_save = t_save, dt = dt, effort = 1)
    expect_equal(dim(sim@effort)[1],
                 length(seq(from = 0, to = t_max, by = t_save)))
    expect_equal(dim(sim@n)[1],
                 length(seq(from = 0, to = t_max, by = t_save)))
    expect_identical(dimnames(sim@effort)[[1]],
                     as.character(seq(from = 0, to = t_max, by = t_save)))
    expect_identical(dimnames(sim@n)[[1]],
                     as.character(seq(from = 0, to = t_max, by = t_save)))
    t_save <- 0.5
    dt <- 0.5
    sim <- project(params, t_max = t_max, t_save = t_save, dt = dt, effort = 1)
    expect_equal(dim(sim@effort)[1], t_max / t_save + 1)
    expect_equal(dim(sim@n)[1], t_max / t_save + 1)
    expect_identical(dimnames(sim@effort)[[1]],
                     as.character(seq(from = 0, to = t_max, by = t_save)))
    expect_identical(dimnames(sim@n)[[1]],
                     as.character(seq(from = 0, to = t_max, by = t_save)))
    # append
    sim <- project(sim, t_max = t_max, t_save = t_save, dt = dt, effort = 1)
    expect_equal(dim(sim@effort)[1], 2 * t_max/t_save + 1)
    expect_equal(dim(sim@n)[1], 2 * t_max/t_save + 1)
    expect_identical(dimnames(sim@effort)[[1]],
                     as.character(seq(from = 0, to = 2 * t_max, by = t_save)))
    expect_identical(dimnames(sim@n)[[1]],
                     as.character(seq(from = 0, to = 2 * t_max, by = t_save)))
    

    # Effort is an effort vector
    effort <- c(Industrial = 1, Pelagic = 0.5, Otter = 0)
    t_max <- 2
    t_save <- 2
    sim <- project(params, t_max = t_max, t_save = t_save, effort = effort)
    expect_identical(names(dimnames(sim@effort)), c("time", "gear"))
    expect_equal(dim(sim@effort)[1], length(seq(from = 0, to = t_max, by = t_save)))
    expect_equal(dim(sim@n)[1], length(seq(from = 0, to = t_max, by = t_save)))
    expect_identical(dimnames(sim@effort)[[1]],
                     as.character(seq(from = 0, to = t_max, by = t_save)))
    expect_identical(dimnames(sim@n)[[1]],
                     as.character(seq(from = 0, to = t_max, by = t_save)))
    dt <- 0.5
    sim <- project(params, t_max = t_max, t_save = t_save, effort = effort)
    expect_equal(dim(sim@effort)[1], length(seq(from = 0, to = t_max, by = t_save)))
    expect_equal(dim(sim@n)[1], length(seq(from = 0, to = t_max, by = t_save)))
    expect_identical(dimnames(sim@effort)[[1]],
                     as.character(seq(from = 0, to = t_max, by = t_save)))
    expect_identical(dimnames(sim@n)[[1]],
                     as.character(seq(from = 0, to = t_max, by = t_save)))
    t_save <- 0.5
    sim <- project(params, t_max = t_max, t_save = t_save, effort = effort)
    expect_equal(dim(sim@effort)[1], length(seq(from = 0, to = t_max, by = t_save)))
    expect_equal(dim(sim@n)[1], length(seq(from = 0, to = t_max, by = t_save)))
    expect_identical(dimnames(sim@effort)[[1]],
                     as.character(seq(from = 0, to = t_max, by = t_save)))
    expect_identical(dimnames(sim@n)[[1]],
                     as.character(seq(from = 0, to = t_max, by = t_save)))

    ## No effort argument but t_start
    sim <- project(params, t_start = 2019, t_max = 2, dt = 1)
    expect_equal(dimnames(sim@n)$time, c("2019", "2020", "2021"))
})

test_that("warns and includes t_max when not a multiple of t_save", {
    # Scalar effort path
    expect_warning(
        sim <- project(params, t_max = 10, t_save = 3, effort = 1,
                       progress_bar = FALSE),
        "`t_max`.*`t_save`"
    )
    t_dims <- as.numeric(dimnames(sim@n)$time)
    expect_equal(tail(t_dims, 1), 10)
    expect_equal(t_dims, c(0, 3, 6, 9, 10))

    # No warning when t_max is exactly a multiple of t_save
    expect_no_warning(
        project(params, t_max = 9, t_save = 3, effort = 1,
                progress_bar = FALSE)
    )
})

test_that("project method selects consumer time stepper", {
    params_small <- newMultispeciesParams(NS_species_params_gears_small[nrow(NS_species_params_gears_small), , drop = FALSE],
                                          info_level = 0)

    sim_default <- project(params_small, t_max = 0.2, t_save = 0.2,
                           dt = 0.1, effort = 1, progress_bar = FALSE)
    sim_euler <- project(params_small, t_max = 0.2, t_save = 0.2,
                         dt = 0.1, effort = 1, progress_bar = FALSE,
                         method = "euler")
    sim_pc <- project(params_small, t_max = 0.2, t_save = 0.2,
                      dt = 0.1, effort = 1, progress_bar = FALSE,
                      method = "predictor_corrector")

    expect_equal(sim_default@n, sim_euler@n)
    expect_true(all(is.finite(sim_pc@n)))
    expect_true(all(sim_pc@n >= 0))
    expect_error(
        project(params_small, t_max = 0.2, t_save = 0.2, dt = 0.1,
                effort = 1, progress_bar = FALSE, method = "bogus"),
        "should be one of"
    )

    expect_warning(
        sim_appended <- project(sim_default, t_max = 0.2, t_save = 0.2,
                                dt = 0.1, effort = 1, progress_bar = FALSE,
                                method = "predictor_corrector"),
        "method"
    )
    expect_s3_class(sim_appended, "MizerSim")
    expect_equal(dim(sim_appended@n)[1], 3)
})

test_that("project.MizerSim defaults dt and method from sim_params", {
    params_small <- newMultispeciesParams(NS_species_params_gears_small[nrow(NS_species_params_gears_small), , drop = FALSE],
                                          info_level = 0)
    sim_pc <- project(params_small, t_max = 0.2, t_save = 0.2, dt = 0.05,
                      effort = 1, progress_bar = FALSE,
                      method = "predictor_corrector")
    # Continuing without supplying dt or method should use stored values
    sim_cont <- project(sim_pc, t_max = 0.2, t_save = 0.2,
                        effort = 1, progress_bar = FALSE)
    expect_equal(getSimParams(sim_cont)$dt, 0.05)
    expect_equal(getSimParams(sim_cont)$method, "predictor_corrector")
    # Mismatched dt should warn
    expect_warning(
        project(sim_pc, t_max = 0.2, t_save = 0.2, dt = 0.1,
                effort = 1, progress_bar = FALSE),
        "dt"
    )
    # No warning when dt and method match
    expect_no_warning(
        project(sim_pc, t_max = 0.2, t_save = 0.2, dt = 0.05,
                method = "predictor_corrector", effort = 1,
                progress_bar = FALSE)
    )
})

# pass in initial species ----
test_that("Can pass in initial species", {
    no_gear <- dim(params@catchability)[1]
    no_sp <- dim(params@catchability)[2]
    max_t_effort <- 10
    effort <- array(abs(rnorm(max_t_effort * no_gear)),
                    dim = c(max_t_effort, no_gear))

    # No time dimnames - fail
    t_max <- 5
    start_year <- 1980
    time_step <- 0.5
    end_year <- start_year + t_max - 1
    time <- seq(from = start_year, to = end_year, by = time_step)
    effort <- array(NA, dim = c(length(time), 4), 
                    dimnames = list(NULL, gear = c("industrial", "pelagic",
                                                   "otter_trawl", "beam_trawl")
                                    )
                    )
    effort[,1] <- seq(from = 0, to = 1, length = nrow(effort))
    effort[,2] <- 0.5
    effort[,3] <- seq(from = 1, to = 0.5, length = nrow(effort))
    effort[,4] <- 0
    expect_error(project(params, effort = effort))
})


# w_min array reference ----
test_that("w_min array reference is working OK", {
    NS_species_params_gears_small$w_min <- 0.001
    NS_species_params_gears_small$w_min[1] <- 1
    params2 <- newMultispeciesParams(NS_species_params_gears_small, inter_small, info_level = 0)
    sim <- project(params2, effort = 1, t_max = 1)
    final <- dim(sim@n)[1]
    expect_equal(sim@n[final, 1, 1:(sim@params@w_min_idx[1] - 1)],
                 rep(0, sim@params@w_min_idx[1] - 1), ignore_attr = TRUE)
})


# Gear checking and sorting ----
test_that("Gear checking and sorting is OK", {
    # Set up trait based model for easy testing ground
    no_sp <- 10
    min_w_max <- 10
    max_w_max <- 1e5
    w_max <- 10^seq(from = log10(min_w_max), to = log10(max_w_max),
                    length = no_sp)
    knife_edges <- w_max * 0.05
    industrial_gears <- w_max <= 500
    other_gears <- w_max > 500
    gear_names <- rep("Industrial", no_sp)
    gear_names[other_gears] <- "Other"
    params_gear <- newTraitParams(no_sp = no_sp, 
                                  min_w_max = min_w_max, 
                                  max_w_max = max_w_max, 
                                  knife_edge_size = knife_edges, 
                                  gear_names = gear_names)
    gear_names <- dimnames(params_gear@catchability)[[1]]
    # Single vector of effort
    sim <- project(params_gear, effort = 0.3, t_max = 2)
    expect_true(all(sim@effort == 0.3))
    # Also check that order of gear names in resulting effort matches catchability
    expect_true(all(dimnames(sim@effort)$gear == gear_names))
    # Effort vector
    # Should give same result
    effort_vec <- c(Other = 1, Industrial = 0)
    effort_vec2 <- c(Industrial = 0, Other = 1)
    sim <- project(params_gear, effort = effort_vec, t_max = 2)
    sim2 <- project(params_gear, effort = effort_vec2, t_max = 2)
    expect_true(all(sim@effort[, "Industrial"] == 0))
    expect_true(all(sim@effort[, "Other"] == 1))
    expect_true(all(sim2@effort[, "Industrial"] == 0))
    expect_true(all(sim2@effort[, "Other"] == 1))
    expect_true(all(dimnames(sim@effort)$gear == gear_names)) 
    expect_true(all(dimnames(sim2@effort)$gear == gear_names)) 
    # Should fail - number of gears wrong
    effort_vec3 <- c(Industrial = 0, Other = 1, Dummy = 0.5)
    expect_error(project(params_gear, effort = effort_vec3, t_max = 2))
    effort_vec4 <- c(Industrial = 0) # Is OK because that gear exists
    expect_error(project(params_gear, effort = effort_vec4, t_max = 2), NA)
    # Should fail - names of gears wrong
    effort_vec5 <- c(Industrial = 0, Dummy = 1)
    expect_error(project(params_gear, effort = effort_vec5, t_max = 2))
    # Array effort
    t_steps <- 3
    effort1 <- array(1, dim = c(t_steps, 2))
    expect_error(project(params_gear, effort = effort1))
    # Different order - should give same result
    effort2 <- array(
      rep(c(1, 0), each = t_steps),
      dim = c(t_steps, 2),
      dimnames = list(
        time = 1:t_steps,
        gear = c("Other", "Industrial")
      )
    )
    effort3 <- array(
      rep(c(0, 1), each = t_steps),
      dim = c(t_steps, 2),
      dimnames = list(
        time = 1:t_steps,
        gear = c("Industrial", "Other")
      )
    )
    sim2 <- project(params_gear, effort = effort2)
    sim3 <- project(params_gear, effort = effort3)
    expect_identical(sim2, sim3)
    # These should all fail - gears incorrectly specified
    effort4 <-
      array(
        rep(c(0, 1, 0.5), each = t_steps),
        dim = c(t_steps, 3),
        dimnames = list(
          time = 1:t_steps,
          gear = c("Industrial", "Other", "Dummy")
        )
      )
    effort5 <- array(
      rep(c(0, 1), each = t_steps),
      dim = c(t_steps, 2),
      dimnames = list(
        time = 1:t_steps,
        gear = c("Industrial", "Dummy")
      )
    )
    effort6 <- array(
      rep(c(1), each = t_steps),
      dim = c(t_steps, 1),
      dimnames = list(time = 1:t_steps, gear = c("Industrial"))
    )
    expect_error(project(params_gear, effort = effort4))
    expect_error(project(params_gear, effort = effort5))
    expect_error(project(params_gear, effort = effort6))
})


# same numerical results as previously ----
test_that("Simulation gives same numerical results as previously",{
  # Snapshot values were recorded with edition 1
  old <- getOption("mizer_defaults_edition")
  on.exit(options(mizer_defaults_edition = old), add = TRUE)
  options(mizer_defaults_edition = 1)
  params <- newMultispeciesParams(NS_species_params_gears_small, inter_small,
                                  n = 2/3, p = 0.7, lambda = 2.8 - 2/3, info_level = 0)
  sim <- project(params, t_max = 1)
  # expect_known_value(sim@n[2, 3, ], "values/projectn")
  # expect_known_value(sim@n_pp[2, ], "values/projectp")
  expect_snapshot(sim@n[2, 3, ])
  expect_snapshot(sim@n_pp[2, ])
  
})

test_that("Final result the same when called with sim or params", {
  params <- NS_params_small
  sim <- project(params, t_max = 1)
  params@initial_n[] <- sim@n[2, , ]
  params@initial_n_pp[] <- sim@n_pp[2, ]
  params@initial_n_other <- sim@n_other[2, ]
  sim1 <- project(params, t_max = 1)
  sim2 <- project(sim, t_max = 1)
  expect_identical(sim1@n[2, 3, ], sim2@n[3, 3, ])
})

# dimnames ----
# This test is motivated by the bug in 
# https://github.com/sizespectrum/mizer/issues/173
test_that("Dimnames on effort have correct names", {
  gear_names <- unique(gear_params(NS_params_small)$gear)
  effort <- array(1, dim = c(3, length(gear_names)), 
                  dimnames = list(1:3,
                                  gear_names))
  sim <- project(NS_params_small, effort, t_max = 0.1)
  expect_identical(names(dimnames(sim@effort)), c("time", "gear"))
})

# t_max and t_save with effort arrays ----
test_that("t_max extends simulation beyond effort array", {
  gear_names <- unique(gear_params(NS_params_small)$gear)
  # Create effort array that goes from year 1 to year 5
  effort <- array(0.5, dim = c(5, length(gear_names)),
                  dimnames = list(time = 1:5, gear = gear_names))
  
  # Extend simulation to year 10 with t_max
  sim <- project(NS_params_small, effort = effort, t_max = 9, dt = 0.1)
  
  # Should save at times 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
  expect_equal(as.numeric(dimnames(sim@n)[[1]]), 1:10)
  expect_equal(dim(sim@n)[1], 10)
  
  # Effort should be extrapolated with last value (0.5)
  expect_true(all(sim@effort[, ] == 0.5))
})

test_that("t_save controls save frequency with effort array", {
  gear_names <- unique(gear_params(NS_params_small)$gear)
  # Create effort array with times 0, 1, 2, 3, 4, 5
  effort <- array(0.5, dim = c(6, length(gear_names)),
                  dimnames = list(time = 0:5, gear = gear_names))
  
  # Use t_save to control output frequency
  sim <- project(NS_params_small, effort = effort, t_save = 2, dt = 0.1)
  
  # Should save at times 0, 2, 4
  expect_equal(as.numeric(dimnames(sim@n)[[1]]), c(0, 2, 4))
  expect_equal(dim(sim@n)[1], 3)
})

test_that("t_max and t_save work together with effort array", {
  gear_names <- unique(gear_params(NS_params_small)$gear)
  # Create effort array from 0 to 3
  effort <- array(0.5, dim = c(4, length(gear_names)),
                  dimnames = list(time = 0:3, gear = gear_names))
  
  # Extend to 6 and save every 2 years
  sim <- project(NS_params_small, effort = effort, t_max = 6, t_save = 2, dt = 0.1)
  
  # Should save at times 0, 2, 4, 6
  expect_equal(as.numeric(dimnames(sim@n)[[1]]), c(0, 2, 4, 6))
  expect_equal(dim(sim@n)[1], 4)
})

test_that("Effort array times used when t_max and t_save not provided", {
  gear_names <- unique(gear_params(NS_params_small)$gear)
  # Create effort array with irregular times
  effort <- array(0.5, dim = c(4, length(gear_names)),
                  dimnames = list(time = c(0, 1, 3, 7), gear = gear_names))
  
  # Without t_max or t_save, should use effort array times
  sim <- project(NS_params_small, effort = effort, dt = 0.1)
  
  expect_equal(as.numeric(dimnames(sim@n)[[1]]), c(0, 1, 3, 7))
  expect_equal(dim(sim@n)[1], 4)
})

test_that("Effort values interpolated correctly", {
  gear_names <- unique(gear_params(NS_params_small)$gear)
  # Create effort array with varying effort
  effort <- array(NA, dim = c(3, length(gear_names)),
                  dimnames = list(time = c(0, 5, 10), gear = gear_names))
  effort[, 1] <- c(0, 0.5, 1.0)  # Industrial
  effort[, 2:3] <- 0.3  # Others constant
  
  # Use t_save to create intermediate time points
  sim <- project(NS_params_small, effort = effort, t_save = 2.5, dt = 0.1)
  
  # Should have times 0, 2.5, 5, 7.5, 10
  expect_equal(as.numeric(dimnames(sim@n)[[1]]), c(0, 2.5, 5, 7.5, 10))
  
  # Check effort interpolation (step function)
  # At time 0, effort is 0
  expect_equal(sim@effort[1, "Industrial"], 0)
  # At time 2.5, effort should be 0 (constant from time 0)
  expect_equal(sim@effort[2, "Industrial"], 0)
  # At time 5, effort is 0.5
  expect_equal(sim@effort[3, "Industrial"], 0.5)
  # At time 7.5, effort should be 0.5 (constant from time 5)
  expect_equal(sim@effort[4, "Industrial"], 0.5)
  # At time 10, effort is 1.0
  expect_equal(sim@effort[5, "Industrial"], 1.0)
})

test_that("Can extend simulation with NA in final effort", {
  # This is the motivating use case from the issue
  gear_names <- unique(gear_params(NS_params_small)$gear)
  effort <- array(0.5, dim = c(3, length(gear_names)),
                  dimnames = list(time = 2017:2019, gear = gear_names))
  # The last year has NA which gets replaced by default (1 for edition >= 2)
  effort[3, ] <- NA
  
  # Run until 2020 without needing to specify effort for 2019
  sim <- project(NS_params_small, effort = effort, t_max = 3, dt = 0.1)
  
  # Should have years 2017, 2018, 2019, 2020
  expect_equal(as.numeric(dimnames(sim@n)[[1]]), 2017:2020)
  expect_equal(dim(sim@n)[1], 4)
  
  # The NA gets replaced by default effort during validation
  # Then extrapolated for 2020
  expect_true(all(!is.na(sim@effort)))
})

test_that("Named effort vectors fill missing gears with the default effort", {
  gear_names <- unique(gear_params(NS_params_small)$gear)
  effort <- c(Industrial = 0.25, Pelagic = 0.4)

  sim <- project(NS_params_small, effort = effort, t_max = 1, progress_bar = FALSE)

  default_effort <- ifelse(defaults_edition(NS_params_small) < 2, 0, 1)
  expect_equal(sim@effort[1, "Industrial"], 0.25)
  expect_equal(sim@effort[1, "Pelagic"], 0.4)
  expect_true(all(sim@effort[1, setdiff(gear_names, c("Industrial", "Pelagic"))] ==
                    default_effort))
})

test_that("NA in an effort array is replaced with the model's default", {
  gear_names <- unique(gear_params(NS_params_small)$gear)
  effort <- array(NA_real_, dim = c(2, length(gear_names)),
                  dimnames = list(time = 0:1, gear = gear_names))
  params <- NS_params_small
  withr::local_options(mizer_defaults_edition = 2)
  expect_true(all(validEffortArray(effort, params) == 0))
  defaults_edition(params) <- 2
  withr::local_options(mizer_defaults_edition = 1)
  expect_true(all(validEffortArray(effort, params) == 1))
})

test_that("t_max less than effort array duration uses effort times", {
  gear_names <- unique(gear_params(NS_params_small)$gear)
  effort <- array(0.5, dim = c(6, length(gear_names)),
                  dimnames = list(time = 0:5, gear = gear_names))
  
  # t_max = 3 should run until time 3, not 5
  sim <- project(NS_params_small, effort = effort, t_max = 3, dt = 0.1)
  
  # Should stop at year 3
  expect_equal(max(as.numeric(dimnames(sim@n)[[1]])), 3)
})

test_that("project does not change the params object", {
    params <- NS_params_small
    params@ext_diffusion[] <- 1
    old_params <- unserialize(serialize(params, NULL))
    sim <- project(params, t_max = 1)
    expect_identical(params, old_params)
})


# predation diffusion ----
test_that("Simulation gives same numerical results with predation diffusion", {
    params_d <- NS_params_small
    params_d@use_predation_diffusion <- TRUE
    sim_d <- project(params_d, t_max = 1)
    expect_snapshot(sim_d@n[2, 3, ])
    expect_snapshot(sim_d@n_pp[2, ])
})

test_that("Final result the same when called with sim or params, predation diffusion on", {
    params_d <- NS_params_small
    params_d@use_predation_diffusion <- TRUE
    sim_d <- project(params_d, t_max = 1)
    params_d@initial_n[] <- sim_d@n[2, , ]
    params_d@initial_n_pp[] <- sim_d@n_pp[2, ]
    params_d@initial_n_other <- sim_d@n_other[2, ]
    sim1 <- project(params_d, t_max = 1)
    sim2 <- project(sim_d, t_max = 1)
    expect_identical(sim1@n[2, 3, ], sim2@n[3, 3, ])
})

test_that("Predation diffusion changes simulation trajectory", {
    params_d <- NS_params_small
    params_d@use_predation_diffusion <- TRUE
    sim_d <- project(params_d, t_max = 1)
    sim_base <- project(NS_params_small, t_max = 1)
    expect_false(identical(sim_d@n, sim_base@n))
})

test_that("callback works during simulation", {
    params <- NS_params_small
    
    # Test callback execution
    counter <- 0
    test_callback <- function(sim, t_idx) {
        counter <<- counter + 1
    }
    sim <- project(params, t_max = 2, callback = test_callback, progress_bar = FALSE)
    expect_s3_class(sim, "MizerSim")
    expect_equal(counter, 3)
})

# project_simple ----

test_that("project_simple matches project for one saved step", {
    params <- NS_params_small
    dt <- 0.1
    t_save <- 1
    steps <- as.integer(t_save / dt)
    effort <- initial_effort(params)
    res_fn <- get(params@resource_dynamics)
    other_fns <- lapply(params@other_dynamics, get)
    rates_fns <- lapply(params@rates_funcs, get)
    start <- list(
        n = params@initial_n,
        n_pp = params@initial_n_pp,
        n_other = params@initial_n_other
    )
    out <- project_simple(
        params,
        n = start$n,
        n_pp = start$n_pp,
        n_other = start$n_other,
        t = 0,
        dt = dt,
        steps = steps,
        effort = effort,
        resource_dynamics_fn = res_fn,
        other_dynamics_fns = other_fns,
        rates_fns = rates_fns
    )
    sim <- project(params, t_max = t_save, dt = dt, t_save = t_save, progress_bar = FALSE)
    expected_n <- params@initial_n
    expected_n[] <- N(sim)[2, , ]
    expect_equal(out$n, expected_n, tolerance = 1e-12)
    expect_equal(out$n_pp, NResource(sim)[2, ], tolerance = 1e-12,
                 ignore_attr = TRUE)
})

test_that("Euler reselects its boundary from current custom diffusion", {
    observed <- new.env(parent = emptyenv())
    assign("boundary_test_growth", function(params, ...) {
        growth <- params@initial_n * 0
        growth[1, ] <- 1
        growth
    }, envir = .GlobalEnv)
    assign("boundary_test_diffusion", function(params, n, n_other, t, ...) {
        top <- support_top_idx(params)[1]
        observed$times <- c(observed$times, t)
        observed$clocks <- c(observed$clocks, n_other$clock)
        observed$tail <- c(observed$tail, n[1, top + 1L])
        diffusion <- n * 0
        # Only Herring diffuses, and only during the middle Euler step.
        if (n_other$clock > 0.0625 && n_other$clock < 0.1875) {
            diffusion[2, top] <- 1
        }
        diffusion
    }, envir = .GlobalEnv)
    assign("boundary_test_clock", function(n_other, component, dt, ...) {
        n_other[[component]] + dt
    }, envir = .GlobalEnv)
    withr::defer(rm(list = c("boundary_test_growth", "boundary_test_diffusion",
                             "boundary_test_clock"), envir = .GlobalEnv))

    for (flux in c("upwind", "van_leer", "centred")) {
        for (bin_average in c(FALSE, TRUE)) {
            p <- NS_params_small
            second_order_w(p) <- c(flux = flux, bin_average = bin_average)
            top <- support_top_idx(p)[1]
            p@initial_n[] <- 0
            p@initial_n[1, top] <- 1
            p <- setComponent(p, "clock", initial_value = 0,
                              dynamics_fun = "boundary_test_clock")
            p <- setRateFunction(p, "EGrowth", "boundary_test_growth")
            p <- setRateFunction(p, "Diffusion", "boundary_test_diffusion")
            expect_true(all(p@ext_diffusion == 0))
            observed$times <- observed$clocks <- observed$tail <- numeric()
            dt <- 0.125
            out <- project_simple(p, dt = dt, steps = 3, method = "euler",
                                  info_level = 0)
            expect_equal(observed$times, c(0, dt, 2 * dt))
            expect_equal(observed$clocks, c(0, dt, 2 * dt))
            expect_equal(out$n_other$clock, 3 * dt)
            expect_equal(observed$tail[1], 0)
            expect_equal(observed$tail[3], 0)
            if (flux == "upwind") {
                # No diffusion lets the legacy step fill Sprat's tail. The
                # next step clears it because another species diffuses; once
                # diffusion stops, the legacy behaviour returns immediately.
                expect_gt(observed$tail[2], 0)
                expect_gt(out$n[1, top + 1L], 0)
            } else {
                # Higher-order fluxes keep the hard boundary throughout,
                # including both steps without diffusion.
                expect_equal(observed$tail[2], 0)
                expect_equal(unname(out$n[1, top + 1L]), 0)
            }
        }
    }
})

test_that("project_simple returns rates from the final update step", {
    params <- NS_params_small
    dt <- 0.1
    steps <- 2
    effort <- initial_effort(params)
    res_fn <- get(params@resource_dynamics)
    other_fns <- lapply(params@other_dynamics, get)
    rates_fns <- lapply(params@rates_funcs, get)
    out <- project_simple(
        params,
        n = params@initial_n,
        n_pp = params@initial_n_pp,
        n_other = params@initial_n_other,
        t = 0,
        dt = dt,
        steps = steps,
        effort = effort,
        resource_dynamics_fn = res_fn,
        other_dynamics_fns = other_fns,
        rates_fns = rates_fns
    )

    step1 <- project_simple(
        params,
        n = params@initial_n,
        n_pp = params@initial_n_pp,
        n_other = params@initial_n_other,
        t = 0,
        dt = dt,
        steps = 1,
        effort = effort,
        resource_dynamics_fn = res_fn,
        other_dynamics_fns = other_fns,
        rates_fns = rates_fns
    )
    expected_rates <- getRates(
        params,
        n = step1$n,
        n_pp = step1$n_pp,
        n_other = step1$n_other,
        effort = effort,
        t = dt
    )

    expect_identical(names(out$rates), names(expected_rates))
    expect_equal(out$rates$encounter, expected_rates$encounter)
    expect_equal(out$rates$rdd, expected_rates$rdd)
})

test_that("project_simple accepts predictor-corrector method", {
    params <- NS_params_small
    effort <- initial_effort(params)
    res_fn <- get(params@resource_dynamics)
    other_fns <- lapply(params@other_dynamics, get)
    rates_fns <- lapply(params@rates_funcs, get)

    out <- project_simple(
        params,
        n = params@initial_n,
        n_pp = params@initial_n_pp,
        n_other = params@initial_n_other,
        t = 0,
        dt = 0.1,
        steps = 1,
        effort = effort,
        resource_dynamics_fn = res_fn,
        other_dynamics_fns = other_fns,
        rates_fns = rates_fns,
        method = "predictor_corrector"
    )

    expect_true(all(is.finite(out$n)))
    expect_true(all(out$n >= 0))
})

# time resampling ----

test_that("Time resampling in project behaves as documented", {
    params <- NS_params_small
    gear_names <- unique(gear_params(params)$gear)

    # Create effort array with irregular times
    times <- c(1, 10, 11)
    effort <- array(1,
        dim = c(length(times), length(gear_names)),
        dimnames = list(time = times, gear = gear_names)
    )

    # Case 1: No t_max or t_save provided
    # Should preserve original times
    sim <- project(params, effort = effort, dt = 0.1)
    expect_equal(as.numeric(dimnames(sim@n)[[1]]), times)

    # Case 2: t_max provided
    # Should resample based on default t_save (which is 1 or inferred)
    # Here inferred t_save from first interval is 10 - 1 = 9.
    # So times would be 1, 10, 19...
    # Wait, the code says:
    # if (length(effort_times) > 1) {
    #     save_freq <- effort_times[2] - effort_times[1]
    # }
    # So save_freq becomes 9.
    # New times: seq(1, 13, by = 9) -> 1, 10.
    # Time 11 is lost.
    sim <- project(params, effort = effort, t_max = 13, dt = 1)
    saved_times <- as.numeric(dimnames(sim@n)[[1]])

    expect_equal(saved_times, c(1, 10))
    expect_false(11 %in% saved_times)

    # Case 3: t_save provided explicitly
    # t_save = 1. Should capture 1, 2, ..., 11, ...
    sim <- project(params, effort = effort, t_max = 13, t_save = 1, dt = 1)
    saved_times <- as.numeric(dimnames(sim@n)[[1]])
    expect_true(11 %in% saved_times)

    # Case 4: t_save provided but does not align with irregular time
    # t_save = 5. Times: 1, 6, 11, 16, 21...
    # Here 11 IS captured by coincidence (1 + 2*5 = 11).
    # Let's try t_save = 4. Times: 1, 5, 9, 13... 11 is lost.
    sim <- project(params, effort = effort, t_max = 13, t_save = 4, dt = 1)
    saved_times <- as.numeric(dimnames(sim@n)[[1]])
    expect_false(11 %in% saved_times)
})

# check_steady ----

# A settled model and a copy that has been knocked off its steady state.
steady_small <- NS_params_steady_small
off_steady_small <- local({
    p <- steady_small
    initialN(p)[1, ] <- initialN(p)[1, ] * 3
    p
})

test_that("project() checks the steady state only when asked", {
    # Off by default, however far off steady state the model is.
    expect_silent(project(off_steady_small, t_max = 0.1, progress_bar = FALSE))
    # On when asked.
    expect_warning(project(off_steady_small, t_max = 0.1, progress_bar = FALSE,
                           check_steady = TRUE),
                   "not at its steady state")
    # And silent on a model that is at its steady state.
    expect_silent(project(steady_small, t_max = 0.1, progress_bar = FALSE,
                          check_steady = TRUE))
})

test_that("project(check_steady) judges at the model's own effort", {
    # Running a scenario at a new effort legitimately starts away from the
    # steady state of that new effort, and must not warn for it.
    expect_silent(project(steady_small, effort = initial_effort(steady_small) + 1,
                          t_max = 0.1, progress_bar = FALSE, check_steady = TRUE))
})

test_that("project.MizerSim checks the steady state only when asked", {
    sim <- project(off_steady_small, t_max = 0.1, progress_bar = FALSE)
    for (append in c(TRUE, FALSE)) {
        expect_silent(project(sim, t_max = 0.1, append = append,
                              progress_bar = FALSE))
        expect_warning(project(sim, t_max = 0.1, append = append,
                               progress_bar = FALSE, check_steady = TRUE),
                       "not at its steady state")
    }
})

test_that("project(check_steady) is silenced by info_level 0", {
    withr::local_options(mizer_info_level = 0)
    expect_silent(project(off_steady_small, t_max = 0.1, progress_bar = FALSE,
                          check_steady = TRUE))
})

# method aliases ----

test_that("normalise_project_method resolves aliases to canonical names", {
    expect_identical(normalise_project_method("second_order"), "tr_bdf2")
    expect_identical(normalise_project_method("predictor-corrector"),
                     "predictor_corrector")
    # Canonical names pass through unchanged.
    for (m in c("euler", "predictor_corrector", "tr_bdf2")) {
        expect_identical(normalise_project_method(m), m)
    }
    # A vector of choices, as supplied by the default of the `method` argument,
    # resolves to its first entry.
    expect_identical(
        normalise_project_method(c("euler", "predictor_corrector", "tr_bdf2")),
        "euler")
    expect_error(normalise_project_method("nonesuch"))
})

test_that("method = 'second_order' projects exactly as 'tr_bdf2' does", {
    params <- NS_params_small
    sim_alias <- project(params, t_max = 1, dt = 0.2, t_save = 1,
                         method = "second_order", progress_bar = FALSE)
    sim_canon <- project(params, t_max = 1, dt = 0.2, t_save = 1,
                         method = "tr_bdf2", progress_bar = FALSE)
    expect_equal(finalN(sim_alias), finalN(sim_canon))
    expect_equal(finalNResource(sim_alias), finalNResource(sim_canon))
    # The alias is resolved before being stored, so the sim reports the
    # canonical name and continuing the run does not warn about a mismatch.
    expect_identical(getSimParams(sim_alias)$method, "tr_bdf2")
    expect_no_warning(project(sim_alias, t_max = 0.2, dt = 0.2,
                              method = "second_order", progress_bar = FALSE))
})

test_that("growth Courant advice is once per call and respects info level", {
    p <- NS_params_small
    second_order_w(p) <- c(flux = "van_leer")
    messages <- character()
    withCallingHandlers({
        project(p, dt = 0.2, t_max = 0.4, t_save = 0.2,
                progress_bar = FALSE, info_level = 1)
    }, warning = function(w) {
        messages <<- c(messages, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    expect_length(messages, 1)
    expect_match(messages, "Growth Courant number")
    expect_match(messages, "size bin")
    expect_match(messages, "current-state suggestion")
    expect_no_warning(project(p, dt = 0.2, t_max = 0.4, t_save = 0.2,
                              progress_bar = FALSE, info_level = 0))
    second_order_w(p) <- c(flux = "upwind")
    expect_no_warning(project(p, dt = 0.2, t_max = 0.4, t_save = 0.2,
                              progress_bar = FALSE))
})

test_that("uniform substeps advance consumers resources components and time", {
    p <- NS_params_small
    second_order_w(p) <- c(flux = "van_leer")
    p <- setComponent(p, "tracked", initial_value = 0,
                      dynamics_fun = "constant_other")
    fns <- projectRateFunctions(p)
    original_rates <- fns$Rates
    g <- p@initial_n * 0
    g[, 2] <- 20 * p@dw[1]
    observed <- new.env(parent = emptyenv())
    observed$times <- numeric()
    fns$Rates <- function(params, t, ...) {
        observed$times <- c(observed$times, t)
        r <- original_rates(params, t = t, ...)
        r$e_growth <- g
        r
    }
    components <- list(tracked = function(params, n, n_pp, n_other, component,
                                          t, dt, ...) {
        n_other[[component]] + dt * (t + sum(n[, 1]) / 1e12 + n_pp[1] / 1e20)
    })
    for (method in c("euler", "tr_bdf2")) {
        observed$times <- numeric()
        adaptive <- project_simple(p, dt = 0.2, steps = 2, method = method,
                                   rates_fns = fns,
                                   other_dynamics_fns = components,
                                   subcycle = TRUE, info_level = 0)
        adaptive_times <- observed$times
        observed$times <- numeric()
        fine <- project_simple(p, dt = 0.05, steps = 8, method = method,
                               rates_fns = fns,
                               other_dynamics_fns = components, info_level = 0)
        expect_equal(adaptive$n, fine$n, tolerance = 1e-12)
        expect_equal(adaptive$n_pp, fine$n_pp, tolerance = 1e-12)
        expect_equal(adaptive$n_other, fine$n_other, tolerance = 1e-12)
        expect_equal(adaptive_times, observed$times, tolerance = 1e-14)
    }
})

test_that("subcycling adapts to current growth and counts actual substeps", {
    p <- NS_params_small
    second_order_w(p) <- c(flux = "van_leer")
    fns <- projectRateFunctions(p)
    original_rates <- fns$Rates
    times <- numeric()
    fns$Rates <- function(params, t, ...) {
        times <<- c(times, t)
        r <- original_rates(params, t = t, ...)
        r$e_growth[] <- 0
        r$e_growth[, 2] <- if (t < 0.04) 20 * p@dw[1] else 0
        r
    }
    # Initially four substeps are forecast, but growth stops after the first.
    expect_no_error(project_simple(p, dt = 0.2, steps = 1, rates_fns = fns,
                                   subcycle = TRUE, max_substeps = 2,
                                   info_level = 0))
    expect_equal(times, c(0, 0.05))
    times <- numeric()
    fns$Rates <- function(params, t, ...) {
        times <<- c(times, t)
        r <- original_rates(params, t = t, ...)
        r$e_growth[] <- 0
        r$e_growth[, 2] <- if (t < 0.04) 20 * p@dw[1] else 80 * p@dw[1]
        r
    }
    project_simple(p, dt = 0.2, steps = 1, rates_fns = fns,
                   subcycle = TRUE, info_level = 0)
    expect_gt(length(times), 4)
    expect_lte(max(diff(c(times, 0.2)) * ifelse(times < 0.04, 20, 80)),
               1 + 1e-12)
    expect_error(project_simple(p, dt = 0.2, steps = 1, rates_fns = fns,
                                subcycle = TRUE, max_substeps = 2,
                                info_level = 0), "max_substeps")
})

test_that("subcycling controls validate and persist on continuation", {
    p <- NS_params_small
    for (bad in list(NA, 1, c(TRUE, FALSE))) {
        expect_error(project(p, subcycle = bad), "subcycle")
    }
    expect_error(project(p, max_courant = 0), "max_courant")
    expect_error(project(p, max_courant = Inf), "max_courant")
    expect_error(project(p, max_substeps = 1.5), "max_substeps")
    callback_times <- numeric()
    callback <- function(sim, t_idx, ...) {
        callback_times <<- c(callback_times, getTimes(sim)[t_idx])
    }
    s <- project(p, dt = 0.1, t_max = 0.2, t_save = 0.1,
                 subcycle = TRUE, max_courant = 0.25, max_substeps = 500,
                 callback = callback, progress_bar = FALSE, info_level = 0)
    expect_equal(callback_times, c(0, 0.1, 0.2))
    continued <- project(s, t_max = 0.2, t_save = 0.1,
                         progress_bar = FALSE, info_level = 0)
    direct <- project(p, dt = 0.1, t_max = 0.4, t_save = 0.1,
                      subcycle = TRUE, max_courant = 0.25, max_substeps = 500,
                      progress_bar = FALSE, info_level = 0)
    expect_equal(continued@n, direct@n, tolerance = 1e-12)
    expect_equal(continued@n_pp, direct@n_pp, tolerance = 1e-12)
    expect_true(getSimParams(continued)$subcycle)
    expect_equal(getSimParams(continued)$max_courant, 0.25)
    expect_equal(getSimParams(continued)$max_substeps, 500)
    old <- s
    old@sim_params[c("subcycle", "max_courant", "max_substeps")] <- NULL
    resumed <- project(old, t_max = 0.1, t_save = 0.1,
                       progress_bar = FALSE, info_level = 0)
    expect_false(getSimParams(resumed)$subcycle)
})

test_that("subcycling removes the frozen-rate Cod checkerboard", {
    for (average in c(FALSE, TRUE)) {
        p <- NS_params_default_small
        second_order_w(p) <- c(flux = "van_leer", bin_average = average)
        rates <- getRates(p)
        fns <- projectRateFunctions(p)
        fns$Rates <- function(...) rates
        residual <- function(dt, steps, subcycle) {
            out <- project_simple(p, dt = dt, steps = steps, rates_fns = fns,
                                  method = "tr_bdf2", subcycle = subcycle,
                                  info_level = 0)
            n <- out$n
            attr(n, "chi") <- NULL
            r <- consumer_residual(p, n, p@initial_n_pp, p@initial_n_other,
                                   p@initial_effort, rdd = rates$rdd,
                                   rates = rates)
            max(abs((r / n)["Cod", 1:25]))
        }
        coarse <- residual(1 / 52, 52, FALSE)
        adaptive <- residual(1 / 52, 52, TRUE)
        fine <- residual(1 / 520, 520, FALSE)
        expect_gt(coarse, 1)
        expect_lt(adaptive, 1e-6)
        expect_lt(fine, 1e-6)
    }
})

test_that("subcycling lands exactly on requested time boundaries", {
    p <- NS_params_small
    second_order_w(p) <- c(flux = "van_leer")
    p <- setComponent(p, "clock", initial_value = 0,
                      dynamics_fun = "constant_other")
    fns <- projectRateFunctions(p)
    original_rates <- fns$Rates
    observed <- new.env(parent = emptyenv())
    observed$times <- numeric()
    fns$Rates <- function(params, t, ...) {
        observed$times <- c(observed$times, t)
        r <- original_rates(params, t = t, ...)
        r$e_growth[] <- 0
        r$e_growth[, 2] <- 100 * p@dw[1]
        r
    }
    components <- list(clock = function(params, n_other, component,
                                        t, dt, ...) {
        n_other[[component]] + dt * as.numeric(t >= 0.1)
    })
    out <- project_simple(p, dt = 0.1, steps = 2, rates_fns = fns,
                          other_dynamics_fns = components,
                          subcycle = TRUE, info_level = 0)
    expect_identical(observed$times[11], 0.1)
    expect_equal(out$n_other$clock, 0.1, tolerance = 1e-14)
})

# Component equations and complete numerical contexts ----
test_that("component derivatives have Euler and Heun temporal order", {
    assign("projection_decay_rate", function(params, n_other, component, ...) {
        -2 * n_other[[component]]
    }, envir = .GlobalEnv)
    withr::defer(rm("projection_decay_rate", envir = .GlobalEnv))
    initial <- c(first = 1, second = 2)
    p <- setComponent(NS_params_small, "pool", initial_value = initial,
                      rate_fun = "projection_decay_rate")
    exact <- initial * exp(-0.8)
    for (method in c("euler", "predictor_corrector", "tr_bdf2")) {
        errors <- vapply(c(0.1, 0.05, 0.025), function(dt) {
            out <- project_simple(p, dt = dt, steps = round(0.4 / dt),
                                  method = method, info_level = 0)
            expect_identical(names(out$n_other$pool), names(initial))
            max(abs(out$n_other$pool - exact))
        }, numeric(1))
        order <- log2(errors[-3] / errors[-1])
        target <- if (method == "euler") 1 else 2
        expect_true(all(order > target - 0.15 & order < target + 0.2),
                    info = paste(method,
                                 paste(round(order, 3), collapse = ", ")))
    }
})

test_that("Heun uses the prediction of every coupled component", {
    oscillator_rate <- function(params, n_other, component, ...) {
        if (component == "x") n_other$y else -n_other$x
    }
    assign("projection_oscillator_rate", oscillator_rate, envir = .GlobalEnv)
    withr::defer(rm("projection_oscillator_rate", envir = .GlobalEnv))
    p <- setComponent(NS_params_small, "x", initial_value = 1,
                      rate_fun = "projection_oscillator_rate")
    p <- setComponent(p, "y", initial_value = 0,
                      rate_fun = "projection_oscillator_rate")
    exact <- c(x = cos(0.4), y = -sin(0.4))
    for (method in c("euler", "predictor_corrector", "tr_bdf2")) {
        errors <- vapply(c(0.1, 0.05, 0.025), function(dt) {
            out <- project_simple(p, dt = dt, steps = round(0.4 / dt),
                                  method = method, info_level = 0)
            max(abs(unlist(out$n_other) - exact))
        }, numeric(1))
        order <- log2(errors[-3] / errors[-1])
        target <- if (method == "euler") 1 else 2
        expect_true(all(order > target - 0.15 & order < target + 0.2),
                    info = paste(method,
                                 paste(round(order, 3), collapse = ", ")))
    }
})

test_that("custom correctors receive all midpoint rates", {
    p <- setComponent(NS_params_small, "tracked", initial_value = 0,
                      dynamics_fun = "constant_other")
    fns <- projectRateFunctions(p)
    original_rates <- fns$Rates
    fns$Rates <- function(params, t, ...) {
        rates <- original_rates(params, t = t, ...)
        rates$f_mort[] <- 1 + t
        rates$pred_rate[] <- 2 + 2 * t
        rates
    }
    components <- list(tracked = function(n_other, rates, dt, component, ...) {
        n_other[[component]] + dt *
            (mean(rates$f_mort) + mean(rates$pred_rate))
    })
    resource <- function(n_pp, rates, dt, ...) {
        n_pp[] <- dt * (mean(rates$f_mort) + mean(rates$pred_rate))
        n_pp
    }
    for (method in c("predictor_corrector", "tr_bdf2")) {
        out <- project_simple(p, dt = 0.2, steps = 1, method = method,
                              rates_fns = fns, other_dynamics_fns = components,
                              resource_dynamics_fn = resource, info_level = 0)
        # Integral of (1 + t) + (2 + 2*t) between 0 and 0.2.
        expect_equal(out$n_other$tracked, 0.66)
        expect_equal(unname(out$n_pp), rep(0.66, length(p$w_full)))
    }
})

test_that("custom steps receive predictions without changing the start state", {
    p <- setComponent(NS_params_small, "tracked", initial_value = 1,
                      dynamics_fun = "constant_other")
    observations <- list()
    resource_observations <- list()
    components <- list(tracked = function(params, n, n_pp, n_other, rates,
                                          t, dt, component,
                                          predicted_state = NULL,
                                          endpoint_rates = NULL) {
        observations[[length(observations) + 1L]] <<-
            list(start = n_other[[component]], t = t, dt = dt,
                 prediction = predicted_state, endpoint = endpoint_rates)
        if (is.null(predicted_state)) {
            n_other[[component]] + dt
        } else {
            predicted_state$n_other[[component]]
        }
    })
    resource <- function(params, n, n_pp, n_other, rates, t, dt,
                         resource_rate, resource_capacity,
                         predicted_state = NULL, endpoint_rates = NULL) {
        resource_observations[[length(resource_observations) + 1L]] <<-
            list(prediction = predicted_state, endpoint = endpoint_rates)
        n_pp
    }
    for (method in c("predictor_corrector", "tr_bdf2")) {
        observations <- list()
        resource_observations <- list()
        out <- project_simple(p, dt = 0.2, steps = 1, method = method,
                              other_dynamics_fns = components,
                              resource_dynamics_fn = resource, info_level = 0)
        expect_length(observations, 2L)
        expect_null(observations[[1]]$prediction)
        expect_null(observations[[1]]$endpoint)
        expect_identical(observations[[2]]$start, 1)
        expect_equal(observations[[2]]$prediction$n_other$tracked, 1.2)
        expect_equal(observations[[2]]$prediction$t, 0.2)
        expect_named(observations[[2]]$prediction,
                     c("n", "n_pp", "n_other", "t"))
        expect_true(all(c("f_mort", "pred_rate", "rdd") %in%
                            names(observations[[2]]$endpoint)))
        expect_length(resource_observations, 2L)
        expect_null(resource_observations[[1]]$prediction)
        expect_null(resource_observations[[1]]$endpoint)
        expect_equal(resource_observations[[2]]$prediction,
                     observations[[2]]$prediction)
        expect_equal(resource_observations[[2]]$endpoint,
                     observations[[2]]$endpoint)
        expect_equal(out$n_other$tracked, 1.2)
    }
})

test_that("legacy strict callbacks keep arbitrary component state", {
    p <- setComponent(NS_params_small, "record", initial_value = list(time = 0),
                      dynamics_fun = "constant_other")
    # Every original argument, but no ellipsis or new optional arguments.
    strict <- function(params, n, n_pp, n_other, rates, t, dt, component) {
        list(time = t + dt, note = "retained")
    }
    for (method in c("euler", "tr_bdf2")) {
        out <- project_simple(p, dt = 0.1, steps = 2, method = method,
                              other_dynamics_fns = list(record = strict),
                              info_level = 0)
        expect_identical(names(out$n_other), "record")
        expect_equal(out$n_other$record,
                     list(time = 0.2, note = "retained"))
    }
})

test_that("custom steps take precedence over a component equation", {
    assign("projection_unused_rate", function(...) {
        stop("the projection must use the custom step")
    }, envir = .GlobalEnv)
    withr::defer(rm("projection_unused_rate", envir = .GlobalEnv))
    p <- setComponent(NS_params_small, "pool", initial_value = 42,
                      dynamics_fun = "constant_other",
                      rate_fun = "projection_unused_rate")
    sim <- project(p, dt = 0.1, t_max = 0.2, t_save = 0.1,
                   method = "tr_bdf2", progress_bar = FALSE, info_level = 0)
    expect_equal(unname(unlist(sim$n_other[, "pool"])), rep(42, 3))
    integration <- sim$sim_params$component_integration
    expect_identical(integration$component, "pool")
    expect_identical(integration$method, "custom")
    expect_true(is.na(integration$order))
})

test_that("equation components use the actual substep length", {
    subcycle_rate <- function(params, n_other, component, ...) {
        -n_other[[component]]
    }
    assign("projection_subcycle_rate", subcycle_rate, envir = .GlobalEnv)
    withr::defer(rm("projection_subcycle_rate", envir = .GlobalEnv))
    p <- setComponent(NS_params_small, "pool", initial_value = 1,
                      rate_fun = "projection_subcycle_rate")
    second_order_w(p) <- c(flux = "van_leer")
    fns <- projectRateFunctions(p)
    original_rates <- fns$Rates
    fns$Rates <- function(params, ...) {
        rates <- original_rates(params, ...)
        rates$e_growth[] <- 0
        rates$e_growth[, 2] <- 20 * params$dw[1]
        rates
    }
    for (method in c("euler", "tr_bdf2")) {
        adaptive <- project_simple(p, dt = 0.2, steps = 2, method = method,
                                   rates_fns = fns, subcycle = TRUE,
                                   info_level = 0)
        fine <- project_simple(p, dt = 0.05, steps = 8, method = method,
                               rates_fns = fns, info_level = 0)
        expect_equal(adaptive$n_other, fine$n_other, tolerance = 1e-14)
    }
})

test_that("simulations record actual component integration methods", {
    p <- setComponent(NS_params_small, "equation", initial_value = 1,
                      rate_fun = "constant_other")
    p <- setComponent(p, "custom", initial_value = 42,
                      dynamics_fun = "constant_other", integration_order = 2)
    for (method in c("euler", "tr_bdf2")) {
        sim <- project(p, dt = 0.1, t_max = 0.1, t_save = 0.1,
                       method = method, progress_bar = FALSE, info_level = 0)
        integration <- sim$sim_params$component_integration
        expect_identical(integration$component, c("equation", "custom"))
        expected <- if (method == "euler") "euler" else "heun"
        expect_identical(integration$method, c(expected, "custom"))
        order <- c(if (method == "euler") 1L else 2L, 2L)
        expect_equal(integration$order, order)
    }
})
