params <- newMultispeciesParams(NS_species_params_gears, inter, info_level = 0)

# time dimension ----
test_that("time dimension is dealt with properly", {

    # Effort is a single numeric
    t_max <- 5
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
    effort <- c(Industrial = 1, Pelagic = 0.5, Beam = 0.3, Otter = 0)
    t_max <- 5
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
    NS_species_params_gears$w_min <- 0.001
    NS_species_params_gears$w_min[1] <- 1
    params2 <- newMultispeciesParams(NS_species_params_gears, inter, info_level = 0)
    sim <- project(params2, effort = 1, t_max = 5)
    expect_equal(sim@n[6, 1, 1:(sim@params@w_min_idx[1] - 1)],
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
  	sim <- project(params_gear, effort = 0.3, t_max = 10)
  	expect_true(all(sim@effort == 0.3))
    # Also check that order of gear names in resulting effort matches catchability
    expect_true(all(dimnames(sim@effort)$gear == gear_names))
    # Effort vector
    # Should give same result
    effort_vec <- c(Other = 1, Industrial = 0)
    effort_vec2 <- c(Industrial = 0, Other = 1)
    sim <- project(params_gear, effort = effort_vec, t_max = 10)
    sim2 <- project(params_gear, effort = effort_vec2, t_max = 10)
    expect_true(all(sim@effort[, "Industrial"] == 0))
    expect_true(all(sim@effort[, "Other"] == 1))
    expect_true(all(sim2@effort[, "Industrial"] == 0))
    expect_true(all(sim2@effort[, "Other"] == 1))
    expect_true(all(dimnames(sim@effort)$gear == gear_names)) 
    expect_true(all(dimnames(sim2@effort)$gear == gear_names)) 
    # Should fail - number of gears wrong
    effort_vec3 <- c(Industrial = 0, Other = 1, Dummy = 0.5)
    expect_error(project(params_gear, effort = effort_vec3, t_max = 10))
    effort_vec4 <- c(Industrial = 0) # Is OK because that gear exists
    expect_error(project(params_gear, effort = effort_vec4, t_max = 10), NA)
    # Should fail - names of gears wrong
    effort_vec5 <- c(Industrial = 0, Dummy = 1)
    expect_error(project(params_gear, effort = effort_vec5, t_max = 10))
    # Array effort
    t_steps <- 10
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
  params <- newMultispeciesParams(NS_species_params_gears, inter,
                                  n = 2/3, p = 0.7, lambda = 2.8 - 2/3, info_level = 0)
  sim <- project(params, t_max = 1)
  # expect_known_value(sim@n[2, 3, ], "values/projectn")
  # expect_known_value(sim@n_pp[2, ], "values/projectp")
  expect_snapshot(sim@n[2, 3, ])
  expect_snapshot(sim@n_pp[2, ])
  
})

test_that("Final result the same when called with sim or params", {
  params <- NS_params
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
  gear_names <- unique(gear_params(NS_params)$gear)
  effort <- array(1, dim = c(3, length(gear_names)), 
                  dimnames = list(1:3,
                                  gear_names))
  sim <- project(NS_params, effort, t_max = 0.1)
  expect_identical(names(dimnames(sim@effort)), c("time", "gear"))
})

# t_max and t_save with effort arrays ----
test_that("t_max extends simulation beyond effort array", {
  gear_names <- unique(gear_params(NS_params)$gear)
  # Create effort array that goes from year 1 to year 5
  effort <- array(0.5, dim = c(5, length(gear_names)),
                  dimnames = list(time = 1:5, gear = gear_names))
  
  # Extend simulation to year 10 with t_max
  sim <- project(NS_params, effort = effort, t_max = 9, dt = 0.1)
  
  # Should save at times 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
  expect_equal(as.numeric(dimnames(sim@n)[[1]]), 1:10)
  expect_equal(dim(sim@n)[1], 10)
  
  # Effort should be extrapolated with last value (0.5)
  expect_true(all(sim@effort[, ] == 0.5))
})

test_that("t_save controls save frequency with effort array", {
  gear_names <- unique(gear_params(NS_params)$gear)
  # Create effort array with times 0, 1, 2, 3, 4, 5
  effort <- array(0.5, dim = c(6, length(gear_names)),
                  dimnames = list(time = 0:5, gear = gear_names))
  
  # Use t_save to control output frequency
  sim <- project(NS_params, effort = effort, t_save = 2, dt = 0.1)
  
  # Should save at times 0, 2, 4
  expect_equal(as.numeric(dimnames(sim@n)[[1]]), c(0, 2, 4))
  expect_equal(dim(sim@n)[1], 3)
})

test_that("t_max and t_save work together with effort array", {
  gear_names <- unique(gear_params(NS_params)$gear)
  # Create effort array from 0 to 3
  effort <- array(0.5, dim = c(4, length(gear_names)),
                  dimnames = list(time = 0:3, gear = gear_names))
  
  # Extend to 6 and save every 2 years
  sim <- project(NS_params, effort = effort, t_max = 6, t_save = 2, dt = 0.1)
  
  # Should save at times 0, 2, 4, 6
  expect_equal(as.numeric(dimnames(sim@n)[[1]]), c(0, 2, 4, 6))
  expect_equal(dim(sim@n)[1], 4)
})

test_that("Effort array times used when t_max and t_save not provided", {
  gear_names <- unique(gear_params(NS_params)$gear)
  # Create effort array with irregular times
  effort <- array(0.5, dim = c(4, length(gear_names)),
                  dimnames = list(time = c(0, 1, 3, 7), gear = gear_names))
  
  # Without t_max or t_save, should use effort array times
  sim <- project(NS_params, effort = effort, dt = 0.1)
  
  expect_equal(as.numeric(dimnames(sim@n)[[1]]), c(0, 1, 3, 7))
  expect_equal(dim(sim@n)[1], 4)
})

test_that("Effort values interpolated correctly", {
  gear_names <- unique(gear_params(NS_params)$gear)
  # Create effort array with varying effort
  effort <- array(NA, dim = c(3, length(gear_names)),
                  dimnames = list(time = c(0, 5, 10), gear = gear_names))
  effort[, 1] <- c(0, 0.5, 1.0)  # Industrial
  effort[, 2:4] <- 0.3  # Others constant
  
  # Use t_save to create intermediate time points
  sim <- project(NS_params, effort = effort, t_save = 2.5, dt = 0.1)
  
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
  gear_names <- unique(gear_params(NS_params)$gear)
  effort <- array(0.5, dim = c(3, length(gear_names)),
                  dimnames = list(time = 2017:2019, gear = gear_names))
  # The last year has NA which gets replaced by default (1 for edition >= 2)
  effort[3, ] <- NA
  
  # Run until 2020 without needing to specify effort for 2019
  sim <- project(NS_params, effort = effort, t_max = 3, dt = 0.1)
  
  # Should have years 2017, 2018, 2019, 2020
  expect_equal(as.numeric(dimnames(sim@n)[[1]]), 2017:2020)
  expect_equal(dim(sim@n)[1], 4)
  
  # The NA gets replaced by default effort during validation
  # Then extrapolated for 2020
  expect_true(all(!is.na(sim@effort)))
})

test_that("t_max less than effort array duration uses effort times", {
  gear_names <- unique(gear_params(NS_params)$gear)
  effort <- array(0.5, dim = c(6, length(gear_names)),
                  dimnames = list(time = 0:5, gear = gear_names))
  
  # t_max = 3 should run until time 3, not 5
  sim <- project(NS_params, effort = effort, t_max = 3, dt = 0.1)
  
  # Should stop at year 3
  expect_equal(max(as.numeric(dimnames(sim@n)[[1]])), 3)
})
