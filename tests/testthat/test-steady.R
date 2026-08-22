# projectUntilSettled ----
# Built once here; tests get a fresh copy via R's copy-on-modify semantics.
small_steady_params <- suppressMessages(
    newTraitParams(no_sp = 2, no_w = 20, max_w_max = 100,
                   min_w = 1e-3, w_pp_cutoff = 5, ks = 4,
                   reproduction_level = 0.25, info_level = 0)
)
test_that("projectUntilSettled() works", {
    params <- small_steady_params
    effort <- params@initial_effort * 1.1
    initialN(params)[1, ] <- initialN(params)[1, ] * 3
    expect_error(projectUntilSettled(params, dt = 1, t_per = 0.5),
                 "t_per must be a positive multiple of dt")
    expect_error(projectUntilSettled(params, t_max = 0.1),
                 "t_max not greater than or equal to t_per")
    expect_message(projectUntilSettled(params, t_max = 0.1, t_per = 0.1,
                                       info_level = 3),
                   "Simulation run did not converge after 0.1 years.")
    # It always returns the trajectory, whatever the run settled on.
    expect_message(sim <- projectUntilSettled(params, t_per = 1, dt = 1,
                                              tol = 1000, effort = effort),
                   "Reached the convergence tolerance")
    expect_s4_class(sim, "MizerSim")
    expect_identical(sim@params@initial_effort, effort)
    # shouldn't take long the second time we run to steady
    expect_message(projectUntilSettled(finalParams(sim), t_per = 1, dt = 1,
                                       tol = 1000),
                   "Reached the convergence tolerance")

    # Alternative distance function
    expect_message(sim <- projectUntilSettled(params, t_per = 1, dt = 1,
                                              distance_func = distanceMaxRelRDI,
                                              tol = 10),
                   "Reached the convergence tolerance")
    # shouldn't take long the second time we run to steady
    expect_message(projectUntilSettled(finalParams(sim), t_per = 1, dt = 1,
                                       distance_func = distanceMaxRelRDI,
                                       tol = 10),
                   "Reached the convergence tolerance")

    # Check extinction
    params@psi[1:2, ] <- 0
    sp1 <- params@species_params$species[1]
    sp2 <- params@species_params$species[2]
    expect_warning(projectUntilSettled(params) |> suppressMessages(),
                   paste0(sp1, ", ", sp2, " are going extinct."))
})

test_that("projectUntilSettled accepts the documented effort forms", {
    params <- NS_params_small
    s1 <- projectUntilSettled(params, effort = 0.5, t_per = 0.1,
                              t_max = 0.1, tol = 10, info_level = 0)
    expect_equal(unname(s1@params@initial_effort),
                 rep(0.5, length(initial_effort(params))))

    named <- initial_effort(params)[1:2]
    named[] <- c(0.2, NA)
    s2 <- projectUntilSettled(params, effort = named, t_per = 0.1,
                              t_max = 0.1, tol = 10, info_level = 0)
    expected <- validEffortVector(named, params)
    expect_equal(s2@params@initial_effort, expected)
})
test_that("projectUntilSettled accepts consumer update method", {
    params <- small_steady_params

    sim <- projectUntilSettled(params, t_per = 1, dt = 1,
                               t_max = 1, tol = 1000,
                               method = "predictor_corrector", info_level = 0)
    expect_s4_class(sim, "MizerSim")
    expect_true(all(is.finite(finalN(sim))))
    expect_true(all(finalN(sim) >= 0))

    expect_error(
        projectUntilSettled(params, t_per = 1, dt = 1, t_max = 1,
                            method = "bogus"),
        "should be one of"
    )
})

# valid_species_arg ----
test_that("valid_species_arg works", {
    # character species argument
    expect_warning(s <- valid_species_arg(NS_params_small, c("non", "sense")),
                   "The following species do not exist: non, sense")
    expect_identical(s, vector(mode = "character"))

    sp1 <- NS_params_small@species_params$species[3]
    sp2 <- NS_params_small@species_params$species[2]
    sp_sprat <- NS_params_small@species_params$species[1]
    sp3 <- NS_params_small@species_params$species[3]

    expect_identical(valid_species_arg(NS_params_small, c(sp1, sp2)),
                     c(sp1, sp2))
    expect_identical(valid_species_arg(NS_params_small, c(sp_sprat, sp2),
                                       return.logical = TRUE),
                     c(TRUE, TRUE, FALSE))
    expect_error(
        valid_species_arg(NS_params_small, c("non", "sense"), error_on_empty = TRUE) |>
            suppressWarnings(),
        "No species have been selected.")
    # numeric species argument
    expect_warning(s <- valid_species_arg(NS_params_small, c(2.5, 4)),
                 "A numeric 'species' argument should only contain the integers 1 to 3")
    expect_identical(s, vector(mode = "character"))
    expect_identical(valid_species_arg(NS_params_small, c(3, 1)),
                     c(sp3, sp_sprat))
    expect_identical(valid_species_arg(NS_params_small, c(1, 3)),
                     c(sp_sprat, sp3))
    expect_identical(valid_species_arg(NS_params_small, c(3, 1),
                                       return.logical = TRUE),
                     c(TRUE, FALSE, TRUE))
    expect_error(
        valid_species_arg(NS_params_small, 4, error_on_empty = TRUE) |>
            suppressWarnings(),
        "No species have been selected.")
    # logical species argument
    expect_error(valid_species_arg(NS_params_small, c(TRUE, FALSE)),
                 "The boolean `species` argument has the wrong length")
    expect_identical(valid_species_arg(NS_params_small,
                                       c(TRUE, FALSE, TRUE)),
                     c(sp_sprat, sp3))
    expect_identical(valid_species_arg(NS_params_small,
                                       c(TRUE, FALSE, TRUE),
                                       return.logical = TRUE),
                     c(TRUE, FALSE, TRUE))
    expect_error(
        valid_species_arg(NS_params_small, rep(FALSE, 3), error_on_empty = TRUE) |>
            suppressWarnings(),
        "No species have been selected.")
    # called with MizerSim object
    sim <- project(NS_params_small, t_max = 1, dt = 1)
    expect_identical(valid_species_arg(sim, sp1),
                     valid_species_arg(NS_params_small, sp1))
    # called without species
    expect_identical(valid_species_arg(NS_params_small),
                     valid_species_arg(NS_params_small, 
                                       NS_params_small@species_params$species))
})
test_that("projectUntilSettled() converges with use_predation_diffusion", {
    params_d <- small_steady_params
    params_d@use_predation_diffusion <- TRUE
    initialN(params_d)[1, ] <- initialN(params_d)[1, ] * 3
    expect_message(
        projectUntilSettled(params_d, t_per = 1, dt = 1, tol = 1000),
        "Reached the convergence tolerance"
    )
})
test_that("valid_gears_arg works", {
    all_gears <- unique(NS_params_small@gear_params$gear)
    expect_identical(valid_gears_arg(NS_params_small), all_gears)
    expect_identical(valid_gears_arg(NS_params_small, all_gears[2:1]), all_gears[2:1])
    expect_warning(
        gears <- valid_gears_arg(NS_params_small, c("nope", all_gears[1])),
        "The following gears do not exist: nope"
    )
    expect_identical(gears, all_gears[1])
    expect_error(valid_gears_arg(NS_params_small, "nope", error_on_empty = TRUE) |>
                     suppressWarnings(),
                 "No gears have been selected.")
})

# constant_other ----
test_that("constant_other returns component value", {
    params <- NS_params_small
    # Create a simple n_other list with test components
    n_other <- list(
        component1 = 100,
        component2 = c(1, 2, 3),
        component3 = matrix(1:6, nrow = 2)
    )
    
    # Test that constant_other returns the correct component
    expect_equal(constant_other(params, n_other, "component1"), 100)
    expect_equal(constant_other(params, n_other, "component2"), c(1, 2, 3))
    expect_equal(constant_other(params, n_other, "component3"), 
                 matrix(1:6, nrow = 2))
})

# distance functions ----
test_that("distance functions implement their documented formulas", {
    params <- NS_params_small
    previous <- list(
        n = initialN(params),
        n_pp = initialNResource(params),
        n_other = list()
    )
    current <- list(
        n = initialN(params) * 1.2,
        n_pp = initialNResource(params),
        n_other = list()
    )

    current_rdi <- getRDI(params, n = current$n, n_pp = current$n_pp,
                          n_other = current$n_other)
    previous_rdi <- getRDI(params, n = previous$n, n_pp = previous$n_pp,
                           n_other = previous$n_other)
    expected_rel <- max(abs((current_rdi - previous_rdi) / previous_rdi))
    dist_rel <- distanceMaxRelRDI(params, current, previous)
    expect_equal(dist_rel, expected_rel)
    expect_true(dist_rel >= 0)
    expect_equal(length(dist_rel), 1)
    expect_equal(distanceMaxRelRDI(params, previous, previous), 0)

    current$n[1, 1] <- 0
    previous$n[1, 2] <- 0
    sel <- current$n > 0 & previous$n > 0
    expected_log <- sum((log(current$n[sel]) - log(previous$n[sel]))^2)
    dist_log <- distanceSSLogN(params, current, previous)
    expect_equal(dist_log, expected_log)
    expect_true(dist_log >= 0)
    expect_equal(length(dist_log), 1)
    expect_equal(distanceSSLogN(params, previous, previous), 0)
})

# limit cycle detection ----

# Limit-cycle detection and the "convergence" attribute --------------------

# Helper functions -----------------------------------------------------------

test_that("find_first_acf_peak returns the lag of the first local maximum", {
    # ac[1] is lag 0; a peak at position k is lag k - 1.
    ac <- c(1, 0.2, -0.3, 0.1, 0.8, 0.4)   # local max at position 5 (lag 4)
    expect_equal(find_first_acf_peak(ac, threshold = 0.5), 4L)
    # No peak above the threshold
    expect_true(is.na(find_first_acf_peak(c(1, 0.4, 0.2, 0.1), threshold = 0.5)))
    # Too short
    expect_true(is.na(find_first_acf_peak(c(1, 0.5), threshold = 0.5)))
})
test_that("amp_window gives the largest per-species relative amplitude", {
    bio <- cbind(c(1, 2, 3),        # (3 - 1) / 2 = 1
                 c(10, 10, 10),     # 0
                 c(4, 6, 8))        # (8 - 4) / 6 = 0.667
    expect_equal(amp_window(bio), 1)
    # A species with zero mean contributes zero, not NaN/Inf
    expect_equal(amp_window(cbind(c(0, 0, 0))), 0)
})
test_that("detect_limit_cycle finds a sustained oscillation", {
    t_save <- 0.25; period <- 5; a <- 0.1
    idx <- 0:199
    bio <- matrix(1000 * (1 + a * sin(2 * pi * idx * t_save / period)), ncol = 1)
    cyc <- detect_limit_cycle(bio, t_save, amplitude_tol = 0.01)
    expect_type(cyc, "list")
    expect_equal(cyc$period, period)
    # For a sinusoid about the mean the relative peak-to-trough amplitude is 2a.
    expect_equal(cyc$amplitude, 2 * a, tolerance = 1e-3)
})
test_that("detect_limit_cycle rejects non-cycles", {
    t_save <- 0.25; period <- 5; a <- 0.1
    idx <- 0:199
    # A decaying oscillation (spiral toward a fixed point) is not a cycle.
    bio_decay <- matrix(1000 * (1 + a * exp(-0.03 * idx) *
                                    sin(2 * pi * idx * t_save / period)), ncol = 1)
    expect_null(detect_limit_cycle(bio_decay, t_save, amplitude_tol = 0.01))
    # A flat series is not a cycle.
    expect_null(detect_limit_cycle(matrix(1000, nrow = 200, ncol = 1),
                                   t_save, amplitude_tol = 0.01))
    # Too few samples.
    expect_null(detect_limit_cycle(matrix(1000, nrow = 10, ncol = 1),
                                   t_save, amplitude_tol = 0.01))
})
test_that("detect_limit_cycle respects amplitude_tol", {
    t_save <- 0.25; period <- 5; a <- 0.03   # ~6% peak-to-trough amplitude
    idx <- 0:199
    bio <- matrix(1000 * (1 + a * sin(2 * pi * idx * t_save / period)), ncol = 1)
    # A 6% cycle counts when the floor is 1% but not when it is 10%.
    expect_type(detect_limit_cycle(bio, t_save, amplitude_tol = 0.01), "list")
    expect_null(detect_limit_cycle(bio, t_save, amplitude_tol = 0.1))
})

# Fixtures for the integration tests -----------------------------------------

cd_params <- suppressMessages(
    newTraitParams(no_sp = 2, no_w = 20, max_w_max = 100,
                   min_w = 1e-3, w_pp_cutoff = 5, ks = 4,
                   reproduction_level = 0.25, info_level = 0)
)
test_that("projectUntilSettled() reports non-convergence", {
    p <- cd_params
    initialN(p)[1, ] <- initialN(p)[1, ] * 3
    sim <- suppressMessages(
        projectUntilSettled(p, t_max = 0.5, t_per = 0.5, dt = 0.1,
                            tol = 1e-12, info_level = 0)
    )
    conv <- attr(sim, "convergence")
    expect_identical(conv$type, "not_converged")
    expect_false(conv$settled)
})
test_that("fine t_save sampling does not change the result", {
    # Sub-blocking the run at the t_save resolution must be numerically
    # identical to stepping a whole t_per block at once.
    args <- list(t_max = 3, t_per = 1.5, dt = 0.5, tol = 1e-12, info_level = 0)
    p_fine   <- suppressMessages(do.call(findSteadyState,
                                         c(list(cd_params, t_save = 0.5), args)))
    p_coarse <- suppressMessages(do.call(findSteadyState,
                                         c(list(cd_params, t_save = 1.5), args)))
    expect_identical(p_fine@initial_n, p_coarse@initial_n)
})
test_that("t_save is validated", {
    expect_error(projectUntilSettled(cd_params, dt = 0.1, t_save = 0.15),
                 "t_save must be a positive multiple of dt")
    expect_error(projectUntilSettled(cd_params, dt = 0.1, t_per = 1,
                                     t_save = 0.3),
                 "t_per must be a positive multiple of t_save")
})
test_that("projectUntilSettled() detects a limit cycle", {
    skip_on_cran()
    # The full North Sea model driven at high fishing effort settles onto a
    # limit cycle rather than a fixed point.
    sim <- suppressMessages(
        projectUntilSettled(NS_params, effort = 2, t_max = 200,
                            t_per = 1.5, dt = 0.1, tol = 1e-8, info_level = 0)
    )
    conv <- attr(sim, "convergence")
    expect_identical(conv$type, "cycle")
    expect_true(conv$settled)
    expect_gt(conv$period, 0)
    expect_gt(conv$amplitude, 0.1)
})