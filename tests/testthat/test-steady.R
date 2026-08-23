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
    expect_error(projectUntilSettled(params, dt = 1, t_check = 0.5),
                 "t_check must be a positive multiple of dt")
    expect_error(projectUntilSettled(params, t_max = 0.1),
                 "t_max not greater than or equal to t_check")
    # The default check interval follows `dt`, so a `dt` that does not divide
    # 1.5 is no longer in conflict with it.
    expect_no_error(suppressMessages(
        projectUntilSettled(params, dt = 0.2, t_max = 3, distance_tol = 1e3,
                            residual_tol = Inf, info_level = 0)))
    expect_message(projectUntilSettled(params, t_max = 0.1, t_check = 0.1,
                                       info_level = 3),
                   "Simulation run did not converge after 0.1 years.")
    # It always returns the trajectory, whatever the run settled on. These runs
    # take a crude `dt = 1` step and stop on an absurdly loose `tol` after a
    # single block, which is the cheapest way to exercise the plumbing. That is
    # not a converged state, and the drift criterion would rightly say so, so it
    # is switched off here rather than the tolerance being pretended tight.
    expect_message(sim <- projectUntilSettled(params, t_check = 1, dt = 1,
                                              distance_tol = 1000, residual_tol = Inf,
                                              effort = effort),
                   "Reached the convergence tolerance")
    expect_s4_class(sim, "MizerSim")
    expect_identical(sim@params@initial_effort, effort)
    # shouldn't take long the second time we run to steady
    expect_message(projectUntilSettled(finalParams(sim), t_check = 1, dt = 1,
                                       distance_tol = 1000, residual_tol = Inf),
                   "Reached the convergence tolerance")

    # Alternative distance function
    expect_message(sim <- projectUntilSettled(params, t_check = 1, dt = 1,
                                              distance_func = distanceMaxRelRDI,
                                              distance_tol = 10, residual_tol = Inf),
                   "Reached the convergence tolerance")
    # shouldn't take long the second time we run to steady
    expect_message(projectUntilSettled(finalParams(sim), t_check = 1, dt = 1,
                                       distance_func = distanceMaxRelRDI,
                                       distance_tol = 10, residual_tol = Inf),
                   "Reached the convergence tolerance")

    # Check extinction
    params@psi[1:2, ] <- 0
    sp1 <- params@species_params$species[1]
    sp2 <- params@species_params$species[2]
    expect_warning(sim_ext <- projectUntilSettled(params) |> suppressMessages(),
                   paste0(sp1, ", ", sp2, " are going extinct."))
    conv_ext <- attr(sim_ext, "convergence")
    expect_identical(conv_ext$termination, "extinction")
    expect_identical(conv_ext$extinct, c(sp1, sp2))
})

test_that("projectUntilSettled accepts the documented effort forms", {
    params <- NS_params_small
    s1 <- projectUntilSettled(params, effort = 0.5, t_check = 0.1,
                              t_max = 0.1, distance_tol = 10, info_level = 0)
    expect_equal(unname(s1@params@initial_effort),
                 rep(0.5, length(initial_effort(params))))

    named <- initial_effort(params)[1:2]
    named[] <- c(0.2, NA)
    s2 <- projectUntilSettled(params, effort = named, t_check = 0.1,
                              t_max = 0.1, distance_tol = 10, info_level = 0)
    expected <- validEffortVector(named, params)
    expect_equal(s2@params@initial_effort, expected)
})
test_that("projectUntilSettled accepts consumer update method", {
    params <- small_steady_params

    sim <- projectUntilSettled(params, t_check = 1, dt = 1,
                               t_max = 1, distance_tol = 1000,
                               method = "predictor_corrector", info_level = 0)
    expect_s4_class(sim, "MizerSim")
    expect_true(all(is.finite(finalN(sim))))
    expect_true(all(finalN(sim) >= 0))

    expect_error(
        projectUntilSettled(params, t_check = 1, dt = 1, t_max = 1,
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
        projectUntilSettled(params_d, t_check = 1, dt = 1, distance_tol = 1000,
                            residual_tol = Inf),
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
    t_sample <- 0.25; period <- 5; a <- 0.1
    idx <- 0:199
    bio <- matrix(1000 * (1 + a * sin(2 * pi * idx * t_sample / period)), ncol = 1)
    cyc <- detect_limit_cycle(bio, t_sample, amplitude_tol = 0.01)
    expect_type(cyc, "list")
    expect_equal(cyc$period, period)
    # For a sinusoid about the mean the relative peak-to-trough amplitude is 2a.
    expect_equal(cyc$amplitude, 2 * a, tolerance = 1e-3)
})
test_that("detect_limit_cycle rejects non-cycles", {
    t_sample <- 0.25; period <- 5; a <- 0.1
    idx <- 0:199
    # A decaying oscillation (spiral toward a fixed point) is not a cycle.
    bio_decay <- matrix(1000 * (1 + a * exp(-0.03 * idx) *
                                    sin(2 * pi * idx * t_sample / period)), ncol = 1)
    expect_null(detect_limit_cycle(bio_decay, t_sample, amplitude_tol = 0.01))
    # A flat series is not a cycle.
    expect_null(detect_limit_cycle(matrix(1000, nrow = 200, ncol = 1),
                                   t_sample, amplitude_tol = 0.01))
    # Too few samples.
    expect_null(detect_limit_cycle(matrix(1000, nrow = 10, ncol = 1),
                                   t_sample, amplitude_tol = 0.01))
})
test_that("detect_limit_cycle respects amplitude_tol", {
    t_sample <- 0.25; period <- 5; a <- 0.03   # ~6% peak-to-trough amplitude
    idx <- 0:199
    bio <- matrix(1000 * (1 + a * sin(2 * pi * idx * t_sample / period)), ncol = 1)
    # A 6% cycle counts when the floor is 1% but not when it is 10%.
    expect_type(detect_limit_cycle(bio, t_sample, amplitude_tol = 0.01), "list")
    expect_null(detect_limit_cycle(bio, t_sample, amplitude_tol = 0.1))
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
        projectUntilSettled(p, t_max = 0.5, t_check = 0.5, dt = 0.1,
                            distance_tol = 1e-12, info_level = 0)
    )
    conv <- attr(sim, "convergence")
    expect_identical(conv$termination, "time_limit")
    expect_false(conv$converged)
    expect_true(is.na(conv$attractor))
    expect_identical(conv$extinct, character(0))
})
test_that("stepping one dt at a time reproduces project()", {
    # The run advances a single time step at a time so that the check and save
    # grids need no common divisor. That must not change the trajectory: an
    # unconverged run has to land exactly where project() lands.
    p <- cd_params
    initialN(p)[1, ] <- initialN(p)[1, ] * 3
    sim_settled <- suppressMessages(
        projectUntilSettled(p, t_max = 3, t_check = 1.5, dt = 0.5,
                            distance_tol = 1e-12, progress_bar = FALSE,
                            info_level = 0)
    )
    sim_project <- project(p, t_max = 3, dt = 0.5, t_save = 1,
                           progress_bar = FALSE)
    # ignore_attr because finalN() carries the model it came from along with
    # the numbers, and the two runs report a different `initial_n` on it.
    expect_equal(finalN(sim_settled), finalN(sim_project), ignore_attr = TRUE)
    expect_equal(finalNResource(sim_settled), finalNResource(sim_project),
                 ignore_attr = TRUE)
})
test_that("t_check and t_save need not be commensurate", {
    # `t_save` is the sim save interval, as in project(), and is independent of
    # the interval at which convergence is checked.
    sim <- suppressMessages(
        projectUntilSettled(cd_params, t_max = 6, t_check = 1.5, dt = 0.5,
                            t_save = 1, distance_tol = 1e-12,
                            progress_bar = FALSE, info_level = 0)
    )
    conv <- attr(sim, "convergence")
    # Saves on the t_save grid, with the state the run settled on appended.
    expect_equal(getTimes(sim), c(0, 1, 2, 3, 4, 5, 6))
    expect_equal(unname(getTimes(sim)[length(getTimes(sim))]), conv$years)
    # A run that stops between two saves still ends on the settled state.
    sim2 <- suppressMessages(
        projectUntilSettled(cd_params, t_max = 6, t_check = 1.5, dt = 0.5,
                            t_save = 2, distance_tol = 1e3, residual_tol = Inf,
                            progress_bar = FALSE, info_level = 0)
    )
    conv2 <- attr(sim2, "convergence")
    expect_equal(conv2$years, 1.5)
    expect_equal(getTimes(sim2), c(0, 1.5))
    expect_equal(finalN(sim2), finalParams(sim2)@initial_n, ignore_attr = TRUE)
})
test_that("t_save is validated", {
    expect_error(projectUntilSettled(cd_params, dt = 0.1, t_save = 0.15),
                 "t_save must be a positive multiple of dt")
})
test_that("projectUntilSettled() detects a limit cycle", {
    skip_on_cran()
    # The full North Sea model driven at high fishing effort settles onto a
    # limit cycle rather than a fixed point.
    sim <- suppressMessages(
        projectUntilSettled(NS_params, effort = 2, t_max = 200,
                            t_check = 1.5, dt = 0.1, distance_tol = 1e-8, info_level = 0)
    )
    conv <- attr(sim, "convergence")
    expect_identical(conv$termination, "cycle_detected")
    expect_true(conv$converged)
    expect_identical(conv$attractor, "limit_cycle")
    expect_gt(conv$period, 0)
    expect_gt(conv$amplitude, 0.1)
})
test_that("projectUntilSettled() passes absolute time to the dynamics", {
    # The run advances one `dt` at a time. Every step has to continue the clock
    # rather than restart it, or a rate or component function that depends on
    # `t` sees the same interval over and over and the result stops matching
    # project().
    assign("mizer_test_recorded_times", numeric(0), envir = .GlobalEnv)
    assign("mizer_test_record_time", function(params, n_other, component, t,
                                              dt, ...) {
        # The steady-state diagnostic differences the dynamics over a very
        # short step; those calls are not part of the trajectory.
        if (dt > 1e-3) {
            .GlobalEnv$mizer_test_recorded_times <-
                c(.GlobalEnv$mizer_test_recorded_times, t)
        }
        n_other[[component]]
    }, envir = .GlobalEnv)
    withr::defer(rm("mizer_test_record_time", "mizer_test_recorded_times",
                    envir = .GlobalEnv))

    p <- setComponent(cd_params, "recorder", initial_value = 1,
                      dynamics_fun = "mizer_test_record_time")

    suppressMessages(
        project(p, t_max = 4.5, dt = 0.1, t_save = 1.5, progress_bar = FALSE)
    )
    from_project <- .GlobalEnv$mizer_test_recorded_times

    .GlobalEnv$mizer_test_recorded_times <- numeric(0)
    suppressMessages(
        projectUntilSettled(p, t_max = 4.5, t_check = 1.5, dt = 0.1,
                            t_save = 1.5, distance_tol = 1e-12, progress_bar = FALSE,
                            info_level = 0)
    )
    from_settled <- .GlobalEnv$mizer_test_recorded_times

    expect_equal(from_settled, from_project)
    expect_false(is.unsorted(from_settled))
    expect_gt(max(from_settled), 4)
})

test_that("a state that is still drifting is not called a fixed point", {
    # A distance function that is blind to the motion that is left used to be
    # enough on its own to stop the run and label the state a fixed point.
    blind <- function(params, current, previous, ...) 0
    p <- cd_params
    initialN(p)[1, ] <- initialN(p)[1, ] * 10

    conv <- attr(suppressMessages(
        projectUntilSettled(p, distance_func = blind, t_max = 3, t_check = 1.5,
                            dt = 0.1, progress_bar = FALSE, info_level = 0)
    ), "convergence")
    expect_identical(conv$termination, "time_limit")
    expect_false(conv$converged)
    expect_true(is.na(conv$attractor))
    expect_gt(conv$residual, steady_residual_tol())

    # The drift is the only thing keeping it from being accepted: allow it and
    # the very same run stops on the distance criterion alone.
    conv <- attr(suppressMessages(
        projectUntilSettled(p, distance_func = blind, residual_tol = Inf,
                            t_max = 3, t_check = 1.5, dt = 0.1,
                            progress_bar = FALSE, info_level = 0)
    ), "convergence")
    expect_identical(conv$termination, "residual_tolerance")
    expect_true(conv$converged)
    # `residual_tol = Inf` says that any drift counts as steady, so the verdict
    # follows the tolerance it was given. That is the user's call to make, and
    # the reported `residual` still says what the drift actually was.
    expect_identical(conv$attractor, "fixed_point")
    expect_gt(conv$residual, steady_residual_tol())
})

test_that("the report says when the drift is what stopped convergence", {
    blind <- function(params, current, previous, ...) 0
    p <- cd_params
    initialN(p)[1, ] <- initialN(p)[1, ] * 10
    expect_message(
        projectUntilSettled(p, distance_func = blind, t_max = 3, t_check = 1.5,
                            dt = 0.1, progress_bar = FALSE, info_level = 3),
        "which is below the distance tolerance, but the biomasses are"
    )
})

test_that("a converged distance function does not mask a limit cycle", {
    skip_on_cran()
    # The same model as the limit-cycle test above, but with a distance function
    # that reports perfect convergence at every block, standing in for a cycle
    # whose period divides `t_check` and is therefore sampled at one phase. The
    # cycle detection works on the finely sampled biomass series, so it sees the
    # oscillation regardless.
    blind <- function(params, current, previous, ...) 0
    sim <- suppressMessages(
        projectUntilSettled(NS_params, effort = 2, distance_func = blind,
                            t_max = 200, t_check = 1.5, dt = 0.1, info_level = 0)
    )
    conv <- attr(sim, "convergence")
    expect_identical(conv$attractor, "limit_cycle")
    expect_gt(conv$period, 0)
    expect_gt(conv$amplitude, 0.1)
})
