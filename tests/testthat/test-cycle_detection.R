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

# The "convergence" attribute ------------------------------------------------

test_that("steady() attaches a 'convergence' attribute for a steady state", {
    p <- suppressWarnings(suppressMessages(steady(cd_params, progress_bar = FALSE)))
    conv <- attr(p, "convergence")
    expect_type(conv, "list")
    expect_named(conv, c("type", "converged", "distance", "years",
                         "period", "amplitude"))
    expect_identical(conv$type, "steady")
    expect_true(conv$converged)
    expect_true(is.na(conv$period))
    expect_true(is.na(conv$amplitude))
})

test_that("the 'convergence' attribute survives return_sim = TRUE", {
    sim <- suppressWarnings(suppressMessages(steady(cd_params, return_sim = TRUE,
                                                    progress_bar = FALSE)))
    expect_s4_class(sim, "MizerSim")
    expect_identical(attr(sim, "convergence")$type, "steady")
})

test_that("projectToSteady() reports non-convergence", {
    p <- cd_params
    initialN(p)[1, ] <- initialN(p)[1, ] * 3
    p <- suppressMessages(
        projectToSteady(p, t_max = 0.5, t_per = 0.5, dt = 0.1,
                        tol = 1e-12, info_level = 0)
    )
    conv <- attr(p, "convergence")
    expect_identical(conv$type, "not_converged")
    expect_false(conv$converged)
})

test_that("fine t_save sampling does not change the result", {
    # Sub-blocking the run at the t_save resolution must be numerically
    # identical to stepping a whole t_per block at once.
    args <- list(t_max = 3, t_per = 1.5, dt = 0.5, tol = 1e-12, info_level = 0)
    p_fine   <- suppressMessages(do.call(projectToSteady,
                                         c(list(cd_params, t_save = 0.5), args)))
    p_coarse <- suppressMessages(do.call(projectToSteady,
                                         c(list(cd_params, t_save = 1.5), args)))
    expect_identical(p_fine@initial_n, p_coarse@initial_n)
})

test_that("t_save is validated", {
    expect_error(projectToSteady(cd_params, dt = 0.1, t_save = 0.15),
                 "t_save must be a positive multiple of dt")
    expect_error(projectToSteady(cd_params, dt = 0.1, t_per = 1, t_save = 0.3),
                 "t_per must be a positive multiple of t_save")
})

test_that("projectToSteady() detects a limit cycle", {
    skip_on_cran()
    # The full North Sea model driven at high fishing effort settles onto a
    # limit cycle rather than a fixed point.
    p <- suppressMessages(
        projectToSteady(NS_params, effort = 2, t_max = 200,
                        t_per = 1.5, dt = 0.1, tol = 1e-8, info_level = 0)
    )
    conv <- attr(p, "convergence")
    expect_identical(conv$type, "cycle")
    expect_true(conv$converged)
    expect_gt(conv$period, 0)
    expect_gt(conv$amplitude, 0.1)
})
