# Tests for the two exported steady-state finders, tuneSteadyState() and
# findSteadyState(), across both of their solvers. The projection engine they
# share is tested in test-steady.R, the Newton machinery in test-steadyNewton.R,
# and the superseded steady()/projectToSteady() wrappers in
# test-backwards_compatibility.R.

# Built once here; tests get a fresh copy via R's copy-on-modify semantics.
tss_params <- suppressMessages(
    newTraitParams(no_sp = 2, no_w = 20, max_w_max = 100,
                   min_w = 1e-3, w_pp_cutoff = 5, ks = 4,
                   reproduction_level = 0.25, info_level = 0)
)

# A steadied model is a good starting guess for the root finder.
delayedAssign("p_steady",
              tuneSteadyState(NS_params_small, t_max = 50,
                              progress_bar = FALSE))

# Deviation of the resource from the semichemostat fixed point
# n_steady = r * cc / (r + mu); zero when the resource is exactly at steady
# state.
resource_deviation <- function(params) {
    mu <- getResourceMort(params)
    n_pp <- params@initial_n_pp
    n_steady <- params@rr_pp * params@cc_pp / (params@rr_pp + mu)
    sel <- n_pp > 0 & is.finite(n_steady)
    max(abs(n_pp[sel] - n_steady[sel]) / n_pp[sel])
}


# tuneSteadyState(solver = "project") ----------------------------------------

test_that("tuneSteadyState works", {
    params <- tss_params
    params@species_params$gamma[2] <- 2000
    params <- setSearchVolume(params)
    p <- tuneSteadyState(params, t_per = 1, t_max = 1, dt = 1, distance_tol = 10) |>
        suppressMessages()
    expect_s4_class(p, "MizerParams")
    expect_snapshot_value(getRDD(p), style = "deparse")
})

test_that("tuneSteadyState accepts consumer update method", {
    params <- tss_params
    params@species_params$gamma[2] <- 2000
    params <- setSearchVolume(params)

    p <- tuneSteadyState(params, t_per = 1, t_max = 1, dt = 1, distance_tol = 10,
                         method = "predictor_corrector") |>
        suppressMessages()

    expect_s4_class(p, "MizerParams")
    expect_true(all(is.finite(p@initial_n)))
    expect_true(all(p@initial_n >= 0))
})

test_that("tuneSteadyState() preserves parameters", {
    params <- tss_params

    params_rdd <- params
    params_rdd@rates_funcs$RDD <- "noRDD"
    p_rdd <- tuneSteadyState(params_rdd, t_per = 1, t_max = 1, dt = 1) |>
        suppressMessages()
    expect_equal(p_rdd@rates_funcs$RDD, "noRDD")

    params_rmax <- params
    species_params(params_rmax)$R_max <- 1.01 * species_params(params_rmax)$R_max
    p_erepro <- tuneSteadyState(params_rmax, t_per = 1, t_max = 1, dt = 1,
                                distance_tol = 10, preserve = "erepro") |>
        suppressMessages()
    expect_equal(p_erepro@species_params$erepro, params_rmax@species_params$erepro)
    p_rl <- tuneSteadyState(params_rmax, t_per = 1, t_max = 1, dt = 1, distance_tol = 10,
                            preserve = "reproduction_level") |> suppressMessages()
    expect_equal(reproduction_level(p_rl), reproduction_level(params_rmax))

    params_erepro <- params
    species_params(params_erepro)$erepro <- 1.01 * species_params(params_erepro)$erepro
    p_rmax <- tuneSteadyState(params_erepro, t_per = 1, t_max = 1, dt = 1,
                              distance_tol = 10, preserve = "R_max") |> suppressMessages()
    expect_equal(p_rmax@species_params$R_max, params_erepro@species_params$R_max)
    expect_false(identical(p_rmax@time_modified, params_erepro@time_modified))
})

test_that("tuneSteadyState() rebalances the resource to a genuine fixed point", {
    # For a normal (non-frozen) resource, tuneSteadyState() rebalances the
    # resource so that it is a genuine fixed point and issues no warning.
    expect_no_warning(
        ps <- suppressMessages(tuneSteadyState(NS_params_small,
                                               progress_bar = FALSE))
    )
    expect_lt(resource_deviation(ps), 1e-8)
    # From this state both the second-order and Euler steppers stay put and
    # agree closely, confirming the resource really is at steady state.
    sim_t <- project(ps, t_max = 5, dt = 0.1, method = "tr_bdf2",
                     progress_bar = FALSE)
    sim_e <- project(ps, t_max = 5, dt = 0.1, method = "euler",
                     progress_bar = FALSE)
    rel_diff_pp <- max(abs(finalNResource(sim_t) - finalNResource(sim_e)) /
                           (finalNResource(sim_e) + 1e-30))
    expect_lt(rel_diff_pp, 1e-3)
})

test_that("tuneSteadyState() rebalances a frozen resource capacity", {
    # Freeze the resource capacity away from its balanced value.
    pf <- suppressMessages(setResource(NS_params_small,
                                       resource_capacity = NS_params_small@cc_pp * 1.3,
                                       balance = FALSE))
    expect_false(is.null(comment(pf@cc_pp)))
    npp_before <- pf@initial_n_pp

    # tuneSteadyState() rebalances the frozen capacity silently, as earlier
    # mizer versions did, without a warning or message about the resource.
    expect_no_warning(ps <- suppressMessages(tuneSteadyState(pf,
                                                             progress_bar = FALSE)))
    # The resource is now a genuine fixed point and no longer frozen ...
    expect_lt(resource_deviation(ps), 1e-8)
    expect_null(comment(ps@cc_pp))
    # ... and the resource abundance itself was preserved (only the capacity
    # was adjusted to balance it).
    expect_equal(unclass(ps@initial_n_pp), unclass(npp_before))
})

test_that("tuneSteadyState() attaches a 'convergence' attribute", {
    p <- suppressWarnings(suppressMessages(
        tuneSteadyState(tss_params, progress_bar = FALSE)))
    conv <- attr(p, "convergence")
    expect_type(conv, "list")
    expect_named(conv, c("termination", "converged", "attractor", "distance",
                         "residual", "years", "period", "amplitude"))
    expect_identical(conv$termination, "residual_tolerance")
    expect_true(conv$converged)
    expect_identical(conv$attractor, "fixed_point")
    expect_true(is.na(conv$period))
    expect_true(is.na(conv$amplitude))
    # The residual measures how far the state actually is from a fixed point,
    # which the distance function only approximates.
    expect_true(is.finite(conv$residual))
    expect_lt(conv$residual, steady_residual_tol())
})


# findSteadyState(solver = "project") ----------------------------------------

test_that("findSteadyState() returns a MizerParams at the settled state", {
    params <- tss_params
    initialN(params)[1, ] <- initialN(params)[1, ] * 3
    effort <- params@initial_effort * 1.1
    # A single crude block, stopped by an absurdly loose distance tolerance;
    # the drift criterion is switched off rather than pretended to be met.
    expect_message(p <- findSteadyState(params, t_per = 1, dt = 1,
                                        distance_tol = 1000,
                                        residual_tol = Inf, effort = effort),
                   "Reached the convergence tolerance")
    expect_s4_class(p, "MizerParams")
    expect_identical(p@initial_effort, effort)
    # It is exactly projectUntilSettled() with the trajectory thrown away.
    sim <- suppressMessages(projectUntilSettled(params, t_per = 1, dt = 1,
                                                distance_tol = 1000,
                                                residual_tol = Inf,
                                                effort = effort))
    expect_equal(p@initial_n, finalN(sim), ignore_attr = TRUE)
    expect_equal(attr(p, "convergence"), attr(sim, "convergence"))
})

test_that("findSteadyState() accepts the documented effort forms", {
    params <- NS_params_small
    p1 <- findSteadyState(params, effort = 0.5, t_per = 0.1,
                          t_max = 0.1, distance_tol = 10, info_level = 0)
    expect_equal(unname(p1@initial_effort), rep(0.5, length(initial_effort(params))))

    named <- initial_effort(params)[1:2]
    named[] <- c(0.2, NA)
    p2 <- findSteadyState(params, effort = named, t_per = 0.1,
                          t_max = 0.1, distance_tol = 10, info_level = 0)
    expected <- validEffortVector(named, params)
    expect_equal(p2@initial_effort, expected)
})

test_that("findSteadyState() reports non-convergence", {
    p <- tss_params
    initialN(p)[1, ] <- initialN(p)[1, ] * 3
    p <- suppressMessages(
        findSteadyState(p, t_max = 0.5, t_per = 0.5, dt = 0.1,
                        distance_tol = 1e-12, info_level = 0)
    )
    conv <- attr(p, "convergence")
    expect_identical(conv$termination, "time_limit")
    expect_false(conv$converged)
    # Half a year from a state knocked well off its steady state: not a fixed
    # point by any measure, so no attractor is claimed.
    expect_true(is.na(conv$attractor))
})

test_that("findSteadyState() changes no parameter", {
    p <- tss_params
    initialN(p)[1, ] <- initialN(p)[1, ] * 3
    # One block and no drift criterion: what is under test is which parameters
    # the two functions touch, not where either of them ends up.
    out <- suppressMessages(findSteadyState(p, t_per = 1, dt = 1, t_max = 1,
                                            distance_tol = 1000,
                                            residual_tol = Inf))
    expect_equal(out@species_params$erepro, p@species_params$erepro)
    expect_equal(out@species_params$R_max, p@species_params$R_max)
    expect_equal(unclass(out@cc_pp), unclass(p@cc_pp))
    # ... whereas tuneSteadyState() adjusts exactly those.
    tuned <- suppressMessages(tuneSteadyState(p, t_per = 1, dt = 1, t_max = 1,
                                              distance_tol = 1000,
                                              residual_tol = Inf))
    expect_false(isTRUE(all.equal(unclass(tuned@cc_pp), unclass(p@cc_pp))))
})


# The Newton solver ----------------------------------------------------------

test_that("tuneSteadyState(solver = 'newton') returns a fixed point", {
    skip_unless_experimental()
    pn <- tuneSteadyState(p_steady, solver = "newton")
    expect_s4_class(pn, "MizerParams")

    # The returned spectra should barely move when projected forward.
    sim <- project(pn, t_max = 1, dt = 0.25, t_save = 1)
    n0 <- pn@initial_n
    n1 <- finalN(sim)
    support <- n0 > 0
    drift <- max((abs(n1 - n0) / n0)[support])
    expect_lt(drift, 1e-3)
})

test_that("findSteadyState(solver = 'newton') returns a fixed point", {
    skip_unless_experimental()
    pn <- findSteadyState(p_steady, solver = "newton")
    expect_s4_class(pn, "MizerParams")
    sim <- project(pn, t_max = 1, dt = 0.25, t_save = 1)
    n0 <- pn@initial_n
    support <- n0 > 0
    drift <- max((abs(finalN(sim) - n0) / n0)[support])
    expect_lt(drift, 1e-3)
})

test_that("the two solvers agree on a stable model", {
    skip_unless_experimental()
    pn <- tuneSteadyState(p_steady, solver = "newton")
    # Compare where the projected spectra are well above zero. The Newton
    # solver finds a tighter fixed point than the time-stepping one reaches, so
    # a few percent difference on the well-resolved bins is expected (and is
    # the point of the exercise).
    big <- p_steady@initial_n > max(p_steady@initial_n) * 1e-6
    rel <- abs(pn@initial_n - p_steady@initial_n) / p_steady@initial_n
    expect_lt(max(rel[big]), 0.1)
})

test_that("tuneSteadyState(solver = 'newton') holds the resource and honours preserve", {
    skip_unless_experimental()
    npp_before <- p_steady@initial_n_pp

    pn_level <- tuneSteadyState(p_steady, solver = "newton",
                                preserve = "reproduction_level")
    expect_equal(reproduction_level(pn_level),
                 reproduction_level(p_steady), tolerance = 1e-5)
    # The resource abundance is what is held; the capacity is what moves.
    expect_equal(unclass(pn_level@initial_n_pp), unclass(npp_before))
    expect_lt(resource_deviation(pn_level), 1e-8)

    pn_rmax <- tuneSteadyState(p_steady, solver = "newton", preserve = "R_max")
    expect_equal(pn_rmax@species_params$R_max,
                 p_steady@species_params$R_max, tolerance = 1e-5)

    pn_erepro <- suppressWarnings(tuneSteadyState(p_steady, solver = "newton",
                                                  preserve = "erepro"))
    expect_equal(pn_erepro@species_params$erepro,
                 p_steady@species_params$erepro, tolerance = 1e-5)
})

test_that("findSteadyState(solver = 'newton') leaves the reproduction parameters alone", {
    skip_unless_experimental()
    pn <- findSteadyState(p_steady, solver = "newton")
    expect_equal(pn@species_params$R_max,
                 p_steady@species_params$R_max, tolerance = 1e-5)
    expect_equal(pn@species_params$erepro,
                 p_steady@species_params$erepro, tolerance = 1e-5)
    expect_equal(unclass(pn@cc_pp), unclass(p_steady@cc_pp))
})

test_that("the Newton solvers attach a well-formed 'convergence' attribute", {
    skip_unless_experimental()
    for (p in list(tuneSteadyState(p_steady, solver = "newton"),
                   findSteadyState(p_steady, solver = "newton"))) {
        conv <- attr(p, "convergence")
        expect_named(conv, c("termination", "converged", "attractor",
                             "distance", "residual", "years", "period",
                             "amplitude"))
        expect_identical(conv$termination, "solver_converged")
        expect_true(conv$converged)
        expect_identical(conv$attractor, "fixed_point")
        expect_true(is.finite(conv$residual))
        expect_true(is.na(conv$distance))
        expect_true(is.na(conv$years))
    }
})

test_that("the `jacobian` argument selects the nleqslv method", {
    skip_unless_experimental()
    pu <- tuneSteadyState(p_steady, solver = "newton", jacobian = "update")
    pr <- tuneSteadyState(p_steady, solver = "newton", jacobian = "recompute")
    # Two routes to the same root.
    expect_equal(unclass(pu@initial_n), unclass(pr@initial_n),
                 tolerance = 1e-4)
    expect_error(tuneSteadyState(p_steady, solver = "newton",
                                 jacobian = "bogus"),
                 "should be one of")
})

test_that("findSteadyState(solver = 'newton') handles extinctions", {
    skip_unless_experimental()
    # Make species 3 (Cod) unviable by setting its reproduction efficiency
    # extremely low.
    p_extinct <- p_steady
    p_extinct@species_params$erepro[3] <- 1e-12
    p_extinct@initial_n[3, ] <- p_steady@initial_n[3, ]

    expect_warning(pn_ext <- findSteadyState(p_extinct, solver = "newton",
                                             extinction_floor = 1e-6),
                   "went extinct and were set to zero")
    lo <- p_extinct@w_min_idx[3]
    expect_equal(pn_ext@initial_n[3, lo], 0)
    conv <- attr(pn_ext, "convergence")
    expect_identical(conv$termination, "extinction")
    # `attractor` is not forced to NA here: the state left behind, with that
    # species at zero, may well be a perfectly good fixed point of what remains,
    # and it says so only if the measured drift agrees. `termination` is what
    # says that a species was lost on the way.
    expect_true(is.na(conv$attractor) ||
                identical(conv$attractor, "fixed_point"))
    expect_true(is.finite(conv$residual))
})

test_that("only the resource-solving Newton branch needs a semichemostat", {
    skip_unless_experimental()
    p_log <- setResource(NS_params_small,
                         resource_dynamics = "resource_logistic")
    # findSteadyState() makes the resource an unknown, so it needs the
    # semichemostat equation ...
    expect_error(findSteadyState(p_log, solver = "newton"), "semichemostat")
    # ... whereas tuneSteadyState() holds the resource fixed and does not.
    pn <- suppressWarnings(tuneSteadyState(p_log, solver = "newton"))
    expect_s4_class(pn, "MizerParams")
    expect_true(all(is.finite(pn@initial_n)))
})

test_that("the Newton solver works with the second-order (van Leer) scheme", {
    skip_unless_experimental()
    p <- NS_params_small
    sow <- second_order_w(p)
    sow$flux <- "van_leer"
    sow$bin_average <- TRUE
    second_order_w(p) <- sow
    ps <- tuneSteadyState(p, t_max = 100, progress_bar = FALSE)

    pn <- tuneSteadyState(ps, solver = "newton")
    expect_s4_class(pn, "MizerParams")
    # The returned state must still be a fixed point of the (second-order)
    # dynamics. The van Leer limiter is only Lipschitz, so the Newton residual
    # cannot be driven to machine precision, but the projected drift is the
    # honest test and stays small.
    sim <- project(pn, t_max = 1, dt = 0.25, t_save = 1)
    n0 <- pn@initial_n
    n1 <- finalN(sim)
    support <- n0 > 0
    drift <- max((abs(n1 - n0) / n0)[support])
    expect_lt(drift, 1e-3)
})

test_that("the Newton solver handles initial guesses that are non-zero above the support", {
    skip_unless_experimental()
    p <- p_steady
    # Fill the trailing tail (which is zero in p_steady) with positive numbers
    # above its actual support. In p_steady, species 1 (Sprat) only grows to
    # its w_max (0.33g); we set non-zero values up to the end of the grid.
    grid_top <- mizer:::support_top_idx(p)
    no_w <- length(p@w)
    idx_zeros <- (grid_top[1] + 1):no_w
    if (length(idx_zeros) > 0) {
        p@initial_n[1, idx_zeros] <- 1e-3
    }

    pn <- tuneSteadyState(p, solver = "newton")
    expect_s4_class(pn, "MizerParams")

    # The tail should be correctly zeroed out in the result
    if (length(idx_zeros) > 0) {
        expect_true(all(pn@initial_n[1, idx_zeros] == 0))
    }
})

test_that("the Newton solver reports the residual it achieved and respects info_level", {
    skip_unless_experimental()
    expect_message(tuneSteadyState(p_steady, solver = "newton"),
                   "change at up to")
    expect_silent(tuneSteadyState(p_steady, solver = "newton", info_level = 0))
})

test_that("verbose = TRUE traces the Newton iterations", {
    skip_unless_experimental()
    out <- capture.output(tuneSteadyState(p_steady, solver = "newton",
                                          verbose = TRUE))
    expect_true(any(grepl("Iteration report", out)))
})

test_that("an unknown solver is rejected", {
    expect_error(tuneSteadyState(tss_params, solver = "bogus"),
                 "should be one of")
    expect_error(findSteadyState(tss_params, solver = "bogus"),
                 "should be one of")
})

# What the analysis actually covers ------------------------------------------

test_that("tuneSteadyState() says when it has held a component fixed", {
    assign("mizer_test_growing_component", function(params, n_other, component,
                                                    dt, ...) {
        n_other[[component]] * (1 + dt)
    }, envir = .GlobalEnv)
    withr::defer(rm("mizer_test_growing_component", envir = .GlobalEnv))

    p <- setComponent(NS_params_small, "grower", initial_value = 1,
                      dynamics_fun = "mizer_test_growing_component")

    # The search holds it fixed and never rebalances it, so the model that comes
    # back is at a fixed point of the fish and the resource only. That has to be
    # said, not assumed — but it is not a refusal: extension models are built
    # this way and still need to be tuned.
    expect_warning(
        p_tuned <- suppressMessages(
            tuneSteadyState(p, t_per = 1, t_max = 1, dt = 1,
                            distance_tol = 10, progress_bar = FALSE)),
        "has dynamics of their own|has\n?\\s*dynamics"
    )
    expect_s4_class(p_tuned, "MizerParams")

    # And the residual it reports does cover the component, so the number the
    # warning points at is one that would notice.
    conv <- attr(p_tuned, "convergence")
    expect_true(is.finite(conv$residual))
    # The component grows at 1/year, and the residual is the largest drift of
    # anything in the model, so it cannot be below that.
    expect_gt(conv$residual, 0.9)
})

test_that("a constant component is nothing to report", {
    p <- setComponent(NS_params_small, "constant", initial_value = 1,
                      dynamics_fun = "constant_other")
    expect_false(mizer:::warn_other_components_fixed(p, "context"))
    expect_true(mizer:::warn_other_components_fixed(
        setComponent(NS_params_small, "moving", initial_value = 1,
                     dynamics_fun = "resource_constant"), "context") |>
        suppressWarnings())
})

test_that("tuneSteadyState() says when the resource could not be rebalanced", {
    assign("mizer_test_resource_dyn", function(params, n_pp, ...) n_pp,
           envir = .GlobalEnv)
    withr::defer(rm("mizer_test_resource_dyn", envir = .GlobalEnv))

    p <- NS_params_small
    p@resource_dynamics <- "mizer_test_resource_dyn"

    # There is no `balance_mizer_test_resource_dyn()`, so setResource() cannot
    # derive a capacity that makes the preserved abundance steady. It used to
    # skip the balancing without a word.
    expect_warning(
        suppressMessages(tuneSteadyState(p, t_per = 1, t_max = 1, dt = 1,
                                         distance_tol = 10,
                                         progress_bar = FALSE)),
        "could not be rebalanced"
    )
})
