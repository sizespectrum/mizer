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
    target <- getRates(params)$rdd
    p <- tuneSteadyState(params, t_check = 1, t_max = 1, dt = 1, distance_tol = 10) |>
        suppressMessages()
    expect_s3_class(p, "MizerParams")
    expect_equal(getRates(p)$rdd, target, ignore_attr = TRUE)
    expect_snapshot_value(getRDD(p), style = "deparse")
})

test_that("tuneSteadyState accepts consumer update method", {
    params <- tss_params
    params@species_params$gamma[2] <- 2000
    params <- setSearchVolume(params)

    p <- tuneSteadyState(params, t_check = 1, t_max = 1, dt = 1, distance_tol = 10,
                         method = "predictor_corrector") |>
        suppressMessages()

    expect_s3_class(p, "MizerParams")
    expect_true(all(is.finite(p@initial_n)))
    expect_true(all(p@initial_n >= 0))
})

test_that("tuneSteadyState() preserves parameters", {
    params <- tss_params

    params_rdd <- setReproduction(params, RDD = "noRDD")
    p_rdd <- tuneSteadyState(params_rdd, t_check = 1, t_max = 1, dt = 1) |>
        suppressMessages()
    expect_equal(p_rdd@rates_funcs$RDD, "noRDD")

    params_rmax <- params
    species_params(params_rmax)$R_max <- 1.01 * species_params(params_rmax)$R_max
    p_erepro <- tuneSteadyState(params_rmax, t_check = 1, t_max = 1, dt = 1,
                                distance_tol = 10, preserve = "erepro") |>
        suppressMessages()
    expect_equal(p_erepro@species_params$erepro, params_rmax@species_params$erepro)
    p_rl <- tuneSteadyState(params_rmax, t_check = 1, t_max = 1, dt = 1, distance_tol = 10,
                            preserve = "reproduction_level") |> suppressMessages()
    expect_equal(reproduction_level(p_rl), reproduction_level(params_rmax))

    params_erepro <- params
    species_params(params_erepro)$erepro <- 1.01 * species_params(params_erepro)$erepro
    p_rmax <- tuneSteadyState(params_erepro, t_check = 1, t_max = 1, dt = 1,
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
                         "residual", "years", "period", "amplitude", "extinct",
                         "assessment", "scope", "search_assessment"))
    expect_identical(conv$termination, "residual_tolerance")
    expect_true(conv$converged)
    expect_identical(conv$attractor, "fixed_point")
    expect_true(is.na(conv$period))
    expect_true(is.na(conv$amplitude))
    expect_identical(conv$extinct, character(0))
    # The residual measures how far the state actually is from a fixed point,
    # which the distance function only approximates.
    expect_true(is.finite(conv$residual))
    expect_lt(conv$residual, steady_residual_tol())
})

test_that("tuneSteadyState() leaves an initially absent species alone", {
    params <- tss_params
    initialN(params)[1, ] <- 0
    # The species has no egg production to balance, which setBevertonHolt()
    # used to refuse with "Some species have no reproduction.", whichever
    # parameter was being preserved.
    for (preserve in c("reproduction_level", "erepro", "R_max")) {
        expect_no_warning(
            p <- tuneSteadyState(params, preserve = preserve,
                                 progress_bar = FALSE) |>
                suppressMessages())
        expect_identical(attr(p, "convergence")$termination,
                         "residual_tolerance")
        expect_true(all(initialN(p)[1, ] == 0))
        expect_equal(p@species_params$erepro[1],
                     params@species_params$erepro[1])
        expect_equal(p@species_params$R_max[1],
                     params@species_params$R_max[1])
        # The species that is present is still tuned to its steady state.
        expect_equal(getRates(p)$rdd[[2]], getRates(params)$rdd[[2]])
        expect_true(attr(p, "convergence")$assessment$verdict)
    }
})


# findSteadyState(solver = "project") ----------------------------------------

test_that("findSteadyState() returns a MizerParams at the settled state", {
    params <- tss_params
    initialN(params)[1, ] <- initialN(params)[1, ] * 3
    effort <- params@initial_effort * 1.1
    # A single crude block, stopped by an absurdly loose distance tolerance;
    # the drift criterion is switched off rather than pretended to be met.
    expect_message(p <- findSteadyState(params, t_check = 1, dt = 1,
                                        distance_tol = 1000,
                                        residual_tol = Inf, effort = effort),
                   "Reached the convergence tolerance")
    expect_s3_class(p, "MizerParams")
    expect_identical(p@initial_effort, effort)
    # It is exactly projectUntilSettled() with the trajectory thrown away.
    sim <- suppressMessages(projectUntilSettled(params, t_check = 1, dt = 1,
                                                distance_tol = 1000,
                                                residual_tol = Inf,
                                                effort = effort))
    expect_equal(p@initial_n, finalN(sim), ignore_attr = TRUE)
    expect_equal(attr(p, "convergence"), attr(sim, "convergence"))
})

test_that("findSteadyState() accepts the documented effort forms", {
    params <- NS_params_small
    p1 <- findSteadyState(params, effort = 0.5, t_check = 0.1,
                          t_max = 0.1, distance_tol = 10, info_level = 0)
    expect_equal(unname(p1@initial_effort), rep(0.5, length(initial_effort(params))))

    named <- initial_effort(params)[1:2]
    named[] <- c(0.2, NA)
    p2 <- findSteadyState(params, effort = named, t_check = 0.1,
                          t_max = 0.1, distance_tol = 10, info_level = 0)
    expected <- validEffortVector(named, params)
    expect_equal(p2@initial_effort, expected)
})

test_that("findSteadyState() reports non-convergence", {
    p <- tss_params
    initialN(p)[1, ] <- initialN(p)[1, ] * 3
    p <- suppressMessages(
        findSteadyState(p, t_max = 0.5, t_check = 0.5, dt = 0.1,
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
    out <- suppressMessages(findSteadyState(p, t_check = 1, dt = 1, t_max = 1,
                                            distance_tol = 1000,
                                            residual_tol = Inf))
    expect_equal(out@species_params$erepro, p@species_params$erepro)
    expect_equal(out@species_params$R_max, p@species_params$R_max)
    expect_equal(unclass(out@cc_pp), unclass(p@cc_pp))
    # ... whereas tuneSteadyState() adjusts exactly those.
    tuned <- suppressMessages(tuneSteadyState(p, t_check = 1, dt = 1, t_max = 1,
                                              distance_tol = 1000,
                                              residual_tol = Inf))
    # Compare size by size: `all.equal()` takes the mean relative difference,
    # which the capacity at the smallest sizes dominates by many orders of
    # magnitude.
    pos <- p@cc_pp > 0
    expect_gt(max(abs(tuned@cc_pp[pos] / p@cc_pp[pos] - 1)), 0.01)
})


# The Newton solver ----------------------------------------------------------

test_that("tuneSteadyState(solver = 'newton') returns a fixed point", {
    skip_unless_experimental()
    pn <- tuneSteadyState(p_steady, solver = "newton")
    expect_s3_class(pn, "MizerParams")

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
    expect_s3_class(pn, "MizerParams")
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

    feasible <- setBevertonHolt(p_steady, reproduction_level = 0.25)
    pn_erepro <- tuneSteadyState(feasible, solver = "newton",
                                 preserve = "erepro")
    expect_equal(pn_erepro@species_params$erepro,
                 feasible@species_params$erepro, tolerance = 1e-5)
})

test_that("tuneSteadyState(solver = 'newton') leaves an absent species alone", {
    skip_unless_experimental()
    params <- tss_params
    initialN(params)[1, ] <- 0
    expect_warning(
        pn <- tuneSteadyState(params, solver = "newton") |> suppressMessages(),
        "already absent and were held at zero")
    expect_identical(attr(pn, "convergence")$termination, "solver_converged")
    expect_true(all(initialN(pn)[1, ] == 0))
    expect_equal(pn@species_params$erepro[1], params@species_params$erepro[1])
    expect_equal(pn@species_params$R_max[1], params@species_params$R_max[1])
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
        expect_named(conv[names(conv) != "search_assessment"],
                     c("termination", "converged", "attractor", "distance",
                       "residual", "years", "period", "amplitude", "extinct",
                       "assessment", "scope"))
        expect_identical(conv$termination, "solver_converged")
        expect_true(conv$converged)
        expect_identical(conv$attractor, "fixed_point")
        expect_true(is.finite(conv$residual))
        expect_true(is.na(conv$distance))
        expect_true(is.na(conv$years))
        expect_identical(conv$extinct, character(0))
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
    expect_identical(conv$extinct, "Cod")
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
    expect_s3_class(pn, "MizerParams")
    expect_true(all(is.finite(pn@initial_n)))
})

test_that("the Newton solver works with the second-order (van Leer) scheme", {
    withr::local_options(mizer_info_level = 0)
    skip_unless_experimental()
    p <- NS_params_small
    sow <- second_order_w(p)
    sow$flux <- "van_leer"
    sow$bin_average <- TRUE
    second_order_w(p) <- sow
    ps <- tuneSteadyState(p, t_max = 100, progress_bar = FALSE)

    pn <- tuneSteadyState(ps, solver = "newton")
    expect_s3_class(pn, "MizerParams")
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
    expect_s3_class(pn, "MizerParams")

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

test_that("convergence support alone does not authorise tuning", {
    assign("tuning_growing_component", function(params, n_other, component,
                                                dt, ...) {
        n_other[[component]] * (1 + dt)
    }, envir = .GlobalEnv)
    assign("tuning_component_assessment", function(...) {
        data.frame(measure = "growth", value = 1, tolerance = 0.05,
                   units = "1/year")
    }, envir = .GlobalEnv)
    withr::defer(rm("tuning_growing_component", "tuning_component_assessment",
                    envir = .GlobalEnv))
    p <- setComponent(NS_params_small, "grower", initial_value = 1,
                      dynamics_fun = "tuning_growing_component",
                      convergence_fun = "tuning_component_assessment")
    expect_error(tuneSteadyState(p, t_check = 0.1, t_max = 0.1, dt = 0.1,
                                 distance_tol = 10, progress_bar = FALSE,
                                 info_level = 0),
                 "stationary tuning support for component:grower")
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

test_that("tuning preflights unsupported custom resource convergence", {
    assign("mizer_test_resource_dyn", function(params, n_pp, ...) n_pp,
           envir = .GlobalEnv)
    withr::defer(rm("mizer_test_resource_dyn", envir = .GlobalEnv))

    p <- NS_params_small
    p@resource_dynamics <- "mizer_test_resource_dyn"

    expect_error(
        tuneSteadyState(p, t_check = 1, t_max = 1, dt = 1,
                        distance_tol = 10, progress_bar = FALSE),
        "Missing convergence support for resource"
    )
})

test_that("tuneSteadyState() does not complain about a constant resource", {
    # `resource_constant()` hands back the abundance it was given, so the
    # preserved abundance is steady whatever the capacity is and there is
    # nothing to rebalance and nothing to report.
    p <- NS_params_small
    p <- setResource(p, resource_dynamics = "resource_constant")
    expect_no_warning(
        suppressMessages(tuneSteadyState(p, t_check = 1, t_max = 1, dt = 1,
                                         distance_tol = 10,
                                         progress_bar = FALSE))
    )
})

test_that("projection calibration forwards Courant subcycling controls", {
    # Fail after one step, before calibration can hide an ignored control.
    expect_error(tuneSteadyState(NS_params_small, solver = "project",
                                 dt = 0.1, t_check = 0.1, t_max = 0.1,
                                 subcycle = TRUE, max_courant = 0.001,
                                 max_substeps = 1, progress_bar = FALSE,
                                 info_level = 0), "max_substeps")
})

test_that("tuning preflights final convergence support", {
    p <- setComponent(NS_params_small, "unmeasured", initial_value = 1,
                      dynamics_fun = "constant_other")
    for (solver in c("project", "newton")) {
        expect_error(tuneSteadyState(p, solver = solver, info_level = 0),
                     "Missing convergence support for component:unmeasured")
    }
})

test_that("Newton requires explicit exclusions for dynamic components", {
    assign("newton_component_check", function(...) {
        data.frame(measure = "change", value = 0, tolerance = 0,
                   units = "1/year")
    }, envir = .GlobalEnv)
    withr::defer(rm("newton_component_check", envir = .GlobalEnv))
    p <- setComponent(NS_params_small, "input", initial_value = 1,
                      dynamics_fun = "constant_other",
                      convergence_fun = "newton_component_check")
    for (finder in list(findSteadyState)) {
        expect_error(finder(p, solver = "newton", info_level = 0),
                     "Newton solver cannot solve dynamic components")
        expect_error(finder(p, solver = "newton", scope = "conditional",
                            fixed = "consumers", info_level = 0),
                     "cannot hold consumers fixed")
    }
    skip_unless_experimental()
    p <- setComponent(p_steady, "input", initial_value = 1,
                      dynamics_fun = "constant_other",
                      convergence_fun = "newton_component_check")
    out <- findSteadyState(p, solver = "newton", scope = "conditional",
                           fixed = "component:input", info_level = 0)
    conv <- attr(out, "convergence")
    expect_identical(conv$scope, "conditional")
    expect_identical(conv$assessment$fixed, "component:input")
    expect_equal(out$initial_n_other$input, 1)
})

# Declared tuning strategies and constrained effective rates ----

test_that("tuning balances held components and evolves selected alternatives", {
    assign("strategy_pool_rate", function(params, n_other, component, ...) {
        params$other_params$supply - n_other[[component]]
    }, envir = .GlobalEnv)
    assign("strategy_pool_check", function(params, derivative, ...) {
        data.frame(measure = "change", value = derivative, tolerance = 1e-8,
                   units = "pool/year")
    }, envir = .GlobalEnv)
    assign("strategy_pool_equations", function(params, n_other, ...) {
        list(input = params$other_params$supply, loss = n_other$pool)
    }, envir = .GlobalEnv)
    assign("strategy_pool_targets", function(params, n_other, ...) {
        list(input = n_other$pool, loss = n_other$pool)
    }, envir = .GlobalEnv)
    assign("strategy_pool_balance", function(params, n_other, ...) {
        params$other_params$supply <- n_other$pool
        params
    }, envir = .GlobalEnv)
    withr::defer(rm("strategy_pool_rate", "strategy_pool_check",
                    "strategy_pool_equations", "strategy_pool_targets",
                    "strategy_pool_balance", envir = .GlobalEnv))
    p <- p_steady
    p$other_params$supply <- 0
    hold <- list(action = "balance",
                 hold = list(states = "component:pool", rates = character()),
                 adjust = list(supply = list(path = c("other_params", "supply"),
                                             lower = 0)),
                 target_fun = "strategy_pool_targets",
                 balance_fun = "strategy_pool_balance")
    evolve <- list(action = "none", evolve = "component:pool",
                   solvers = "project")
    p <- setComponent(
        p, "pool", initial_value = 2, rate_fun = "strategy_pool_rate",
        convergence_fun = "strategy_pool_check", stationary = TRUE,
        equation_rates_fun = "strategy_pool_equations",
        tuning = list(default = "hold", strategies = list(hold = hold,
                                                          evolve = evolve))
    )
    solvers <- if (testing_experimental()) c("project", "newton") else "project"
    for (solver in solvers) {
        out <- tuneSteadyState(p, solver = solver, info_level = 0)
        expect_equal(out$initial_n_other$pool, 2)
        expect_equal(out$other_params$supply, 2)
        conv <- attr(out, "convergence")
        expect_identical(conv$search_assessment$scope, "conditional")
        expect_true("component:pool" %in% conv$search_assessment$fixed)
        expect_identical(conv$search_assessment$rate_constraints,
                         "consumers/recruitment")
        measures <- conv$assessment$measures
        expect_true(all(measures$passed[measures$owner == "component:pool"]))
        expect_identical(out$other_dynamics, p$other_dynamics)
    }
    evolved <- tuneSteadyState(p, strategies = c("component:pool" = "evolve"),
                               t_max = 0.1, t_check = 0.1, dt = 0.1,
                               residual_tol = Inf, distance_tol = 1e10,
                               info_level = 0)
    expect_lt(evolved$initial_n_other$pool, 2)
    expect_equal(evolved$other_params$supply, 0)
    expect_false("component:pool" %in%
                     attr(evolved, "convergence")$search_assessment$fixed)
    expect_error(tuneSteadyState(p, solver = "newton",
                                 strategies = c("component:pool" = "evolve"),
                                 info_level = 0), "[Nn]ewton|solver")
    expect_error(tuneSteadyState(p, scope = "conditional",
                                 fixed = "component:pool",
                                 strategies = c("component:pool" = "evolve"),
                                 info_level = 0), "compatible strategy")
})

test_that("recruitment-only search constraints are reported as conditional", {
    p <- p_steady
    p$capabilities$processes$resource$tuning$strategies$evolve <-
        list(action = "none", evolve = "resource")
    out <- tuneSteadyState(p, strategies = c(resource = "evolve"),
                           t_max = 0.1, t_check = 0.1, dt = 0.1,
                           residual_tol = Inf, distance_tol = 1e10,
                           info_level = 0)
    conv <- attr(out, "convergence")
    expect_identical(conv$search_assessment$scope, "conditional")
    expect_identical(conv$search_assessment$fixed, character())
    expect_identical(conv$search_assessment$rate_constraints,
                     "consumers/recruitment")
    expect_identical(conv$scope, "whole")
    expect_equal(out$cc_pp, p$cc_pp, ignore_attr = TRUE)
})

test_that("tuning freezes recruitment after an immigration Rates wrapper", {
    assign("tuning_immigration_rates", function(params, ...) {
        rates <- mizerRates(params, ...)
        rates$rdd <- rates$rdd + params$other_params$immigration
        rates
    }, envir = .GlobalEnv)
    assign("tuning_immigration_balance", function(params, rates, rdd_target,
                                                  preserved, ...) {
        balance_BevertonHoltRDD(
            params, rates = rates,
            rdd_target = rdd_target - params$other_params$immigration,
            preserve = "erepro", preserved = preserved
        )
    }, envir = .GlobalEnv)
    withr::defer(rm("tuning_immigration_rates", "tuning_immigration_balance",
                    envir = .GlobalEnv))
    p <- p_steady
    declarations <- p$capabilities$processes
    p$other_params$immigration <- 0.2 * getRates(p)$rdd
    p$rates_funcs$Rates <- "tuning_immigration_rates"
    # This provider explicitly declares that immigration only changes the
    # final recruitment, and supplies the matching inverse target transform.
    binding <- tuning_provider_binding(p)
    for (owner in names(declarations)) {
        declarations[[owner]]$tuning_binding <- binding
    }
    strategy <- declarations$reproduction$tuning$strategies$preserve_erepro
    strategy$balance_fun <- "tuning_immigration_balance"
    declarations$reproduction$tuning <- list(
        default = "immigration", strategies = list(immigration = strategy)
    )
    p$capabilities$processes <- declarations
    target <- getRates(p)$rdd
    reference <- p
    reference$rates_funcs$Rates <- "mizerRates"
    reference$rates_funcs$RDD <- "constantRDD"
    reference$species_params$constant_reproduction <- target
    reference$resource_dynamics <- "resource_constant"
    expected <- project_simple(reference, dt = 0.1, steps = 1, info_level = 0)
    out <- tuneSteadyState(p, t_max = 0.1, t_check = 0.1, dt = 0.1,
                           residual_tol = Inf, distance_tol = 1e10,
                           info_level = 0)
    expect_equal(out$initial_n, expected$n, ignore_attr = TRUE)
    expect_equal(getRates(out)$rdd, target, ignore_attr = TRUE)
    expect_identical(out$rates_funcs, p$rates_funcs)
    search <- attr(out, "convergence")$search_assessment
    expect_identical(search$rate_constraints, "consumers/recruitment")
    skip_unless_experimental()
    expected <- newton_steady_state(reference, effort = p$initial_effort,
                                    rdd_const = target, resource = "fixed")
    out <- tuneSteadyState(p, solver = "newton", info_level = 0)
    expect_equal(out$initial_n, expected$n, ignore_attr = TRUE)
    expect_equal(getRates(out)$rdd, target, ignore_attr = TRUE)
    expect_identical(attr(out, "convergence")$search_assessment$scope,
                     "conditional")
})

test_that("Failed tuning assessments cannot mutate caller reference data", {
    rlang::local_bindings(
        failing_reference_assessment = function(params, ...) {
            params$other_params$audit$value <- 2
            stop("deliberate assessment failure")
        }, .env = globalenv()
    )
    p <- NS_params_small
    p$other_params$audit <- new.env(parent = emptyenv())
    p$other_params$audit$value <- 1
    p <- setProcess(p, "consumers",
                    convergence_fun = "failing_reference_assessment")
    expect_error(tuneSteadyState(p, t_max = 0.1, dt = 0.1, t_check = 0.1,
                                 progress_bar = FALSE, info_level = 0),
                 "deliberate assessment failure")
    expect_identical(p$other_params$audit$value, 1)
})
