# projectToSteady ----
# Built once here; tests get a fresh copy via R's copy-on-modify semantics.
small_steady_params <- suppressMessages(
    newTraitParams(no_sp = 2, no_w = 20, max_w_max = 100,
                   min_w = 1e-3, w_pp_cutoff = 5, ks = 4,
                   reproduction_level = 0.25, info_level = 0)
)

test_that("projectToSteady() works", {
    params <- small_steady_params
    effort <- params@initial_effort * 1.1
    initialN(params)[1, ] <- initialN(params)[1, ] * 3
    expect_error(projectToSteady(params, dt = 1, t_per = 0.5),
                 "t_per must be a positive multiple of dt")
    expect_error(projectToSteady(params, t_max = 0.1),
                 "t_max not greater than or equal to t_per")
    expect_message(projectToSteady(params, t_max = 0.1, t_per = 0.1, info_level = 3),
                   "Simulation run did not converge after 0.1 years.")
    expect_message(paramsc <- projectToSteady(params, t_per = 1, dt = 1,
                                              tol = 1000, effort = effort),
                   "Convergence was achieved")
    expect_s4_class(paramsc, "MizerParams")
    expect_identical(paramsc@initial_effort, effort)
    # shouldn't take long the second time we run to steady
    expect_message(projectToSteady(paramsc, t_per = 1, dt = 1, tol = 1000),
                   "Convergence was achieved")
    
    # return sim
    expect_message(sim <- projectToSteady(params, return_sim = TRUE,
                                          t_per = 1, dt = 1, tol = 1000),
                   "Convergence was achieved")
    expect_s4_class(sim, "MizerSim")
    
    # Alternative distance function
    expect_message(paramsc <- projectToSteady(params, t_per = 1, dt = 1,
                                              distance_func = distanceMaxRelRDI,
                                              tol = 10),
                   "Convergence was achieved")
    # shouldn't take long the second time we run to steady
    expect_message(projectToSteady(paramsc, t_per = 1, dt = 1,
                                   distance_func = distanceMaxRelRDI,
                                   tol = 10),
                   "Convergence was achieved")
    
    # Check extinction
    params@psi[1:2, ] <- 0
    sp1 <- params@species_params$species[1]
    sp2 <- params@species_params$species[2]
    expect_warning(projectToSteady(params) |> suppressMessages(),
                   paste0(sp1, ", ", sp2, " are going extinct."))
})

test_that("projectToSteady accepts the documented effort forms", {
    params <- NS_params_small
    p1 <- projectToSteady(params, effort = 0.5, t_per = 0.1,
                          t_max = 0.1, tol = 10, info_level = 0)
    expect_equal(unname(p1@initial_effort), rep(0.5, length(initial_effort(params))))

    named <- initial_effort(params)[1:2]
    named[] <- c(0.2, NA)
    p2 <- projectToSteady(params, effort = named, t_per = 0.1,
                          t_max = 0.1, tol = 10, info_level = 0)
    expected <- validEffortVector(named, params)
    expect_equal(p2@initial_effort, expected)
})

test_that("projectToSteady accepts consumer update method", {
    params <- small_steady_params

    pc <- projectToSteady(params, t_per = 1, dt = 1,
                          t_max = 1, tol = 1000,
                          method = "predictor_corrector", info_level = 0)
    expect_s4_class(pc, "MizerParams")
    expect_true(all(is.finite(pc@initial_n)))
    expect_true(all(pc@initial_n >= 0))

    expect_error(
        projectToSteady(params, t_per = 1, dt = 1, t_max = 1,
                        method = "bogus"),
        "should be one of"
    )
})

# steady ----
test_that("steady works", {
    params <- small_steady_params
    params@species_params$gamma[2] <- 2000
    params <- setSearchVolume(params)
    sim <- steady(params, t_per = 1, t_max = 1, dt = 1, tol = 10,
                  return_sim = TRUE) |>
        suppressMessages()
    expect_s4_class(sim, "MizerSim")
    expect_snapshot_value(getRDD(sim@params), style = "deparse")
})

test_that("steady accepts consumer update method", {
    params <- small_steady_params
    params@species_params$gamma[2] <- 2000
    params <- setSearchVolume(params)

    p <- steady(params, t_per = 1, t_max = 1, dt = 1, tol = 10,
                method = "predictor_corrector") |>
        suppressMessages()

    expect_s4_class(p, "MizerParams")
    expect_true(all(is.finite(p@initial_n)))
    expect_true(all(p@initial_n >= 0))
})

test_that("steady() preserves parameters", {
    params <- small_steady_params

    params_rdd <- params
    params_rdd@rates_funcs$RDD <- "noRDD"
    p_rdd <- steady(params_rdd, t_per = 1, t_max = 1, dt = 1) |> suppressMessages()
    expect_equal(p_rdd@rates_funcs$RDD, "noRDD")

    params_rmax <- params
    species_params(params_rmax)$R_max <- 1.01 * species_params(params_rmax)$R_max
    p_erepro <- steady(params_rmax, t_per = 1, t_max = 1, dt = 1, tol = 10,
                       preserve = "erepro") |> suppressMessages()
    expect_equal(p_erepro@species_params$erepro, params_rmax@species_params$erepro)
    p_rl <- steady(params_rmax, t_per = 1, t_max = 1, dt = 1, tol = 10,
                   preserve = "reproduction_level") |> suppressMessages()
    expect_equal(getReproductionLevel(p_rl), getReproductionLevel(params_rmax))

    params_erepro <- params
    species_params(params_erepro)$erepro <- 1.01 * species_params(params_erepro)$erepro
    p_rmax <- steady(params_erepro, t_per = 1, t_max = 1, dt = 1, tol = 10,
                     preserve = "R_max") |> suppressMessages()
    expect_equal(p_rmax@species_params$R_max, params_erepro@species_params$R_max)
    expect_false(identical(p_rmax@time_modified, params_erepro@time_modified))
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

test_that("projectToSteady() converges with use_predation_diffusion", {
    params_d <- small_steady_params
    params_d@use_predation_diffusion <- TRUE
    initialN(params_d)[1, ] <- initialN(params_d)[1, ] * 3
    expect_message(
        projectToSteady(params_d, t_per = 1, dt = 1, tol = 1000),
        "Convergence was achieved"
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

test_that("steady() rebalances the resource to a genuine fixed point", {
    # For a normal (non-frozen) resource, steady() rebalances the resource so
    # that it is a genuine fixed point and issues no warning.
    expect_no_warning(
        ps <- suppressMessages(steady(NS_params_small, progress_bar = FALSE))
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

test_that("steady() rebalances a frozen resource capacity", {
    # Freeze the resource capacity away from its balanced value.
    pf <- suppressMessages(setResource(NS_params_small,
                                       resource_capacity = NS_params_small@cc_pp * 1.3,
                                       balance = FALSE))
    expect_false(is.null(comment(pf@cc_pp)))
    npp_before <- pf@initial_n_pp

    # steady() rebalances the frozen capacity silently, as earlier mizer
    # versions did, without a warning or message about the resource.
    expect_no_warning(ps <- suppressMessages(steady(pf, progress_bar = FALSE)))
    # The resource is now a genuine fixed point and no longer frozen ...
    expect_lt(resource_deviation(ps), 1e-8)
    expect_null(comment(ps@cc_pp))
    # ... and the resource abundance itself was preserved (only the capacity
    # was adjusted to balance it).
    expect_equal(unclass(ps@initial_n_pp), unclass(npp_before))
})
