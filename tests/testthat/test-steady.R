# projectToSteady ----
small_steady_params <- function() {
    suppressMessages(newTraitParams(no_sp = 2, no_w = 20, max_w_max = 100,
                                    min_w = 1e-3, w_pp_cutoff = 5, ks = 4,
                                    reproduction_level = 0.25,
                                    info_level = 0))
}

test_that("projectToSteady() works", {
    params <- small_steady_params()
    effort <- params@initial_effort * 1.1
    initialN(params)[1, ] <- initialN(params)[1, ] * 3
    expect_error(projectToSteady(params, dt = 1, t_per = 0.5),
                 "t_per must be a positive multiple of dt")
    expect_error(projectToSteady(params, t_max = 0.1),
                 "t_max not greater than or equal to t_per")
    expect_message(projectToSteady(params, t_max = 0.1, t_per = 0.1),
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
    params <- NS_params
    p1 <- suppressMessages(projectToSteady(params, effort = 0.5, t_per = 0.1,
                                           t_max = 0.1, tol = 10))
    expect_equal(unname(p1@initial_effort), rep(0.5, length(initial_effort(params))))

    named <- initial_effort(params)[1:2]
    named[] <- c(0.2, NA)
    p2 <- suppressMessages(projectToSteady(params, effort = named, t_per = 0.1,
                                           t_max = 0.1, tol = 10))
    expected <- validEffortVector(named, params)
    expect_equal(p2@initial_effort, expected)
})

test_that("projectToSteady accepts consumer update method", {
    params <- small_steady_params()

    pc <- suppressMessages(projectToSteady(params, t_per = 1, dt = 1,
                                           t_max = 1, tol = 1000,
                                           method = "predictor_corrector"))
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
    params <- small_steady_params()
    params@species_params$gamma[2] <- 2000
    params <- setSearchVolume(params)
    p <- steady(params, t_per = 1, t_max = 1, dt = 1, tol = 10) |>
        suppressMessages()
    expect_snapshot_value(getRDD(p), style = "deparse")
    # and works the same when returning sim
    sim <- steady(params, t_per = 1, t_max = 1, dt = 1, tol = 10,
                  return_sim = TRUE) |>
        suppressMessages()
    expect_s4_class(sim, "MizerSim")
    expect_snapshot_value(getRDD(sim@params), style = "deparse")
})

test_that("steady accepts consumer update method", {
    params <- small_steady_params()
    params@species_params$gamma[2] <- 2000
    params <- setSearchVolume(params)

    p <- steady(params, t_per = 1, t_max = 1, dt = 1, tol = 10,
                method = "predictor_corrector") |>
        suppressMessages()

    expect_s4_class(p, "MizerParams")
    expect_true(all(is.finite(p@initial_n)))
    expect_true(all(p@initial_n >= 0))
})

test_that("steady() preserves reproduction function", {
    params <- small_steady_params()
    params@rates_funcs$RDD <- "noRDD"
    p2 <- steady(params, t_per = 1, t_max = 1, dt = 1) |>
        suppressMessages()
    expect_equal(params@rates_funcs$RDD, "noRDD")
})
test_that("steady() preserves erepro", {
    params <- small_steady_params()
    species_params(params)$R_max <- 1.01 * species_params(params)$R_max
    p2 <- steady(params, t_per = 1, t_max = 1, dt = 1, tol = 10,
                 preserve = "erepro") |>
        suppressMessages()
    expect_equal(p2@species_params$erepro, params@species_params$erepro)
})
test_that("steady() preserves reproduction level", {
    params <- small_steady_params()
    species_params(params)$R_max <- 1.01 * species_params(params)$R_max
    p2 <- steady(params, t_per = 1, t_max = 1, dt = 1, tol = 10,
                 preserve = "reproduction_level") |>
        suppressMessages()
    expect_equal(getReproductionLevel(p2), getReproductionLevel(params))
})
test_that("steady() preserves R_max", {
    params <- small_steady_params()
    species_params(params)$erepro <- 1.01 * species_params(params)$erepro
    p2 <- steady(params, t_per = 1, t_max = 1, dt = 1, tol = 10,
                 preserve = "R_max") |>
        suppressMessages()
    expect_equal(p2@species_params$R_max, params@species_params$R_max)
    # Test that steady updates time_modified
    expect_false(identical(p2@time_modified, params@time_modified))
})

# valid_species_arg ----
test_that("valid_species_arg works", {
    # character species argument
    expect_warning(s <- valid_species_arg(NS_params, c("non", "sense")),
                   "The following species do not exist: non, sense")
    expect_identical(s, vector(mode = "character"))

    sp1 <- NS_params@species_params$species[3]
    sp2 <- NS_params@species_params$species[2]
    sp_sprat <- NS_params@species_params$species[1]
    sp3 <- NS_params@species_params$species[3]

    expect_identical(valid_species_arg(NS_params, c(sp1, sp2)),
                     c(sp1, sp2))
    expect_identical(valid_species_arg(NS_params, c(sp_sprat, sp2),
                                       return.logical = TRUE),
                     c(TRUE, TRUE, FALSE))
    expect_error(
        valid_species_arg(NS_params, c("non", "sense"), error_on_empty = TRUE) |>
            suppressWarnings(),
        "No species have been selected.")
    # numeric species argument
    expect_warning(s <- valid_species_arg(NS_params, c(2.5, 4)),
                 "A numeric 'species' argument should only contain the integers 1 to 3")
    expect_identical(s, vector(mode = "character"))
    expect_identical(valid_species_arg(NS_params, c(3, 1)),
                     c(sp3, sp_sprat))
    expect_identical(valid_species_arg(NS_params, c(1, 3)),
                     c(sp_sprat, sp3))
    expect_identical(valid_species_arg(NS_params, c(3, 1),
                                       return.logical = TRUE),
                     c(TRUE, FALSE, TRUE))
    expect_error(
        valid_species_arg(NS_params, 4, error_on_empty = TRUE) |>
            suppressWarnings(),
        "No species have been selected.")
    # logical species argument
    expect_error(valid_species_arg(NS_params, c(TRUE, FALSE)),
                 "The boolean `species` argument has the wrong length")
    expect_identical(valid_species_arg(NS_params,
                                       c(TRUE, FALSE, TRUE)),
                     c(sp_sprat, sp3))
    expect_identical(valid_species_arg(NS_params,
                                       c(TRUE, FALSE, TRUE),
                                       return.logical = TRUE),
                     c(TRUE, FALSE, TRUE))
    expect_error(
        valid_species_arg(NS_params, rep(FALSE, 3), error_on_empty = TRUE) |>
            suppressWarnings(),
        "No species have been selected.")
    # called with MizerSim object
    sim <- project(NS_params, t_max = 1, dt = 1)
    expect_identical(valid_species_arg(sim, sp1),
                     valid_species_arg(NS_params, sp1))
    # called without species
    expect_identical(valid_species_arg(NS_params),
                     valid_species_arg(NS_params, 
                                       NS_params@species_params$species))
})

test_that("projectToSteady() converges with use_predation_diffusion", {
    params_d <- small_steady_params()
    params_d@use_predation_diffusion <- TRUE
    initialN(params_d)[1, ] <- initialN(params_d)[1, ] * 3
    expect_message(
        projectToSteady(params_d, t_per = 1, dt = 1, tol = 1000),
        "Convergence was achieved"
    )
})

test_that("valid_gears_arg works", {
    all_gears <- unique(NS_params@gear_params$gear)
    expect_identical(valid_gears_arg(NS_params), all_gears)
    expect_identical(valid_gears_arg(NS_params, all_gears[2:1]), all_gears[2:1])
    expect_warning(
        gears <- valid_gears_arg(NS_params, c("nope", all_gears[1])),
        "The following gears do not exist: nope"
    )
    expect_identical(gears, all_gears[1])
    expect_error(valid_gears_arg(NS_params, "nope", error_on_empty = TRUE) |>
                     suppressWarnings(),
                 "No gears have been selected.")
})

# constant_other ----
test_that("constant_other returns component value", {
    params <- NS_params
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
test_that("distanceMaxRelRDI calculates max relative RDI change", {
    params <- NS_params
    # Create two states with different abundances
    current <- list(
        n = initialN(params) * 1.1,  # 10% increase
        n_pp = initialNResource(params),
        n_other = list()
    )
    previous <- list(
        n = initialN(params),
        n_pp = initialNResource(params),
        n_other = list()
    )
    
    # Distance should be non-negative
    distance <- distanceMaxRelRDI(params, current, previous)
    expect_true(distance >= 0)
    expect_true(is.numeric(distance))
    expect_equal(length(distance), 1)
    
    # Distance should be 0 when states are identical
    distance_zero <- distanceMaxRelRDI(params, previous, previous)
    expect_equal(distance_zero, 0)
})

test_that("distance functions implement their documented formulas", {
    params <- NS_params
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
    expect_equal(distanceMaxRelRDI(params, current, previous), expected_rel)

    current$n[1, 1] <- 0
    previous$n[1, 2] <- 0
    sel <- current$n > 0 & previous$n > 0
    expected_log <- sum((log(current$n[sel]) - log(previous$n[sel]))^2)
    expect_equal(distanceSSLogN(params, current, previous), expected_log)
})

test_that("distanceSSLogN calculates sum of squared log differences", {
    params <- NS_params
    # Create two states with different abundances
    current <- list(
        n = initialN(params) * 1.1,  # 10% increase
        n_pp = initialNResource(params),
        n_other = list()
    )
    previous <- list(
        n = initialN(params),
        n_pp = initialNResource(params),
        n_other = list()
    )
    
    # Distance should be non-negative
    distance <- distanceSSLogN(params, current, previous)
    expect_true(distance >= 0)
    expect_true(is.numeric(distance))
    expect_equal(length(distance), 1)
    
    # Distance should be 0 when states are identical
    distance_zero <- distanceSSLogN(params, previous, previous)
    expect_equal(distance_zero, 0)
    
    # Distance should handle zero abundances gracefully
    current_with_zeros <- current
    current_with_zeros$n[1, 1] <- 0
    expect_error(distanceSSLogN(params, current_with_zeros, previous), NA)
})
