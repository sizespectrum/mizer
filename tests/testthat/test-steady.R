# projectToSteady ----
test_that("projectToSteady() works", {
    params <- NS_params
    effort <- params@initial_effort * 1.1
    initialN(params)[1, ] <- initialN(params)[1, ] * 3
    expect_error(projectToSteady(NS_params, dt = 1, t_per = 0.5),
                 "t_per must be a positive multiple of dt")
    expect_error(projectToSteady(NS_params, t_max = 0.1),
                 "t_max not greater than or equal to t_per")
    expect_message(projectToSteady(params, t_max = 0.1, t_per = 0.1),
                   "Simulation run did not converge after 0.1 years.")
    expect_message(paramsc <- projectToSteady(params, tol = 10, effort = effort),
                   "Convergence was achieved in 4.5 years")
    expect_s4_class(paramsc, "MizerParams")
    expect_identical(paramsc@initial_effort, effort)
    # shouldn't take long the second time we run to steady
    expect_message(projectToSteady(paramsc, tol = 10),
                   "Convergence was achieved in 1.5 years")
    
    # return sim
    expect_message(sim <- projectToSteady(params, return_sim = TRUE, tol = 10),
                   "Convergence was achieved in 4.5 years")
    expect_s4_class(sim, "MizerSim")
    
    # Alternative distance function
    expect_message(paramsc <- projectToSteady(params,
                                              distance_func = distanceMaxRelRDI,
                                              tol = 0.1),
                   "Convergence was achieved in 9 years.")
    # shouldn't take long the second time we run to steady
    expect_message(projectToSteady(paramsc,
                                   distance_func = distanceMaxRelRDI,
                                   tol = 0.1),
                   "Convergence was achieved in 1.5 years")
    
    # Check extinction
    params@psi[5:6, ] <- 0
    expect_warning(projectToSteady(params),
                   "Dab, Whiting are going extinct.")
})

# steady ----
test_that("steady works", {
    expect_message(params <- newTraitParams(no_sp = 4, no_w = 30,
                                            n = 2/3, lambda = 2 + 3/4 - 2/3,
                                            max_w_inf = 1e3, min_w = 1e-4,
                                            w_pp_cutoff = 10, ks = 4,
                                            reproduction_level = 0),
                   "Increased no_w to 36")
    params@species_params$gamma[2] <- 2000
    params <- setSearchVolume(params)
    p <- steady(params, t_per = 2)
    expect_known_value(getRDD(p), "values/steady")
    # and works the same when returning sim
    sim <- steady(params, t_per = 2, return_sim = TRUE)
    expect_is(sim, "MizerSim")
    expect_known_value(getRDD(sim@params), "values/steady")
})

test_that("steady() preserves erepro", {
    params <- NS_params
    species_params(params)$R_max <- 1.01 * species_params(params)$R_max
    p2 <- steady(params, t_per = 1, preserve = "erepro")
    expect_equal(p2@species_params$erepro, params@species_params$erepro)
})
test_that("steady() preserves reproduction level", {
    params <- NS_params
    species_params(params)$R_max <- 1.01 * species_params(params)$R_max
    p2 <- steady(params, t_per = 1, preserve = "reproduction_level")
    expect_equal(getReproductionLevel(p2), getReproductionLevel(params))
})
test_that("steady() preserves R_max", {
    params <- NS_params
    species_params(params)$erepro <- 1.01 * species_params(params)$erepro
    expect_warning(p2 <- steady(params, t_per = 1, preserve = "R_max"),
                   "The following species require an unrealistic reproductive")
    expect_equal(p2@species_params$R_max, params@species_params$R_max)
})

# valid_species_arg ----
test_that("valid_species_arg works", {
    expect_warning(s <- valid_species_arg(NS_params, c("non", "sense")),
                   "The following species do not exist: non, sense")
    expect_identical(s, vector(mode = "character"))

    expect_identical(valid_species_arg(NS_params, c("Cod", "Sandeel")),
                     c("Cod", "Sandeel"))
    expect_identical(valid_species_arg(NS_params, c("Sprat", "Sandeel"),
                                       return.logical = TRUE),
                     c(TRUE, TRUE, rep(FALSE, 10)))
    # numeric species argument
    expect_warning(s <- valid_species_arg(NS_params, c(2.5, 13)),
                 "A numeric 'species' argument should only contain the integers 1 to 12")
    expect_identical(s, vector(mode = "character"))
    expect_identical(valid_species_arg(NS_params, c(3, 1)),
                     c("N.pout", "Sprat"))
    expect_identical(valid_species_arg(NS_params, c(1, 3)),
                     c("Sprat", "N.pout"))
    expect_identical(valid_species_arg(NS_params, c(3, 1),
                                       return.logical = TRUE),
                     c(TRUE, FALSE, TRUE, rep(FALSE, 9)))
    # logical species argument
    expect_error(valid_species_arg(NS_params, c(TRUE, FALSE)),
                 "The boolean `species` argument has the wrong length")
    expect_identical(valid_species_arg(NS_params, 
                                       c(TRUE, FALSE, TRUE, rep(FALSE, 9))),
                     c("Sprat", "N.pout"))
    expect_identical(valid_species_arg(NS_params, 
                                       c(TRUE, FALSE, TRUE, rep(FALSE, 9)),
                                       return.logical = TRUE),
                     c(TRUE, FALSE, TRUE, rep(FALSE, 9)))
    # called with MizerSim object
    sim <- project(NS_params, t_max = 1, dt = 1)
    expect_identical(valid_species_arg(sim, "Cod"),
                     valid_species_arg(NS_params, "Cod"))
    # called without species
    expect_identical(valid_species_arg(NS_params),
                     valid_species_arg(NS_params, 
                                       NS_params@species_params$species))
})
