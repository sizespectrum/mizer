test_that("projectToSteady() works", {
    params <- NS_params
    initialN(params)[1, ] <- initialN(params)[1, ] * 3
    expect_error(projectToSteady(NS_params, dt = 1, t_per = 0.5),
                 "t_per must be a positive multiple of dt")
    expect_error(projectToSteady(NS_params, t_max = 0.1),
                 "t_max not greater than or equal to t_per")
    expect_message(projectToSteady(params, t_max = 0.1, t_per = 0.1),
                   "Simulation run did not converge after 0.1 years.")
    expect_message(paramsc <- projectToSteady(params, tol = 10),
                   "Convergence was achieved in 4.5 years")
    expect_s4_class(paramsc, "MizerParams")
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
    expect_message(params <- newTraitParams(no_sp = 4, no_w = 30, R_factor = Inf,
                                            n = 2/3, lambda = 2 + 3/4 - 2/3,
                                            max_w_inf = 1e3, min_w = 1e-4,
                                            w_pp_cutoff = 10, ks = 4),
                   "Increased no_w to 36")
    params@species_params$gamma[2] <- 2000
    params <- setSearchVolume(params)
    p <- steady(params, t_per = 2)
    expect_known_value(getRDI(p), "values/steady")
    # and works the same when returning sim
    sim <- steady(params, t_per = 2, return_sim = TRUE)
    expect_is(sim, "MizerSim")
    expect_known_value(getRDI(sim@params), "values/steady")
})

# retune_erepro ----
test_that("retune_erepro works", {
    params <- NS_params
    params1 <- retune_erepro(params, species = "Cod")
    expect_equal(which(params1@species_params$erepro != 
                               params@species_params$erepro), 11)
    expect_identical(params1@rates_funcs$RDD, "noRDD")
})

# valid_species_arg ----
test_that("valid_species_arg works", {
    expect_warning(valid_species_arg(NS_params, c("Cod", "non", "sense")),
                   "The following species do not exist: non, sense")
    suppressWarnings(
        expect_error(valid_species_arg(NS_params, c("non", "sense")),
                   "The species argument matches none of the species in the params object")
    )
    expect_identical(valid_species_arg(NS_params, c("Cod", "Sandeel")),
                     c("Cod", "Sandeel"))
    expect_identical(valid_species_arg(NS_params, c("Sprat", "Sandeel"),
                                       return.logical = TRUE),
                     c(TRUE, TRUE, rep(FALSE, 10)))
    # numeric species argument
    expect_warning(valid_species_arg(NS_params, c(2.5, 3)),
                 "A numeric 'species' argument should only contain the integers 1 to 12")
    suppressWarnings(
        expect_error(valid_species_arg(NS_params, c(2.5, 13)),
                     "None of the numbers in the species argument are valid species indices.")
    )
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
})
