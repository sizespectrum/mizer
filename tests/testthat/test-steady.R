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
