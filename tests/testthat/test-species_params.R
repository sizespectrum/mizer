test_that("set_species_param_default sets default correctly", {
    params <- NS_params_small
    no_sp <- nrow(params@species_params)
    
    # Add comments to test that they are preserved
    comment(params@species_params) <- "top"
    comment(params@species_params$w_max) <- "test"
    # creates new column correctly
    expect_condition(set_species_param_default(params, "hype", 2, "hi"),
                   "hi", class = "info_about_default")
    p2 <- set_species_param_default(params, "hype", 2, "hi")
    expect_identical(p2@species_params$hype, rep(2, no_sp))
    expect_identical(comment(p2@species_params$w_max), "test")
    expect_identical(comment(p2@species_params), "top")
    expect_message(sp2 <- set_species_param_default(params@species_params, "hype", 3), NA)
    expect_identical(sp2$hype, rep(3, no_sp))
    # does not change existing colunn
    p2 <- set_species_param_default(params, "species", "a")
    expect_identical(p2, params)
    # changes NA's correctly
    sp1 <- params@species_params$species[1]
    params@species_params$species[1] <- NA
    params <- set_species_param_default(params, "species", sp1)
    expect_identical(p2, params)
    # Should throw errors
    expect_error(set_species_param_default(params, 1, "a"),
                 "parname is not a string")
})



test_that("default for gamma is correct", {
    params <- NS_params_small
    # check that missing h is o.k.
    params@species_params$alpha <- 0.1
    species_params <- params@species_params
    gamma_default <- get_gamma_default(params)
    # Compare to the analytic result
    lm2 <- params@resource_params$lambda - 2
    ae <- sqrt(2 * pi) * species_params$sigma * species_params$beta^lm2 *
        exp(lm2^2 * species_params$sigma^2 / 2) *
        # The factor on the following lines takes into account the cutoff
        # of the integral at 0 and at beta + 3 sigma
        (pnorm(3 - lm2 * species_params$sigma) + 
             pnorm(log(species_params$beta)/species_params$sigma + 
                       lm2 * species_params$sigma) - 1)
    if (!"h" %in% names(params@species_params) || 
        any(is.na(species_params$h))) {
        species_params$h <- get_h_default(params)
    }
    gamma_analytic <- (species_params$h / (params@resource_params$kappa * ae)) * 
        (species_params$f0 / (1 - species_params$f0))
    # The analytic formula is only approximate; tolerance depends on species params
    expect_equal(gamma_default / gamma_analytic,
                 rep(1, length(gamma_default)),
                 tolerance = 0.5)
})


test_that("Setting species params works", {
    params <- newMultispeciesParams(NS_species_params_small, info_level = 0)
    # changing h changes intake_max
    h_old <- params@species_params$h[[1]]
    intake_max_old <- params@intake_max[1, 1]
    species_params(params)$h[[1]] <- 1
    expect_identical(params@species_params$h[[1]], 1)
    expect_equal(params@intake_max[1, 1], 0.01)
    # setting to NA leads to recalculation of defaults
    species_params(params)$h[[1]] <- NA
    expect_identical(params@species_params$h[[1]], h_old)
    expect_identical(params@intake_max[1, 1], intake_max_old)
    
    # changing k_vb does not immediately change anything
    species_params(params)$k_vb[[1]] <- 2 * species_params(params)$k_vb[[1]]
    expect_identical(params@intake_max[1, 1], intake_max_old)
    # but clearing the default on h will lead to change
    species_params(params)$h[[1]] <- NA
    expect_equal(params@species_params$h[[1]], 2 * h_old)
    
    # increasing f0 increases gamma
    gamma_old <- species_params(params)$gamma[[1]]
    species_params(params)$f0 <- max(getFeedingLevel(params)) + 0.1
    species_params(params)$gamma <- NA
    expect_gt(species_params(params)$gamma[[1]], gamma_old)
    
    # increasing fc increases ks
    ks_old <- species_params(params)$ks[[1]]
    species_params(params)$fc <- max(getCriticalFeedingLevel(params)) + 0.1
    species_params(params)$ks <- NA
    expect_gt(species_params(params)$ks[[1]], ks_old)
    
    # changing w_min changes w_min_idx
    species_params(params)$w_min[[1]] <- 1
    expect_identical(params@w_min_idx[[1]], 40)
    
    # given species params are not affected
    beta <- params@given_species_params$beta
    species_params(params)$beta <- 1
    expect_identical(params@given_species_params$beta, beta)
})

test_that("Error if species names don't match", {
    sp <- NS_species_params_small
    sp$species[2] <- "not a species"
    expect_error(species_params(NS_params_small) <- sp,
                 "The species names in the new species parameter data frame do not match")
    expect_error(given_species_params(NS_params_small) <- sp,
                 "The species names in the new species parameter data frame do not match")
})


test_that("set_species_params_from_length works", {
    sp <- data.frame(species = 1:2, a = 0.01, b = 3)
    # Does nothing if no length
    expect_identical(set_species_param_from_length(sp, "w_mat", "l_mat"),
                     sp)
    # Converts as expected
    sp$l_mat <- c(1, 2)
    sp2 <- set_species_param_from_length(sp, "w_mat", "l_mat")
    expect_identical(sp2$w_mat, c(0.01, 0.08))
    # Can deal with NAs
    sp2$w_mat[2] <- NA
    sp2 <- set_species_param_from_length(sp2, "w_mat", "l_mat")
    expect_identical(sp2$w_mat, c(0.01, 0.08))
    # negative or zero lengths give error
    sp2$l_mat[2] <- 0
    expect_error(set_species_param_from_length(sp2, "w_mat", "l_mat"),
                 "All lengths should be positive and non-zero.")
})

test_that("`given_species_params<-()` gives correct warnings", {
    params <- NS_params_small
    
    no_sp <- nrow(params@species_params)
    expect_warning(given_species_params(params)$f0 <- 1)
    expect_warning(given_species_params(params)$fc <- 1)
    expect_warning(given_species_params(params)$age_mat <- 1)
    expect_warning(given_species_params(params)$catchability <- 2)
    expect_warning(given_species_params(params)$yield_observed <- 1)
    
    # No warning if NA
    params@given_species_params$gamma[-1] <- NA
    expect_warning(given_species_params(params)$f0 <- c(NA, rep(2, no_sp - 1)),
                   NA)
    
})

test_that("`given_species_params<-()` triggers recalculation", {
    params <- NS_params_small
    params@given_species_params$gamma <- NULL
    gamma <- params@species_params$gamma
    given_species_params(params)$f0 <- 0.1
    expect_gt(sum(gamma - params@species_params$gamma), 0)
})

test_that("`given_species_params<-()` can remove columns", {
    params <- NS_params_small
    given_species_params(params)$gamma <- NULL
    expect_false("gamma" %in% names(params@given_species_params))
    expect_true("gamma" %in% names(params@species_params))
    expect_error(given_species_params(params)$species <- NULL)
})

test_that("calculated_species_params returns only non-given values", {
    params <- NS_params_small

    calculated <- calculated_species_params(params)
    expect_false("species" %in% names(calculated))
    expect_true(all(vapply(calculated, function(col) !all(is.na(col)), logical(1))))

    given_species_params(params)$gamma <- NULL
    calculated <- calculated_species_params(params)
    expect_true("gamma" %in% names(calculated))
    expect_identical(calculated$gamma, species_params(params)$gamma)
})

test_that("species_params and given_species_params accessors return stored tables", {
    expect_identical(species_params(NS_params_small), NS_params_small@species_params)
    expect_identical(given_species_params(NS_params_small), NS_params_small@given_species_params)
})

test_that("species_params setter validates and recalculates", {
    params <- NS_params_small
    sp <- species_params(params)
    sp$w_min[1] <- 1
    species_params(params) <- sp
    expect_identical(species_params(params)$w_min[1], 1)
    idx <- unname(params@w_min_idx[1])
    expect_lte(w(params)[idx], 1)
    if (idx < length(w(params))) {
        expect_gt(w(params)[idx + 1], 1)
    }
})

test_that("given_species_params setter can add new explicit columns", {
    params <- NS_params_small
    sp <- given_species_params(params)
    sp$custom <- seq_len(nrow(sp))
    given_species_params(params) <- sp
    expect_identical(given_species_params(params)$custom, seq_len(nrow(sp)))
    expect_identical(species_params(params)$custom, seq_len(nrow(sp)))
})

test_that("set_species_param_default converts factors to character and fills NAs only", {
    sp <- species_params(NS_params_small)
    sp$dummy <- factor(rep("blue", nrow(sp)))
    sp$dummy[1] <- NA
    sp2 <- set_species_param_default(sp, "dummy", "black")
    expect_true(is.character(sp2$dummy))
    expect_identical(sp2$dummy[1], "black")
    expect_identical(sp2$dummy[-1], rep("blue", nrow(sp) - 1))
})

test_that("get_h_default, get_f0_default and get_ks_default follow documented defaults", {
    params <- NS_params_small

    sp <- species_params(params)
    sp$h <- rep(NA_real_, nrow(sp))
    sp$age_mat <- rep(NA_real_, nrow(sp))
    sp$k_vb <- rep(NA_real_, nrow(sp))
    h <- get_h_default(sp)
    expect_identical(h, rep(30, nrow(sp)))

    params2 <- params
    params2@species_params$f0[] <- 0.6
    params2@species_params$gamma[] <- NA
    expect_identical(get_f0_default(params2), rep(0.6, nrow(species_params(params2))))

    params3 <- params
    params3@species_params$ks[] <- NA
    params3@species_params$fc <- rep(0.2, nrow(species_params(params3)))
    expected_ks <- with(
        species_params(params3),
        fc * alpha * h * w_mat^(n - p)
    )
    expect_equal(get_ks_default(params3), expected_ks)
})
