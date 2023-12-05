local_edition(3)

test_that("completeSpeciesParams() works", {
    species_params <- NS_species_params
    
    # test w_mat
    sp <- species_params
    sp$w_mat[1] <- NA
    expect_message(sp <- completeSpeciesParams(sp), NA)
    expect_equal(sp$w_mat[1], sp$w_max[1] / 4)
    sp$w_mat[2:4] <- 100
    expect_warning(sp <- completeSpeciesParams(sp),
                   "For the species Sandeel, N.pout the value")
    expect_equal(sp$w_mat[2], sp$w_max[2] / 4)
    
    # test w_mat25
    sp <- species_params
    sp$w_mat25 <- c(NA, 1:11)
    expect_warning(completeSpeciesParams(sp), NA)
    sp$w_mat25[2:5] <- 21
    expect_warning(sp <- completeSpeciesParams(sp),
                   "For the species Sandeel, Dab the value")
    expect_true(is.na(sp$w_mat25[[2]]))
    expect_identical(sp$w_mat25[[3]], 21)
    
    # test w_min
    sp <- species_params
    sp$w_min <- c(NA, 1:11)
    expect_warning(completeSpeciesParams(sp), NA)
    sp$w_min[2:5] <- 21
    expect_warning(sp <- completeSpeciesParams(sp),
                   "For the species Sandeel, Dab the value")
    expect_identical(sp$w_min[[2]], 0.001)
    expect_identical(sp$w_min[[3]], 21)
    
    # test misspelling detection
    sp$wmat <- 1
    expect_warning(completeSpeciesParams(sp),
                   "very close to standard parameter names") %>%
        suppressWarnings()
    
    # minimal species_params
    sp <- data.frame(species = c(2, 1),
                     w_max = c(100, 1000),
                     stringsAsFactors = FALSE)
    expect_s3_class(sp <- completeSpeciesParams(sp), "data.frame")
    expect_equal(sp$w_mat, sp$w_max / 4)
    expect_equal(sp$alpha, c(0.6, 0.6))
    expect_equal(sp$interaction_resource, c(1, 1))
    expect_identical(rownames(sp), c("2", "1"))
})

test_that("completeSpeciesParams converts from length to weight", {
    sp <- data.frame(species = 1:2,
                     l_max = 1:2,
                     a = 0.01, b = 3)
    sp2 <- completeSpeciesParams(sp)
    expect_identical(sp2$w_max, c(0.01, 0.08))
})


test_that("set_species_param_default sets default correctly", {
    params <- NS_params
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
    params <- NS_params
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
    # TODO: reduce the tolerance below
    expect_equal(gamma_default / gamma_analytic, 
                 rep(1, length(gamma_default)),
                 tolerance = 0.1)
})


test_that("Setting species params works", {
    params <- newMultispeciesParams(NS_species_params, info_level = 0)
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
