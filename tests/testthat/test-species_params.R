
## Initialise ----
params <- NS_params
no_sp <- nrow(params@species_params)

## validSpeciesParams ----
test_that("validSpeciesParams() works", {
    species_params <- NS_species_params
    # test w_mat25
    species_params$w_mat25 <- c(NA, 1:11)
    expect_message(validSpeciesParams(species_params), NA)
    species_params$w_mat25[2:5] <- 21
    expect_message(validSpeciesParams(species_params),
                   "For the species Sandeel, Dab the value")
    # minimal species_params
    sp <- data.frame(species = c("species1", "species2"),
                     w_inf = c(100, 1000),
                     stringsAsFactors = FALSE)
    expect_s3_class(sp <- validSpeciesParams(sp), "data.frame")
    expect_equal(sp$w_mat, sp$w_inf / 4)
    expect_equal(sp$alpha, c(0.6, 0.6))
    expect_equal(sp$interaction_resource, c(1, 1))
    expect_identical(rownames(sp), c("species1", "species2"))
})

## set_species_param_default ----
test_that("set_species_param_default sets default correctly", {
    # creates new column correctly
    expect_message(p2 <- set_species_param_default(params, "hype", 2, "hi"),
                   "hi")
    expect_identical(p2@species_params$hype, rep(2, no_sp))
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


# Test default values ----
test_that("default for gamma is correct", {
    params <- NS_params
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
    expect_equal(gamma_default/ gamma_analytic, 
                 rep(1, length(gamma_default)),
                 tolerance = 0.1)
})

# species_params<-() ----
test_that("Setting species params works", {
    params <- newMultispeciesParams(NS_species_params)
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