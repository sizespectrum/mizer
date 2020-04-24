
## Initialise ----
params <- NS_params
no_sp <- nrow(params@species_params)

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

# validSpeciesParams
test_that("validSpeciesParams() works", {
    species_params <- NS_species_params
    # test w_mat25
    species_params$w_mat25 <- c(NA, 1:11)
    expect_message(validSpeciesParams(species_params), NA)
    species_params$w_mat25[2:5] <- 21
    expect_message(validSpeciesParams(species_params),
                 "For the species Sandeel, Dab the value")
})