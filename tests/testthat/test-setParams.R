context("Setting parameters")
## Initialise ----
data(NS_species_params_gears)
data(NS_species_params)
data(inter)
no_sp <- nrow(NS_species_params)
params <- MizerParams(NS_species_params, inter)

resource_dynamics <-
    list("detritus" = function(params, n, n_pp, B, rates, dt, ...) B["detritus"],
         "carrion" = function(params, n, n_pp, B, rates, dt, ...) B["carrion"])
NS_species_params$rho_detritus <- no_sp:1
NS_species_params$rho_carrion <- 1:no_sp
params_res <- MizerParams(NS_species_params, inter,
                          resource_dynamics = resource_dynamics)

## set_species_param_default ----
test_that("set_species_param_default sets default correctly", {
    # creates new column correctly
    expect_message(p2 <- set_species_param_default(params, "hype", 2, "hi"),
                   "hi")
    expect_identical(p2@species_params$hype, rep(2, no_sp))
    expect_message(sp2 <- set_species_param_default(params@species_params, "hype", 2), NA)
    expect_identical(sp2$hype, rep(2, no_sp))
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

## get_phi ----
test_that("get_phi works", {
    NS_species_params$pred_kernel_type <- "box"
    NS_species_params$ppmr_min <- 2
    NS_species_params$ppmr_max <- 4
    phi <- get_phi(NS_species_params, 1:5)
    expect_identical(phi[1, ], phi[2, ])
    expect_identical(phi[1, 1], 0)
    expect_identical(phi[1, 2], 1)
    expect_identical(phi[1, 5], 0)
})


## setInteraction ----
test_that("setInteraction works", {
    expect_identical(setInteraction(params, interaction = params@interaction),
                     params)
    inter <- matrix(1/2, nrow = no_sp, ncol = no_sp)
    p2 <- setInteraction(params, inter)
    expect_equivalent(p2@interaction, inter)
    inter[1, 1] <- 2
    expect_error(setInteraction(params, inter),
                 "Values in the interaction matrix must be between 0 and 1")
    expect_error(setInteraction(params, inter[1:(no_sp - 1), ]),
                 "interaction matrix is not of the right dimensions")
    inter[1, 1] <- 0
    dimnames(inter) <- list(sp = params@species_params$species,
                            sp = params@species_params$species)
    expect_message(setInteraction(params, inter),
                 "Your interaction matrix has dimensions called: `sp, sp`. I expected 'predator, prey'")
    dimnames(inter) <- list(predator = rev(params@species_params$species),
                            prey = params@species_params$species)
    expect_message(setInteraction(params, inter),
                   "Dimnames of interaction matrix do not match")
    params@species_params$interaction_p <- -1
    expect_error(setInteraction(params),
                 "Values in the plantkon interaction vector should be between 0 and 1")
})

## setPredKernel ----
test_that("setPredKernel works", {
    expect_identical(setPredKernel(params), params)
    expect_identical(setPredKernel(params, pred_kernel = NULL), 
                     params)
    params@species_params$pred_kernel_type <- "box"
    params@species_params$ppmr_min <- 2
    expect_error(setPredKernel(params), 
                 "missing from the parameter dataframe: ppmr_max")
    params@species_params$ppmr_max <- 4
    p2 <- setPredKernel(params)
    pred_kernel <- getPredKernel(params)
    expect_error(setPredKernel(params, pred_kernel[1:2, ]),
                 "incorrect number of dimensions")
    expect_error(setPredKernel(params, pred_kernel - 1),
                 "pred_kernel >= 0 are not true")
    p2 <- setPredKernel(params, pred_kernel)
    expect_equal(p2@ft_pred_kernel_e, array())
    expect_equal(p2@ft_pred_kernel_p, array())
    expect_equivalent(p2@pred_kernel, pred_kernel)
    expect_identical(p2@pred_kernel, getPredKernel(p2))
})


## setResourceEncounter ----
test_that("setResourceEncounter works", {
    species_params <- NS_species_params
    species_params$rho_detritus <- 1:no_sp
    species_params$rho_carrion <- no_sp:1
    params <- MizerParams(species_params, resource_dynamics = resource_dynamics)
    expect_equal(params@rho[2, 1, 1], 2 * params@w[1]^params@n)
    expect_equal(params@rho[2, 2, 1], (no_sp - 1) * params@w[1]^params@n)
})

## setParams ----
test_that("setParams can leave params unchanged", {
    expect_equal(setParams(params_res), params_res)
})