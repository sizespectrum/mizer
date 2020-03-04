context("Setting parameters")
## Initialise ----
no_sp <- nrow(NS_species_params)
params <- MizerParams(NS_species_params, inter)


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
                 "Values in the plankton interaction vector should be between 0 and 1")
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
test_that("Comment works on pred kernel", {
    pred_kernel <- getPredKernel(params)
    comment(pred_kernel) <- "test"
    params_c <- setPredKernel(params, pred_kernel = pred_kernel)
    expect_identical(comment(params_c@pred_kernel), "test")
    expect_message(setPredKernel(params_c),
                   "has been commented")
})

## setSearchVolume ----
test_that("setSearchVolume works", {
    expect_identical(setSearchVolume(params, params@search_vol), params)
    params@species_params$gamma <- 2 * params@species_params$gamma
    p2 <- setSearchVolume(params)
    expect_identical(2 * params@search_vol, p2@search_vol)
})
test_that("Comment works on search volume", {
    comment(params@search_vol) <- "test"
    params <- setSearchVolume(params, search_vol = params@search_vol)
    expect_identical(comment(params@search_vol), "test")
    expect_message(setSearchVolume(params),
                   "has been commented")
})

## setIntakeMax ----
test_that("ssetIntakeMax works", {
    expect_identical(setIntakeMax(params, params@intake_max), params)
    params@species_params$h <- 2 * params@species_params$h
    p2 <- setIntakeMax(params)
    expect_identical(2 * params@intake_max, p2@intake_max)
})
test_that("Comment works on intake_max", {
    comment(params@intake_max) <- "test"
    params <- setIntakeMax(params, intake_max = params@intake_max)
    expect_identical(comment(params@intake_max), "test")
    expect_message(setIntakeMax(params),
                   "has been commented")
})

## setMetab ----
test_that("setMetab works", {
    expect_identical(setMetab(params, params@metab), params)
    params@species_params$ks <- 2 * params@species_params$ks
    p2 <- setMetab(params)
    expect_identical(2 * params@metab, p2@metab)
})
test_that("Comment works on metab", {
    comment(params@metab) <- "test"
    params <- setMetab(params, metab = params@metab)
    expect_identical(comment(params@metab), "test")
    expect_message(setMetab(params),
                   "has been commented")
})

## setBMort ----
test_that("setBMort works", {
    expect_identical(setBMort(params, params@mu_b), params)
    params@species_params$z0 <- 2 * params@species_params$z0
    p2 <- setBMort(params)
    expect_identical(2 * params@mu_b, p2@mu_b)
})
test_that("Comment works on mu_b", {
    comment(params@mu_b) <- "test"
    params <- setBMort(params, mu_b = params@mu_b)
    expect_identical(comment(params@mu_b), "test")
    expect_message(setBMort(params),
                   "has been commented")
})

## setReproduction ----
test_that("setReproduction works", {
    expect_equal(setReproduction(params), params)
    maturity <- array(1, dim = c(no_sp, length(params@w)))
    p2 <- setReproduction(params, maturity = maturity)
    expect_equal(p2, setReproduction(p2, maturity = maturity,
                                     repro_prop = p2@psi))
    expect_equal(params, setReproduction(params, repro_prop = p2@psi))
    expect_error(setReproduction(params, srr = "sum"),
                 "Arguments of srr function must be 'rdi' and 'species_params'")
    params@species_params$erepro[1] <- NA
    p2 <- setReproduction(params, srr = "srrSheperd")
    expect_equal(p2@species_params$erepro[1], 1)
    p2@species_params$sheperd_b <- 0
    expect_error(getRDD(p2),
                 "The species_params dataframe must contain columns sheperd_b and sheperd_c.")
    p2 <- setReproduction(params, srr = "srrRicker")
    expect_error(getRDD(p2),
                 "The ricker_b column is missing in species_params")
    p2@species_params$ricker_b <- 0
    expect_equal(getRDI(p2), getRDD(p2))
})
test_that("Comment works on maturity", {
    comment(params@maturity) <- "test"
    params <- setReproduction(params, maturity = params@maturity)
    expect_identical(comment(params@maturity), "test")
    expect_message(setReproduction(params),
                   "maturity ogive has been commented")
})
test_that("Comment works on psi", {
    repro_prop <- params@psi
    comment(repro_prop) <- "test"
    params <- setReproduction(params, repro_prop = repro_prop)
    expect_identical(comment(params@psi), "test")
    expect_message(setReproduction(params),
                   "has been commented")
})

## setParams ----
test_that("setParams can leave params unchanged", {
    expect_equal(setParams(params), params)
})

## upgradeParams ----
test_that("upgradeParams leaves new params unchanged", {
    expect_identical(upgradeParams(params), params)
})
test_that("upgradeParams preserves comments", {
    comment(params) <- "test"
    for (slot in (slotNames(params))) {
        comment(slot(params, slot)) <- slot
    }
    expect_identical(upgradeParams(params), params)
})
