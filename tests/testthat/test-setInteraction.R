params <- NS_params
no_sp <- nrow(params@species_params)

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
    params@species_params$interaction_resource <- -1
    expect_error(setInteraction(params),
                 "Values in the resource interaction vector should be between 0 and 1")
})
