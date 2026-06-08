test_that("setInteraction works", {
    params <- NS_params_small
    no_sp <- nrow(params@species_params)
    
    expect_unchanged(setInteraction(params, interaction = params@interaction),
                     params)
    inter_small <- matrix(1/2, nrow = no_sp, ncol = no_sp)
    p2 <- setInteraction(params, inter_small)
    expect_equal(p2@interaction, inter_small, ignore_attr = TRUE)
    expect_error(setInteraction(params, inter_small[1:(no_sp - 1), ]),
                 "interaction matrix is not of the right dimensions")
    
    intera <- inter_small
    intera[1, 1] <- "a"
    expect_error(setInteraction(params, intera),
                 "The entries of the interaction matrix should be numeric.")
    intera <- inter_small
    intera[1, 1] <- -1
    expect_error(setInteraction(params, intera),
                 "All entries in the interaction matrix must be non-negative.")
    
    dimnames(inter_small) <- list(sp = params@species_params$species,
                            sp = params@species_params$species)
    expect_message(setInteraction(params, inter_small),
                   "Your interaction matrix has dimensions called: `sp, sp`. I expected 'predator, prey'")
    dimnames(inter_small) <- list(predator = rev(params@species_params$species),
                            prey = params@species_params$species)
    expect_message(setInteraction(params, inter_small),
                   "Dimnames of interaction matrix do not match")
    
    # If user only specifies column names, these are also used as rownames
    rownames(inter_small) <- as.character(1:no_sp)
    p2 <- setInteraction(params, inter_small)
    expect_identical(colnames(params@interaction), rownames(params@interaction))
    
    params@species_params$interaction_resource <- -1
    expect_error(setInteraction(params),
                 "Values in the resource interaction vector must be non-negative.")
    
})


test_that("getInteraction works", {
    params <- NS_params_small
    p <- setInteraction(params, interaction = params@interaction)
    expect_unchanged(params, p)
    expect_warning(expect_identical(getInteraction(params),
                                    interaction_matrix(params)),
                   "deprecated")
    expect_identical(interaction_matrix(params), params@interaction)
})

test_that("interaction_matrix setter works", {
    params <- NS_params_small
    new <- params@interaction
    new[1, 2] <- new[1, 2] / 2
    interaction_matrix(params) <- new
    expect_identical(interaction_matrix(params), new)
})
