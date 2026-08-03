test_that("markBackground marks the requested species in a MizerParams", {
    params <- markBackground(NS_params_small, species = "Sprat")
    expect_identical(unname(species_params(params)$is_background),
                     c(TRUE, FALSE, FALSE))
    # the other species parameters are untouched
    expect_identical(species_params(params)$species,
                     species_params(NS_params_small)$species)
})

test_that("markBackground with no species argument marks all species", {
    params <- markBackground(NS_params_small)
    expect_true(all(species_params(params)$is_background))
})

test_that("markBackground works on a MizerSim", {
    sim <- markBackground(NS_sim_small, species = "Cod")
    expect_identical(unname(species_params(sim@params)$is_background),
                     c(FALSE, FALSE, TRUE))
})

test_that("markBackground rejects objects of the wrong class", {
    expect_error(markBackground(1:3),
                 "The `object` argument must be of type MizerParams or MizerSim.")
})

test_that("removeBackgroundSpecies removes exactly the marked species", {
    params <- markBackground(NS_params_small, species = c("Sprat", "Cod"))
    reduced <- removeBackgroundSpecies(params)
    expect_identical(species_params(reduced)$species, "Herring")
    expect_identical(nrow(reduced@initial_n), 1L)
})

test_that("removeBackgroundSpecies is a no-op when nothing is marked", {
    reduced <- removeBackgroundSpecies(NS_params_small)
    expect_identical(species_params(reduced)$species,
                     species_params(NS_params_small)$species)
})
