# newMultispeciesParams ----
# * Dimensions are correct ----
test_that("constructor with species_params and interaction signature gives the right dimensions", {
    expect_message(
        params <- newMultispeciesParams(NS_species_params_small, inter_small, info_level = 3),
        "Because"
    )
    # expect_that(params, is_a("MizerParams")) # deprecated, trying to find alternative
    expect_equal(class(params)[1], "MizerParams") # alternative?
    expect_equal(dim(params@psi)[1], nrow(NS_species_params_small))
    expect_equal(dimnames(params@psi)$sp, NS_species_params_small$species)
    params_gears <- newMultispeciesParams(NS_species_params_gears_small, inter_small, info_level = 0)
    expect_equal(unique(dimnames(params_gears@selectivity)$gear),
                 unique(params_gears@species_params$gear))
    # pass in other arguments
    params_gears <- newMultispeciesParams(NS_species_params_gears_small,
                                          inter_small, no_w = 50, info_level = 0)
    expect_length(params_gears@w, 50)
    expect_equal(dimnames(params_gears@selectivity)$gear,
                 unique(NS_species_params_gears_small$gear))
})

test_that("constructor with only species_params signature gives the right dimensions", {
    params <- newMultispeciesParams(NS_species_params_small, info_level = 0)
    expect_true(all(params@interaction == 1))
    expect_equal(dim(params@interaction), c(dim(params@psi)[1],
                                            dim(params@psi)[1]))
})

# * w_min_idx is correct ----
test_that("w_min_idx is being set correctly", {
    # default - no w_min in params data so set to first size
    params <- newMultispeciesParams(NS_species_params_gears_small, inter_small, info_level = 0)
    expect_true(all(params@species_params$w_min == params@w[1]))
    expect_true(all(params@w_min_idx == 1))
    # Set w_min to be the min by hand
    NS_species_params_gears_small$w_min <- 0.001
    params <- newMultispeciesParams(NS_species_params_gears_small, inter_small, info_level = 0)
    expect_true(all(params@w_min_idx == 1))
    # Change w_min of one of the species
    NS_species_params_gears_small$w_min <- 0.001
    NS_species_params_gears_small$w_min[2] <- 10
    params <- newMultispeciesParams(NS_species_params_gears_small, inter_small, info_level = 0)
    expect_true(all(params@w_min_idx[c(1, 3)] == 1))
    expect_equal(as.integer(params@w_min_idx[2]), max(which(params@w <= 10)))
})

test_that("Errors are reported", {
    expect_error(newMultispeciesParams(NS_species_params_small, min_w_pp = 1,
                                       info_level = 0),
                 "min_w_pp must be larger than min_w")
})

test_that("Sets given_species_params", {
    # Calling `given_species_params<-()` should not make a change
    sp <- data.frame(species = "sp1", w_max = 1000)
    params <- newMultispeciesParams(sp, info_level = 0)
    p2 <- params
    given_species_params(p2) <- given_species_params(p2)
    expect_unchanged(p2, params)
})

test_that("newMultispeciesParams sets initial resource spectrum and cutoff", {
    params <- newMultispeciesParams(NS_species_params_small,
                                    kappa = 10,
                                    lambda = 2,
                                    w_pp_cutoff = 1,
                                    info_level = 0)
    expect_equal(initialNResource(params)[w_full(params) < 1],
                 10 * w_full(params)[w_full(params) < 1] ^ (-2),
                 ignore_attr = TRUE)
    expect_true(all(initialNResource(params)[w_full(params) >= 1] == 0))
})

# setParams ----
test_that("setParams can leave params unchanged", {
    params <- setParams(NS_params_small, info_level = 0)
    expect_unchanged(setParams(params, info_level = 0), params)
    params@species_params$h <- NA
    expect_message(setParams(params, info_level = 3), "Because")
})

test_that("setParams handles change in w_max", {
    params <- NS_params_small
    # Check that ft_mask is recalculated correctly
    params@species_params$w_max[1] <- 1000
    params <- setParams(params)
    expect_equal(sum(params@ft_mask[1, ]), 38)

    # Check that warning is given if w_max is too large
    params@species_params$w_max[1] <- max(params@w) + 10
    params@species_params$w_repro_max[1] <- max(params@w) + 10
    expect_warning(setParams(params),
                 "The maximum weight of a species is larger than")
})

test_that("setParams reapplies line colours and linetypes from species_params", {
    params <- NS_params_small
    params@species_params$linecolour <- rep("#123456", nrow(species_params(params)))
    params@species_params$linetype <- rep("dashed", nrow(species_params(params)))
    params2 <- setParams(params)
    sp <- species_params(params2)$species
    expect_true(all(getColours(params2)[sp] == "#123456"))
    expect_true(all(getLinetypes(params2)[sp] == "dashed"))
    expect_identical(names(getColours(params2))[seq_along(sp)], sp)
    expect_identical(names(getLinetypes(params2))[seq_along(sp)], sp)
})
