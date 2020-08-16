# newMultispeciesParams ----
# * Dimensions are correct ----
test_that("constructor with species_params and interaction signature gives the right dimensions", {
    params <- newMultispeciesParams(NS_species_params, inter)
    expect_that(params, is_a("MizerParams"))
    expect_equal(dim(params@psi)[1], nrow(NS_species_params))
    expect_equal(dimnames(params@psi)$sp, as.character(NS_species_params$species))
    params_gears <- newMultispeciesParams(NS_species_params_gears, inter)  
    expect_equal(unique(dimnames(params_gears@selectivity)$gear), 
                 as.character(unique(params_gears@species_params$gear)))
    # pass in other arguments
    params_gears <- newMultispeciesParams(NS_species_params_gears, inter, no_w = 50)  
    expect_length(params_gears@w, 50)
    expect_equal(dimnames(params_gears@selectivity)$gear,
                 unique(NS_species_params_gears$gear))
})

test_that("constructor with only species_params signature gives the right dimensions", {
    params <- newMultispeciesParams(NS_species_params)  
    expect_true(all(params@interaction == 1))
    expect_equal(dim(params@interaction), c(dim(params@psi)[1],
                                            dim(params@psi)[1]))
})

# * w_min_idx is correct ----
test_that("w_min_idx is being set correctly", {
    # default - no w_min in params data so set to first size
    params <- newMultispeciesParams(NS_species_params_gears, inter)
    expect_true(all(params@species_params$w_min == params@w[1]))
    expect_true(all(params@w_min_idx == 1))
    # Set w_min to be the min by hand
    NS_species_params_gears$w_min <- 0.001
    params <- newMultispeciesParams(NS_species_params_gears, inter)
    expect_true(all(params@w_min_idx == 1))
    # Change w_min of one of the species
    NS_species_params_gears$w_min <- 0.001
    NS_species_params_gears$w_min[7] <- 10
    params <- newMultispeciesParams(NS_species_params_gears, inter)
    expect_true(all(params@w_min_idx[c(1:6, 8:12)] == 1))
    expect_equal(params@w_min_idx[7], max(which(params@w <= 10)), 
                 check.names = FALSE)
})

## setParams ----
test_that("setParams can leave params unchanged", {
    expect_equal(setParams(NS_params), NS_params)
})