context("MizerParams constructor dimension checks")
no_sp <- nrow(NS_species_params)
params <- newMultispeciesParams(NS_species_params, inter)

# test dimensions ----
test_that("basic constructor sets dimensions properly", {
    species_params <- NS_species_params[c(6, 10, 11), ]
    species_names <- species_params$species
    no_sp <- 3
    min_w <- 0.1
    max_w <- 5000
    no_w <- 200
    min_w_pp <- 1e-8
    expect_error(emptyParams(species_params, min_w = min_w, max_w = max_w,
                             no_w = no_w, min_w_pp = min_w_pp),
                 "Some of your species have an maximum size larger than max_w: Cod")
    max_w <- 40000
    test_params <- 
        emptyParams(species_params, min_w = min_w, max_w = max_w,
                    no_w = no_w, min_w_pp = min_w_pp)
    # Lengths of sizes OK?
    expect_length(test_params@w, no_w)
    expect_length(test_params@dw, no_w)
    no_w_full <- length(test_params@w_full)
    
    # Check that that log of w_full is evenly spaced
    expect_equal(max(diff(log(test_params@w_full))), 
                 min(diff(log(test_params@w_full))))
    # values of sizes OK?
    expect_equal(test_params@w[1], min_w)
    expect_equal(test_params@w[no_w], max_w)
    expect_equal(test_params@dw[1], test_params@w[2] - test_params@w[1])
    expect_equal(test_params@w_full[1], min_w_pp)
    # Test that first weight entry after resource spectrum equals smallest 
    # fish weight 
    expect_equal(test_params@w_full[1 + no_w_full - no_w], test_params@w[1])
    # Dimensions of array slots
    expect_equal(dim(test_params@psi), c(no_sp, no_w))
    expect_equal(dim(test_params@intake_max), c(no_sp, no_w))
    expect_equal(dim(test_params@search_vol), c(no_sp,no_w))
    expect_equal(dim(test_params@metab), c(no_sp, no_w))
    expect_equal(dim(test_params@ft_pred_kernel_e), c(no_sp, no_w_full))
    expect_equal(dim(test_params@catchability), c(1, no_sp))
    expect_equal(dim(test_params@selectivity), c(1, no_sp, no_w))
    expect_equal(dim(test_params@interaction), c(no_sp, no_sp))
    # lengths of the other slots
    expect_length(test_params@rr_pp, no_w_full) 
    expect_length(test_params@cc_pp, no_w_full) 
    # Final check to make sure that the gears are being treated properly
    gear_names <- c("Trawl", "Pelagic")
    species_params$gear <- c("Trawl", "Pelagic", "Trawl")
    test_params_gears <-
        emptyParams(species_params, min_w = min_w, max_w = max_w,
                    no_w = no_w, min_w_pp = min_w_pp)
    expect_equal(dim(test_params_gears@catchability), 
                 c(length(gear_names),no_sp))
    expect_equal(dim(test_params_gears@selectivity), 
                 c(length(gear_names),no_sp, no_w))
    # dimnames of species and gears - just do a couple because the validity 
    # check should ensure the consistency of the others
    expect_equal(dimnames(test_params_gears@psi)$sp, species_names)
    expect_equal(dimnames(test_params_gears@catchability)$gear, gear_names)
})

test_that("constructor with species_params and interaction signature gives the right dimensions", {
    expect_that(params, is_a("MizerParams"))
    expect_equal(dim(params@psi)[1], nrow(NS_species_params))
    expect_equal(dimnames(params@psi)$sp, as.character(NS_species_params$species))
    expect_equal(dimnames(params@selectivity)$gear, "knife_edge_gear")
    params_gears <- newMultispeciesParams(NS_species_params_gears, inter)  
    expect_equal(unique(dimnames(params_gears@selectivity)$gear), 
                as.character(unique(params_gears@species_params$gear)))
    # pass in other arguments
    params_gears <- newMultispeciesParams(NS_species_params_gears, inter, no_w = 50)  
    expect_length(params_gears@w, 50)
})

test_that("constructor with only species_params signature gives the right dimensions", {
    params <- newMultispeciesParams(NS_species_params)  
    expect_true(all(params@interaction == 1))
    expect_equal(dim(params@interaction), c(dim(params@psi)[1],
                                                 dim(params@psi)[1]))
})

# w_min_idx is correct ----
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


test_that("Slots are allowed to have comments", {
    params <- NS_params
    comment(params) <- "All slots are given comments"
    for (slot in (slotNames(params))) {
        comment(slot(params, slot)) <- slot
    }
    expect_error(validObject(params), NA)
    expect_error(project(params, t_max = 0.1), NA)
})
