context("MizerParams constructor dimension checks")

test_that("basic constructor sets dimensions properly", {
    expect_is(emptyParams(1), "MizerParams")
    no_sp <- 3
    min_w <- 0.1
    max_w <- 5000
    no_w <- 200
    min_w_pp <- 1e-8
    species_names <- c("Cod","Haddock","Whiting")
    test_params <- 
        emptyParams(no_sp, min_w = min_w, max_w = max_w, no_w = no_w, 
                    min_w_pp = min_w_pp, species_names = species_names)
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
    # Test that first weight entry after plankton spectrum equals smallest 
    # fish weight 
    expect_equal(test_params@w_full[1 + no_w_full - no_w], test_params@w[1])
    # Dimensions of array slots
    expect_equal(dim(test_params@psi), c(no_sp, no_w))
    expect_equal(dim(test_params@intake_max), c(no_sp, no_w))
    expect_equal(dim(test_params@search_vol), c(no_sp,no_w))
    expect_equal(dim(test_params@metab), c(no_sp, no_w))
    expect_equal(dim(test_params@ft_pred_kernel_e), c(no_sp, no_w_full))
    expect_equal(dim(test_params@catchability), c(no_sp, no_sp))
    expect_equal(dim(test_params@selectivity), c(no_sp, no_sp, no_w))
    expect_equal(dim(test_params@interaction), c(no_sp, no_sp))
    # lengths of the other slots
    expect_length(test_params@rr_pp, no_w_full) 
    expect_length(test_params@cc_pp, no_w_full) 
    # Final check to make sure that the gears are being treated properly
    gear_names <- c("Trawl", "Pelagic")
    test_params_gears <-
        emptyParams(no_sp, min_w = min_w, max_w = max_w,
                    no_w = no_w, min_w_pp = min_w_pp,
                    species_names = species_names, gear_names = gear_names)
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
    data(NS_species_params_gears)
    data(NS_species_params)
    data(inter)
    test_params <- MizerParams(NS_species_params, inter) # seems fine
    expect_that(test_params, is_a("MizerParams"))
    expect_equal(dim(test_params@psi)[1], nrow(NS_species_params))
    expect_equal(dimnames(test_params@psi)$sp, as.character(NS_species_params$species))
    expect_equal(dimnames(test_params@selectivity)$gear, dimnames(test_params@selectivity)$sp)
    test_params_gears <- MizerParams(NS_species_params_gears, inter)  
    expect_equal(unique(dimnames(test_params_gears@selectivity)$gear), 
                as.character(unique(test_params_gears@species_params$gear)))
    # pass in other arguments
    test_params_gears <- MizerParams(NS_species_params_gears, inter, no_w = 50)  
    expect_length(test_params_gears@w, 50)
})

test_that("constructor with only species_params signature gives the right dimensions", {
    data(NS_species_params_gears)
    data(NS_species_params)
    test_params <- MizerParams(NS_species_params)  
    expect_true(all(test_params@interaction == 1))
    expect_equal(dim(test_params@interaction), c(dim(test_params@psi)[1],
                                                 dim(test_params@psi)[1]))
})


test_that("w_min_idx is being set correctly", {
    data(NS_species_params_gears)
    data(inter)
    # default - no w_min in params data so set to first size
    params <- MizerParams(NS_species_params_gears, inter)
    expect_true(all(params@species_params$w_min == params@w[1]))
    expect_true(all(params@w_min_idx == 1))
    # Set w_min to be the min by hand
    NS_species_params_gears$w_min <- 0.001
    params <- MizerParams(NS_species_params_gears, inter)
    expect_true(all(params@w_min_idx == 1))
    # Change w_min of one of the species
    NS_species_params_gears$w_min <- 0.001
    NS_species_params_gears$w_min[7] <- 10
    params <- MizerParams(NS_species_params_gears, inter)
    expect_true(all(params@w_min_idx[c(1:6, 8:12)] == 1))
    expect_equal(params@w_min_idx[7], max(which(params@w <= 10)), 
                 check.names = FALSE)
    expect_error(MizerParams(NS_species_params_gears, inter, min_w = 1))
})
