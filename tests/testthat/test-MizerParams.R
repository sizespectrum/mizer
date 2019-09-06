context("MizerParams constructor dimension checks")
data(NS_species_params_gears)
data(NS_species_params)
data(inter)
no_sp <- nrow(NS_species_params)
params <- set_multispecies_model(NS_species_params, inter)

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
    # Test that first weight entry after plankton spectrum equals smallest 
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
    params_gears <- set_multispecies_model(NS_species_params_gears, inter)  
    expect_equal(unique(dimnames(params_gears@selectivity)$gear), 
                as.character(unique(params_gears@species_params$gear)))
    # pass in other arguments
    params_gears <- set_multispecies_model(NS_species_params_gears, inter, no_w = 50)  
    expect_length(params_gears@w, 50)
})

test_that("constructor with only species_params signature gives the right dimensions", {
    params <- set_multispecies_model(NS_species_params)  
    expect_true(all(params@interaction == 1))
    expect_equal(dim(params@interaction), c(dim(params@psi)[1],
                                                 dim(params@psi)[1]))
})

# w_min_idx is correct ----
test_that("w_min_idx is being set correctly", {
    # default - no w_min in params data so set to first size
    params <- set_multispecies_model(NS_species_params_gears, inter)
    expect_true(all(params@species_params$w_min == params@w[1]))
    expect_true(all(params@w_min_idx == 1))
    # Set w_min to be the min by hand
    NS_species_params_gears$w_min <- 0.001
    params <- set_multispecies_model(NS_species_params_gears, inter)
    expect_true(all(params@w_min_idx == 1))
    # Change w_min of one of the species
    NS_species_params_gears$w_min <- 0.001
    NS_species_params_gears$w_min[7] <- 10
    params <- set_multispecies_model(NS_species_params_gears, inter)
    expect_true(all(params@w_min_idx[c(1:6, 8:12)] == 1))
    expect_equal(params@w_min_idx[7], max(which(params@w <= 10)), 
                 check.names = FALSE)
})

# min_w_pp is correct ----
test_that("min_w_pp is being set correctly", {
    sp <- params@species_params
    sp$pred_kernel_type = "box"
    sp$ppmr_min <- 2
    sp$ppmr_max <- 4
    params <- set_multispecies_model(sp)
    min_w_feeding <- min(params@species_params$w_min / 4)
    expect_gte(min_w_feeding, params@w_full[1])
    expect_lte(min_w_feeding, params@w_full[3])
    # A single species can make a difference
    sp$ppmr_max[1] <- 100
    params <- set_multispecies_model(sp)
    expect_gte(params@species_params$w_min[1] / 100, params@w_full[1])
    expect_lte(params@species_params$w_min[1] / 100, params@w_full[2])
    # but only if it feeds on plankton
    sp$interaction_p <- 1
    sp$interaction_p[1] <- 0
    params <- set_multispecies_model(sp)
    expect_lte(min_w_feeding, params@w_full[3])
    # respect explicitly set min_w_pp
    expect_error(set_multispecies_model(sp, min_w_pp = 1),
                 "min_w_pp not less than or equal to min_w")
    expect_message(set_multispecies_model(sp, min_w_pp = 0.001),
                   "feeding kernels that extend below")
    params <- set_multispecies_model(sp, min_w_pp = 0.001)
    expect_identical(params@w_full[1], 0.001)
    # if none of the species feed on plankton, min_w_pp = min_w
    sp$interaction_p <- 0
    params <- set_multispecies_model(sp)
    expect_identical(params@w[1], params@w_full[1])
})