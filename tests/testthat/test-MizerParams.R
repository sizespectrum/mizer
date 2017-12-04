context("MizerParams constructor dimension checks")

test_that("basic constructor sets dimensions properly",{
    expect_that(MizerParams(1), is_a("MizerParams"))
    no_sp <- 3
    min_w <- 0.1
    max_w <- 5000
    no_w <- 200
    min_w_pp <- 1e-8
    species_names <- c("Cod","Haddock","Whiting")
    test_params <- MizerParams(no_sp, min_w = min_w, max_w = max_w, no_w = no_w, min_w_pp = min_w_pp, species_names = species_names)
    # Lengths of sizes OK?
    expect_that(length(test_params@w), equals(no_w))
    expect_that(length(test_params@dw), equals(no_w))
    
    # Check that that log of w_full is evenly spaced
    expect_that(max(diff(log(test_params@w_full))), equals(min(diff(log(test_params@w_full)))))
    # values of sizes OK?
    expect_that(test_params@w[1], equals(min_w))
    expect_that(test_params@w[length(test_params@w)], equals(max_w))
    expect_that(test_params@dw[1], equals(test_params@w[2]-test_params@w[1]))
    expect_that(test_params@dw[length(test_params@dw)], equals(test_params@dw[length(test_params@dw - 1)]))
    expect_that(test_params@w_full[1], equals(min_w_pp))
    # Test that first weight entry after plankton spectrum equals smallest fish weight 
    expect_that(test_params@w_full[1+length(test_params@w_full)-length(test_params@w)], equals(test_params@w[1]))
    # Dimensions of array slots
    expect_that(dim(test_params@psi), equals(c(no_sp,no_w)))
    expect_that(dim(test_params@intake_max), equals(c(no_sp,no_w)))
    expect_that(dim(test_params@search_vol), equals(c(no_sp,no_w)))
    expect_that(dim(test_params@activity), equals(c(no_sp,no_w)))
    expect_that(dim(test_params@std_metab), equals(c(no_sp,no_w)))
    expect_that(dim(test_params@ft_pred_kernel_e), equals(c(no_sp,length(test_params@w_full))))
    expect_that(dim(test_params@catchability), equals(c(no_sp,no_sp)))
    expect_that(dim(test_params@selectivity), equals(c(no_sp,no_sp, no_w)))
    expect_that(dim(test_params@interaction), equals(c(no_sp,no_sp)))
    # lengths of the other slots
    expect_that(length(test_params@rr_pp), equals(length(test_params@w_full))) 
    expect_that(length(test_params@cc_pp), equals(length(test_params@w_full))) 
    # Final check to make sure that the gears are being treated properly
    gear_names <- c("Trawl","Pelagic")
    test_params_gears <- MizerParams(no_sp, min_w = min_w, max_w = max_w, no_w = no_w, min_w_pp = min_w_pp, species_names = species_names, gear_names = gear_names)
    expect_that(dim(test_params_gears@catchability), equals(c(length(gear_names),no_sp)))
    expect_that(dim(test_params_gears@selectivity), equals(c(length(gear_names),no_sp, no_w)))
    # dimnames of species and gears - just do a couple because the validity check should ensure the consistency of the others
    expect_that(dimnames(test_params_gears@psi)$sp, equals(species_names))
    expect_that(dimnames(test_params_gears@catchability)$gear, equals(gear_names))
})

test_that("constructor with species_params and interaction signature gives the right dimensions",{
    data(NS_species_params_gears)
    data(NS_species_params)
    data(inter)
    test_params <- MizerParams(4) # is fine  
    test_params <- MizerParams(NS_species_params) # is fine
    test_params <- MizerParams(NS_species_params, inter) # seems fine
    expect_that(test_params, is_a("MizerParams"))
    expect_that(dim(test_params@psi)[1], equals(nrow(NS_species_params)))
    expect_that(dimnames(test_params@psi)$sp, equals(as.character(NS_species_params$species)))
    expect_that(dimnames(test_params@selectivity)$gear, equals(dimnames(test_params@selectivity)$sp))
    test_params_gears <- MizerParams(NS_species_params_gears, inter)  
    expect_that(unique(dimnames(test_params_gears@selectivity)$gear), equals(as.character(unique(test_params_gears@species_params$gear))))
    # pass in other arguments
    test_params_gears <- MizerParams(NS_species_params_gears, inter, no_w = 50)  
    expect_that(length(test_params_gears@w), equals(50))
})

test_that("constructor with only species_params signature gives the right dimensions",{
    data(NS_species_params_gears)
    data(NS_species_params)
    test_params <- MizerParams(NS_species_params)  
    expect_that(all(test_params@interaction==1),is_true()) 
    expect_that(dim(test_params@interaction), equals(c(dim(test_params@psi)[1],dim(test_params@psi)[1])))
})


test_that("w_min_idx is being set correctly",{
    data(NS_species_params_gears)
    data(inter)
    # default - no w_min in params data so set to first size
    params <- MizerParams(NS_species_params_gears, inter)
    expect_that(all(params@species_params$w_min == params@w[1]), is_true())
    expect_that(all(params@species_params$w_min_idx == 1), is_true())
    # Set w_min to be the min by hand
    NS_species_params_gears$w_min <- 0.001
    params <- MizerParams(NS_species_params_gears, inter)
    expect_that(all(params@species_params$w_min_idx == 1), is_true())
    # Change w_min of one of the species
    NS_species_params_gears$w_min <- 0.001
    NS_species_params_gears$w_min[7] <- 10
    params <- MizerParams(NS_species_params_gears, inter)
    expect_that(all(params@species_params$w_min_idx[c(1:6,8:12)] == 1), is_true())
    expect_that(params@species_params$w_min_idx[7], equals(max(which(params@w <= 10))))
    expect_error(MizerParams(NS_species_params_gears,inter, min_w = 1))
})

test_that("h and gamma and z0 are being set properly if not passed in",{

})
