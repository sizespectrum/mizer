# emptyParams ----
# * test dimensions ----
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
    expect_equal(dim(test_params@catchability), c(0, no_sp))
    expect_equal(dim(test_params@selectivity), c(0, no_sp, no_w))
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

# validMizerParams ----
test_that("Slots are allowed to have comments", {
    params <- NS_params
    comment(params) <- "All slots are given comments"
    for (slot in (slotNames(params))) {
        comment(slot(params, slot)) <- slot
    }
    expect_error(validObject(params), NA)
    expect_error(project(params, t_max = 0.1), NA)
})

# setColours, getColours ----
test_that("setColours and getColours works", {
    params <- NS_params
    no_col <- length(getColours(params))
    # set new entry
    params <- setColours(params, list("test" = "orange"))
    expect_equal(length(getColours(params)), no_col + 1)
    expect_identical(getColours(params)[["test"]], "orange")
    # overwrite existing and set new
    params <- setColours(params, list("test" = "blue", test2 = "orange"))
    expect_equal(length(getColours(params)), no_col + 2)
    expect_identical(getColours(params)[["test"]], "blue")
    expect_identical(getColours(params)[["test2"]], "orange")
})

# setLinetypes, getLinetypes ----
test_that("setLinetypes and getLinetypes works", {
    params <- NS_params
    no_types <- length(getLinetypes(params))
    # set new entry
    params <- setLinetypes(params, list("test" = "dashed"))
    expect_equal(length(getLinetypes(params)), no_types + 1)
    expect_identical(getLinetypes(params)[["test"]], "dashed")
    # overwrite existing and set new
    params <- setLinetypes(params, list("test" = "dotted", test2 = "dashed"))
    expect_equal(length(getLinetypes(params)), no_types + 2)
    expect_identical(getLinetypes(params)[["test"]], "dotted")
    expect_identical(getLinetypes(params)[["test2"]], "dashed")
})

# size bins ----
test_that("w, w_full, dw, dw_full work", {
    params <- NS_params
    expect_identical(w(params), params@w)
    expect_identical(w_full(params), params@w_full)
    expect_identical(dw(params), params@dw)
    expect_identical(dw_full(params), params@dw_full)
})

# validParams ----
test_that("validParams works", {
    simc.0.4 <- readRDS("assets/simc.0.4.rds")
    expect_warning(p <- validParams(simc.0.4@params),
                   "You need to upgrade your MizerParams object")
    expect_true(validObject(p))
    simc.1.0 <- readRDS("assets/simc.1.0.rds")
    expect_warning(p <- validParams(simc.1.0@params),
                   "You need to upgrade your MizerParams object")
    expect_true(validObject(p))
})