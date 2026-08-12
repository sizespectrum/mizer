# emptyParams ----
# * test dimensions ----
test_that("basic constructor sets dimensions properly", {
    species_params <- NS_species_params_small[c(1, 2, 3), ]
    species_names <- species_params$species
    no_sp <- 3
    min_w <- 0.1
    max_w <- 5000
    no_w <- 200
    min_w_pp <- 1e-8
    expect_error(emptyParams(species_params, min_w = min_w, max_w = max_w,
                             no_w = no_w, min_w_pp = min_w_pp),
                 paste0("Some of your species have an maximum size larger than max_w: ", species_params$species[3]))
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
    expect_equal(dim(test_params@search_vol), c(no_sp, no_w))
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
                 c(length(gear_names), no_sp))
    expect_equal(dim(test_params_gears@selectivity), 
                 c(length(gear_names), no_sp, no_w))
    # dimnames of species and gears - just do a couple because the validity 
    # check should ensure the consistency of the others
    expect_equal(dimnames(test_params_gears@psi)$sp, species_names)
    expect_equal(dimnames(test_params_gears@catchability)$gear, gear_names)
})

test_that("emptyParams validates min_w_pp against min_w", {
    species_params <- NS_species_params_small[c(1, 2, 3), ]
    expect_error(emptyParams(species_params,
                             min_w = 0.1,
                             max_w = 40000,
                             no_w = 200,
                             min_w_pp = 0.1),
                 "min_w_pp must be larger than min_w")
})

# validMizerParams ----
test_that("Slots are allowed to have comments", {
    params <- NS_params_small
    comment(params) <- "All slots are given comments"
    for (slot in (slotNames(params))) {
        comment(slot(params, slot)) <- slot
    }
    expect_error(validObject(params), NA)
})


# size bins ----
test_that("w, w_full, dw, dw_full work", {
    params <- NS_params_small
    expect_identical(w(params), params@w)
    expect_identical(w_full(params), params@w_full)
    expect_identical(dw(params), params@dw)
    expect_identical(dw_full(params), params@dw_full)
})

# validParams ----
test_that("validParams works", {
    simc.0.4 <- readRDS("assets/simc.0.4.rds")
    (p <- validParams(simc.0.4@params)) |> 
        expect_message() |>
        expect_warning("Your MizerParams object was created with an earlier")
    expect_true(validObject(p))
    simc.1.0 <- readRDS("assets/simc.1.0.rds")
    (p <- validParams(simc.1.0@params)) |> 
        expect_message() |>
        expect_warning("Your MizerParams object was created with an earlier")
    expect_true(validObject(p))
})

test_that("validParams checks w_min and w_max", {
    # w_max
    params <- NS_params_small
    params@species_params$w_max[1] <- 1e6
    expect_warning(validParams(params), "The maximum weight of a species is larger")
    # w_min
    params <- NS_params_small
    params@species_params$w_min[1:3] <- 1e-6
    expect_warning(params <- validParams(params), "smaller than the minimum")
    expect_true(validObject(params))
    expect_equal(params@w_min_idx[1:3], rep(1, 3), ignore_attr = TRUE)
})

test_that("validParams skips the repair only for an unchanged object", {
    clear_validated_params()
    params <- NS_params_small

    # The first validation does the full work and records the fingerprint.
    expect_false(is_validated(validation_key(params)))
    params <- validParams(params)
    expect_true(is_validated(validation_key(params)))

    # Every change to a slot that the repair or the validity checks depend on
    # gives a new fingerprint, so the full validation runs again.
    p <- params
    p@species_params$w_mat[1] <- p@species_params$w_mat[1] / 2
    expect_false(is_validated(validation_key(p)))
    p <- params
    p@given_species_params$new_column <- 1
    expect_false(is_validated(validation_key(p)))
    p <- params
    p@gear_params$catchability[1] <- 0.5
    expect_false(is_validated(validation_key(p)))
    p <- params
    p@w_min_idx[1] <- p@w_min_idx[1] + 1
    expect_false(is_validated(validation_key(p)))
    p <- params
    p@ft_mask[1, 1] <- !p@ft_mask[1, 1]
    expect_false(is_validated(validation_key(p)))
    p <- params
    p@second_order_w$bin_average <- TRUE
    expect_false(is_validated(validation_key(p)))
    p <- params
    p@mizer_version <- "99.0.0"
    expect_false(is_validated(validation_key(p)))

    # Changing the shape of a rate array is caught, so that the structural
    # checks in validObject() cannot be evaded by the fast path.
    p <- params
    p@psi <- p@psi[, 1:10]
    expect_false(is_validated(validation_key(p)))
    p <- params
    dimnames(p@interaction)[[1]][1] <- "Wrong"
    expect_false(is_validated(validation_key(p)))

    # Changing values that the repair does not depend on does not.
    p <- params
    p@psi[1, 1] <- 0
    expect_true(is_validated(validation_key(p)))
    p@initial_effort[1] <- 100
    expect_true(is_validated(validation_key(p)))
})

test_that("validParams repairs and validates an object that is not yet known", {
    clear_validated_params()
    params <- NS_params_small
    # A model whose w_min_idx and ft_mask have been made inconsistent is
    # invalid, but validParams() recalculates them.
    params@w_min_idx[] <- 1
    params@ft_mask[] <- FALSE
    expect_error(validObject(params))
    expect_true(validObject(validParams(params)))

    # The repair is not skipped for an object that merely resembles one that
    # has been validated before.
    valid <- validParams(NS_params_small)
    params <- valid
    params@ft_mask[] <- FALSE
    expect_error(validObject(params))
    expect_true(validObject(validParams(params)))
})

test_that("validParams checks the array values even on the fast path", {
    clear_validated_params()
    params <- validParams(NS_params_small)
    expect_true(is_validated(validation_key(params)))
    # A bad value does not change the shape of the array and hence not the
    # fingerprint, so it has to be caught by the unconditional check.
    params@metab[1, 1] <- NaN
    expect_true(is_validated(validation_key(params)))
    expect_error(validParams(params), "metab must not contain non-finite values")
})

test_that("validParams leaves an already valid object untouched", {
    clear_validated_params()
    params <- validParams(NS_params_small)
    expect_identical(validParams(params), params)
    clear_validated_params()
    # The result must not depend on whether the fingerprint is known.
    expect_identical(validParams(params), params)
})

test_that("the record of validated params objects is bounded", {
    clear_validated_params()
    for (i in 1:3) record_validated(paste0("key", i), max_size = 3)
    expect_length(ls(validated_params), 3)
    record_validated("key4", max_size = 3)
    expect_identical(ls(validated_params), "key4")
    clear_validated_params()
})

test_that("validParams rejects non-finite rate arrays and allows infinite intake_max", {
    params <- NS_params_small
    params@search_vol[1, 1] <- NaN
    expect_error(validParams(params), "search_vol must not contain non-finite values")

    params <- NS_params_small
    params@intake_max[1, 1] <- Inf
    expect_no_error(validParams(params))

    params <- NS_params_small
    params@intake_max[1, 1] <- NaN
    expect_error(validParams(params),
                 "intake_max must contain finite or infinite numeric values only")
})
