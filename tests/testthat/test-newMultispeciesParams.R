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

test_that("w_min survives a given_species_params<- round-trip", {
    sp <- data.frame(species = c("a", "b"), w_inf = c(100, 200))
    # min_w smaller than default 0.001
    params <- newMultispeciesParams(sp, min_w = 1e-4, info_level = 0)
    expect_true("w_min" %in% names(given_species_params(params)))
    expect_equal(species_params(params)$w_min[[1]], 1e-4)
    given_species_params(params)$alpha <- 0.6
    expect_equal(species_params(params)$w_min[[1]], 1e-4,
                 label = "w_min after round-trip (min_w < 0.001)")

    # min_w larger than default 0.001
    params2 <- newMultispeciesParams(sp, min_w = 0.01, info_level = 0)
    expect_equal(species_params(params2)$w_min[[1]], 0.01)
    expect_no_warning(given_species_params(params2)$alpha <- 0.6)
    expect_equal(species_params(params2)$w_min[[1]], 0.01,
                 label = "w_min after round-trip (min_w > 0.001)")
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

test_that("setParams rejects arguments it does not use", {
    params <- NS_params_small
    # A typo used to be silently ignored
    expect_error(setParams(params, metabolic = 99),
                 "does not have an argument `metabolic`")
    expect_error(setParams(params, metabolic = 99),
                 "Check for a typo")
    # Arguments belonging to other setters point at those setters
    expect_error(setParams(params, resource_rate = 5),
                 "Use `setResource\\(\\)` to set `resource_rate`")
    expect_error(setParams(params, kappa = 5),
                 "Use `setResource\\(\\)` to set `kappa`")
    expect_error(setParams(params, gear_params = gear_params(params)),
                 "Use `gear_params<-\\(\\)` to set `gear_params`")
    # Several at once are all reported
    expect_error(setParams(params, resource_rate = 5, metabolic = 99),
                 "does not have arguments `resource_rate`, `metabolic`")
    # Unnamed arguments beyond `interaction` and `info_level` are rejected
    expect_error(setParams(params, NULL, 0, 1),
                 "must be named")
    # Setter arguments are still recognised and warn when they are ignored.
    params@given_species_params$z0 <- params@species_params$z0
    expect_warning(setParams(params, z0pre = 0.5,
                             RDD = "BevertonHoltRDD", info_level = 0),
                   NA)
    expect_warning(setParams(params, z0pre = 0.5,
                             RDD = "BevertonHoltRDD", info_level = 1),
                   "already present in `given_species_params` for every species")
})

test_that("constructor records z0 only when z0 arguments are explicit", {
    sp <- data.frame(species = c("a", "b"), w_inf = c(100, 1000))

    defaults <- suppressMessages(newMultispeciesParams(sp, no_w = 20))
    expect_false("z0" %in% names(given_species_params(defaults)))

    expect_no_warning(params <- suppressMessages(
        newMultispeciesParams(sp, z0pre = 2, z0exp = -0.5, no_w = 20)))

    expected <- 2 * sp$w_inf^(-0.5)
    expect_equal(species_params(params)$z0, expected, ignore_attr = TRUE)
    expect_equal(given_species_params(params)$z0, expected,
                 ignore_attr = TRUE)
    expect_equal(params@mu_b,
                 outer(expected, rep(1, length(params@w))),
                 ignore_attr = TRUE)

    # Explicit z0 values still win, while only missing entries use the
    # construction arguments.
    sp$z0 <- c(0.1, NA)
    params <- suppressMessages(newMultispeciesParams(
        sp, z0pre = 2, z0exp = -0.5, no_w = 20))
    expect_equal(species_params(params)$z0, c(0.1, expected[[2]]),
                 ignore_attr = TRUE)
    expect_equal(given_species_params(params)$z0, c(0.1, expected[[2]]),
                 ignore_attr = TRUE)

    # Explicit default arguments are still given arguments.
    params <- suppressMessages(newMultispeciesParams(
        sp[c("species", "w_inf")], z0pre = 0.6, no_w = 20))
    expect_equal(given_species_params(params)$z0,
                 species_params(params)$z0, ignore_attr = TRUE)
})

test_that("setParams passes reset on to all the setters", {
    params <- NS_params_small
    # Freeze two of the arrays
    params <- setSearchVolume(params, search_vol = params@search_vol * 2)
    params <- setMetabolicRate(params, metab = params@metab * 2)
    expect_false(is.null(comment(params@search_vol)))
    expect_false(is.null(comment(params@metab)))

    # Without reset they stay frozen at their custom values
    unchanged <- setParams(params, info_level = 0)
    expect_false(different(unchanged@search_vol, params@search_vol))
    expect_false(different(unchanged@metab, params@metab))

    # With reset they are recalculated from the species parameters
    thawed <- setParams(params, reset = TRUE, info_level = 0)
    expect_null(comment(thawed@search_vol))
    expect_null(comment(thawed@metab))
    expect_false(different(thawed@search_vol,
                           setSearchVolume(params, reset = TRUE)@search_vol))

    expect_error(setParams(params, reset = "yes"), "not a flag")
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
