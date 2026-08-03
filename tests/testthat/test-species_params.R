test_that("set_species_param_default sets default correctly", {
    params <- NS_params_small
    no_sp <- nrow(params@species_params)

    # Add comments to test that they are preserved
    comment(params@species_params) <- "top"
    comment(params@species_params$w_max) <- "test"
    # creates new column correctly
    expect_condition(set_species_param_default(params, "hype", 2, "hi"),
                   "hi", class = "info_about_default")
    p2 <- set_species_param_default(params, "hype", 2, "hi")
    expect_identical(p2@species_params$hype, rep(2, no_sp), ignore_attr = TRUE)
    expect_identical(comment(p2@species_params$w_max), "test")
    expect_identical(comment(p2@species_params), "top")
    expect_message(sp2 <- set_species_param_default(params@species_params, "hype", 3), NA)
    expect_identical(sp2$hype, rep(3, no_sp), ignore_attr = TRUE)
    # does not change existing colunn
    p2 <- set_species_param_default(params, "species", "a")
    expect_identical(p2, params)
    # changes NA's correctly
    sp1 <- params@species_params$species[1]
    params@species_params$species[1] <- NA
    params <- set_species_param_default(params, "species", sp1)
    expect_identical(p2, params)
    # Should throw errors
    expect_error(set_species_param_default(params, 1, "a"),
                 "parname is not a string")
})



test_that("default for gamma is correct", {
    params <- NS_params_small
    # check that missing h is o.k.
    params@species_params$alpha <- 0.1
    species_params <- params@species_params
    gamma_default <- get_gamma_default(params)
    # Compare to the analytic result
    lm2 <- params@resource_params$lambda - 2
    ae <- sqrt(2 * pi) * species_params$sigma * species_params$beta^lm2 *
        exp(lm2^2 * species_params$sigma^2 / 2) *
        # The factor on the following lines takes into account the cutoff
        # of the integral at 0 and at beta + 3 sigma
        (pnorm(3 - lm2 * species_params$sigma) +
             pnorm(log(species_params$beta)/species_params$sigma +
                       lm2 * species_params$sigma) - 1)
    if (!"h" %in% names(params@species_params) ||
        any(is.na(species_params$h))) {
        species_params$h <- get_h_default(params)
    }
    gamma_analytic <- (species_params$h / (params@resource_params$kappa * ae)) *
        (species_params$f0 / (1 - species_params$f0))
    # The analytic formula is only approximate; tolerance depends on species params
    expect_equal(unname(gamma_default / gamma_analytic),
                 rep(1, length(gamma_default)),
                 tolerance = 0.5)
})


test_that("Setting species params works", {
    params <- newMultispeciesParams(NS_species_params_small, info_level = 0)
    # changing h changes intake_max
    h_old <- params@species_params$h[[1]]
    intake_max_old <- params@intake_max[1, 1]
    species_params(params)$h[[1]] <- 1
    expect_identical(params@species_params$h[[1]], 1)
    expect_equal(params@intake_max[1, 1], 0.01)
    # setting to NA leads to recalculation of defaults
    species_params(params)$h[[1]] <- NA
    expect_identical(params@species_params$h[[1]], h_old)
    expect_identical(params@intake_max[1, 1], intake_max_old)

    # changing k_vb changes h because h is missing from given_species_params
    species_params(params)$k_vb[[1]] <- 2 * species_params(params)$k_vb[[1]]
    expect_equal(params@species_params$h[[1]], 2 * h_old)

    # increasing f0 increases gamma
    gamma_old <- species_params(params)$gamma[[1]]
    species_params(params)$f0 <- max(getFeedingLevel(params)) + 0.1
    species_params(params)$gamma <- NA
    expect_gt(species_params(params)$gamma[[1]], gamma_old)

    # increasing fc increases ks
    ks_old <- species_params(params)$ks[[1]]
    species_params(params)$fc <- max(getCriticalFeedingLevel(params)) + 0.1
    species_params(params)$ks <- NA
    expect_gt(species_params(params)$ks[[1]], ks_old)

    # changing w_min changes w_min_idx
    species_params(params)$w_min[[1]] <- 1
    expect_identical(params@w_min_idx[[1]], 40)

    # given species params are updated by species_params<-
    species_params(params)$beta <- 1
    expect_identical(unname(params@given_species_params$beta), rep(1, 3))
})

test_that("Error if species names don't match", {
    sp <- NS_species_params_small
    sp$species[2] <- "not a species"
    expect_error(species_params(NS_params_small) <- sp,
                 "The species names in the new species parameter data frame do not match")
    expect_error(given_species_params(NS_params_small) <- sp,
                 "The species names in the new species parameter data frame do not match")
})


test_that("set_species_params_from_length works", {
    sp <- data.frame(species = 1:2, a = 0.01, b = 3)
    # Does nothing if no length
    expect_identical(set_species_param_from_length(sp, "w_mat", "l_mat"),
                     sp)
    # Converts as expected
    sp$l_mat <- c(1, 2)
    sp2 <- set_species_param_from_length(sp, "w_mat", "l_mat")
    expect_identical(sp2$w_mat, c(0.01, 0.08))
    # Can deal with NAs
    sp2$w_mat[2] <- NA
    sp2 <- set_species_param_from_length(sp2, "w_mat", "l_mat")
    expect_identical(sp2$w_mat, c(0.01, 0.08))
    # negative or zero lengths give error
    sp2$l_mat[2] <- 0
    expect_error(set_species_param_from_length(sp2, "w_mat", "l_mat"),
                 "All lengths should be positive and non-zero.")
})

test_that("`given_species_params<-()` gives correct warnings", {
    params <- NS_params_small

    no_sp <- nrow(params@species_params)
    expect_warning(given_species_params(params)$f0 <- 1)
    expect_warning(given_species_params(params)$fc <- 1)
    expect_warning(given_species_params(params)$age_mat <- 1)
    expect_warning(given_species_params(params)$catchability <- 2)
    expect_warning(given_species_params(params)$yield_observed <- 1)

    # No warning if NA
    params@given_species_params$gamma[-1] <- NA
    expect_warning(given_species_params(params)$f0 <- c(NA, rep(2, no_sp - 1)),
                   NA)

})

test_that("`given_species_params<-()` triggers recalculation", {
    params <- NS_params_small
    params@given_species_params$gamma <- NULL
    gamma <- params@species_params$gamma
    given_species_params(params)$f0 <- 0.1
    expect_gt(sum(gamma - params@species_params$gamma), 0)
})

test_that("`given_species_params<-()` can remove columns", {
    params <- NS_params_small
    given_species_params(params)$gamma <- NULL
    expect_false("gamma" %in% names(params@given_species_params))
    expect_true("gamma" %in% names(params@species_params))
    expect_error(given_species_params(params)$species <- NULL)
})

test_that("calculated_species_params returns only non-given values", {
    params <- NS_params_small

    calculated <- calculated_species_params(params)
    expect_true("species" %in% names(calculated))
    expect_identical(names(calculated)[1], "species")
    expect_true(all(vapply(calculated[, -1, drop = FALSE], function(col) !all(is.na(col)), logical(1))))

    given_species_params(params)$gamma <- NULL
    calculated <- calculated_species_params(params)
    expect_true("gamma" %in% names(calculated))
    expect_identical(calculated$gamma, species_params(params)$gamma)
})

test_that("species_params and given_species_params accessors return stored tables", {
    expect_identical(species_params(NS_params_small), NS_params_small@species_params)
    expect_identical(given_species_params(NS_params_small), NS_params_small@given_species_params)
})

test_that("species_params setter validates and recalculates", {
    params <- NS_params_small
    sp <- species_params(params)
    sp$w_min[1] <- 1
    species_params(params) <- sp
    expect_identical(species_params(params)$w_min[1], 1, ignore_attr = TRUE)
    idx <- unname(params@w_min_idx[1])
    expect_lte(w(params)[idx], 1)
    if (idx < length(w(params))) {
        expect_gt(w(params)[idx + 1], 1)
    }
})

test_that("species_params setter re-derives the calculated species parameters", {
    # The calculated species parameters that a rate setter owns must not be
    # carried over from the previous species parameter table, or a value
    # calculated earlier would look like a given one and stop responding to the
    # parameters it is derived from.
    params <- NS_params_small
    for (par in c("h", "gamma", "ks")) {
        params@given_species_params[[par]] <- NULL
    }
    old_h <- species_params(params)$h

    sp <- species_params(params)
    # Scaled down by a factor that stays above the calculated `w_mat25`
    # (`w_mat / 3^(1/10)`), so that `w_mat25` is not invalidated as well.
    sp$w_mat <- sp$w_mat * 0.95
    changed <- params
    species_params(changed) <- sp

    # `h` is derived from `w_mat`, so it has to follow, and `gamma` and `ks`
    # with it.
    expect_false(isTRUE(all.equal(species_params(changed)$h, old_h,
                                  check.attributes = FALSE)))

    # `given_species_params<-()` makes the same change, so both routes must
    # arrive at the same model.
    reference <- params
    given_sp <- given_species_params(reference)
    given_sp$w_mat <- given_sp$w_mat * 0.95
    given_species_params(reference) <- given_sp
    for (par in c("h", "gamma", "ks", "w_mat25")) {
        expect_equal(species_params(changed)[[par]],
                     species_params(reference)[[par]],
                     ignore_attr = TRUE, info = par)
    }
    expect_equal(changed@search_vol, reference@search_vol, ignore_attr = TRUE)
    expect_equal(changed@metab, reference@metab, ignore_attr = TRUE)
})

test_that("species_params setter keeps a given value of a calculated parameter", {
    # `NS_params_small` has `gamma` among the given species parameters, so it
    # must survive a change to a parameter that `gamma` is derived from.
    params <- NS_params_small
    old_gamma <- species_params(params)$gamma
    sp <- species_params(params)
    sp$w_mat <- sp$w_mat * 0.95
    species_params(params) <- sp
    expect_equal(species_params(params)$gamma, old_gamma, ignore_attr = TRUE)
})

test_that("species_params setter preserves columns mizer does not calculate", {
    # Columns that mizer knows nothing about are not rebuilt from the given
    # species parameters, so they have to be carried over.
    params <- NS_params_small
    params@species_params$my_own_param <- seq_len(nrow(species_params(params)))
    sp <- species_params(params)
    sp$w_mat <- sp$w_mat * 0.95
    species_params(params) <- sp
    expect_equal(species_params(params)$my_own_param,
                 seq_len(nrow(species_params(params))), ignore_attr = TRUE)
})

test_that("species_params setter handles list and matrix columns", {
    # The old-vs-new diff in `species_params<-()` must not use `==` on columns
    # where it is undefined (list, S4) or does not reduce to one logical per
    # species (matrix). It should fall back to a per-species `identical()`.
    params <- NS_params_small
    no_sp <- nrow(species_params(params))

    # List column, round-tripped unchanged: no error, value preserved.
    sp <- species_params(params)
    sp$groups <- lapply(seq_len(no_sp), function(i) letters[i])
    expect_error(species_params(params) <- sp, NA)
    expect_identical(species_params(params)$groups[[2]], "b")

    # Changing one element of the list column is detected and recorded.
    sp <- species_params(params)
    sp$groups[[1]] <- c("a", "z")
    species_params(params) <- sp
    expect_identical(species_params(params)$groups[[1]], c("a", "z"))
    expect_identical(given_species_params(params)$groups[[1]], c("a", "z"))

    # Matrix column, round-tripped unchanged: no error (a bare `==` would return
    # a matrix rather than one logical per species).
    params <- NS_params_small
    sp <- species_params(params)
    sp$mat <- matrix(seq_len(2 * no_sp), nrow = no_sp)
    expect_error(species_params(params) <- sp, NA)
    expect_equal(species_params(params)$mat, matrix(seq_len(2 * no_sp), nrow = no_sp),
                 ignore_attr = TRUE)
})

test_that("species_params setter with recalculate = FALSE records but does not recalculate", {
    params <- NS_params_small
    sp <- species_params(params)
    new_ks <- sp$ks * 2
    sp$ks <- new_ks

    quick <- params
    species_params(quick, recalculate = FALSE) <- sp

    # The new value is stored and recorded among the given parameters
    expect_equal(species_params(quick)$ks, new_ks, ignore_attr = TRUE)
    expect_equal(given_species_params(quick)$ks, new_ks, ignore_attr = TRUE)
    expect_s3_class(species_params(quick), "species_params")

    # ... but the rate it determines is left untouched, whereas the default
    # `recalculate = TRUE` would have doubled it.
    expect_identical(quick@metab, params@metab)
    recalculated <- params
    species_params(recalculated) <- sp
    expect_false(isTRUE(all.equal(recalculated@metab, params@metab)))
})

test_that("species_params setter with recalculate = FALSE does not re-derive calculated params", {
    params <- NS_params_small
    sp <- species_params(params)
    # `w_mat25` is a calculated parameter. Clearing it asks for it to be
    # derived afresh, which only the recalculation does.
    sp$w_mat25 <- NA

    quick <- params
    species_params(quick, recalculate = FALSE) <- sp
    expect_true(all(is.na(species_params(quick)$w_mat25)))

    recalculated <- params
    species_params(recalculated) <- sp
    expect_false(any(is.na(species_params(recalculated)$w_mat25)))
})

test_that("species_params setter with recalculate = FALSE fills in no defaults", {
    # Filling in a default would record it as a given value, which would turn a
    # parameter the model had left calculated into one that no longer responds
    # to the parameters it is derived from.
    params <- NS_params_small
    sp <- species_params(params)
    sp$a <- NULL
    sp$b <- NULL
    params@species_params <- sp
    params@given_species_params$a <- NULL
    params@given_species_params$b <- NULL

    sp$w_min <- sp$w_min * 2
    species_params(params, recalculate = FALSE) <- sp
    expect_false("a" %in% names(species_params(params)))
    expect_false("b" %in% names(given_species_params(params)))
    expect_setequal(names(species_params(params)), names(sp))
})

test_that("parameters set with recalculate = FALSE survive a later recalculation", {
    params <- NS_params_small
    sp <- species_params(params)
    new_ks <- sp$ks * 2
    sp$ks <- new_ks
    species_params(params, recalculate = FALSE) <- sp
    # The caller sets the matching rate itself
    params@metab[] <- params@metab * 2

    # Any later parameter change recalculates the rates, but from the recorded
    # `ks` rather than from a default.
    species_params(params)$w_mat <- species_params(params)$w_mat
    expect_equal(species_params(params)$ks, new_ks, ignore_attr = TRUE)
})

test_that("species_params setter with recalculate = FALSE still validates", {
    params <- NS_params_small
    sp <- species_params(params)
    sp$species <- rev(sp$species)
    expect_error(species_params(params, recalculate = FALSE) <- sp,
                 "species names in the new species parameter data frame do not match")
})

# length and weight parameters -------------------------------------------

# A model in which the maturity size is given as a length
length_based_params <- function() {
    params <- NS_params_small
    sp <- species_params(params)
    sp$l_mat <- w2l(sp$w_mat, sp)
    suppressMessages(species_params(params) <- sp)
    params
}

test_that("setting a weight parameter is not converted away by its length", {
    params <- length_based_params()
    sp <- species_params(params)
    new_w_mat <- sp$w_mat[[1]] * 1.5

    species_params(params)$w_mat[1] <- new_w_mat
    expect_equal(species_params(params)$w_mat[[1]], new_w_mat,
                 ignore_attr = TRUE)
    # and the length has followed the new weight
    expect_equal(species_params(params)$l_mat[[1]],
                 w2l(new_w_mat, species_params(params))[[1]],
                 ignore_attr = TRUE)
    # the other species are untouched
    expect_equal(species_params(params)$w_mat[-1], sp$w_mat[-1],
                 ignore_attr = TRUE)
})

test_that("all three assignment forms keep length and weight consistent", {
    params <- length_based_params()
    new_w_mat <- species_params(params)$w_mat[[1]] * 1.5

    for (assign in list(function(p) {species_params(p)$w_mat[1] <- new_w_mat; p},
                        function(p) {species_params(p)[1, "w_mat"] <- new_w_mat; p},
                        function(p) {species_params(p)[[1, "w_mat"]] <- new_w_mat; p})) {
        out <- species_params(assign(params))
        expect_equal(out$w_mat[[1]], new_w_mat, ignore_attr = TRUE)
        expect_equal(out$l_mat[[1]], w2l(new_w_mat, out)[[1]],
                     ignore_attr = TRUE)
    }
})

test_that("editing a species parameter data frame on its own does not validate", {
    # The checks and conversions happen when the data frame is validated, which
    # includes when it is assigned into a model, not on every assignment to it.
    sp <- species_params(length_based_params())
    l_mat_before <- sp$l_mat[[1]]

    sp$w_mat[1] <- sp$w_mat[[1]] * 1.5
    expect_equal(sp$l_mat[[1]], l_mat_before, ignore_attr = TRUE)
    expect_no_warning(sp$wmin <- 0.1)

    # ... and the class survives all three assignment forms
    expect_s3_class(sp, "species_params")
    sp[1, "w_mat"] <- sp$w_mat[[1]]
    expect_s3_class(sp, "species_params")
    sp[[1, "w_mat"]] <- sp$w_mat[[1]]
    expect_s3_class(sp, "species_params")
})

test_that("setting a length parameter still determines the weight", {
    params <- length_based_params()
    sp <- species_params(params)
    new_l_mat <- sp$l_mat[[1]] * 1.05

    # The weight follows the length, and that is not reported as an
    # inconsistency because the length is what the caller just set.
    expect_warning(species_params(params)$l_mat[1] <- new_l_mat, NA)
    expect_equal(species_params(params)$l_mat[[1]], new_l_mat,
                 ignore_attr = TRUE)
    expect_equal(species_params(params)$w_mat[[1]],
                 l2w(new_l_mat, species_params(params))[[1]],
                 ignore_attr = TRUE)
})

test_that("length to weight conversion is per species", {
    # A species whose length is not known must keep its weight even when
    # another species has inconsistent values.
    sp <- species_params(length_based_params())
    sp <- as.data.frame(sp)
    sp$l_mat[2] <- NA
    w_mat_2 <- sp$w_mat[[2]]
    sp$w_mat[1] <- sp$w_mat[[1]] * 1.1  # now inconsistent with its length

    out <- suppressWarnings(validSpeciesParams(sp))
    # The species whose length is unknown keeps its weight
    expect_equal(out$w_mat[[2]], w_mat_2, ignore_attr = TRUE)
    # The other species keeps its weight too, and its length now matches
    expect_equal(out$w_mat[[1]], sp$w_mat[[1]], ignore_attr = TRUE)
    expect_equal(out$l_mat[[1]], w2l(sp$w_mat, sp)[[1]], ignore_attr = TRUE)
})

test_that("a missing weight is still taken from the length", {
    sp <- as.data.frame(species_params(length_based_params()))
    sp$w_mat[2] <- NA
    out <- expect_warning(validSpeciesParams(sp), NA)
    expect_equal(out$w_mat[[2]], l2w(sp$l_mat, sp)[[2]], ignore_attr = TRUE)
})

test_that("when both are given at once the weight wins, with a warning", {
    # Both supplied together in one data frame, disagreeing: the weight is
    # authoritative and the length is set to match.
    sp <- data.frame(species = c("a", "b"), l_max = c(10, 20),
                     w_max = c(1, 1000), a = 0.01, b = 3)
    expect_warning(validSpeciesParams(sp),
                   paste0("the value of `l_max` is not consistent with the ",
                          "value of `w_max`.*a, b"))
    out <- suppressWarnings(validSpeciesParams(sp))
    expect_equal(out$w_max, c(1, 1000), ignore_attr = TRUE)
    expect_equal(out$l_max, w2l(c(1, 1000), sp), ignore_attr = TRUE)

    # No warning when they agree
    expect_warning(validSpeciesParams(as.data.frame(
        species_params(length_based_params()))), NA)
})

test_that("setting weight and length together lets the weight win", {
    params <- length_based_params()
    # Edit a plain data frame so that both changes arrive in the same
    # assignment rather than being reconciled one at a time.
    sp <- as.data.frame(species_params(params))
    new_w_mat <- sp$w_mat[[1]] * 1.5
    sp$w_mat[1] <- new_w_mat
    sp$l_mat[1] <- sp$l_mat[[1]] * 1.1

    suppressMessages(species_params(params) <- sp)
    expect_equal(species_params(params)$w_mat[[1]], new_w_mat,
                 ignore_attr = TRUE)
    expect_equal(species_params(params)$l_mat[[1]],
                 w2l(new_w_mat, species_params(params))[[1]],
                 ignore_attr = TRUE)
})

test_that("of two successive changes to the model the later one wins", {
    # Each of these assignments goes through `species_params<-()`, so the model
    # sees them one at a time and the rule can tell which came last.
    params <- length_based_params()

    # weight first, then length: the length is the later information
    p <- params
    species_params(p)$w_mat[1] <- species_params(p)$w_mat[[1]] * 1.5
    species_params(p)$l_mat[1] <- species_params(p)$l_mat[[1]] * 1.05
    sp <- species_params(p)
    expect_equal(sp$w_mat[[1]], l2w(sp$l_mat, sp)[[1]], ignore_attr = TRUE)

    # length first, then weight: now the weight is
    p <- params
    new_w_mat <- species_params(p)$w_mat[[1]] * 1.5
    species_params(p)$l_mat[1] <- species_params(p)$l_mat[[1]] * 1.05
    species_params(p)$w_mat[1] <- new_w_mat
    sp <- species_params(p)
    expect_equal(sp$w_mat[[1]], new_w_mat, ignore_attr = TRUE)
    expect_equal(sp$l_mat[[1]], w2l(new_w_mat, sp)[[1]], ignore_attr = TRUE)
})

test_that("given_species_params setter can add new explicit columns", {
    params <- NS_params_small
    sp <- given_species_params(params)
    sp$custom <- seq_len(nrow(sp))
    given_species_params(params) <- sp
    expect_equal(given_species_params(params)$custom, seq_len(nrow(sp)),
                 ignore_attr = TRUE)
    expect_equal(species_params(params)$custom, seq_len(nrow(sp)),
                 ignore_attr = TRUE)
})

test_that("set_species_param_default converts factors to character and fills NAs only", {
    sp <- species_params(NS_params_small)
    sp$dummy <- factor(rep("blue", nrow(sp)))
    sp$dummy[1] <- NA
    sp2 <- set_species_param_default(sp, "dummy", "black")
    expect_true(is.character(sp2$dummy))
    expect_identical(sp2$dummy[1], "black", ignore_attr = TRUE)
    expect_identical(sp2$dummy[-1], rep("blue", nrow(sp) - 1),
                     ignore_attr = TRUE)
})

test_that("get_h_default, get_f0_default and get_ks_default follow documented defaults", {
    params <- NS_params_small

    sp <- species_params(params)
    sp$h <- rep(NA_real_, nrow(sp))
    sp$age_mat <- rep(NA_real_, nrow(sp))
    sp$k_vb <- rep(NA_real_, nrow(sp))
    h <- get_h_default(sp)
    expect_identical(unname(h), rep(30, nrow(sp)))

    params2 <- params
    params2@species_params$f0[] <- 0.6
    params2@species_params$gamma[] <- NA
    expect_identical(unname(get_f0_default(params2)), rep(0.6, nrow(species_params(params2))))

    params3 <- params
    params3@species_params$ks[] <- NA
    params3@species_params$fc <- rep(0.2, nrow(species_params(params3)))
    expected_ks <- with(
        species_params(params3),
        fc * alpha * h * w_mat^(n - p)
    )
    expect_equal(unname(get_ks_default(params3)), unname(expected_ks))
})

test_that("species_params S3 class properties work", {
    params <- NS_params_small
    sim <- NS_sim_small

    # Test getters return species_params objects
    expect_true(is.species_params(species_params(params)))
    expect_true(is.species_params(given_species_params(params)))
    expect_true(is.species_params(species_params(sim)))

    # Test given_species_params specific class
    expect_true(is.given_species_params(given_species_params(params)))
    expect_false(is.given_species_params(species_params(params)))
    expect_true(is.given_species_params(given_species_params(sim)))

    # Test constructor on data frame
    df <- data.frame(species = c("Sprat", "Herring"), w_inf = c(10, 100))
    sp_df <- species_params(df)
    expect_true(is.species_params(sp_df))
    expect_identical(class(sp_df)[1], "species_params")

    given_df <- given_species_params(df)
    expect_true(is.given_species_params(given_df))
    expect_identical(class(given_df)[1:2], c("given_species_params", "species_params"))

    # Test constructor on already S3 species_params object
    expect_identical(species_params(sp_df), sp_df)
    expect_identical(given_species_params(given_df), given_df)

    # Test class preservation on subsetting and modifications
    expect_true(is.species_params(sp_df[1, ]))
    expect_true(is.species_params(sp_df[, 1, drop = FALSE]))

    expect_true(is.given_species_params(given_df[1, ]))
    expect_true(is.given_species_params(given_df[, 1, drop = FALSE]))

    sp_df$w_inf[1] <- 12
    expect_true(is.species_params(sp_df))

    given_df$w_inf[1] <- 12
    expect_true(is.given_species_params(given_df))

    sp_df[1, "w_inf"] <- 15
    expect_true(is.species_params(sp_df))

    given_df[1, "w_inf"] <- 15
    expect_true(is.given_species_params(given_df))

    sp_df[[1, "w_inf"]] <- 18
    expect_true(is.species_params(sp_df))

    given_df[[1, "w_inf"]] <- 18
    expect_true(is.given_species_params(given_df))
})

test_that("Validation and conversions work", {
    # These happen when a species parameter data frame is validated, which is
    # what `species_params()` does to a plain data frame and what
    # `species_params<-()` does to the table it is given.

    # 1. Misspelling warnings
    df <- data.frame(species = c("Sprat", "Herring"), w_inf = c(10, 100))
    df$wmin <- 0.1
    # (the misspelling check runs at both validation stages, hence the
    # `suppressWarnings()` around the expectation to swallow the duplicate)
    suppressWarnings(
        expect_warning(species_params(df),
                       "very close to standard parameter names"))
    # Fuzzy typo of a recognised name is flagged with a suggestion ...
    df2 <- data.frame(species = c("Sprat", "Herring"), w_inf = c(10, 100),
                      w_maxx = 100)
    suppressWarnings(expect_warning(species_params(df2), "did you mean `w_max`"))
    # ... but a genuine custom column is not
    df3 <- data.frame(species = c("Sprat", "Herring"), w_inf = c(10, 100),
                      my_note = "x")
    expect_no_warning(species_params(df3))

    # 2. Length-to-weight conversion
    dfl <- data.frame(species = "Sprat", a = 0.01, b = 3, l_mat = 10)
    sp <- species_params(dfl)
    # Check that w_mat was automatically calculated (0.01 * 10^3 = 10)
    expect_equal(sp$w_mat, 10, ignore_attr = TRUE)

    # On a standalone table there is nothing to say which of the pair is the
    # more recent, so both count as given at once and the weight wins.
    sp$l_mat <- 20
    suppressWarnings(expect_warning(out <- species_params(sp),
                                    "is not consistent with"))
    expect_equal(out$w_mat, 10, ignore_attr = TRUE)
    expect_equal(out$l_mat, 10, ignore_attr = TRUE)

    # In a model the change can be attributed, so an edited length wins there
    params <- suppressMessages(newMultispeciesParams(
        data.frame(species = "Sprat", w_inf = 100, a = 0.01, b = 3, l_mat = 10)))
    species_params(params)$l_mat <- 20
    expect_equal(species_params(params)$w_mat[[1]], 80, ignore_attr = TRUE)

    # 3. Consistency checks
    sp <- suppressWarnings(
        species_params(data.frame(species = "Sprat", w_inf = 50, w_mat = 80)))
    expect_equal(sp$w_mat, 12.5, ignore_attr = TRUE)  # corrected to w_inf / 4
    sp$w_mat <- 60
    suppressWarnings(
        expect_warning(species_params(sp),
                       "the value for `w_mat` is not smaller than that of `w_inf`"))
})

test_that("Print and Summary methods work", {
    df <- data.frame(species = c("Sprat", "Herring"), w_inf = c(10, 100), extra_param = c(1, 2))
    sp <- species_params(df)
    given <- given_species_params(df)
    gp <- gear_params(data.frame(gear = "g", species = "Sprat", catchability = 0.5, extra = 1))

    expect_output(print(sp), "An object of class \"species_params\"")
    expect_output(print(given), "An object of class \"given_species_params\"")
    expect_output(print(gp), "An object of class \"gear_params\"")
    expect_output(summary(sp), "Summary of species_params")
})

test_that("$ on species_params returns named vectors for non-character columns", {
    params <- NS_params_small
    sp <- species_params(params)
    sp_names <- sp$species  # character column

    # Numeric columns get species names
    expect_named(sp$w_mat, sp_names)

    # Character column stays unnamed (avoids self-referential names)
    expect_null(names(sp$species))

    # Works via accessor too
    expect_named(species_params(params)$w_mat, sp_names)
    expect_named(given_species_params(params)$w_inf, sp_names)

    # gear_params columns get row names
    gp <- gear_params(params)
    expect_named(gp$catchability, rownames(gp))
})

test_that("get_h_default accepts MizerParams, species_params and data.frame", {
    params <- NS_params_small
    sp <- species_params(params)
    df <- as.data.frame(sp)

    h_params <- get_h_default(params)
    h_sp <- get_h_default(sp)
    h_df <- get_h_default(df)

    expect_equal(h_params, h_sp)
    expect_equal(h_params, h_df)
})

# record_given_species_params ---------------------------------------------

test_that("record_given_species_params records changed values", {
    params <- NS_params_small
    sp_before <- species_params(params)
    value <- sp_before
    value$gamma <- value$gamma * 2

    given <- record_given_species_params(given_species_params(params),
                                         value, sp_before)

    expect_equal(given$gamma, sp_before$gamma * 2, ignore_attr = TRUE)
})

test_that("record_given_species_params records nothing when nothing changed", {
    params <- NS_params_small
    sp <- species_params(params)
    given <- given_species_params(params)

    expect_identical(record_given_species_params(given, sp, sp), given)
})

test_that("record_given_species_params only records the species that changed", {
    params <- NS_params_small
    sp_before <- species_params(params)
    value <- sp_before
    value$gamma[2] <- value$gamma[2] * 2

    given_before <- given_species_params(params)
    given <- record_given_species_params(given_before, value, sp_before)

    expect_equal(given$gamma[2], sp_before$gamma[2] * 2, ignore_attr = TRUE)
    expect_equal(given$gamma[-2], given_before$gamma[-2], ignore_attr = TRUE)
})

test_that("record_given_species_params treats NA as a value", {
    params <- NS_params_small
    sp_before <- species_params(params)
    sp_before$my_par <- NA_real_
    value <- sp_before
    value$my_par[1] <- 3.5

    given <- record_given_species_params(given_species_params(params),
                                         value, sp_before)

    expect_equal(unname(given$my_par[1]), 3.5)
    expect_true(all(is.na(given$my_par[-1])))
})

test_that("record_given_species_params records a wholly new column", {
    params <- NS_params_small
    sp_before <- species_params(params)
    value <- sp_before
    value$my_par <- seq_len(nrow(value)) + 0.5

    given <- record_given_species_params(given_species_params(params),
                                         value, sp_before)

    expect_equal(given$my_par, seq_len(nrow(value)) + 0.5, ignore_attr = TRUE)
})

test_that("recorded parameters survive a recalculation", {
    params <- NS_params_small
    sp_before <- species_params(params)
    new_ks <- sp_before$ks * 2

    # Without recording, the recalculation triggered by any species parameter
    # change derives ks afresh from the given parameters.
    unprotected <- params
    unprotected@species_params$ks <- new_ks
    species_params(unprotected)$w_mat <- species_params(unprotected)$w_mat
    expect_false(isTRUE(all.equal(species_params(unprotected)$ks, new_ks)))

    protected <- params
    protected@species_params$ks <- new_ks
    protected@given_species_params <-
        record_given_species_params(given_species_params(protected),
                                    species_params(protected), sp_before)
    species_params(protected)$w_mat <- species_params(protected)$w_mat
    expect_equal(species_params(protected)$ks, new_ks)
})

test_that("record_given_species_params checks the number of species", {
    params <- NS_params_small
    sp <- species_params(params)
    given <- given_species_params(params)

    expect_error(record_given_species_params(given, sp, sp[-1, ]),
                 "must all have one row per species")
    expect_error(record_given_species_params(given, sp[-1, ], sp),
                 "must all have one row per species")
})

test_that("record_given_species_params preserves the class of `given`", {
    params <- NS_params_small
    sp_before <- species_params(params)
    given_before <- given_species_params(params)

    # Both with and without a wholly new column, which used to be added with
    # `cbind()` and thereby strip the class.
    value <- sp_before
    value$my_par <- seq_len(nrow(value)) + 0.5
    for (v in list(sp_before, value)) {
        given <- record_given_species_params(given_before, v, sp_before)
        expect_s3_class(given, "given_species_params")
        expect_s3_class(given, "species_params")
        expect_identical(rownames(given), rownames(given_before))
        # `$` names the entries after the species, as for species_params()
        expect_identical(names(given$w_mat), rownames(given))
    }
})
