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

test_that("Reactive validation and conversions work", {
    # 1. Misspelling warnings
    df <- data.frame(species = c("Sprat", "Herring"), w_inf = c(10, 100))
    sp <- species_params(df)
    expect_warning(sp$wmin <- 0.1, "very close to standard parameter names")

    # 2. Length-to-weight conversion
    df2 <- data.frame(species = "Sprat", a = 0.01, b = 3, l_mat = 10)
    sp2 <- species_params(df2)
    # Check that w_mat was automatically calculated (0.01 * 10^3 = 10)
    expect_equal(sp2$w_mat, 10, ignore_attr = TRUE)

    # Check conversion updates when l_mat is edited
    sp2$l_mat <- 20
    expect_equal(sp2$w_mat, 80, ignore_attr = TRUE)

    # 3. Consistency checks
    # Remove l_mat so the length-to-weight conversion does not override w_mat
    sp2$l_mat <- NULL
    # Setting w_inf < current w_mat (80) should warn and auto-correct w_mat
    expect_warning(sp2$w_inf <- 50, "the value for `w_mat` is not smaller than that of `w_inf`")
    # After auto-correction, w_mat should now be w_inf/4 = 12.5, so setting
    # w_mat <- 60 (> w_inf=50) should warn again
    expect_warning(sp2$w_mat <- 60, "the value for `w_mat` is not smaller than that of `w_inf`")
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

test_that("get_h_default S3 methods work", {
    params <- NS_params_small
    sp <- species_params(params)
    df <- as.data.frame(sp)

    h_params <- get_h_default(params)
    h_sp <- get_h_default(sp)
    h_df <- get_h_default(df)

    expect_equal(h_params, h_sp)
    expect_equal(h_params, h_df)
})
