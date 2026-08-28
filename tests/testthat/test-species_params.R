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



test_that("the length-weight defaults report themselves", {
    sp <- data.frame(species = c("a", "b"), w_max = c(10, 100))
    expect_condition(species_params(sp), "using a = 0.01",
                     class = "info_about_default")
    expect_condition(species_params(sp), "isometric default b = 3",
                     class = "info_about_default")
    # but say nothing when the user has supplied the parameters
    sp$a <- c(0.01, 0.02)
    sp$b <- c(3, 3.1)
    msgs <- capture_messages(with_info_level(species_params(sp), info_level = 3))
    expect_false(any(grepl("w = a l^b", msgs, fixed = TRUE)))
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


test_that("gamma and f0 defaults ignore a frozen search volume (#488)", {
    params <- NS_params_small
    frozen <- params
    sv <- search_vol(frozen)
    sv[] <- sv * 10
    search_vol(frozen) <- sv
    expect_false(is.null(comment(frozen@search_vol)))

    # `get_gamma_default()` measures the available energy with a search volume
    # coefficient of 1, so the user's array must not enter the calculation.
    params@species_params$gamma[] <- NA
    frozen@species_params$gamma[] <- NA
    # `gamma` is of order 1e-11, far below the default tolerance, so the
    # comparison has to be of the ratio.
    expect_equal(unname(get_gamma_default(frozen) / get_gamma_default(params)),
                 rep(1, nrow(species_params(params))))

    # `get_f0_default()` is the inverse and must use the search volume implied
    # by `gamma`, not the frozen one.
    expect_equal(get_f0_default(NS_params_small),
                 suppressMessages(get_f0_default(
                     setSearchVolume(NS_params_small, search_vol = sv))))

    # The frozen array itself is left untouched by the default calculations.
    expect_equal(search_vol(frozen), sv, ignore_attr = TRUE)
})

test_that("gamma and f0 defaults ignore extension dispatch (#577)", {
    ext <- paste0("mizerTestGammaExt", Sys.getpid())
    chain <- setNames(NA_character_, ext)
    # An extension that scales the search volume, the way therMizer scales it
    # with its temperature scalar.
    registerS3method(
        "projectEncounter", ext,
        function(params, ...) {
            scalar <- params@other_params[["test_scalar"]]
            if (is.null(scalar)) scalar <- 0.5
            params@search_vol <- params@search_vol * scalar
            NextMethod()
        },
        envir = asNamespace("mizer")
    )

    params <- NS_params_small
    ext_params <- params
    ext_params@extensions <- chain
    ext_params <- coerceToExtensionClass(ext_params)

    # The extension really does change the encounter rate ...
    expect_equal(getEncounter(ext_params), getEncounter(params) / 2,
                 ignore_attr = TRUE)
    # ... but it must not change the defaults, which measure the available
    # energy with mizer's own encounter rate.
    expect_equal(get_f0_default(ext_params), get_f0_default(params))
    params@species_params$gamma[] <- NA
    ext_params@species_params$gamma[] <- NA
    # `gamma` is of order 1e-11, far below the default tolerance, so the
    # comparison has to be of the ratio.
    expect_equal(unname(get_gamma_default(ext_params) /
                            get_gamma_default(params)),
                 rep(1, nrow(species_params(params))))

    # A scalar of zero for one species used to turn this into an error
    zero <- ext_params
    zero@other_params$test_scalar <- c(1, 0, 1)
    expect_equal(unname(get_gamma_default(zero) / get_gamma_default(params)),
                 rep(1, nrow(species_params(params))))
})

test_that("a rebuild on an extension object leaves gamma alone (#577)", {
    ext <- paste0("mizerTestGammaRebuild", Sys.getpid())
    chain <- setNames(NA_character_, ext)
    registerS3method(
        "projectEncounter", ext,
        function(params, ...) {
            params@search_vol <- params@search_vol / 2
            NextMethod()
        },
        envir = asNamespace("mizer")
    )

    # Hand `gamma` back to mizer so that it is recalculated on every rebuild
    params <- NS_params_small
    suppressMessages(given_species_params(params)$gamma <- NULL)
    params@extensions <- chain
    params <- coerceToExtensionClass(params)
    gamma <- species_params(params)$gamma

    # Adding a column of the extension's own triggers a rebuild but changes
    # nothing that the encounter rate depends on
    suppressMessages(
        species_params(params)$my_extension_note <-
            seq_len(nrow(species_params(params)))
    )
    expect_equal(unname(species_params(params)$gamma / gamma),
                 rep(1, length(gamma)))
})

test_that("gamma and f0 defaults ignore a custom encounter function", {
    e <- globalenv()
    e$test_half_encounter <- function(params, n, n_pp, n_other, t = 0, ...) {
        mizerEncounter(params, n = n, n_pp = n_pp, n_other = n_other,
                       t = t, ...) / 2
    }
    withr::defer(rm("test_half_encounter", envir = e))

    params <- NS_params_small
    custom <- setRateFunction(params, "Encounter", "test_half_encounter")
    expect_equal(get_f0_default(custom), get_f0_default(params))

    params@species_params$gamma[] <- NA
    custom@species_params$gamma[] <- NA
    expect_equal(unname(get_gamma_default(custom) / get_gamma_default(params)),
                 rep(1, nrow(species_params(params))))
})

test_that("gamma and f0 defaults ignore additive encounters (#586)", {
    e <- globalenv()
    e$test_component_dynamics_586 <- function(params, n_other, component, ...) {
        n_other[[component]]
    }
    e$test_component_encounter_586 <- function(params, ...) {
        params@intake_max / 4
    }
    withr::defer(rm(list = c("test_component_dynamics_586",
                             "test_component_encounter_586"), envir = e))

    for (bin_average in c(FALSE, TRUE)) {
        params <- NS_params_small
        suppressMessages(
            second_order_w(params) <- c(bin_average = bin_average)
        )
        no_sp <- nrow(species_params(params))
        base_encounter <- getEncounter(params)
        base_f0 <- get_f0_default(params)
        target_f0 <- species_params(params)$f0

        with_ext <- params
        with_ext@ext_encounter[] <- with_ext@ext_encounter +
            with_ext@intake_max / 2
        expect_gt(sum(getEncounter(with_ext) - base_encounter), 0)
        expect_equal(get_f0_default(with_ext), base_f0)

        with_component <- setComponent(
            params,
            component = "test_food_586",
            initial_value = 0,
            dynamics_fun = "test_component_dynamics_586",
            encounter_fun = "test_component_encounter_586"
        )
        expect_gt(sum(getEncounter(with_component) - base_encounter), 0)
        expect_equal(get_f0_default(with_component), base_f0)

        no_gamma <- params
        no_gamma@species_params$gamma[] <- NA
        ext_no_gamma <- with_ext
        ext_no_gamma@species_params$gamma[] <- NA
        component_no_gamma <- with_component
        component_no_gamma@species_params$gamma[] <- NA
        gamma <- suppressMessages(get_gamma_default(no_gamma))
        expect_equal(
            unname(suppressMessages(get_gamma_default(ext_no_gamma)) / gamma),
            rep(1, no_sp), tolerance = 1e-14
        )
        expect_equal(
            unname(suppressMessages(get_gamma_default(component_no_gamma)) /
                       gamma),
            rep(1, no_sp), tolerance = 1e-14
        )

        with_both <- setComponent(
            with_ext,
            component = "test_food_586",
            initial_value = 0,
            dynamics_fun = "test_component_dynamics_586",
            encounter_fun = "test_component_encounter_586"
        )
        with_both@species_params$gamma[] <- NA
        with_both@species_params$gamma <-
            suppressMessages(get_gamma_default(with_both))
        expect_equal(get_f0_default(with_both), target_f0)
    }
})

test_that("a gamma that cannot be defaulted names the species (#586)", {
    withr::local_options(list(mizer_defaults_edition = 2))
    params <- NS_params_small
    sp <- species_params(params)$species
    # A species that eats nothing in the reference state has no `gamma` that
    # can give it the target feeding level. Under edition 2 the zero
    # interaction is respected, so the available energy comes out as zero.
    params@species_params$interaction_resource <- c(0, 0, 1)
    params@species_params$gamma[] <- NA
    # An external encounter used to stand in for the missing predation
    # encounter and produce a `gamma` off by orders of magnitude
    params@ext_encounter[] <- params@intake_max
    expect_error(suppressMessages(get_gamma_default(params)),
                 paste0("Could not calculate a default `gamma` for the ",
                        "following species: ", sp[[1]], ", ", sp[[2]], "\\. ",
                        "The available energy measured for them .* was 0, 0\\."),
                 fixed = FALSE)
})

test_that("invalid f0 is rejected even when gamma is given (#517)", {
    sp <- data.frame(species = "a", w_inf = 100)
    params <- newMultispeciesParams(sp, info_level = 0)
    expect_false("gamma" %in% names(given_species_params(params)))

    expect_error(given_species_params(params)$f0 <- 1,
                 "`f0` must be finite and in the interval")
    expect_error(given_species_params(params)$f0 <- 2,
                 "`f0` must be finite and in the interval")
    expect_error(given_species_params(params)$f0 <- Inf,
                 "`f0` must be finite and in the interval")
    expect_error(given_species_params(params)$f0 <- NaN,
                 "`f0` must be finite and in the interval")

    sp$gamma <- 1
    params <- newMultispeciesParams(sp, info_level = 0)
    expect_true("gamma" %in% names(given_species_params(params)))
    expect_error(given_species_params(params)$f0 <- 1,
                 "`f0` must be finite and in the interval")
    sp_invalid <- data.frame(species = "a", w_inf = 100, gamma = 1, f0 = 1)
    expect_error(given_species_params(sp_invalid),
                 "`f0` must be finite and in the interval")

    params@species_params$f0 <- 1
    expect_error(get_gamma_default(params),
                 "`f0` must be finite and in the interval")
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
    expect_warning(given_species_params(params)$f0 <- 0.5)
    expect_warning(given_species_params(params)$fc <- 1)
    expect_warning(given_species_params(params)$age_mat <- 1)
    expect_warning(given_species_params(params)$catchability <- 2)
    expect_warning(given_species_params(params)$yield_observed <- 1)

    # No warning if NA
    params@given_species_params$gamma[-1] <- NA
    f0 <- c(NA, rep(0.2, no_sp - 1))
    expect_warning(given_species_params(params)$f0 <- f0, NA)

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

test_that("`species_params<-()` removes a column the table drops (#578)", {
    params <- NS_params_small
    species_params(params)$my_col <- 1
    expect_true("my_col" %in% names(given_species_params(params)))

    sp <- species_params(params)
    sp$my_col <- NULL
    expect_message(species_params(params) <- sp,
                   "removed the species parameter column `my_col`",
                   fixed = TRUE)

    # Gone from the model, from the record of what the user gave, and from the
    # values mizer claims to have calculated.
    expect_false("my_col" %in% names(species_params(params)))
    expect_false("my_col" %in% names(given_species_params(params)))
    expect_false("my_col" %in% names(calculated_species_params(params)))
})

test_that("the `$` and `[[` idioms remove a species parameter column (#578)", {
    params <- NS_params_small
    species_params(params)$my_col <- 1
    suppressMessages(species_params(params)$my_col <- NULL)
    expect_false("my_col" %in% names(species_params(params)))
    expect_false("my_col" %in% names(given_species_params(params)))

    species_params(params)$my_col <- 1
    suppressMessages(species_params(params)[["my_col"]] <- NULL)
    expect_false("my_col" %in% names(species_params(params)))
    expect_false("my_col" %in% names(given_species_params(params)))
})

test_that("`given_species_params<-()` removes what mizer cannot calculate (#578)", {
    params <- NS_params_small
    species_params(params)$my_col <- 1
    suppressMessages(given_species_params(params)$my_col <- NULL)
    # Nothing recalculates a column mizer knows nothing about, so leaving it in
    # the species parameters would report the user's own value as one that
    # mizer had produced.
    expect_false("my_col" %in% names(species_params(params)))
    expect_false("my_col" %in% names(calculated_species_params(params)))
})

test_that("withdrawing a calculated parameter hands it back to mizer (#578)", {
    params <- NS_params_small
    sp <- species_params(params)
    sp$gamma <- NULL
    suppressMessages(species_params(params) <- sp)

    # `gamma` is one that mizer knows how to calculate, so it comes straight
    # back -- as a calculated value that follows `f0` again.
    expect_true("gamma" %in% names(species_params(params)))
    expect_false("gamma" %in% names(given_species_params(params)))
    expect_true("gamma" %in% names(calculated_species_params(params)))

    gamma <- species_params(params)$gamma
    suppressMessages(given_species_params(params)$f0 <- 0.1)
    expect_gt(sum(gamma - species_params(params)$gamma), 0)
})

test_that("withdrawing a parameter warns when its rate is frozen (#578)", {
    params <- NS_params_small
    search_vol(params) <- search_vol(params)
    withr::local_options(mizer_info_level = 1)

    expect_warning(given_species_params(params)$gamma <- NULL,
                   paste0("Your change to the species parameter `gamma`.*",
                          "search volume"))
})

test_that("a column is withdrawn without recalculation too (#578)", {
    params <- NS_params_small
    species_params(params)$my_col <- 1
    sp <- species_params(params)
    sp$my_col <- NULL
    metab_before <- params@metab

    suppressMessages(species_params(params, recalculate = FALSE) <- sp)

    expect_false("my_col" %in% names(species_params(params)))
    expect_false("my_col" %in% names(given_species_params(params)))
    # The two tables agree even though nothing was rebuilt.
    expect_equal(params@metab, metab_before, ignore_attr = TRUE)
})

test_that("withdrawing a required species parameter errors (#578)", {
    params <- NS_params_small
    sp <- species_params(params)
    sp$w_inf <- NULL
    sp$w_max <- NULL
    sp$w_repro_max <- NULL
    expect_error(species_params(params) <- sp,
                 "You need to specify the asymptotic size")
})

test_that("the removal of a column is reported at info level 3 (#578)", {
    params <- NS_params_small
    species_params(params)$my_col <- 1
    species_params(params)$my_other_col <- 2
    sp <- species_params(params)
    sp$my_col <- NULL
    sp$my_other_col <- NULL
    expect_message(species_params(params) <- sp,
                   "columns `my_col`, `my_other_col`", fixed = TRUE)

    # Removing a column is what the user asked for, not something that went
    # differently from how they asked, so it is not one of the level-1 reports.
    params <- NS_params_small
    species_params(params)$my_col <- 1
    sp <- species_params(params)
    sp$my_col <- NULL
    withr::local_options(mizer_info_level = 1)
    expect_message(species_params(params) <- sp, NA)
})

test_that("a failed replacement does not report a removal (#578)", {
    params <- NS_params_small
    species_params(params)$my_col <- 1
    sp <- species_params(params)
    sp$my_col <- NULL
    sp$pred_kernel_type <- "definitely_missing"
    messages <- character()

    expect_error(
        withCallingHandlers(
            species_params(params) <- sp,
            message = function(cnd) {
                messages <<- c(messages, conditionMessage(cnd))
                invokeRestart("muffleMessage")
            }),
        "pred_kernel_func is not a function")
    expect_length(messages, 0)
    expect_true("my_col" %in% names(species_params(params)))
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

test_that("species_params setter does not record a default it filled in (#496)", {
    # A model whose species parameters do not carry the length-weight
    # parameters, as `NS_params` does not. `validSpeciesParams()` fills in
    # defaults for them, and those must not be mistaken for user input.
    params <- NS_params_small
    params@species_params$a <- NULL
    params@species_params$b <- NULL
    expect_false(any(c("a", "b") %in% names(given_species_params(params))))

    species_params(params)$beta <- 150

    given <- given_species_params(params)
    # The change the user made is recorded ...
    expect_equal(given$beta, rep(150, nrow(given)), ignore_attr = TRUE)
    # ... but the defaults mizer filled in are not, so that they keep
    # responding to later changes instead of being frozen as given values.
    expect_false(any(c("a", "b") %in% names(given)) &&
                     any(!is.na(given[["a"]])))
    expect_true(all(c("a", "b") %in% names(species_params(params))))
})

test_that("species_params setter still records a column the user adds (#496)", {
    # The counterpart of the test above: a column that is genuinely new to the
    # model, rather than one mizer filled in, is user input and is recorded.
    params <- NS_params_small
    no_sp <- nrow(species_params(params))
    species_params(params)$my_own_param <- seq_len(no_sp)
    expect_equal(given_species_params(params)$my_own_param, seq_len(no_sp),
                 ignore_attr = TRUE)
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

test_that("species_params setter is quiet when a frozen rate blocks the change (#496)", {
    params <- NS_params_small
    # Freeze the metabolic rate at its current value
    metab(params) <- metab(params)
    before <- params@metab

    sp <- species_params(params)
    sp$ks <- sp$ks * 2
    # The diagnostics belong to `given_species_params<-()`; this setter is the
    # quiet one, so that scripts written against earlier versions of mizer keep
    # running clean.
    expect_silent(species_params(params) <- sp)
    # The table records the change but the model does not.
    expect_equal(species_params(params)$ks, sp$ks, ignore_attr = TRUE)
    expect_equal(params@metab, before, ignore_attr = TRUE)
})

test_that("species_params setter is quiet about an overridden parameter (#496)", {
    params <- NS_params_small
    expect_false(any(is.na(given_species_params(params)$gamma)))
    expect_silent(species_params(params)$f0 <- 0.5)
    # while the given species parameter setter says so
    params <- NS_params_small
    expect_warning(given_species_params(params)$f0 <- 0.5,
                   "will not lead to a re-calculation of `gamma`")
})

test_that("given_species_params setter warns when a frozen rate blocks the change (#489)", {
    params <- NS_params_small
    search_vol(params) <- search_vol(params)
    before <- params@search_vol

    gsp <- given_species_params(params)
    gsp$gamma <- gsp$gamma * 2
    expect_warning(given_species_params(params) <- gsp,
                   "Your change to the species parameter `gamma`.*search volume")
    expect_equal(params@search_vol, before, ignore_attr = TRUE)
})

test_that("given_species_params setter ignores a new all-NA column (#524)", {
    params <- NS_params_small
    ext_mort(params) <- ext_mort(params)
    expect_false("z0" %in% names(given_species_params(params)))

    # Nothing acquires a value, so none of the three diagnostics fires.
    expect_silent(given_species_params(params)$z0 <- NA)
    expect_false("catchability" %in% names(given_species_params(params)))
    expect_silent(given_species_params(params)$catchability <- NA)
})

test_that("given_species_params setter treats clearing to NA as a change (#524)", {
    params <- NS_params_small
    search_vol(params) <- search_vol(params)
    before <- params@search_vol
    expect_false(any(is.na(given_species_params(params)$gamma)))

    # Handing `gamma` back to mizer's calculation is an instruction that the
    # frozen search volume cannot carry out, so the user is told.
    expect_warning(given_species_params(params)$gamma <- NA,
                   "Your change to the species parameter `gamma`.*search volume")
    expect_equal(params@search_vol, before, ignore_attr = TRUE)
})

test_that("the given_species_params diagnostics agree about clearing to NA (#524)", {
    params <- NS_params_small
    expect_warning(given_species_params(params)$catchability <- 2,
                   "you should use `gear_params\\(\\)<-`")
    # Clearing it again is just as much a change, so it is reported again.
    expect_warning(given_species_params(params)$catchability <- NA,
                   "you should use `gear_params\\(\\)<-`")

    # Only a value that is there can be overruled by `gamma`, so clearing `f0`
    # is a change that `signal_ignored_changes()` has nothing to say about.
    params <- NS_params_small
    expect_false(any(is.na(given_species_params(params)$gamma)))
    expect_false(any(is.na(given_species_params(params)$f0)))
    expect_warning(given_species_params(params)$f0 <- 0.5,
                   "will not lead to a re-calculation of `gamma`")
    expect_silent(given_species_params(params)$f0 <- NA)
})

test_that("species_params setter is quiet when no frozen rate is in the way", {
    params <- NS_params_small
    sp <- species_params(params)
    sp$ks <- sp$ks * 2
    expect_silent(species_params(params) <- sp)
    expect_false(isTRUE(all.equal(params@metab, metab(NS_params_small),
                                  check.attributes = FALSE)))
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

test_that("a column restored during validation is not reported as removed", {
    params <- length_based_params()
    # With `a` and `b` explicitly given, validating the given length also puts
    # its matching weight back into the given table.
    given_species_params(params)$a <- species_params(params)$a
    given_species_params(params)$b <- species_params(params)$b

    expect_message(species_params(params)$w_mat <- NULL, NA)
    expect_true("w_mat" %in% names(given_species_params(params)))
})

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

test_that("the two setters apply the length/weight rule the same way", {
    # The documentation promises that the only difference between the two
    # setters is the warnings, so a change made through either has to arrive at
    # the same model.
    params <- length_based_params()

    edits <- list(
        # only the length changes, so it determines the weight
        l_mat = function(sp) {
            sp$l_mat <- sp$l_mat * 1.2
            sp
        },
        # only the weight changes, so the length follows it
        w_mat = function(sp) {
            sp$w_mat <- sp$w_mat * 0.9
            sp
        },
        # both change at once, so the weight wins
        both = function(sp) {
            sp$l_mat <- sp$l_mat * 1.2
            sp$w_mat <- sp$w_mat * 1.1
            sp
        },
        # neither changes but the conversion does, so the weight wins again
        a = function(sp) {
            sp$a <- rep(0.02, nrow(sp))
            sp
        }
    )

    for (edit in names(edits)) {
        via_sp <- params
        suppressWarnings(species_params(via_sp) <-
                             edits[[edit]](species_params(via_sp)))
        via_given <- params
        suppressWarnings(given_species_params(via_given) <-
                             edits[[edit]](given_species_params(via_given)))
        for (par in c("w_mat", "l_mat", "h", "gamma", "ks", "w_mat25")) {
            expect_equal(species_params(via_sp)[[par]],
                         species_params(via_given)[[par]],
                         ignore_attr = TRUE, info = paste(edit, par))
        }
        expect_equal(via_sp@search_vol, via_given@search_vol,
                     ignore_attr = TRUE, info = edit)
        expect_equal(via_sp@psi, via_given@psi,
                     ignore_attr = TRUE, info = edit)
    }
})

test_that("given_species_params setter keeps length and weight consistent", {
    # The given species parameters do not hold `a` and `b`, but the rule still
    # has to be applied to them, or the model would report the same
    # inconsistency on every later change.
    params <- length_based_params()
    expect_false(any(c("a", "b") %in% names(given_species_params(params))))

    given <- given_species_params(params)
    given$l_mat <- given$l_mat * 1.2
    expect_warning(given_species_params(params) <- given, NA)
    expect_equal(given_species_params(params)$w_mat,
                 l2w(given$l_mat, species_params(params)),
                 ignore_attr = TRUE)
    # `a` and `b` were not given and must not have been recorded as given
    expect_false(any(c("a", "b") %in% names(given_species_params(params))))

    # A later, unrelated change finds nothing left to complain about
    given <- given_species_params(params)
    given$beta <- given$beta * 1.1
    expect_warning(given_species_params(params) <- given, NA)
})

test_that("given_species_params setter preserves columns mizer does not calculate", {
    # Columns that live only in the `@species_params` slot are not rebuilt from
    # the given species parameters, so they have to be carried over, just as
    # `species_params<-()` carries them over.
    params <- NS_params_small
    params@species_params$my_own_param <- seq_len(nrow(species_params(params)))
    given <- given_species_params(params)
    given$beta <- given$beta * 1.1
    given_species_params(params) <- given
    expect_equal(species_params(params)$my_own_param,
                 seq_len(nrow(species_params(params))), ignore_attr = TRUE)
})

test_that("given_species_params can protect an unchanged calculated value", {
    params <- NS_params_small
    params@given_species_params$q <- NULL
    current_q <- species_params(params)$q

    # Make a deliberately inconsistent, unfrozen array. A rebuild would repair
    # it, so retaining it demonstrates that this provenance-only change does
    # not call setParams().
    params@search_vol[] <- 2 * params@search_vol
    search_vol_before <- params@search_vol

    given_species_params(params)$q <- current_q

    expect_equal(given_species_params(params)$q, current_q,
                 ignore_attr = TRUE)
    expect_false("q" %in% names(calculated_species_params(params)))
    expect_equal(params@search_vol, search_vol_before, ignore_attr = TRUE)

    # The promotion has an effect on future recalculations: q no longer follows
    # a change in the resource exponent.
    resource_params(params)$lambda <- resource_params(params)$lambda + 0.1
    expect_equal(species_params(params)$q, current_q, ignore_attr = TRUE)
})

test_that("the full species_params table can be declared given without rebuilding", {
    params <- NS_params_small
    params@search_vol[] <- 2 * params@search_vol
    search_vol_before <- params@search_vol
    sp <- species_params(params)

    given_species_params(params) <- sp

    given <- given_species_params(params)
    for (col in names(sp)) {
        if (is.atomic(sp[[col]]) && is.null(dim(sp[[col]]))) {
            supplied <- !is.na(sp[[col]])
            expect_equal(given[[col]][supplied], sp[[col]][supplied],
                         ignore_attr = TRUE, info = col)
        }
    }
    expect_equal(params@search_vol, search_vol_before, ignore_attr = TRUE)
})

test_that("leaf and custom species parameters do not rebuild the model", {
    params <- NS_params_small
    params@search_vol[] <- 2 * params@search_vol
    search_vol_before <- params@search_vol
    observed <- seq_len(nrow(species_params(params))) * 100

    species_params(params)$biomass_observed <- observed
    expect_equal(species_params(params)$biomass_observed, observed,
                 ignore_attr = TRUE)
    expect_equal(params@search_vol, search_vol_before, ignore_attr = TRUE)

    species_params(params)$my_unused_parameter <- observed / 10
    expect_equal(species_params(params)$my_unused_parameter, observed / 10,
                 ignore_attr = TRUE)
    expect_equal(params@search_vol, search_vol_before, ignore_attr = TRUE)

    # The authoritative setter also merges a changed leaf value into the full
    # table without rebuilding unrelated arrays.
    given_species_params(params)$biomass_observed <- observed * 2
    expect_equal(species_params(params)$biomass_observed, observed * 2,
                 ignore_attr = TRUE)
    expect_equal(params@search_vol, search_vol_before, ignore_attr = TRUE)
})

test_that("dependent changes and demotions still rebuild the model", {
    params <- NS_params_small
    params@given_species_params$q <- NULL
    params@search_vol[] <- 2 * params@search_vol
    search_vol_before <- params@search_vol

    species_params(params)$q <- species_params(params)$q + 0.1
    expect_false(isTRUE(all.equal(params@search_vol, search_vol_before,
                                  check.attributes = FALSE)))

    # A demotion must derive the replacement even when the old value happens
    # to be the current species-parameter value.
    params <- NS_params_small
    old_gamma <- species_params(params)$gamma
    given <- given_species_params(params)
    given$f0 <- rep(0.2, nrow(given))
    given$gamma <- NA
    suppressWarnings(given_species_params(params) <- given)
    expect_true(any(species_params(params)$gamma != old_gamma))
})

test_that("custom predation-kernel arguments require recalculation", {
    fun_name <- "dependency_test_pred_kernel"
    assign(fun_name, function(ppmr, custom_shape) ppmr^custom_shape,
           envir = .GlobalEnv)
    on.exit(rm(list = fun_name, envir = .GlobalEnv), add = TRUE)

    params <- NS_params_small
    sp <- species_params(params)
    sp$pred_kernel_type <- "dependency_test"
    sp$custom_shape <- rep(-1, nrow(sp))
    species_params(params) <- sp

    required <- recalculation_species_params(params, species_params(params))
    expect_in("custom_shape", required)
})

test_that("extension objects treat unknown species parameters conservatively", {
    params <- NS_params_small
    class(params) <- c("SpeciesParamsDependencyTestParams", "MizerParams")
    sp <- species_params(params)
    sp$extension_parameter <- seq_len(nrow(sp))

    required <- recalculation_species_params(params, sp)
    expect_in("extension_parameter", required)
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

    # Error message when h = Inf
    params_inf <- params
    params_inf@species_params$gamma[] <- NA
    params_inf@species_params$h[1] <- Inf
    expect_error(get_gamma_default(params_inf),
                 "Cannot calculate default `gamma` for species with `h = Inf`")
    params_inf2 <- params
    params_inf2@species_params$ks[] <- NA
    params_inf2@species_params$h[1] <- Inf
    expect_error(get_ks_default(params_inf2),
                 "Cannot calculate default `ks` for species with `h = Inf`")

    # get_f0_default under edition 2 matches get_gamma_default with interaction_resource != 1
    params_ed2 <- params
    species_params(params_ed2)$interaction_resource <- 0.5
    params_ed2@species_params$gamma[] <- NA
    defaults_edition(2)
    on.exit(defaults_edition(1), add = TRUE)
    gamma_ed2 <- get_gamma_default(params_ed2)
    params_ed2@species_params$gamma <- gamma_ed2
    f0_ed2 <- get_f0_default(params_ed2)
    expect_equal(unname(f0_ed2), rep(0.6, nrow(species_params(params_ed2))))

    # get_h_default informs when n != p
    sp_np <- species_params(params)
    sp_np$h[] <- NA
    sp_np$n[1] <- 0.7
    sp_np$p[1] <- 0.8
    expect_message(with_info_level(get_h_default(sp_np), info_level = 1),
                   "Because you have n != p")
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
    expect_warning(species_params(df),
                   "very close to standard parameter names")
    # Fuzzy typo of a recognised name is flagged with a suggestion ...
    df2 <- data.frame(species = c("Sprat", "Herring"), w_inf = c(10, 100),
                      w_maxx = 100)
    expect_warning(species_params(df2), "did you mean `w_max`")
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

test_that("an invalid w_mat25 really is replaced by its default (#580)", {
    params <- NS_params_small
    warnings <- capture_warnings({
        species_params(params)$w_mat25 <- species_params(params)$w_mat
    })
    expect_length(warnings, 1)
    expect_match(warnings,
                 "marking it as missing so that its default will be used",
                 fixed = TRUE)
    # The message promises the default, so the stored value has to be it.
    expect_false(any(is.na(species_params(params)$w_mat25)))
    expect_equal(species_params(params)$w_mat25,
                 species_params(params)$w_mat / (3 ^ (1 / 10)),
                 ignore_attr = TRUE)
})

test_that("an invalid w_mat25 is not restored from its length (#580)", {
    # With `l_mat25` present, clearing `w_mat25` alone is not enough: the
    # length/weight conversion fills the missing weight straight back in from
    # the length, the rejected value returns to within rounding of `w_mat`, and
    # the maturity ogive degenerates into a step function.
    sp <- data.frame(species = c("A", "B"), w_max = c(100, 1000),
                     w_mat = c(10, 100), a = 0.01, b = 3)
    sp$l_mat25 <- (sp$w_mat / sp$a) ^ (1 / sp$b)
    sp$w_mat25 <- sp$w_mat
    valid <- suppressWarnings(given_species_params(sp))
    expect_true(all(is.na(valid$w_mat25)))
    expect_true(all(is.na(valid$l_mat25)))

    params <- suppressMessages(suppressWarnings(
        newMultispeciesParams(sp, no_w = 20)))
    expect_equal(species_params(params)$w_mat25,
                 species_params(params)$w_mat / (3 ^ (1 / 10)),
                 ignore_attr = TRUE)
    # A genuine sigmoid, not the two-valued step function that an almost-equal
    # w_mat25 produces.
    expect_gt(length(unique(as.vector(params@maturity))), 2)
})

test_that("a misspelled column is reported once, when it appears (#581)", {
    params <- NS_params_small
    warnings <- capture_warnings(species_params(params)$wmat <- 100)
    expect_length(warnings, 1)
    expect_match(warnings, "did you mean `w_mat`", fixed = TRUE)

    # `params` still carries the `wmat` column from above; changing an
    # unrelated parameter must not report it again. That repetition is the
    # second half of #581.
    expect_no_warning(species_params(params)$b <- 3.1)

    # The same for the other setter, which keeps its own baseline.
    params <- NS_params_small
    warnings <- capture_warnings(given_species_params(params)$wmat <- 100)
    expect_length(warnings, 1)
    expect_match(warnings, "did you mean `w_mat`", fixed = TRUE)
    expect_no_warning(given_species_params(params)$b <- 3.1)
})

test_that("the misspelling check is not skipped for a classed table", {
    # A table that has already been through `species_params()` keeps its class,
    # and a typo added afterwards still has to be reported.
    sp <- species_params(data.frame(species = c("A", "B"),
                                    w_max = c(100, 1000)))
    sp$wmat <- 1
    expect_true(is.species_params(sp))
    expect_warning(species_params(sp), "did you mean `w_mat`")
    expect_warning(validSpeciesParams(sp), "did you mean `w_mat`")
    expect_warning(suppressMessages(newMultispeciesParams(sp, no_w = 20)),
                   "did you mean `w_mat`")
    # ... and mizer can still ask for it to be skipped.
    expect_no_warning(species_params(sp, check_misspellings = FALSE))
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

test_that("$ on species_params does not partially match column names", {
    # Without length-weight parameters, `$a` used to return `alpha` and `$b`
    # used to return `beta`, with the species names attached
    sp <- species_params(NS_params_small)
    sp <- sp[, setdiff(names(sp), c("a", "b"))]
    expect_s3_class(sp, "species_params")
    expect_true(all(c("alpha", "beta") %in% names(sp)))

    expect_warning(a <- sp$a, "There is no `a` column")
    expect_null(a)
    expect_warning(b <- sp$b, "There is no `b` column")
    expect_null(b)
    # The warning points at the column that used to be returned
    expect_warning(sp$a, "`alpha` column")

    # A name matching nothing at all stays silent, so that testing for the
    # presence of a parameter remains an option
    expect_silent(out <- sp$no_such_parameter)
    expect_null(out)

    # An ambiguous prefix was already NULL before, so it stays silent too
    expect_silent(out <- sp$w_)
    expect_null(out)

    # Exact matches are unaffected
    expect_identical(unname(sp$alpha), sp[["alpha"]])

    # Reads and writes now agree about which column `$b` means
    sp$b <- 3
    expect_identical(unname(sp$b), rep(3, nrow(sp)))

    # given_species_params inherits the behaviour
    given <- given_species_params(NS_params_small)
    given <- given[, setdiff(names(given), c("a", "b"))]
    expect_s3_class(given, "given_species_params")
    expect_true("beta" %in% names(given))
    expect_warning(expect_null(given$b), "There is no `b` column")
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

# changed_entries ----------------------------------------------------------

test_that("changed_entries compares NA as a value", {
    expect_identical(changed_entries(c(1, NA, 3), c(1, NA, 4), 3),
                     c(FALSE, FALSE, TRUE))
    # A value cleared to NA is a change, an NA that stays NA is not.
    expect_identical(changed_entries(c(NA, NA), c(1, NA), 2), c(TRUE, FALSE))
    # As is a value where there was none.
    expect_identical(changed_entries(c(1, NA), c(NA, NA), 2), c(TRUE, FALSE))
})

test_that("changed_entries treats a column that did not exist as all NA", {
    expect_identical(changed_entries(c(1, 2), NULL, 2), c(TRUE, TRUE))
    expect_identical(changed_entries(c(1, NA), NULL, 2), c(TRUE, FALSE))
    expect_identical(changed_entries(c(NA, NA), NULL, 2), c(FALSE, FALSE))
})

test_that("changed_entries handles list and matrix columns", {
    expect_identical(changed_entries(list(c("a", "b"), "c"),
                                     list(c("a", "b"), "d"), 2),
                     c(FALSE, TRUE))
    expect_identical(changed_entries(matrix(1:4, nrow = 2),
                                     matrix(c(1L, 9L, 3L, 4L), nrow = 2), 2),
                     c(FALSE, TRUE))
    # A list column is new, so every species that has an entry has changed
    expect_identical(changed_entries(list("a", "b"), NULL, 2), c(TRUE, TRUE))
})

test_that("changed_entries says nothing about a column of the wrong length", {
    expect_identical(changed_entries(1, c(1, 2), 2), c(FALSE, FALSE))
    expect_identical(changed_entries(NULL, c(1, 2), 2), c(FALSE, FALSE))
    # Nothing to compare against, so everything counts as changed
    expect_identical(changed_entries(c(1, 2), c(1, 2, 3), 2), c(TRUE, TRUE))
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


# reconcileSpeciesParams --------------------------------------------------

test_that("reconcileSpeciesParams() leaves a consistent model alone", {
    params <- NS_params_small
    expect_identical(reconcileSpeciesParams(params), params)
})

test_that("reconcileSpeciesParams() records only the entries a recalculation would change", {
    params <- NS_params_small
    # `w_mat25` is calculated by `setReproduction()` and is not among the given
    # species parameters of this model, so a value written straight into the
    # slot would be undone by the next recalculation.
    new_value <- unname(species_params(params)$w_mat25[2]) / 2
    params@species_params$w_mat25[2] <- new_value

    reconciled <- reconcileSpeciesParams(params, info_level = 0)

    given <- given_species_params(reconciled)
    expect_equal(unname(given$w_mat25[2]), new_value)
    # The species whose value a recalculation reproduces stay calculated
    expect_true(all(is.na(given$w_mat25[-2])))
})

test_that("reconcileSpeciesParams() changes neither the species parameters nor the rates", {
    params <- NS_params_small
    params@species_params$w_mat25[2] <- species_params(params)$w_mat25[2] / 2

    reconciled <- reconcileSpeciesParams(params, info_level = 0)

    expect_identical(species_params(reconciled), species_params(params))
    expect_identical(reconciled@psi, params@psi)
    expect_identical(reconciled@maturity, params@maturity)
    expect_identical(reconciled@search_vol, params@search_vol)
})

test_that("reconcileSpeciesParams() protects the value against a recalculation", {
    params <- NS_params_small
    new_value <- unname(species_params(params)$w_mat25[2]) / 2
    params@species_params$w_mat25[2] <- new_value

    # Without the reconciliation the next recalculation undoes the change
    unprotected <- params
    species_params(unprotected)$w_mat <- species_params(unprotected)$w_mat
    expect_false(isTRUE(all.equal(species_params(unprotected)$w_mat25[2],
                                  new_value)))

    protected <- reconcileSpeciesParams(params, info_level = 0)
    species_params(protected)$w_mat <- species_params(protected)$w_mat
    expect_equal(unname(species_params(protected)$w_mat25[2]), new_value)
})

test_that("reconcileSpeciesParams() reports what it recorded", {
    params <- NS_params_small
    params@species_params$w_mat25[2] <- species_params(params)$w_mat25[2] / 2

    expect_message(reconcileSpeciesParams(params),
                   "`w_mat25` holds a value that a recalculation would not reproduce")
    expect_silent(reconcileSpeciesParams(params, info_level = 0))
    # A consistent model gives no report at all
    expect_silent(reconcileSpeciesParams(NS_params_small))
})

test_that("reconcileSpeciesParams() keeps the class of the given species parameters", {
    params <- NS_params_small
    params@species_params$w_mat25[2] <- species_params(params)$w_mat25[2] / 2

    given <- given_species_params(reconcileSpeciesParams(params,
                                                         info_level = 0))
    expect_s3_class(given, "given_species_params")
    expect_s3_class(given, "species_params")
    expect_identical(rownames(given),
                     rownames(given_species_params(params)))
    # `a` and `b` are borrowed while the length/weight rules are applied and
    # must not be left behind as given values.
    expect_false(any(c("a", "b") %in% names(given)))
})

test_that("reconcileSpeciesParams() overwrites an out-of-date given value", {
    params <- NS_params_small
    # `gamma` is among the given species parameters of this model, so the entry
    # is already there and only its value has to follow.
    new_value <- unname(species_params(params)$gamma[1]) * 2
    params@species_params$gamma[1] <- new_value

    reconciled <- reconcileSpeciesParams(params, info_level = 0)

    expect_equal(unname(given_species_params(reconciled)$gamma[1]), new_value)
    expect_equal(unname(given_species_params(reconciled)$gamma[-1]),
                 unname(given_species_params(params)$gamma[-1]))
})

test_that("reconcileSpeciesParams() makes the species parameters a fixed point", {
    # `h` is not among the given species parameters of a model built this way,
    # and mizer derives `gamma` and `ks` from it, so a single pass of recording
    # would leave those two to move at the next recalculation.
    plain <- suppressMessages(
        newMultispeciesParams(NS_species_params_small, no_w = 20,
                              info_level = 0))
    edited <- plain
    edited@species_params$h[2] <- 2 * species_params(plain)$h[2]

    # Recalculating the unrepaired model moves `h` back
    expect_true(length(recalculation_moves(edited)) > 0)

    reconciled <- reconcileSpeciesParams(edited, info_level = 0)

    # The species parameters are untouched ...
    expect_identical(species_params(reconciled), species_params(edited))
    # ... and no amount of recalculation moves them
    expect_identical(recalculation_moves(reconciled, times = 3), character())
    # The parameters derived from `h` had to be recorded as well to get there
    given <- given_species_params(reconciled)
    expect_false(is.na(given$h[2]))
    expect_false(is.na(given$gamma[2]))
    expect_false(is.na(given$ks[2]))
})

test_that("reconcileSpeciesParams() is idempotent", {
    plain <- suppressMessages(
        newMultispeciesParams(NS_species_params_small, no_w = 20,
                              info_level = 0))
    plain@species_params$h[2] <- 2 * species_params(plain)$h[2]

    once <- reconcileSpeciesParams(plain, info_level = 0)
    expect_identical(reconcileSpeciesParams(once, info_level = 0), once)
    # A second call has nothing to report either
    expect_silent(reconcileSpeciesParams(once))
})
