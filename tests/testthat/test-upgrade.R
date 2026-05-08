params <- NS_params
sim <- project(params, t_max = 0.1, t_save = 0.1)

## upgradeParams ----
test_that("upgradeParams leaves new params unchanged", {
    expect_identical(upgradeParams(params), params)
})
test_that("upgradeParams preserves comments", {
    comment(params) <- "test"
    for (slot in (slotNames(params))) {
        comment(slot(params, slot)) <- slot
    }
    expect_identical(upgradeParams(params), params)
})

test_that("upgradeParams from pre-2.4 preserves all slot values", {
    params <- NS_params
    # Version "2.4.0" triggers the full newMultispeciesParams() rebuild path
    params@mizer_version <- "2.4.0"

    # Change non-numeric slots to non-default values so they can be checked
    params@resource_dynamics <- "resource_constant"
    params@linecolour[] <- paste0("#", formatC(seq_along(params@linecolour) * 11111,
                                               width = 6, flag = "0"))
    params@linetype[] <- rep("dashed", length(params@linetype))
    params@use_predation_diffusion <- TRUE
    params@metadata[["test_key"]] <- "test_value"
    params@resource_params[["test_key"]] <- 42
    params@rates_funcs[["TestFunc"]] <- "identity"
    params@initial_n_other[["test_comp"]] <- 42
    params@other_dynamics[["test_comp"]] <- "identity"
    params@other_params[["test_comp"]] <- list(x = 1.23)
    params@other_encounter[["test_comp"]] <- "identity"
    params@other_mort[["test_comp"]] <- "identity"
    params@given_species_params[["test_col"]] <-
        seq_len(nrow(params@given_species_params))

    # Randomise all numeric array and vector slots while preserving
    # dims, dimnames, and comments. We keep species_params, given_species_params,
    # and gear_params at their NS_params values so the internal
    # newMultispeciesParams() call succeeds.
    # w, dw, w_full, dw_full define the grid and must stay consistent.
    # w_min_idx is legitimately recalculated by validParams() from species_params.
    set.seed(42)
    skip_randomise <- c("w", "dw", "w_full", "dw_full", "w_min_idx",
                        "maturity", "psi")
    for (s in slotNames(params)) {
        if (s %in% skip_randomise) next
        sl <- slot(params, s)
        if (is.numeric(sl) && length(sl) > 0) {
            cmt <- comment(sl)
            slot(params, s)[] <- runif(length(sl), 0.1, 2)
            comment(slot(params, s)) <- cmt
        }
    }
    # maturity must stay in [0, 1]; psi must not exceed maturity
    params@maturity[] <- runif(length(params@maturity), 0.1, 0.9)
    params@psi[] <- params@maturity * runif(length(params@psi), 0, 0.9)

    params_before <- params
    params_after <- suppressMessages(upgradeParams(params))

    # species_params may gain additional derived columns from newMultispeciesParams;
    # check only that the original columns are preserved unchanged.
    original_cols <- names(params_before@species_params)
    expect_identical(params_after@species_params[, original_cols],
                     params_before@species_params,
                     label = "params@species_params (original columns)")

    slots_to_check <- setdiff(slotNames(params),
                               c("mizer_version", "time_modified", "species_params"))
    for (s in slots_to_check) {
        expect_identical(slot(params_after, s), slot(params_before, s),
                         label = paste0("params@", s))
    }
})

## upgradeSim ----
test_that("upgradeSim leaves new sim unchanged", {
    expect_identical(upgradeSim(sim), sim)
})
test_that("upgradeSim preserves comments", {
    comment(sim) <- "test"
    for (slot in (slotNames(sim))) {
        comment(slot(sim, slot)) <- slot
    }
    expect_identical(upgradeSim(sim), sim)
})

test_that("Object from version 0.4 can be upgraded", {
    simc.0.4 <- readRDS("assets/simc.0.4.rds")
    (sim <- validSim(simc.0.4)) |>
        expect_message("Initial effort has been set to 0") |>
        expect_warning("Your MizerSim object was created with an earlier")
    
    expect_true(validObject(sim))
})
test_that("Object from version 1.0 can be upgraded", {
    simc.1.0 <- readRDS("assets/simc.1.0.rds")
    (sim <- validSim(simc.1.0)) |>
        expect_message("Initial effort has been set to 0") |>
        expect_warning("Your MizerSim object was created with an earlier")
    expect_true(validObject(sim))
})
test_that("r_max is renamed", {
    params <- NS_params
    params@species_params$r_max <- params@species_params$R_max
    params@mizer_version <- "1.1"
    expect_message(upgradeParams(params),
                   "The 'r_max' column has been renamed to 'R_max'.")
})

test_that("Some functions work with params from earlier versions", {
    params.0.4 <- readRDS("assets/simc.0.4.rds")@params
    expect_warning(getEGrowth(params.0.4) |> expect_message(),
                   "Your MizerParams object was created with an earlier")
    expect_warning(plotFeedingLevel(params.0.4) |> expect_message(),
                   "Your MizerParams object was created with an earlier")
    expect_warning(project(params.0.4, t_max = 0.1) |> expect_message(),
                   "Your MizerParams object was created with an earlier")
    # renaming of resource dynamics functions
    slot(params.0.4, "srr", check = FALSE) <- "srrNone"
    (p4 <- suppressWarnings(upgradeParams(params.0.4))) |> expect_message()
    expect_identical(p4@rates_funcs$RDD, "noRDD")
})

test_that("upgradeParams updates `time_modified`", {
    p <- NS_params
    p@mizer_version <- "2.0.0"
    p2 <- suppressMessages(upgradeParams(p))
    expect_false(identical(p2@time_modified, p@time_modified))
})

test_that("upgradeParams preserves subclass information", {
    if (!methods::isClass("TestMizerParams")) {
        methods::setClass(
            "TestMizerParams",
            contains = "MizerParams",
            slots = c(extra = "character"),
            prototype = list(extra = "default")
        )
    }
    p <- as(NS_params, "TestMizerParams")
    p@extra <- "custom"
    p@mizer_version <- "2.0.0"

    p2 <- suppressMessages(upgradeParams(p))

    expect_s4_class(p2, "TestMizerParams")
    expect_identical(p2@extra, "custom")
    expect_true(validObject(p2))
})
