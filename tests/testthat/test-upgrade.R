params <- NS_params_small
sim <- project(params, t_max = 0.1, t_save = 0.1)

## upgradeParams ----
test_that("upgradeParams leaves new params unchanged", {
    expect_identical(upgradeParams(params), params)
})

test_that("legacy component registrations become stored specifications", {
    p <- setComponent(NS_params_small, "pool", initial_value = 42,
                      dynamics_fun = "constant_other",
                      component_params = list(capacity = 80),
                      colour = "red", linetype = "dashed")
    p$capabilities <- NULL
    expected <- legacy_capability_registry(p)
    for (version in c("3.4.0.9003", "2.4.0")) {
        old <- p
        old$mizer_version <- version
        up <- suppressMessages(upgradeParams(old))
        expect_identical(up$capabilities, expected)
        expect_identical(up$other_dynamics, p$other_dynamics)
        expect_identical(up$initial_n_other, p$initial_n_other)
        expect_identical(up$other_params, p$other_params)
        expect_identical(upgradeParams(up), up)
        expect_named(up$capabilities$processes, c("consumers", "resource",
                                                  "reproduction"))
        expect_identical(up$capabilities$components$pool$role, "dynamic")
    }
})

test_that("old-model reconstruction preserves an assembled specification", {
    p <- setComponent(NS_params_small, "pool", initial_value = 42,
                      dynamics_fun = "constant_other", colour = "red")
    stored <- p$capabilities
    expect_identical(stored$components$pool$dynamics_fun, "constant_other")
    p$mizer_version <- "2.4.0"
    up <- suppressMessages(upgradeParams(p))
    expect_identical(up$capabilities, stored)
})

test_that("validSim upgrades its embedded capability registry", {
    old <- NS_sim_small
    old$params$capabilities <- NULL
    old$params$mizer_version <- "3.4.0.9003"
    up <- suppressWarnings(validSim(old))
    expect_identical(up$params$capabilities$schema_version, 7L)
    expect_length(up$params$capabilities$components, 0L)
    expect_identical(up$params$capabilities$processes,
                     default_process_registry(up$params))
    expect_equal(up$n, old$n)
    expect_equal(up$n_pp, old$n_pp)
    expect_identical(upgradeSim(up), up)
})
test_that("upgradeParams preserves comments", {
    comment(params) <- "test"
    for (slot in (slotNames(params))) {
        comment(slot(params, slot)) <- slot
    }
    expect_identical(upgradeParams(params), params)
})

test_that("upgradeParams derives w_inf for objects that predate it", {
    params <- NS_params_small
    # Simulate an object from before w_inf became the required parameter
    params@mizer_version <- "3.0.0.9002"
    w_repro_max <- params@species_params$w_repro_max
    params@species_params$w_inf <- NULL
    params@given_species_params$w_inf <- NULL
    expect_false("w_inf" %in% names(params@species_params))

    up <- suppressWarnings(upgradeParams(params))
    # w_repro_max is present, so w_inf is derived from it (preferred over w_max)
    expect_true("w_inf" %in% names(up@species_params))
    expect_equal(up@species_params$w_inf, w_repro_max)

    # When w_repro_max is also absent, w_inf falls back to w_max
    params2 <- params
    params2@species_params$w_repro_max <- NULL
    params2@given_species_params$w_repro_max <- NULL
    up2 <- suppressWarnings(upgradeParams(params2))
    expect_equal(up2@species_params$w_inf, up2@species_params$w_max)
})

test_that("upgradeParams from pre-2.4 preserves all slot values", {
    params <- NS_params_small
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
    # and gear_params at their NS_params_small values so the internal
    # newMultispeciesParams() call succeeds.
    # w, dw, w_full, dw_full define the grid and must stay consistent.
    # w_min_idx is legitimately recalculated by validParams() from species_params.
    set.seed(42)
    skip_randomise <- c("w", "dw", "w_full", "dw_full", "w_min_idx",
                        "maturity", "psi", "defaults_edition")
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
    skip_on_os("mac")
    for (s in slots_to_check) {
        expect_identical(slot(params_after, s), slot(params_before, s),
                         label = paste0("params@", s))
    }
})

test_that("upgradeParams marks an object without an edition as edition 1", {
    params <- NS_params_small
    # Simulate an object from before the edition was recorded
    params@mizer_version <- "3.4.0.9000"
    params$defaults_edition <- NULL
    # The active edition plays no role
    withr::local_options(mizer_defaults_edition = 2)
    up <- suppressWarnings(upgradeParams(params))
    expect_identical(up@defaults_edition, 1)
    expect_warning(up <- validParams(params),
                   "created with an earlier version of mizer")
    expect_identical(up@defaults_edition, 1)
})

test_that("upgradeParams keeps the recorded edition", {
    params <- NS_params_small
    defaults_edition(params) <- 2
    # Version "2.4.0" triggers the full newMultispeciesParams() rebuild path
    params@mizer_version <- "2.4.0"
    withr::local_options(mizer_defaults_edition = 1)
    up <- suppressMessages(suppressWarnings(upgradeParams(params)))
    expect_identical(up@defaults_edition, 2)
})

test_that("upgradeParams keeps the circular FFT kernels of earlier versions (#618)", {
    # Objects from earlier versions are of edition 1, which is not zero-padded,
    # so their kernels and rates must come through the upgrade unchanged even
    # when edition 3 is active
    params <- kernel_below_grid_params(edition = 1)
    # The rates carry the model, whose time_modified the upgrade changes
    drop_params <- function(x) { attr(x, "params") <- NULL; x }
    encounter <- drop_params(getEncounter(params))
    pred_rate <- drop_params(getPredRate(params))
    params@mizer_version <- "3.4.0"
    params$defaults_edition <- NULL
    withr::local_options(mizer_defaults_edition = 3)
    up <- suppressWarnings(upgradeParams(params))
    expect_identical(up@defaults_edition, 1)
    for (s in c("ft_pred_kernel_e", "ft_pred_kernel_p", "ft_pred_kernel_d")) {
        expect_identical(dim(slot(up, s)), c(2L, length(params@w_full)))
    }
    expect_equal(drop_params(getEncounter(up)), encounter, tolerance = 0)
    expect_equal(drop_params(getPredRate(up)), pred_rate, tolerance = 0)
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
    params <- NS_params_small
    params@species_params$r_max <- params@species_params$R_max
    params@mizer_version <- "1.1"
    expect_message(upgradeParams(params),
                   "The 'r_max' column has been renamed to 'R_max'.")
})

test_that("Some functions work with params from earlier versions", {
    params.0.4 <- readRDS("assets/simc.0.4.rds")@params
    expect_warning(getEGrowth(params.0.4) |> expect_message(),
                   "Your MizerParams object was created with an earlier") |>
        suppressWarnings() |> suppressMessages()
    expect_warning(plotFeedingLevel(params.0.4) |> expect_message(),
                   "Your MizerParams object was created with an earlier") |>
        suppressWarnings()
    expect_warning(project(params.0.4, t_max = 0.1) |> expect_message(),
                   "Your MizerParams object was created with an earlier") |>
        suppressWarnings()
    # renaming of resource dynamics functions
    slot(params.0.4, "srr", check = FALSE) <- "srrNone"
    (p4 <- suppressWarnings(upgradeParams(params.0.4))) |> expect_message()
    expect_identical(p4@rates_funcs$RDD, "noRDD")
})

test_that("upgradeParams updates `time_modified`", {
    p <- NS_params_small
    p@mizer_version <- "2.0.0"
    p2 <- suppressMessages(upgradeParams(p))
    expect_false(identical(p2@time_modified, p@time_modified))
})

test_that("upgradeParams preserves subclass information", {
    p <- NS_params_small
    class(p) <- c("TestMizerParams", "MizerParams")
    p$extra <- "custom"
    p$mizer_version <- "2.0.0"

    p2 <- suppressMessages(upgradeParams(p))

    expect_s3_class(p2, "TestMizerParams")
    expect_true(inherits(p2, "MizerParams"))
    expect_identical(p2$extra, "custom")
    expect_true(validObject(p2))
})

test_that("schema-one definitions upgrade without live resolution", {
    p <- NS_params_small
    old <- list(dynamics_fun = "unloaded_legacy_step",
                encounter_fun = NULL, mort_fun = NULL,
                colour = "orange", linetype = "solid")
    p$capabilities <- list(schema_version = 1L, components = list(pool = old),
                           processes = list())
    p$other_dynamics <- list(pool = "unloaded_legacy_step")
    p$initial_n_other <- list(pool = 1)
    p$mizer_version <- "3.4.0.9004"
    upgraded <- upgradeParams(p)
    expect_identical(upgraded$capabilities$schema_version, 7L)
    expect_identical(upgraded$capabilities$components$pool$dynamics_fun,
                     "unloaded_legacy_step")
    expect_null(upgraded$capabilities$components$pool$rate_fun)
    expect_null(upgraded$capabilities$components$pool$integration_order)
    expect_null(upgraded$capabilities$components$pool$convergence_fun)
    expect_identical(upgraded$capabilities$components$pool$role, "dynamic")
    expect_named(upgraded$capabilities$processes, c("consumers", "resource",
                                                    "reproduction"))
    expect_identical(upgradeParams(upgraded), upgraded)
})

test_that("schema-two upgrades distinguish known and unspecified fixed states", {
    p <- NS_params_small
    p$resource_dynamics <- "resource_constant"
    p$other_dynamics <- list(pool = "constant_other")
    p$initial_n_other <- list(pool = 1)
    p$capabilities <- list(
        schema_version = 2L,
        components = list(pool = list(dynamics_fun = "constant_other",
                                      rate_fun = NULL, encounter_fun = NULL,
                                      mort_fun = NULL, colour = "grey",
                                      linetype = "solid", integration_order = 1L)),
        processes = list()
    )
    p$mizer_version <- "3.4.0.9005"
    upgraded <- upgradeParams(p)
    expect_identical(upgraded$capabilities$schema_version, 7L)
    expect_identical(upgraded$capabilities$components$pool$integration_order, 1L)
    expect_identical(upgraded$capabilities$components$pool$role, "dynamic")
    expect_null(upgraded$capabilities$components$pool$convergence_fun)
    expect_identical(assessment_capability(upgraded, "resource")$role,
                     "prescribed")
    expect_identical(p$capabilities$schema_version, 2L)
})

test_that("Validating setters preserve extension upgrade ordering and stamps", {
    table <- get(".__S3MethodsTable__.", asNamespace("utils"))
    names <- c("upgrade.stats", "upgrade.utils")
    old <- mget(names, table, ifnotfound = list(NULL))
    withr::defer({
        for (name in names) {
            if (is.null(old[[name]])) rm(list = name, envir = table) else
                assign(name, old[[name]], envir = table)
        }
    })
    calls <- character()
    inner <- function(p) {
        calls <<- c(calls, "inner:start")
        restored <- unserialize(serialize(p, NULL))
        expect_false(extension_upgrade_active(restored))
        expect_false(extension_upgrade_active(clone_params_runtime(restored)))
        copied <- clone_params_runtime(p)
        expect_true(extension_upgrade_active(copied))
        expect_identical(attr(copied, "mizer_extension_upgrade"),
                         attr(p, "mizer_extension_upgrade"))
        expect_identical(unname(extensionVersions(p$extensions)["utils"]),
                         "0.0.0")
        p <- setComponent(p, "upgraded", initial_value = 1,
                          dynamics_fun = "constant_other")
        calls <<- c(calls, "inner:end")
        p
    }
    outer <- function(p) {
        calls <<- c(calls, "outer")
        expect_identical(unname(extensionVersions(p$extensions)["utils"]),
                         as.character(utils::packageVersion("utils")))
        expect_identical(unname(extensionVersions(p$extensions)["stats"]),
                         "0.0.0")
        setProcess(p, "consumers", autonomous = TRUE)
    }
    registerS3method("upgrade", "utils", inner, envir = asNamespace("utils"))
    registerS3method("upgrade", "stats", outer, envir = asNamespace("utils"))
    p <- NS_params_small
    p$extensions <- list(stats = c(requirement = NA, version = "0.0.0"),
                         utils = c(requirement = NA, version = "0.0.0"))
    q <- validParams(p)
    expect_identical(calls, c("inner:start", "inner:end", "outer"))
    expect_identical(q$initial_n_other$upgraded, 1)
    expect_false(extension_needs_upgrading(q))
    expect_true(extension_needs_upgrading(p))
    expect_null(attr(q, "mizer_extension_upgrade"))
    expect_length(extension_upgrade_runtime$active, 0)
})

test_that("Failed upgrades clear runtime guards and remain retryable", {
    table <- get(".__S3MethodsTable__.", asNamespace("utils"))
    name <- "upgrade.stats"
    old <- get0(name, table, inherits = FALSE)
    withr::defer({
        if (is.null(old)) rm(list = name, envir = table) else
            assign(name, old, envir = table)
    })
    captured <- NULL
    fail <- function(p) {
        captured <<- p
        p$other_params$audit$value <- 2
        p <- setProcess(p, "consumers", autonomous = TRUE)
        stop("migration failed")
    }
    registerS3method("upgrade", "stats", fail, envir = asNamespace("utils"))
    p <- NS_params_small
    p$extensions <- list(stats = c(requirement = NA, version = "0.0.0"))
    p$other_params$audit <- new.env(parent = emptyenv())
    p$other_params$audit$value <- 1
    expect_error(validParams(p), "migration failed")
    expect_equal(p$other_params$audit$value, 1)
    expect_equal(captured$other_params$audit$value, 2)
    expect_length(extension_upgrade_runtime$active, 0)
    expect_false(extension_upgrade_active(captured))
    expect_null(attr(p, "mizer_extension_upgrade"))
    registerS3method("upgrade", "stats", function(p) p,
                     envir = asNamespace("utils"))
    q <- validParams(captured)
    expect_equal(q$other_params$audit$value, 2)
    expect_equal(p$other_params$audit$value, 1)
    expect_false(extension_needs_upgrading(q))
    expect_null(attr(q, "mizer_extension_upgrade"))
    expect_true(extension_needs_upgrading(p))
})

test_that("Independent models can upgrade inside another migration", {
    table <- get(".__S3MethodsTable__.", asNamespace("utils"))
    name <- "upgrade.stats"
    old <- get0(name, table, inherits = FALSE)
    withr::defer({
        if (is.null(old)) rm(list = name, envir = table) else
            assign(name, old, envir = table)
    })
    p <- NS_params_small
    p$extensions <- list(stats = c(requirement = NA, version = "0.0.0"))
    calls <- 0L
    nested <- NULL
    upgrade <- function(candidate) {
        calls <<- calls + 1L
        if (calls == 1L) nested <<- validParams(p)
        validParams(candidate)
    }
    registerS3method("upgrade", "stats", upgrade,
                     envir = asNamespace("utils"))
    q <- validParams(p)
    expect_identical(calls, 2L)
    expect_false(extension_needs_upgrading(q))
    expect_false(extension_needs_upgrading(nested))
    expect_null(attr(q, "mizer_extension_upgrade"))
    expect_null(attr(nested, "mizer_extension_upgrade"))
    expect_length(extension_upgrade_runtime$active, 0)
})

test_that("Scaling inside an upgrade retains the active guard", {
    table <- get(".__S3MethodsTable__.", asNamespace("utils"))
    name <- "upgrade.stats"
    old <- get0(name, table, inherits = FALSE)
    withr::defer({
        if (is.null(old)) rm(list = name, envir = table) else
            assign(name, old, envir = table)
    })
    hook_calls <- 0L
    rlang::local_bindings(
        scale_upgrade_pool = function(params, before, component, details, ...) {
            hook_calls <<- hook_calls + 1L
            expect_true(extension_upgrade_active(params))
            expect_true(extension_upgrade_active(before))
            params$initial_n_other[[component]] <-
                before$initial_n_other[[component]] * details$factor
            params
        }, .env = globalenv()
    )
    calls <- 0L
    upgrade <- function(candidate) {
        calls <<- calls + 1L
        if (calls > 1L) stop("Scaling restarted the migration")
        candidate <- setComponent(candidate, "stock", initial_value = 2,
            dynamics_fun = "constant_other", transformations = list(
                scale_abundance = list(action = "update",
                                       transform_fun = "scale_upgrade_pool"),
                scale_rates = list(action = "none")
            )
        )
        for (owner in c("consumers", "resource", "reproduction")) {
            candidate <- setProcess(candidate, owner, transformations = list(
                scale_abundance = list(action = "none"),
                scale_rates = list(action = "none")
            ))
        }
        scaleRates(scaleModel(candidate, 2), 3)
    }
    registerS3method("upgrade", "stats", upgrade,
                     envir = asNamespace("utils"))
    p <- NS_params_small
    p$extensions <- list(stats = c(requirement = NA, version = "0.0.0"))
    q <- validParams(p)
    expect_identical(calls, 1L)
    expect_identical(hook_calls, 1L)
    expect_equal(q$initial_n, 2 * p$initial_n)
    expect_equal(q$initial_n_pp, 2 * p$initial_n_pp)
    expect_equal(q$metab, 3 * p$metab, ignore_attr = TRUE)
    expect_equal(q$initial_n_other$stock, 4)
    expect_null(p$initial_n_other$stock)
    expect_false(extension_needs_upgrading(q))
    expect_true(extension_needs_upgrading(p))
    expect_null(attr(q, "mizer_extension_upgrade"))
    expect_length(extension_upgrade_runtime$active, 0)
})

test_that("Grid reconstruction inside an upgrade retains the active guard", {
    table <- get(".__S3MethodsTable__.", asNamespace("utils"))
    name <- "upgrade.stats"
    old <- get0(name, table, inherits = FALSE)
    withr::defer({
        if (is.null(old)) rm(list = name, envir = table) else
            assign(name, old, envir = table)
    })
    calls <- 0L
    upgrade <- function(candidate) {
        calls <<- calls + 1L
        if (calls > 1L) stop("Grid reconstruction restarted the migration")
        for (owner in c("consumers", "resource", "reproduction")) {
            candidate <- setProcess(candidate, owner, transformations = list(
                size_grid = list(action = "none")
            ))
        }
        step <- candidate$w[[2L]] / candidate$w[[1L]]
        adjustSizeGrid(candidate, new_max_w = max(candidate$w) * step)
    }
    registerS3method("upgrade", "stats", upgrade,
                     envir = asNamespace("utils"))
    p <- NS_params_small
    p$extensions <- list(stats = c(requirement = NA, version = "0.0.0"))
    q <- validParams(p)
    expect_identical(calls, 1L)
    expect_length(q$w, length(p$w) + 1L)
    expect_equal(q$initial_n[, seq_along(p$w)], p$initial_n)
    expect_false(extension_needs_upgrading(q))
    expect_true(extension_needs_upgrading(p))
    expect_null(attr(q, "mizer_extension_upgrade"))
    expect_length(extension_upgrade_runtime$active, 0)
})
