test_that("saveParams/readParams round-trip", {
    params <- NS_params_small
    tmp <- tempfile(fileext = ".rds")
    on.exit(unlink(tmp), add = TRUE)
    expect_invisible(saveParams(params, tmp))
    params2 <- readParams(tmp)
    expect_s3_class(params2, "MizerParams")
    expect_identical(dim(params2@initial_n), dim(params@initial_n))
    expect_identical(dimnames(params2@initial_n), dimnames(params@initial_n))
    expect_identical(params2@species_params$species, params@species_params$species)
})

test_that("saving checks callbacks nested in tuning declarations", {
    rlang::local_bindings(
        saved_equation_rates = function(params, ...) list(),
        saved_balance = function(params, ...) params,
        saved_preserved = function(params, ...) 1,
        .env = globalenv()
    )
    strategy <- list(action = "balance", balance_fun = "saved_balance",
                     preserve = list(quantity = list(fun = "saved_preserved")))
    p <- setProcess(
        NS_params_small, "reproduction", stationary = TRUE,
        equation_rates_fun = "saved_equation_rates",
        tuning = list(default = "custom", strategies = list(custom = strategy))
    )
    tmp <- tempfile(fileext = ".rds")
    withr::defer(unlink(tmp))
    expect_warning(saveParams(p, tmp),
                   "saved_equation_rates.*saved_balance.*saved_preserved")
    restored <- readRDS(tmp)
    expect_equal(restored$capabilities, p$capabilities)
})

test_that("saved transformation declarations retain callback behaviour", {
    rlang::local_bindings(
        saved_pool_scale = function(params, before, component, details, ...) {
            state <- before$initial_n_other[[component]]
            state$stock <- state$stock * details$factor
            params$initial_n_other[[component]] <- state
            params
        },
        .env = globalenv()
    )
    p <- setComponent(
        NS_params_small, "pool",
        initial_value = list(stock = 2, temperature = 15, label = "inshore"),
        dynamics_fun = "constant_other",
        transformations = list(
            scale_abundance = list(
                action = "update", transform_fun = "saved_pool_scale"
            )
        )
    )
    tmp <- tempfile(fileext = ".rds")
    withr::defer(unlink(tmp))
    expect_warning(saveParams(p, tmp), "saved_pool_scale")
    restored <- readParams(tmp)
    expect_equal(restored$capabilities, p$capabilities)
    scaled <- scaleModel(restored, 3)
    expect_equal(scaled$initial_n_other$pool,
                 list(stock = 6, temperature = 15, label = "inshore"))
    expect_equal(p$initial_n_other$pool,
                 list(stock = 2, temperature = 15, label = "inshore"))
})

test_that("stored component definitions survive params and sim round trips", {
    definition <- list(dynamics_fun = "constant_other", colour = "red",
                       linetype = "dashed")
    p <- setComponent(NS_params_small, "pool", initial_value = 42,
                      definition = definition,
                      component_params = list(capacity = 80))
    stored <- p$capabilities
    definition$colour <- "blue"
    tmp <- tempfile(fileext = ".rds")
    withr::defer(unlink(tmp))

    saveParams(p, tmp)
    restored <- readParams(tmp)
    expect_identical(restored$capabilities, stored)
    expect_identical(restored$initial_n_other, p$initial_n_other)
    expect_identical(restored$other_params, p$other_params)

    sim <- project(restored, t_max = 0.1, t_save = 0.1,
                   progress_bar = FALSE)
    saveSim(sim, tmp)
    restored_sim <- readSim(tmp)
    expect_identical(restored_sim$params$capabilities, stored)
    expect_identical(restored_sim$n_other, sim$n_other)
})

test_that("saveParams preserves extension classes", {
    ext_a <- paste0("mizerTestParamsSaveA", Sys.getpid())
    chain <- setNames(NA_character_, ext_a)

    params <- NS_params_small
    params@extensions <- chain
    params <- coerceToExtensionClass(params)

    tmp <- tempfile(fileext = ".rds")
    withr::defer(unlink(tmp))
    saveParams(params, tmp)

    saved <- readRDS(tmp)
    expect_identical(class(saved), c(ext_a, "MizerParams"))
    expect_identical(saved@extensions, chain)

    params2 <- readParams(tmp)
    expect_identical(class(params2), class(params))
})

test_that("readParams repairs extension classes on legacy files", {
    ext_a <- paste0("mizerTestParamsReadA", Sys.getpid())
    chain <- setNames(NA_character_, ext_a)

    params <- NS_params_small
    params@extensions <- chain
    class(params) <- "MizerParams"

    tmp <- tempfile(fileext = ".rds")
    withr::defer(unlink(tmp))
    saveRDS(params, tmp)

    params2 <- readParams(tmp)
    expect_identical(class(params2), c(ext_a, "MizerParams"))
})

test_that("saveParams reports missing extension packages by name", {
    params <- NS_params_small
    params@extensions <- c(definitelyMissingPkg = "github::owner/repo")
    tmp <- tempfile(fileext = ".rds")
    on.exit(unlink(tmp), add = TRUE)

    expect_error(saveParams(params, tmp),
                 "Some required extension packages are not installed: definitelyMissingPkg")
})

test_that("saveSim/readSim round-trip", {
    sim <- project(NS_params_small, t_max = 0.1, t_save = 0.1)
    tmp <- tempfile(fileext = ".rds")
    on.exit(unlink(tmp), add = TRUE)

    expect_invisible(saveSim(sim, tmp))
    sim2 <- readSim(tmp)

    expect_s3_class(sim2, "MizerSim")
    expect_s3_class(sim2@params, "MizerParams")
    expect_identical(dim(sim2@n), dim(sim@n))
    expect_identical(dimnames(sim2@n), dimnames(sim@n))
    expect_identical(dim(sim2@n_pp), dim(sim@n_pp))
    expect_identical(dimnames(sim2@effort), dimnames(sim@effort))
})

test_that("saveSim preserves extension classes", {
    ext_a <- paste0("mizerTestSimReadA", Sys.getpid())
    chain <- setNames(NA_character_, ext_a)

    params <- NS_params_small
    params@extensions <- chain
    params <- coerceToExtensionClass(params)
    sim <- project(params, t_max = 0.1, t_save = 0.1)

    tmp <- tempfile(fileext = ".rds")
    withr::defer(unlink(tmp))
    saveSim(sim, tmp)

    saved <- readRDS(tmp)
    expect_identical(class(saved), c(simExtensionClass(ext_a), "MizerSim"))
    expect_identical(class(saved@params), c(ext_a, "MizerParams"))

    sim2 <- readSim(tmp)
    expect_identical(class(sim2), class(sim))
    expect_identical(class(sim2@params), class(params))
})

test_that("readSim repairs extension classes on legacy files", {
    ext_a <- paste0("mizerTestSimReadA", Sys.getpid())
    chain <- setNames(NA_character_, ext_a)

    params <- NS_params_small
    params@extensions <- chain
    params <- coerceToExtensionClass(params)
    sim <- project(params, t_max = 0.1, t_save = 0.1)
    class(sim) <- "MizerSim"
    class(sim@params) <- "MizerParams"

    tmp <- tempfile(fileext = ".rds")
    withr::defer(unlink(tmp))
    saveRDS(sim, tmp)

    sim2 <- readSim(tmp)
    expect_identical(class(sim2), c(simExtensionClass(ext_a), "MizerSim"))
    expect_identical(class(sim2@params), c(ext_a, "MizerParams"))
})

test_that("readParams reconciles the species parameters", {
    params <- NS_params_small
    # A value written straight into the slot, as an old model or code that
    # bypasses `species_params<-()` may hold it.
    new_value <- unname(species_params(params)$w_mat25[2]) / 2
    params@species_params$w_mat25[2] <- new_value
    expect_true(is.null(given_species_params(params)$w_mat25))

    tmp <- tempfile(fileext = ".rds")
    on.exit(unlink(tmp), add = TRUE)
    saveParams(params, tmp)
    params2 <- suppressMessages(readParams(tmp))

    expect_equal(unname(given_species_params(params2)$w_mat25[2]), new_value)
    expect_true(all(is.na(given_species_params(params2)$w_mat25[-2])))
    # The value now survives a recalculation
    species_params(params2)$w_mat <- species_params(params2)$w_mat
    expect_equal(unname(species_params(params2)$w_mat25[2]), new_value)
})

test_that("readParams leaves a consistent model unchanged", {
    params <- NS_params_small
    tmp <- tempfile(fileext = ".rds")
    on.exit(unlink(tmp), add = TRUE)
    saveParams(params, tmp)

    expect_silent(params2 <- readParams(tmp))
    expect_identical(given_species_params(params2),
                     given_species_params(params))
})


test_that("Equation callbacks survive saving with their stored identifiers", {
    p <- setComponent(NS_params_small, "pool", initial_value = c(a = 1, b = 2),
                      rate_fun = "base::identity",
                      component_params = list(capacity = 80))
    tmp <- tempfile(fileext = ".rds")
    withr::defer(unlink(tmp))

    # The internal stepping bridge and exported equation are package callbacks.
    expect_false(is_custom("component_rate_dynamics", packages = "mizer"))
    expect_silent(saveParams(p, tmp))
    restored <- readParams(tmp)
    expect_identical(restored$capabilities, p$capabilities)
    expect_identical(restored$initial_n_other, p$initial_n_other)
    expect_identical(restored$other_params, p$other_params)
    expect_identical(getComponent(restored, "pool")$rate_fun, "base::identity")
    expect_identical(prepare_component_rates(restored)$pool, base::identity)
})

test_that("Saving reports custom component equations", {
    rlang::local_bindings(
        saved_component_equation = function(params, n_other, component, ...) {
            -n_other[[component]]
        }, .env = .GlobalEnv
    )
    p <- setComponent(NS_params_small, "pool", initial_value = 2,
                      rate_fun = "saved_component_equation")
    tmp <- tempfile(fileext = ".rds")
    withr::defer(unlink(tmp))
    expect_warning(saveParams(p, tmp),
                   "Your model is using the functions saved_component_equation")
})

test_that("Saving accounts for component and process assessments", {
    rlang::local_bindings(
        saved_assessment = function(...) {
            data.frame(measure = "change", value = 0,
                       tolerance = 1e-6, units = "1/year")
        }, .env = .GlobalEnv
    )
    p <- setComponent(NS_params_small, "pool", initial_value = 2,
                      dynamics_fun = "constant_other",
                      convergence_fun = "saved_assessment")
    expect_warning(checkCustomFunctions(p), "saved_assessment")
    p <- removeComponent(p, "pool")
    p <- setProcess(p, "resource", convergence_fun = "saved_assessment")
    expect_warning(checkCustomFunctions(p), "saved_assessment")

    p <- setProcess(p, "resource", convergence_fun = "base::identity")
    tmp <- tempfile(fileext = ".rds")
    withr::defer(unlink(tmp))
    expect_silent(saveParams(p, tmp))
    restored <- readParams(tmp)
    expect_identical(getProcess(restored), getProcess(p))
})
