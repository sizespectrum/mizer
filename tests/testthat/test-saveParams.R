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
