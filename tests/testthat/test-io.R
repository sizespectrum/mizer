test_that("saveParams/readParams round-trip", {
    params <- NS_params
    tmp <- tempfile(fileext = ".rds")
    on.exit(unlink(tmp), add = TRUE)
    expect_invisible(saveParams(params, tmp))
    params2 <- readParams(tmp)
    expect_s4_class(params2, "MizerParams")
    expect_identical(dim(params2@initial_n), dim(params@initial_n))
    expect_identical(dimnames(params2@initial_n), dimnames(params@initial_n))
    expect_identical(params2@species_params$species, params@species_params$species)
})

test_that("saveParams reports missing extension packages by name", {
    params <- NS_params
    params@extensions <- c(definitelyMissingPkg = "github::owner/repo")
    tmp <- tempfile(fileext = ".rds")
    on.exit(unlink(tmp), add = TRUE)

    expect_error(saveParams(params, tmp),
                 "Some required extension packages are not installed: definitelyMissingPkg")
})

test_that("saveSim/readSim round-trip", {
    sim <- project(NS_params, t_max = 0.1, t_save = 0.1)
    tmp <- tempfile(fileext = ".rds")
    on.exit(unlink(tmp), add = TRUE)

    expect_invisible(saveSim(sim, tmp))
    sim2 <- readSim(tmp)

    expect_s4_class(sim2, "MizerSim")
    expect_s4_class(sim2@params, "MizerParams")
    expect_identical(dim(sim2@n), dim(sim@n))
    expect_identical(dimnames(sim2@n), dimnames(sim@n))
    expect_identical(dim(sim2@n_pp), dim(sim@n_pp))
    expect_identical(dimnames(sim2@effort), dimnames(sim@effort))
})

test_that("saveSim stores base classes and readSim restores extension classes", {
    clearExtensionChain()
    withr::defer(clearExtensionChain())

    ext_a <- paste0("mizerTestSimReadA", Sys.getpid())
    chain <- setNames(NA_character_, ext_a)
    registerExtensions(chain)

    params <- NS_params
    params@extensions <- chain
    params <- coerceToExtensionClass(params)
    sim <- project(params, t_max = 0.1, t_save = 0.1)

    tmp <- tempfile(fileext = ".rds")
    withr::defer(unlink(tmp))
    saveSim(sim, tmp)

    saved <- readRDS(tmp)
    expect_s4_class(saved, "MizerSim")
    expect_s4_class(saved@params, "MizerParams")
    expect_identical(saved@params@extensions, chain)

    clearExtensionChain()
    sim2 <- readSim(tmp)

    expect_identical(getRegisteredExtensions(), chain)
    expect_s4_class(sim2, simExtensionClass(ext_a))
    expect_s4_class(sim2@params, ext_a)
})
