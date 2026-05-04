test_that("registerExtensions creates marker class chains in dispatch order", {
    resetMizerSession()
    withr::defer(resetMizerSession())

    ext_a <- paste0("mizerTestA", Sys.getpid())
    ext_b <- paste0("mizerTestB", Sys.getpid())
    chain <- setNames(c(NA_character_, NA_character_), c(ext_b, ext_a))

    expect_invisible(registerExtensions(chain))
    expect_identical(getRegisteredExtensions(), chain)
    expect_true(methods::extends(ext_b, ext_a))
    expect_true(methods::extends(ext_a, "MizerParams"))
    expect_true(methods::extends(simExtensionClass(ext_b), simExtensionClass(ext_a)))
    expect_true(methods::extends(simExtensionClass(ext_a), "MizerSim"))

    params <- NS_params
    params@extensions <- chain
    params <- coerceToExtensionClass(params)
    expect_s4_class(params, ext_b)
    expect_true(is(params, ext_a))
})

test_that("registerExtensions accepts suffixes and prepended superchains", {
    resetMizerSession()
    withr::defer(resetMizerSession())

    ext_a <- paste0("mizerTestSuffixA", Sys.getpid())
    ext_b <- paste0("mizerTestSuffixB", Sys.getpid())
    ext_c <- paste0("mizerTestSuffixC", Sys.getpid())
    inner <- setNames(NA_character_, ext_a)
    full <- setNames(rep(NA_character_, 2), c(ext_b, ext_a))
    incompatible <- setNames(rep(NA_character_, 2), c(ext_c, ext_a))

    registerExtensions(inner)
    expect_identical(getRegisteredExtensions(), inner)

    registerExtensions(full)
    expect_identical(getRegisteredExtensions(), full)

    registerExtensions(inner)
    expect_identical(getRegisteredExtensions(), full)

    expect_error(registerExtensions(incompatible),
                 "different extension chain is already active")
})

test_that("coercion uses the object's own suffix chain", {
    resetMizerSession()
    withr::defer(resetMizerSession())

    ext_a <- paste0("mizerTestCoerceA", Sys.getpid())
    ext_b <- paste0("mizerTestCoerceB", Sys.getpid())
    chain <- setNames(rep(NA_character_, 2), c(ext_b, ext_a))
    suffix <- setNames(NA_character_, ext_a)
    registerExtensions(chain)

    params <- NS_params
    params@extensions <- suffix
    params <- coerceToExtensionClass(params)
    expect_s4_class(params, ext_a)
    expect_false(is(params, ext_b))

    sim <- MizerSim(params, t_dimnames = 0)
    expect_s4_class(sim, simExtensionClass(ext_a))
    expect_false(is(sim, simExtensionClass(ext_b)))
})

test_that("S3 dispatch follows registered extension order", {
    resetMizerSession()
    withr::defer(resetMizerSession())

    ext_a <- paste0("mizerTestDispatchA", Sys.getpid())
    ext_b <- paste0("mizerTestDispatchB", Sys.getpid())
    chain <- setNames(rep(NA_character_, 2), c(ext_b, ext_a))
    registerExtensions(chain)

    extensionChainTestDispatch <- function(params) {
        UseMethod("extensionChainTestDispatch")
    }
    method_names <- c(
        "extensionChainTestDispatch.MizerParams",
        paste0("extensionChainTestDispatch.", ext_a),
        paste0("extensionChainTestDispatch.", ext_b)
    )
    assign(method_names[1], function(params) "base", envir = .GlobalEnv)
    assign(method_names[2], function(params) paste("A", NextMethod()),
           envir = .GlobalEnv)
    assign(method_names[3], function(params) paste("B", NextMethod()),
           envir = .GlobalEnv)
    withr::defer(rm(list = method_names, envir = .GlobalEnv))

    params <- NS_params
    params@extensions <- chain
    params <- coerceToExtensionClass(params)
    expect_equal(extensionChainTestDispatch(params), "B A base")
})

test_that("base objects remain valid in extension sessions", {
    resetMizerSession()
    withr::defer(resetMizerSession())

    ext_a <- paste0("mizerTestBaseA", Sys.getpid())
    registerExtensions(setNames(NA_character_, ext_a))

    params <- validParams(NS_params)
    expect_s4_class(params, "MizerParams")

    sim <- MizerSim(params, t_dimnames = 0)
    expect_s4_class(sim, "MizerSim")
})

test_that("classless extensions remain metadata only", {
    resetMizerSession()
    withr::defer(resetMizerSession())

    chain <- c(stats = "0.0")
    registerExtensions(chain)

    params <- NS_params
    params@extensions <- chain
    params <- coerceToExtensionClass(params)

    expect_s4_class(params, "MizerParams")
    expect_false(usesExtensionDispatch(params))
    expect_identical(params@extensions, chain)

    sim <- MizerSim(params, t_dimnames = 0)
    expect_s4_class(sim, "MizerSim")
    expect_identical(sim@params@extensions, chain)
})

test_that("readParams registers and coerces saved extension objects", {
    resetMizerSession()
    withr::defer(resetMizerSession())

    ext_a <- paste0("mizerTestReadA", Sys.getpid())
    chain <- setNames(NA_character_, ext_a)
    registerExtensions(chain)

    params <- NS_params
    params@extensions <- chain
    params <- coerceToExtensionClass(params)

    tmp <- tempfile(fileext = ".rds")
    withr::defer(unlink(tmp))
    saveParams(params, tmp)

    saved <- readRDS(tmp)
    expect_s4_class(saved, "MizerParams")
    expect_identical(saved@extensions, chain)

    resetMizerSession()
    params2 <- readParams(tmp)
    expect_identical(getRegisteredExtensions(), chain)
    expect_s4_class(params2, ext_a)
})
