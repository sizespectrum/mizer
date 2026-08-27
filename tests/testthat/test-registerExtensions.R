test_that("registerExtensions creates marker class chains in dispatch order", {
    clearExtensionChain()
    withr::defer(clearExtensionChain())

    ext_a <- paste0("mizerTestA", Sys.getpid())
    ext_b <- paste0("mizerTestB", Sys.getpid())
    chain <- setNames(c(NA_character_, NA_character_), c(ext_b, ext_a))

    expect_invisible(registerExtensions(chain))
    expect_identical(getRegisteredExtensions(), chain)
    expect_true(methods::extends(ext_b, ext_a))
    expect_true(methods::extends(ext_a, "MizerParams"))
    expect_true(methods::extends(simExtensionClass(ext_b), simExtensionClass(ext_a)))
    expect_true(methods::extends(simExtensionClass(ext_a), "MizerSim"))

    params <- NS_params_small
    params@extensions <- chain
    params <- coerceToExtensionClass(params)
    expect_s4_class(params, ext_b)
    expect_true(is(params, ext_a))
})

test_that("registerExtensions accepts suffixes and prepended superchains", {
    clearExtensionChain()
    withr::defer(clearExtensionChain())

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

    # Registering a suffix also repairs the active superchain. Removing the
    # inner class is what a reload of the inner package does, and it also
    # strips `ext_a` from the `contains` list of `ext_b`, so the outer class
    # has to be rebuilt as well.
    expect_true(methods::removeClass(ext_a))
    expect_true(methods::removeClass(simExtensionClass(ext_a)))
    registerExtensions(inner)
    expect_identical(getRegisteredExtensions(), full)
    expect_true(methods::extends(ext_b, ext_a))
    expect_true(methods::extends(ext_a, "MizerParams"))
    expect_true(methods::extends(simExtensionClass(ext_b),
                                 simExtensionClass(ext_a)))
    expect_true(methods::extends(simExtensionClass(ext_a), "MizerSim"))

    expect_error(registerExtensions(incompatible),
                 "different extension chain is already active")
})

test_that("coercion uses the object's own suffix chain", {
    clearExtensionChain()
    withr::defer(clearExtensionChain())

    ext_a <- paste0("mizerTestCoerceA", Sys.getpid())
    ext_b <- paste0("mizerTestCoerceB", Sys.getpid())
    chain <- setNames(rep(NA_character_, 2), c(ext_b, ext_a))
    suffix <- setNames(NA_character_, ext_a)
    registerExtensions(chain)

    params <- NS_params_small
    params@extensions <- suffix
    params <- coerceToExtensionClass(params)
    expect_s4_class(params, ext_a)
    expect_false(is(params, ext_b))

    sim <- MizerSim(params, t_dimnames = 0)
    expect_s4_class(sim, simExtensionClass(ext_a))
    expect_false(is(sim, simExtensionClass(ext_b)))
})

test_that("S3 dispatch follows registered extension order", {
    clearExtensionChain()
    withr::defer(clearExtensionChain())

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

    params <- NS_params_small
    params@extensions <- chain
    params <- coerceToExtensionClass(params)
    expect_equal(extensionChainTestDispatch(params), "B A base")
})

test_that("base objects remain valid in extension sessions", {
    clearExtensionChain()
    withr::defer(clearExtensionChain())

    ext_a <- paste0("mizerTestBaseA", Sys.getpid())
    registerExtensions(setNames(NA_character_, ext_a))

    params <- validParams(NS_params_small)
    expect_s4_class(params, "MizerParams")

    sim <- MizerSim(params, t_dimnames = 0)
    expect_s4_class(sim, "MizerSim")
})

test_that("classless extensions remain metadata only", {
    clearExtensionChain()
    withr::defer(clearExtensionChain())

    chain <- c(stats = "0.0")
    registerExtensions(chain)

    params <- NS_params_small
    params@extensions <- chain
    params <- coerceToExtensionClass(params)

    expect_s4_class(params, "MizerParams")
    expect_false(usesExtensionDispatch(params))
    expect_identical(params@extensions, chain)

    sim <- MizerSim(params, t_dimnames = 0)
    expect_s4_class(sim, "MizerSim")
    expect_identical(sim@params@extensions, chain)
})

test_that("registerExtension prepends to chain in load order", {
    clearExtensionChain()
    withr::defer(clearExtensionChain())

    ext_a <- paste0("mizerTestIncA", Sys.getpid())
    ext_b <- paste0("mizerTestIncB", Sys.getpid())

    # A loads first (as a dependency of B would), then B prepends itself
    registerExtension(ext_a)
    expect_identical(getRegisteredExtensions(), setNames(NA_character_, ext_a))

    registerExtension(ext_b)
    expected <- setNames(rep(NA_character_, 2), c(ext_b, ext_a))
    expect_identical(getRegisteredExtensions(), expected)

    # B is outermost: extends A extends MizerParams
    expect_true(methods::extends(ext_b, ext_a))
    expect_true(methods::extends(ext_a, "MizerParams"))
    expect_true(methods::extends(simExtensionClass(ext_b),
                                 simExtensionClass(ext_a)))
    expect_true(methods::extends(simExtensionClass(ext_a), "MizerSim"))
})

test_that("registerExtension is idempotent", {
    clearExtensionChain()
    withr::defer(clearExtensionChain())

    ext_a <- paste0("mizerTestIdemA", Sys.getpid())
    ext_b <- paste0("mizerTestIdemB", Sys.getpid())

    registerExtension(ext_a)
    registerExtension(ext_b)
    full_chain <- getRegisteredExtensions()

    # Calling again for either extension leaves the chain unchanged
    registerExtension(ext_a)
    expect_identical(getRegisteredExtensions(), full_chain)

    registerExtension(ext_b)
    expect_identical(getRegisteredExtensions(), full_chain)

    # Idempotent registration also repairs marker classes lost during reload,
    # whether the reloaded package is the outermost one or an inner one.
    expect_true(methods::removeClass(ext_b))
    expect_true(methods::removeClass(simExtensionClass(ext_b)))
    registerExtension(ext_b)
    expect_identical(getRegisteredExtensions(), full_chain)
    expect_true(methods::extends(ext_b, ext_a))
    expect_true(methods::extends(simExtensionClass(ext_b),
                                 simExtensionClass(ext_a)))

    expect_true(methods::removeClass(ext_a))
    expect_true(methods::removeClass(simExtensionClass(ext_a)))
    registerExtension(ext_a)
    expect_identical(getRegisteredExtensions(), full_chain)
    expect_true(methods::extends(ext_b, ext_a))
    expect_true(methods::extends(ext_a, "MizerParams"))
    expect_true(methods::extends(simExtensionClass(ext_b),
                                 simExtensionClass(ext_a)))
    expect_true(methods::extends(simExtensionClass(ext_a), "MizerSim"))

    params <- NS_params_small
    params@extensions <- full_chain
    expect_s4_class(coerceToExtensionClass(params), ext_b)
})

test_that("marker class repair rebuilds a middle class of a longer chain", {
    clearExtensionChain()
    withr::defer(clearExtensionChain())

    ext_a <- paste0("mizerTestMidA", Sys.getpid())
    ext_b <- paste0("mizerTestMidB", Sys.getpid())
    ext_c <- paste0("mizerTestMidC", Sys.getpid())

    registerExtension(ext_a)
    registerExtension(ext_b)
    registerExtension(ext_c)
    full_chain <- getRegisteredExtensions()

    expect_true(methods::removeClass(ext_b))
    expect_true(methods::removeClass(simExtensionClass(ext_b)))
    registerExtension(ext_b)

    expect_identical(getRegisteredExtensions(), full_chain)
    expect_true(methods::extends(ext_c, ext_b))
    expect_true(methods::extends(ext_b, ext_a))
    expect_true(methods::extends(ext_a, "MizerParams"))
    expect_true(methods::extends(simExtensionClass(ext_c),
                                 simExtensionClass(ext_b)))
    expect_true(methods::extends(simExtensionClass(ext_b),
                                 simExtensionClass(ext_a)))

    params <- NS_params_small
    params@extensions <- full_chain
    expect_s4_class(coerceToExtensionClass(params), ext_c)
})

test_that("marker class repair leaves an intact chain untouched", {
    clearExtensionChain()
    withr::defer(clearExtensionChain())

    ext_a <- paste0("mizerTestIntactA", Sys.getpid())
    ext_b <- paste0("mizerTestIntactB", Sys.getpid())

    registerExtension(ext_a)
    registerExtension(ext_b)

    # An intact chain is reported as needing no work and is not rebuilt.
    expect_true(extensionClassesIntact(getRegisteredExtensions()))
    definition <- methods::getClassDef(ext_b)
    expect_false(repairExtensionClasses(getRegisteredExtensions()))
    expect_identical(methods::getClassDef(ext_b), definition)

    # A broken chain is detected and rebuilt.
    expect_true(methods::removeClass(ext_a))
    expect_false(extensionClassesIntact(getRegisteredExtensions()))
    expect_true(repairExtensionClasses(getRegisteredExtensions()))
    expect_true(extensionClassesIntact(getRegisteredExtensions()))
})

test_that("dynamic marker classes survive a cleanEx global wipe", {
    clearExtensionChain()
    withr::defer(clearExtensionChain())

    ext_a <- paste0("mizerTestCleanExA", Sys.getpid())
    chain <- setNames(NA_character_, ext_a)
    registerExtensions(chain)

    classes <- c(ext_a, simExtensionClass(ext_a))
    metadata_names <- vapply(classes, methods::classMetaName, character(1))
    class_environment <- extensionClassEnvironment(create = FALSE)
    expect_true(all(vapply(metadata_names, exists, logical(1),
                           envir = class_environment, inherits = FALSE)))
    expect_false(any(vapply(metadata_names, exists, logical(1),
                            envir = .GlobalEnv, inherits = FALSE)))

    # Remove exactly the bindings that cleanEx() would have removed under the
    # old implementation. The persistent bindings must remain resolvable from
    # inside mizer's namespace without another registration call.
    global_bindings <- metadata_names[vapply(
        metadata_names, exists, logical(1), envir = .GlobalEnv,
        inherits = FALSE
    )]
    if (length(global_bindings) > 0) {
        rm(list = global_bindings, envir = .GlobalEnv)
    }

    params <- NS_params_small
    params@extensions <- chain
    expect_s4_class(coerceToExtensionClass(params), ext_a)
    expect_true(extensionClassesIntact(chain))
})

test_that("marker class repair detects a missing metadata binding", {
    clearExtensionChain()
    withr::defer(clearExtensionChain())

    ext_a <- paste0("mizerTestBindingA", Sys.getpid())
    chain <- setNames(NA_character_, ext_a)
    registerExtensions(chain)

    class_environment <- extensionClassEnvironment(create = FALSE)
    metadata_name <- methods::classMetaName(ext_a)
    rm(list = metadata_name, envir = class_environment)

    # rm() leaves methods' class cache behind, reproducing the state that made
    # isClass() an insufficient integrity check in #587.
    expect_true(methods::isClass(ext_a))
    expect_false(markerClassPresent(ext_a))
    expect_false(extensionClassesIntact(chain))

    expect_true(repairExtensionClasses(chain))
    expect_true(markerClassPresent(ext_a))
    expect_true(extensionClassesIntact(chain))
})

test_that("marker class repair covers extensions with version requirements", {
    clearExtensionChain()
    withr::defer(clearExtensionChain())

    ext_a <- paste0("mizerTestVerA", Sys.getpid())
    ext_b <- paste0("mizerTestVerB", Sys.getpid())
    chain <- setNames(c("1.2.0", "0.4.1"), c(ext_b, ext_a))

    # Stand in for two installed extension packages that register S3 dispatch
    # methods for their marker classes. With a non-NA requirement,
    # `dispatchExtensions()` can no longer fall back on `isClass()` once the
    # class has gone, so keeping the extension in the chain rests entirely on
    # `providesDispatchMethods()`.
    local_mocked_bindings(
        ensureExtensionNamespaces = function(extensions, install = FALSE) {
            invisible(TRUE)
        },
        providesDispatchMethods = function(name) {
            name %in% c(ext_a, ext_b, simExtensionClass(ext_a),
                        simExtensionClass(ext_b))
        }
    )

    registerExtensions(chain)
    expect_true(methods::extends(ext_b, ext_a))

    expect_true(methods::removeClass(ext_a))
    expect_true(methods::removeClass(simExtensionClass(ext_a)))
    expect_false(methods::isClass(ext_a))

    registerExtension(ext_a, "0.4.1")
    expect_identical(getRegisteredExtensions(), chain)
    expect_true(methods::extends(ext_a, "MizerParams"))
    expect_true(methods::extends(ext_b, ext_a))
    expect_true(methods::extends(simExtensionClass(ext_a), "MizerSim"))
    expect_true(methods::extends(simExtensionClass(ext_b),
                                 simExtensionClass(ext_a)))
})

test_that("registerExtension does not install or check other extensions", {
    clearExtensionChain()
    withr::defer(clearExtensionChain())

    ext_a <- paste0("mizerTestNoInstallA", Sys.getpid())
    ext_b <- paste0("mizerTestNoInstallB", Sys.getpid())

    registerExtension(ext_a)
    registerExtension(ext_b)
    full_chain <- getRegisteredExtensions()

    # A sibling extension that is not installed must not turn a repeated
    # registration into an error, because that error would abort the `.onLoad`
    # of the package doing the re-registration.
    .mizerSession$extensions <- c(full_chain,
                                  setNames("9.9.9", "mizerNotInstalledExt"))
    withr::defer(.mizerSession$extensions <- character())

    expect_no_error(registerExtension(ext_b))
    expect_true(methods::extends(ext_b, ext_a))
})

test_that("registerExtension coerces objects to correct subclass", {
    clearExtensionChain()
    withr::defer(clearExtensionChain())

    ext_a <- paste0("mizerTestCoerceIncA", Sys.getpid())
    ext_b <- paste0("mizerTestCoerceIncB", Sys.getpid())

    registerExtension(ext_a)
    registerExtension(ext_b)
    full_chain <- getRegisteredExtensions()

    # Object carrying only A's chain coerces to ext_a, not ext_b
    params_a <- NS_params_small
    params_a@extensions <- setNames(NA_character_, ext_a)
    params_a <- coerceToExtensionClass(params_a)
    expect_s4_class(params_a, ext_a)
    expect_false(is(params_a, ext_b))

    # Object carrying the full chain coerces to ext_b
    params_b <- NS_params_small
    params_b@extensions <- full_chain
    params_b <- coerceToExtensionClass(params_b)
    expect_s4_class(params_b, ext_b)
    expect_true(is(params_b, ext_a))
})

test_that("readParams registers and coerces saved extension objects", {
    clearExtensionChain()
    withr::defer(clearExtensionChain())

    ext_a <- paste0("mizerTestReadA", Sys.getpid())
    chain <- setNames(NA_character_, ext_a)
    registerExtensions(chain)

    params <- NS_params_small
    params@extensions <- chain
    params <- coerceToExtensionClass(params)

    tmp <- tempfile(fileext = ".rds")
    withr::defer(unlink(tmp))
    saveParams(params, tmp)

    saved <- readRDS(tmp)
    expect_s4_class(saved, "MizerParams")
    expect_identical(saved@extensions, chain)

    clearExtensionChain()
    params2 <- readParams(tmp)
    expect_identical(getRegisteredExtensions(), chain)
    expect_s4_class(params2, ext_a)
})

# extension dispatch ----

test_that("getEncounter dispatches through extension chain", {
    clearExtensionChain()
    withr::defer(clearExtensionChain())

    ext <- paste0("mizerTestEncounter", Sys.getpid())
    chain <- setNames(NA_character_, ext)
    registerExtensions(chain)

    method <- function(params, n, n_pp, n_other, t = 0, ...) {
        NextMethod() + 1
    }
    registerS3method(
        "projectEncounter", ext, method,
        envir = asNamespace("mizer")
    )

    params <- NS_params_small
    params@extensions <- chain
    params <- coerceToExtensionClass(params)

    base <- projectEncounter.MizerParams(
        params,
        n = initialN(params),
        n_pp = initialNResource(params),
        n_other = initialNOther(params),
        t = 0
    )

    expect_equal(getEncounter(params), base + 1, ignore_attr = TRUE)
    expect_equal(projectEncounter(params, n = initialN(params),
                                  n_pp = initialNResource(params),
                                  n_other = initialNOther(params),
                                  t = 0),
                 base + 1, ignore_attr = TRUE)
    expect_equal(getRates(params)$encounter, base + 1, ignore_attr = TRUE)
})

test_that("getEncounter composes a two-extension chain outermost-first", {
    clearExtensionChain()
    withr::defer(clearExtensionChain())

    suffix <- Sys.getpid()
    extA <- paste0("mizerTestChainA", suffix)
    extB <- paste0("mizerTestChainB", suffix)
    # The chain is stored outermost-first, so extB dispatches before extA.
    chain <- setNames(c(NA_character_, NA_character_), c(extB, extA))
    registerExtensions(chain)

    order <- new.env(parent = emptyenv())
    order$log <- character()
    registerS3method(
        "projectEncounter", extA,
        function(params, ...) {
            order$log <- c(order$log, "extA")
            NextMethod() + 1
        },
        envir = asNamespace("mizer")
    )
    registerS3method(
        "projectEncounter", extB,
        function(params, ...) {
            order$log <- c(order$log, "extB")
            NextMethod() + 10
        },
        envir = asNamespace("mizer")
    )

    params <- NS_params_small
    params@extensions <- chain
    params <- coerceToExtensionClass(params)
    # The object is coerced to the outermost class, extending extA and MizerParams.
    expect_identical(class(params)[[1]], extB)
    expect_true(is(params, extA))

    base <- projectEncounter.MizerParams(
        params,
        n = initialN(params),
        n_pp = initialNResource(params),
        n_other = initialNOther(params),
        t = 0
    )

    result <- getEncounter(params)
    # Outermost extB runs first, then extA, then the base method.
    expect_identical(order$log, c("extB", "extA"))
    expect_equal(result, base + 11, ignore_attr = TRUE)
})

test_that("projectEncounter base method is mizerEncounter", {
    expect_identical(projectEncounter.MizerParams, mizerEncounter)
    expect_identical(projectFeedingLevel.MizerParams, mizerFeedingLevel)
    expect_identical(projectEReproAndGrowth.MizerParams, mizerEReproAndGrowth)
    expect_identical(projectERepro.MizerParams, mizerERepro)
    expect_identical(projectEGrowth.MizerParams, mizerEGrowth)
    expect_identical(projectDiffusion.MizerParams, mizerDiffusion)
    expect_identical(projectPredRate.MizerParams, mizerPredRate)
    expect_identical(projectPredMort.MizerParams, mizerPredMort)
    expect_identical(projectFMort.MizerParams, mizerFMort)
    expect_identical(projectMort.MizerParams, mizerMort)
    expect_identical(projectRDI.MizerParams, mizerRDI)
    expect_identical(projectResourceMort.MizerParams, mizerResourceMort)
})

test_that("getRates dispatches through all projection hooks", {
    clearExtensionChain()
    withr::defer(clearExtensionChain())

    ext <- paste0("mizerTestRates", Sys.getpid())
    chain <- setNames(NA_character_, ext)
    registerExtensions(chain)

    hooks <- c("projectEncounter", "projectFeedingLevel",
               "projectEReproAndGrowth", "projectERepro", "projectEGrowth",
               "projectDiffusion", "projectPredRate", "projectPredMort",
               "projectFMort", "projectMort", "projectRDI",
               "projectRDD", "projectResourceMort")
    called <- new.env(parent = emptyenv())
    for (hook in hooks) {
        local({
            hook_name <- hook
            registerS3method(
                hook_name, ext,
                function(params, ...) {
                    called[[hook_name]] <- TRUE
                    NextMethod()
                },
                envir = asNamespace("mizer")
            )
        })
    }

    params <- NS_params_small
    params@extensions <- chain
    params <- coerceToExtensionClass(params)

    getRates(params)

    expect_true(all(vapply(hooks, function(hook) {
        isTRUE(called[[hook]])
    }, logical(1))))
})

test_that("getEncounter honours rates_funcs for base objects", {
    params <- NS_params_small
    params@rates_funcs$Encounter <- "constant_encounter_for_dispatch_test"

    assign("constant_encounter_for_dispatch_test",
           function(params, n, n_pp, n_other, t = 0, ...) {
        params@initial_n * 0 + t
    }, envir = .GlobalEnv)
    withr::defer(rm("constant_encounter_for_dispatch_test", envir = .GlobalEnv))

    expect_equal(getEncounter(params, t = 2), params@initial_n * 0 + 2,
                 ignore_attr = TRUE)
})

test_that("setRateFunction is honoured when extension dispatch is active", {
    clearExtensionChain()
    withr::defer(clearExtensionChain())

    ext <- paste0("mizerTestSetRateFunc", Sys.getpid())
    chain <- setNames(NA_character_, ext)
    registerExtensions(chain)

    params <- NS_params_small
    params@extensions <- chain
    params <- coerceToExtensionClass(params)

    assign("constant_encounter_for_dispatch_test",
           function(params, n, n_pp, n_other, t = 0, ...) {
               params@initial_n * 0 + 42
           }, envir = .GlobalEnv)
    withr::defer(rm("constant_encounter_for_dispatch_test", envir = .GlobalEnv))

    params <- setRateFunction(params, "Encounter", "constant_encounter_for_dispatch_test")

    expect_equal(getRates(params)$encounter, params@initial_n * 0 + 42,
                 ignore_attr = TRUE)
})

test_that("classless extensions do not trigger project dispatch", {
    clearExtensionChain()
    withr::defer(clearExtensionChain())

    chain <- c(stats = "0.0")
    registerExtensions(chain)

    params <- NS_params_small
    params@extensions <- chain
    params@rates_funcs$Rates <- "classless_rates_for_dispatch_test"
    params@rates_funcs$Encounter <- "classless_encounter_for_dispatch_test"

    assign("classless_rates_for_dispatch_test",
           function(params, n, n_pp, n_other, t = 0, effort, rates_fns, ...) {
        list(classless = TRUE)
    }, envir = .GlobalEnv)
    assign("classless_encounter_for_dispatch_test",
           function(params, n, n_pp, n_other, t = 0, ...) {
        params@initial_n * 0 + t
    }, envir = .GlobalEnv)
    withr::defer(rm("classless_rates_for_dispatch_test",
                    "classless_encounter_for_dispatch_test",
                    envir = .GlobalEnv))

    expect_identical(getRates(params)$classless, TRUE)
    expect_equal(getEncounter(params, t = 3), params@initial_n * 0 + 3,
                 ignore_attr = TRUE)
})

test_that("dispatch extensions are detected from their registered S3 methods", {
    # An installed extension participates in dispatch if its package registers
    # S3 methods for its marker class, even with a non-NA requirement and no
    # statically defined S4 class. This lets packages omit the static setClass
    # and let mizer build the class hierarchy dynamically so they can be chained.
    ns <- asNamespace("mizer")
    orig <- getNamespaceInfo(ns, "S3methods")
    withr::defer(setNamespaceInfo(ns, "S3methods", orig))
    registerS3method("format", "mizer", function(x, ...) x, envir = ns)

    expect_true(providesDispatchMethods("mizer"))
    expect_true("mizer" %in%
                    names(dispatchExtensions(c(mizer = "1.0.0"))))

    # A package that registers no methods for its own name is not a dispatcher,
    # nor is an unloaded package name.
    expect_false(providesDispatchMethods("methods"))
    expect_false(providesDispatchMethods("nonexistent_pkg_xyz_123"))
    expect_false("stats" %in% names(dispatchExtensions(c(stats = "0.0"))))
})

# extension versioning ----

test_that("upgrade is dispatched as an S3 generic with our methods", {
    methods <- as.character(utils::.S3methods("upgrade"))
    expect_true("upgrade.MizerParams" %in% methods)
    expect_true("upgrade.MizerSim" %in% methods)
    # The base R method must still be reachable.
    expect_false(is.null(getS3method("upgrade", "packageStatus", optional = TRUE)))
})

test_that("extensionRequirements / extensionVersions read both slot forms", {
    # legacy character form: requirements only, versions all NA
    char <- c(extA = "owner/repo", extB = NA_character_)
    expect_identical(extensionRequirements(char), char)
    expect_identical(extensionVersions(char),
                     setNames(c(NA_character_, NA_character_), c("extA", "extB")))

    # versioned list form
    lst <- makeExtensions(c(extA = "owner/repo"), c(extA = "1.2.3"))
    expect_identical(extensionRequirements(lst), c(extA = "owner/repo"))
    expect_identical(extensionVersions(lst), c(extA = "1.2.3"))

    # empty
    expect_identical(extensionRequirements(character()), character())
})

test_that("recordExtension stamps without disturbing other entries", {
    p <- NS_params
    p@extensions <- list(other = c(requirement = NA_character_,
                                   version = NA_character_),
                         extA = c(requirement = "owner/repo",
                                  version = NA_character_))
    before <- p@extensions

    # version = NULL on an existing entry is a no-op on the whole slot
    p2 <- recordExtension(p, "extA")
    expect_identical(p2@extensions, before)

    # stamping updates only the named entry
    p3 <- recordExtension(p, "extA", version = "2.0.0")
    expect_identical(p3@extensions$other, before$other)
    expect_identical(unname(p3@extensions$extA[["version"]]), "2.0.0")
})

test_that("recordExtension prepends genuinely new entries", {
    p <- NS_params
    p@extensions <- list(extA = c(requirement = "owner/repo",
                                  version = NA_character_))

    # unversioned new entry goes to the front
    p2 <- recordExtension(p, "extB")
    expect_identical(names(p2@extensions), c("extB", "extA"))

    # versioned new entry also goes to the front, existing entry untouched
    p3 <- recordExtension(p, "extB", version = "1.0.0")
    expect_identical(names(p3@extensions), c("extB", "extA"))
    expect_identical(unname(p3@extensions$extB[["version"]]), "1.0.0")
    expect_identical(p3@extensions$extA, p@extensions$extA)

    # legacy character slot: new entry prepended, slot stays a character vector
    p@extensions <- c(extA = "owner/repo")
    p4 <- recordExtension(p, "extB")
    expect_false(is.list(p4@extensions))
    expect_identical(names(p4@extensions), c("extB", "extA"))
})

test_that("recordExtension can stamp onto an empty slot", {
    p <- NS_params
    expect_length(p@extensions, 0)
    p <- recordExtension(p, "extA", version = "1.0.0")
    expect_true(is.list(p@extensions))
    expect_identical(unname(p@extensions$extA[["requirement"]]), NA_character_)
    expect_identical(unname(p@extensions$extA[["version"]]), "1.0.0")
})

test_that("registerExtensions accepts the versioned list form", {
    clearExtensionChain()
    withr::defer(clearExtensionChain())
    ext <- paste0("mizerTestVer", Sys.getpid())
    lst <- makeExtensions(setNames(NA_character_, ext), character())
    expect_invisible(registerExtensions(lst))
    expect_identical(getRegisteredExtensions(), setNames(NA_character_, ext))
})

test_that("extension_needs_upgrading reacts to missing/stale stamps", {
    p <- NS_params
    # No extensions -> never needs upgrading on extension grounds
    expect_false(extension_needs_upgrading(p))

    # An extension whose package is not installed is ignored
    p@extensions <- makeExtensions(c(definitelyMissingPkg = "owner/repo"),
                                   character())
    expect_false(extension_needs_upgrading(p))

    # mizer itself, with a stale stamp, does need upgrading; with the installed
    # version it does not. (mizer is, of course, installed.)
    p@extensions <- makeExtensions(c(mizer = "x"), c(mizer = "0.0.1"))
    expect_true(extension_needs_upgrading(p))
    p@extensions <- makeExtensions(
        c(mizer = "x"),
        c(mizer = as.character(utils::packageVersion("mizer"))))
    expect_false(extension_needs_upgrading(p))
})
