test_that("coerceToExtensionClass sets S3 classes in dispatch order", {
    ext_a <- paste0("mizerTestA", Sys.getpid())
    ext_b <- paste0("mizerTestB", Sys.getpid())
    chain <- setNames(c(NA_character_, NA_character_), c(ext_b, ext_a))

    params <- NS_params_small
    params@extensions <- chain
    params <- coerceToExtensionClass(params)
    expect_s3_class(params, ext_b)
    expect_true(inherits(params, ext_a))
    expect_true(inherits(params, "MizerParams"))
    expect_identical(class(params), c(ext_b, ext_a, "MizerParams"))

    sim <- MizerSim(params, t_dimnames = 0)
    expect_s3_class(sim, simExtensionClass(ext_b))
    expect_true(inherits(sim, simExtensionClass(ext_a)))
    expect_true(inherits(sim, "MizerSim"))
    expect_identical(class(sim), c(simExtensionClass(ext_b), simExtensionClass(ext_a), "MizerSim"))
})

test_that("S3 dispatch follows extension order and NextMethod", {
    ext_a <- paste0("mizerTestDispatchA", Sys.getpid())
    ext_b <- paste0("mizerTestDispatchB", Sys.getpid())
    chain <- setNames(rep(NA_character_, 2), c(ext_b, ext_a))

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

test_that("base objects remain valid MizerParams", {
    params <- validParams(NS_params_small)
    expect_s3_class(params, "MizerParams")
    expect_identical(class(params), "MizerParams")

    sim <- MizerSim(params, t_dimnames = 0)
    expect_s3_class(sim, "MizerSim")
    expect_identical(class(sim), "MizerSim")
})

test_that("readParams restores extension class on saved objects", {
    ext_a <- paste0("mizerTestReadA", Sys.getpid())
    chain <- setNames(NA_character_, ext_a)

    params <- NS_params_small
    params@extensions <- chain
    params <- coerceToExtensionClass(params)

    tmp <- tempfile(fileext = ".rds")
    withr::defer(unlink(tmp))
    saveParams(params, tmp)

    saved <- readRDS(tmp)
    expect_s3_class(saved, "MizerParams")
    expect_identical(saved@extensions, chain)

    params2 <- readParams(tmp)
    expect_s3_class(params2, ext_a)
    expect_identical(class(params2), c(ext_a, "MizerParams"))
})

# extension dispatch ----

test_that("getEncounter dispatches through extension class", {
    ext <- paste0("mizerTestEncounter", Sys.getpid())
    chain <- setNames(NA_character_, ext)

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
    suffix <- Sys.getpid()
    extA <- paste0("mizerTestChainA", suffix)
    extB <- paste0("mizerTestChainB", suffix)
    chain <- setNames(c(NA_character_, NA_character_), c(extB, extA))

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
    expect_identical(class(params)[[1]], extB)
    expect_true(inherits(params, extA))

    base <- projectEncounter.MizerParams(
        params,
        n = initialN(params),
        n_pp = initialNResource(params),
        n_other = initialNOther(params),
        t = 0
    )

    result <- getEncounter(params)
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
    ext <- paste0("mizerTestRates", Sys.getpid())
    chain <- setNames(NA_character_, ext)

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
    ext <- paste0("mizerTestSetRateFunc", Sys.getpid())
    chain <- setNames(NA_character_, ext)

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
    p <- validParams(NS_params)
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
    p <- validParams(NS_params)
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
    p <- validParams(NS_params)
    expect_length(p@extensions, 0)
    p <- recordExtension(p, "extA", version = "1.0.0")
    expect_true(is.list(p@extensions))
    expect_identical(unname(p@extensions$extA[["requirement"]]), NA_character_)
    expect_identical(unname(p@extensions$extA[["version"]]), "1.0.0")
})

test_that("extension_needs_upgrading reacts to missing/stale stamps", {
    p <- validParams(NS_params)
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
