test_that("getEncounter dispatches through extension chain", {
    resetMizerSession()
    withr::defer(resetMizerSession())

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

    params <- NS_params
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
    resetMizerSession()
    withr::defer(resetMizerSession())

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

    params <- NS_params
    params@extensions <- chain
    params <- coerceToExtensionClass(params)

    getRates(params)

    expect_true(all(vapply(hooks, function(hook) {
        isTRUE(called[[hook]])
    }, logical(1))))
})

test_that("getEncounter honours rates_funcs for base objects", {
    params <- NS_params
    params@rates_funcs$Encounter <- "constant_encounter_for_dispatch_test"

    assign("constant_encounter_for_dispatch_test",
           function(params, n, n_pp, n_other, t = 0, ...) {
        params@initial_n * 0 + t
    }, envir = .GlobalEnv)
    withr::defer(rm("constant_encounter_for_dispatch_test", envir = .GlobalEnv))

    expect_equal(getEncounter(params, t = 2), params@initial_n * 0 + 2,
                 ignore_attr = TRUE)
})

test_that("classless extensions do not trigger project dispatch", {
    resetMizerSession()
    withr::defer(resetMizerSession())

    chain <- c(stats = "0.0")
    registerExtensions(chain)

    params <- NS_params
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
