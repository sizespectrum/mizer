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
