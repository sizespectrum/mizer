test_that("Whole-model assessments record complete built-in coverage", {
    report <- getSteadyAssessment(NS_params_small)
    expect_s3_class(report, "MizerSteadyAssessment")
    expect_identical(report$scope, "whole")
    expect_identical(report$criterion, "state")
    expect_identical(report$fixed, character())
    expect_setequal(report$required, c("consumers", "resource"))
    expect_setequal(report$covered, report$required)
    expect_identical(report$missing, character())
    expect_equal(nrow(report$unavailable), 0)
    expect_equal(nrow(report$measures),
                 nrow(NS_params_small$species_params) + 1)
    expect_identical(report$verdict, all(report$measures$passed))
    expect_true(all(report$measures$units == "g/year"))
    expect_equal(attr(report, "drift_report")$drift,
                 steady_biomass_drift(NS_params_small))
})

test_that("Unsupported components are unknown and are never silently stepped", {
    rlang::local_bindings(
        unsupported_step = function(...) stop("do not step"),
        unsupported_rate = function(...) stop("do not differentiate"),
        .env = .GlobalEnv
    )
    p <- setComponent(NS_params_small, "pool", initial_value = 2,
                      dynamics_fun = "unsupported_step",
                      rate_fun = "unsupported_rate")
    report <- getSteadyAssessment(p, tol = 1e6)
    expect_identical(report$verdict, NA)
    expect_identical(report$missing, "component:pool")
    expect_identical(report$unavailable$capability, "convergence")
    expect_match(report$unavailable$reason, "No convergence callback")
    expect_error(getSteadyResidual(p), "do not differentiate")
    failing <- getSteadyAssessment(p, tol = 0, atol = 0)
    expect_identical(failing$verdict, FALSE)
    expect_identical(failing$missing, "component:pool")
})

test_that("Conditional scope is explicit and retains its fixed-owner IDs", {
    p <- setComponent(NS_params_small, "resource", initial_value = 2,
                      dynamics_fun = "constant_other")
    report <- getSteadyAssessment(p, tol = 1e6, scope = "conditional",
                                  fixed = "component:resource")
    expect_true(report$verdict)
    expect_identical(report$fixed, "component:resource")
    expect_identical(report$required, c("consumers", "resource"))
    expect_error(getSteadyAssessment(p, fixed = "component:resource"),
                 "conditional")
    expect_error(getSteadyAssessment(p, scope = "conditional"), "nonempty")
    expect_error(getSteadyAssessment(p, scope = "conditional", fixed = "pool"),
                 "canonical owner")
    expect_error(getSteadyAssessment(p, scope = "conditional",
                                     fixed = c("resource", "resource")),
                 "unique canonical")
})

test_that("Prescribed owners are recorded separately from held-fixed owners", {
    p <- setComponent(NS_params_small, "input", initial_value = list(value = 2),
                      dynamics_fun = "constant_other", role = "prescribed")
    resource_dynamics(p) <- "resource_constant"
    report <- getSteadyAssessment(p, tol = 1e6)
    expect_true(report$verdict)
    expect_setequal(report$prescribed, c("resource", "component:input"))
    expect_identical(report$fixed, character())
    expect_identical(report$required, "consumers")

    p <- setProcess(p, "consumers", role = "prescribed", convergence_fun = NULL)
    report <- getSteadyAssessment(p)
    expect_true(report$verdict)
    expect_identical(report$required, character())
    expect_identical(report$covered, character())
    expect_equal(nrow(report$measures), 0)
    expect_equal(assessment_ratio(report), 0)
})

test_that("State norms prevent cancellation under both quadrature schemes", {
    for (second in c(FALSE, TRUE)) {
        p <- NS_params_small
        p$second_order_w[["bin_average"]] <- second
        derivative <- p$initial_n * 0
        weight <- biomass_bin_weight(p)
        derivative[, 2] <- 1 / weight[[2]]
        derivative[, 3] <- -1 / weight[[3]]
        state <- builtin_consumer_assessment(
            p, n = p$initial_n,
            derivative = derivative, tol = 0, atol = 1, criterion = "state"
        )
        biomass <- builtin_consumer_assessment(
            p, n = p$initial_n,
            derivative = derivative, tol = 0, atol = 1, criterion = "biomass"
        )
        expect_equal(state$value, rep(2, nrow(p$initial_n)))
        expect_equal(biomass$value, rep(0, nrow(p$initial_n)),
                     tolerance = 1e-14)
        expect_true(all(state$value > state$tolerance))
        expect_true(all(biomass$value <= biomass$tolerance))
        scaled <- builtin_consumer_assessment(
            p, n = p$initial_n,
            derivative = derivative, tol = 0.2, atol = 1, criterion = "state"
        )
        expect_equal(scaled$tolerance,
                     1 + 0.2 * as.numeric(sizeIntegral(p, weighting = p$w)))
    }
    for (second in c(FALSE, TRUE)) {
        p <- NS_params_small
        p$second_order_w[["bin_average"]] <- second
        derivative <- p$initial_n_pp * 0
        weight <- bin_average_weight(p$w_full, p) * p$dw_full
        derivative[2:3] <- c(1, -1) / weight[2:3]
        state <- builtin_resource_assessment(
            p,
            n_pp = p$initial_n_pp, derivative = derivative,
            tol = 0, atol = 1, criterion = "state"
        )
        biomass <- builtin_resource_assessment(
            p,
            n_pp = p$initial_n_pp, derivative = derivative,
            tol = 0, atol = 1, criterion = "biomass"
        )
        expect_equal(state$value, 2)
        expect_equal(biomass$value, 0)
    }
})

test_that("Callbacks receive one shared rate evaluation and actual time", {
    calls <- 0L
    rlang::local_bindings(
        assessment_rates = function(params, ...) {
            calls <<- calls + 1L
            mizerRates(params, ...)
        },
        assessment_equation = function(params, n_other, component, t, ...) {
            t * n_other[[component]]
        },
        assessment_measure = function(p, n_other, component, rates, t,
                                      derivative) {
            stopifnot(inherits(p, "MizerParams"), !is.null(rates$f_mort))
            stopifnot(identical(derivative, t * n_other[[component]]))
            data.frame(measure = "signed_change", value = -sum(derivative),
                       tolerance = 5, units = "units/year")
        }, .env = .GlobalEnv
    )
    p <- setComponent(NS_params_small, "pool", initial_value = c(a = 1, b = 2),
                      rate_fun = "assessment_equation",
                      convergence_fun = "assessment_measure")
    p$rates_funcs$Rates <- "assessment_rates"
    report <- getSteadyAssessment(p, tol = 1e6, t = 2)
    expect_identical(calls, 1L)
    expect_equal(report$time, 2)
    row <- report$measures[report$measures$owner == "component:pool", ]
    expect_equal(row$value, -6)
    expect_equal(row$tolerance, 5)
    expect_false(row$passed)
    expect_false(report$verdict)
    expect_error(getSteadyAssessment(p, t = NA_real_), "finite number")
})

test_that("Extension tolerances are independent and may be overridden by ID", {
    rlang::local_bindings(
        own_tolerances = function(...) {
            data.frame(measure = c("rate", "imbalance"), value = c(-2, 0),
                       tolerance = c(1, 0), units = c("units/year", "units"))
        }, .env = .GlobalEnv
    )
    p <- setComponent(NS_params_small, "pool", initial_value = 0,
                      dynamics_fun = "constant_other",
                      convergence_fun = "own_tolerances")
    report <- getSteadyAssessment(p, tol = 1e6)
    expect_false(report$verdict)
    expect_equal(assessment_ratio(report), 2)
    report <- getSteadyAssessment(p, tol = 1e6,
                                  tolerances = c("component:pool/rate" = 2))
    expect_true(report$verdict)
    expect_equal(assessment_ratio(report), 1)
    expect_error(getSteadyAssessment(p, tolerances = c("pool/rate" = 2)),
                 "Unknown assessment tolerance")
    expect_error(
        getSteadyAssessment(p, tolerances = c("component:pool/rate" = -1)),
        "nonnegative"
    )
    expect_error(getSteadyAssessment(p, tol = Inf), "finite nonnegative")
})

test_that("Malformed callback results and programming errors propagate", {
    rlang::local_bindings(
        broken_assessment = function(...) stop("intentional failure"),
        malformed_assessment = function(...) {
            data.frame(measure = "x", value = NA_real_, tolerance = 0,
                       units = "units")
        }, .env = .GlobalEnv
    )
    p <- setComponent(NS_params_small, "pool", initial_value = 2,
                      dynamics_fun = "constant_other",
                      convergence_fun = "broken_assessment")
    expect_error(getSteadyAssessment(p),
                 "assessment for 'component:pool' failed: intentional")
    p <- setComponent(p, "pool", convergence_fun = "malformed_assessment")
    expect_error(getSteadyAssessment(p), "component:pool.*finite numeric")
    p$capabilities$components$pool$convergence_fun <- "missing_assessment"
    expect_error(getSteadyAssessment(p), "component:pool.*not found")
})

test_that("A custom resource assessment need not probe its step", {
    rlang::local_bindings(
        unprobed_resource = function(...) stop("do not step resource"),
        resource_measure = function(p, process, derivative, rates, t, dt) {
            stopifnot(inherits(p, "MizerParams"), process == "resource",
                      is.null(derivative), !is.null(rates$f_mort), dt > 0)
            data.frame(measure = "time", value = t, tolerance = 1,
                       units = "year")
        }, .env = .GlobalEnv
    )
    p <- NS_params_small
    resource_dynamics(p) <- "unprobed_resource"
    missing <- getSteadyAssessment(p, tol = 1e6)
    expect_identical(missing$missing, "resource")
    p <- setProcess(p, "resource", convergence_fun = "resource_measure")
    report <- getSteadyAssessment(p, tol = 1e6, t = 2)
    expect_false(report$verdict)
    expect_identical(report$missing, character())
})

test_that("Prepared assessments accept candidate state", {
    p <- NS_params_small
    prepared <- prepare_steady_assessment(p, tol = 1e6)
    n <- p$initial_n * 2
    report <- steady_assessment(p, prepared = prepared, n = n, t = 4)
    candidate <- p
    candidate$initial_n <- n
    direct <- getSteadyAssessment(candidate, tol = 1e6, t = 4)
    expect_equal(report, direct)
    expect_equal(p$initial_n, NS_params_small$initial_n)
    resource_dynamics(p) <- "resource_constant"
    expect_error(steady_assessment(p, prepared = prepared), "providers changed")
})


test_that("Prepared equations and built-in resource rates are reused", {
    rlang::local_bindings(
        reusable_equation = function(params, n_other, component, ...) {
            n_other[[component]] * 0
        },
        reusable_assessment = function(derivative, ...) {
            stopifnot(!is.null(derivative))
            data.frame(measure = "derivative", value = 0,
                       tolerance = 0, units = "units/year")
        }, .env = .GlobalEnv
    )
    p <- setComponent(NS_params_small, "pool", initial_value = 2,
                      rate_fun = "reusable_equation",
                      convergence_fun = "reusable_assessment")
    p <- setProcess(p, "resource", convergence_fun = "reusable_assessment")
    prepared <- prepare_steady_assessment(p, tol = 1e6)
    rm("reusable_equation", envir = .GlobalEnv)
    report <- steady_assessment(p, prepared = prepared)
    expect_true(report$verdict)
    expect_identical(report$missing, character())
    expect_error(getSteadyAssessment(p),
                 "Equation for 'component:pool'.*not found")
})

test_that("Quantitative measure tables reject malformed columns", {
    valid <- data.frame(measure = "rate", value = -1, tolerance = 2,
                        units = "1/year")
    invalid <- list(valid[FALSE, ], rbind(valid, valid))
    bad <- valid
    bad$value <- NA_real_
    invalid[[3]] <- bad
    bad$value <- matrix(c(1, 2), nrow = 1)
    invalid[[4]] <- bad
    bad <- valid
    bad$tolerance <- -1
    invalid[[5]] <- bad
    bad <- valid
    bad$units <- ""
    invalid[[6]] <- bad
    for (value in invalid) {
        expect_error(checked_assessment_measures(value, "component:test"),
                     "component:test.*must return")
    }
    expect_true(checked_assessment_measures(valid, "component:test")$passed)
})

test_that("Custom projection providers do not receive unrelated derivatives", {
    rlang::local_bindings(
        project.assessment_test_subclass = function(params, ...) params,
        external_equation_assessment = function(params, derivative, ...) {
            stopifnot(is.null(derivative))
            data.frame(measure = "own_equation", value = 0,
                       tolerance = 0, units = "units/year")
        }, .env = .GlobalEnv
    )
    p <- NS_params_small
    p$extensions <- c(assessment_test_subclass = NA_character_)
    p <- coerceToExtensionClass(p)
    p <- setProcess(p, "consumers",
                    convergence_fun = "external_equation_assessment")
    p <- setProcess(p, "resource",
                    convergence_fun = "external_equation_assessment")
    report <- getSteadyAssessment(p)
    expect_true(report$verdict)
    expect_identical(report$covered, c("consumers", "resource"))
    expect_true(is.na(attr(report, "drift_report")$drift))
    p <- setProcess(p, "consumers",
                    convergence_fun = "builtin_consumer_assessment")
    expect_error(getSteadyAssessment(p), "No built-in consumer derivative")
})
