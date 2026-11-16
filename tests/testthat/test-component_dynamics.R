test_that("Exported callbacks resolve without attached packages", {
    expect_identical(resolve_model_callback("base::identity"), base::identity)
    expect_identical(resolve_model_callback("constant_other"), constant_other)
    expect_error(resolve_model_callback("base::missing_mizer_callback"),
                 "was not found")
    expect_error(
        resolve_model_callback("unavailableMizerTestPkg::step", "Pool"),
        "Pool requires package"
    )
    expect_error(resolve_model_callback("base:::identity"), "exported_function")
    expect_error(resolve_model_callback("base::"), "exported_function")
    expect_error(resolve_model_callback(NA_character_), "function name")
})

test_that("Prediction arguments are optional without hiding callback errors", {
    strict <- function(params, dt) params + dt
    expect_equal(invoke_model_step(strict, list(params = 2, dt = 3),
                                   predicted_state = list(n = 9)), 5)
    flexible <- function(params, predicted_state = NULL,
                         endpoint_rates = NULL) {
        list(predicted_state, endpoint_rates)
    }
    expect_identical(invoke_model_step(flexible, list(2)), list(NULL, NULL))
    prediction <- list(n = 9)
    rates <- list(mort = 2)
    expect_identical(
        invoke_model_step(flexible, list(params = 2), prediction, rates),
        list(prediction, rates)
    )
    broken <- function(params, ...) stop("callback failed")
    expect_error(invoke_model_step(broken, list(params = 2)), "callback failed")
})

test_that("Numeric component equations retain their state shape", {
    p <- NS_params_small
    state <- matrix(1:4, 2, dimnames = list(c("a", "b"), c("small", "large")))
    derivative <- function(params, n_other, component, ...) {
        -n_other[[component]]
    }
    evaluate <- function(fun, value = state) {
        checked_component_rate(fun, p, p$initial_n, p$initial_n_pp,
                               list(pool = value), list(), 0, "pool")
    }
    expect_identical(evaluate(derivative), -state)
    expect_error(evaluate(function(...) 1), "dimensions and names")
    expect_error(evaluate(function(...) unname(-state)), "dimensions and names")
    expect_error(evaluate(function(...) state * NA_real_), "finite numeric")
    expect_error(evaluate(derivative, list(value = 1)), "finite numeric state")
})

test_that("A temporary step override takes precedence over a stored equation", {
    assign("component_decay_probe", function(params, n_other, component, ...) {
        -n_other[[component]]
    }, .GlobalEnv)
    withr::defer(rm("component_decay_probe", envir = .GlobalEnv))
    p <- setComponent(NS_params_small, "pool", 2,
                      rate_fun = "component_decay_probe")
    prepared <- prepare_component_dynamics(p)[["pool"]]
    expect_true(is.function(attr(prepared, "component_rate")))
    expect_identical(component_integrators(p, "euler")$order, 1L)
    expect_identical(component_integrators(p, "tr_bdf2")$method, "heun")
    p$other_dynamics$pool <- "constant_other"
    expect_null(prepare_component_rates(p)[["pool"]])
    expect_null(attr(prepare_component_dynamics(p)[["pool"]], "component_rate"))
    expect_identical(component_integrators(p, "tr_bdf2")$method, "custom")
    expect_true(is.na(component_integrators(p, "tr_bdf2")$order))
})

test_that("An empty integration report retains its columns", {
    report <- component_integrators(NS_params_small, "euler")
    expect_named(report, c("component", "method", "order", "callback"))
    expect_equal(nrow(report), 0)
})

test_that("Step adapters retain the positional first model argument", {
    legacy <- function(p, dt, ...) p + dt
    expect_equal(invoke_model_step(legacy, list(params = 2, dt = 3)), 5)
    strict <- function(object, dt) object + dt
    expect_equal(invoke_model_step(strict, list(params = 2, dt = 3)), 5)

    rlang::local_bindings(
        positional_component_step = function(p, n_other, component, dt, ...) {
            stopifnot(inherits(p, "MizerParams"))
            n_other[[component]] + dt
        },
        positional_resource_step = function(p, n_pp, ...) {
            stopifnot(inherits(p, "MizerParams"))
            n_pp
        }, .env = .GlobalEnv
    )
    p <- setComponent(NS_params_small, "pool", initial_value = 2,
                      dynamics_fun = "positional_component_step")
    p@resource_dynamics <- "positional_resource_step"
    sim <- project(p, t_max = 0.1, t_save = 0.1, dt = 0.1,
                   progress_bar = FALSE)
    expect_equal(finalNOther(sim)$pool, 2.1)
    expect_equal(as.numeric(finalNResource(sim)), as.numeric(p@initial_n_pp))
    expect_equal(steady_rates(p)$other$pool, 0.5)
})

test_that("Equation and step programming errors retain provider context", {
    bad <- function(...) stop("deliberate callback error")
    p <- NS_params_small
    expect_error(checked_component_rate(bad, p, p$initial_n, p$initial_n_pp,
                                        list(pool = 1), list(), 0, "pool"),
                 "Equation for component 'pool' failed: deliberate callback")
    expect_error(invoke_model_step(bad, list(params = p, component = "pool")),
                 "Step for component 'pool' failed: deliberate callback")
    expect_error(invoke_model_step(bad, list(params = p)),
                 "Resource step failed: deliberate callback")
})
