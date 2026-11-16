test_that("Structural support includes absent species and zero resource bins", {
    p <- NS_params_small
    p$initial_n[1, ] <- 0
    p$initial_n_pp[] <- 0
    p$rr_pp[] <- 0
    for (scheme in c(FALSE, TRUE)) {
        p$second_order_w <- list(flux = if (scheme) "van_leer" else "upwind",
                                 bin_average = scheme)
        state <- prepare_state_capabilities(p)
        expected <- col(p$initial_n) >= p$w_min_idx &
            col(p$initial_n) <= support_top_idx(p)
        dimnames(expected) <- dimnames(p$initial_n)
        expect_equal(state$consumers$structural_mask, expected)
        expect_true(any(state$consumers$dynamic_mask[1, ]))
        expect_true(all(state$resource$dynamic_mask))
        expect_true(state$consumers$autonomous)
        expect_true(state$resource$autonomous)
    }
})

test_that("Partial prescribed masks have exact numeric state geometry", {
    rlang::local_bindings(
        state_test_mask = function(params, component) {
            params$other_params[[component]]$prescribed
        },
        state_bad_mask = function(...) TRUE,
        .env = globalenv()
    )
    mask <- c(a = FALSE, b = TRUE)
    p <- setComponent(NS_params_small, "pool", c(a = 1, b = 2),
                      dynamics_fun = "constant_other",
                      prescribed_mask_fun = "state_test_mask",
                      component_params = list(prescribed = mask),
                      autonomous = TRUE)
    plan <- prepare_state_capabilities(p)[["component:pool"]]
    expect_identical(plan$prescribed_mask, c(a = FALSE, b = TRUE))
    expect_identical(plan$dynamic_mask, c(a = TRUE, b = FALSE))
    expect_true(plan$mask_current)
    expect_true(plan$autonomous)
    bad <- setComponent(p, "pool", prescribed_mask_fun = "state_bad_mask")
    expect_error(prepare_state_capabilities(bad), "exactly match")
    bad$capabilities$components$pool$prescribed_mask_fun <- "missing_mask"
    expect_no_error(prepare_state_capabilities(bad, fixed = "component:pool"))
    expect_error(prepare_state_capabilities(bad), "missing_mask")
    nonnumeric <- setComponent(p, "pool", initial_value = list(stock = 1))
    expect_error(prepare_state_capabilities(nonnumeric), "numeric state")
    prescribed <- setComponent(nonnumeric, "pool", role = "prescribed",
                               prescribed_mask_fun = NULL)
    plan <- prepare_state_capabilities(prescribed)[["component:pool"]]
    expect_true(plan$prescribed_mask)
})

test_that("Temporal declarations are independent of instantaneous zeros", {
    rlang::local_bindings(
        state_forced_rhs = function(n_other, component, t, ...) {
            sin(t) * n_other[[component]]
        },
        state_forced_mort = function(params, t, ...) {
            params$ext_mort * sin(t)
        }, .env = globalenv()
    )
    p <- setComponent(NS_params_small, "pool", 0,
                      rate_fun = "state_forced_rhs", stationary = TRUE,
                      autonomous = FALSE)
    spec <- prepare_state_capabilities(p)[["component:pool"]]
    expect_false(spec$autonomous)
    expect_match(spec$temporal_reason, "nonautonomous")
    p <- setComponent(p, "pool", mort_fun = "state_forced_mort",
                      autonomous = FALSE)
    fixed <- prepare_state_capabilities(p, fixed = "component:pool")
    expect_true(fixed[["component:pool"]]$temporal_required)
    expect_false(fixed[["component:pool"]]$autonomous)
    expect_false(any(fixed[["component:pool"]]$dynamic_mask))
    unset <- setComponent(p, "pool", autonomous = NULL)
    unset_state <- prepare_state_capabilities(unset)[["component:pool"]]
    expect_null(unset_state$autonomous)
})

test_that("Exact resource equations use effective rates without time probes", {
    p <- NS_params_small
    n_pp <- p$initial_n_pp
    rates <- list(resource_mort = n_pp * 0 + 0.2)
    for (law in c("semichemostat", "logistic", "constant")) {
        p <- setResource(p, resource_dynamics = paste0("resource_", law),
                         balance = FALSE)
        spec <- prepare_state_capabilities(p)$resource
        actual <- evaluate_resource_equation(p, spec, p$initial_n, n_pp,
                                             p$initial_n_other, rates, 3)
        expected <- switch(
            law,
            semichemostat = p$rr_pp * (p$cc_pp - n_pp) - 0.2 * n_pp,
            logistic = {
                value <- n_pp * 0
                positive <- n_pp != 0
                growth <- p$rr_pp[positive] *
                    (1 - n_pp[positive] / p$cc_pp[positive])
                value[positive] <- (growth - 0.2) * n_pp[positive]
                value
            },
            constant = n_pp * 0
        )
        expect_equal(actual, expected, ignore_attr = TRUE)
        expect_true(spec$equation_current)
    }
})

test_that("Custom resource equations require independent current adoption", {
    rlang::local_bindings(
        state_custom_step = function(n_pp, ...) n_pp,
        state_custom_rhs = function(n_pp, t, ...) n_pp * t,
        state_bad_rhs = function(...) 1,
        .env = globalenv()
    )
    p <- setResource(NS_params_small,
                     resource_dynamics = "state_custom_step", balance = FALSE)
    expect_false(prepare_state_capabilities(p)$resource$equation_current)
    q <- setProcess(p, "resource", rate_fun = "state_custom_rhs")
    spec <- prepare_state_capabilities(q)$resource
    expect_true(spec$equation_current)
    expect_null(spec$autonomous)
    expect_equal(evaluate_resource_equation(q, spec, q$initial_n,
                                            q$initial_n_pp, list(), list(), 2),
                 q$initial_n_pp * 2)
    bad <- setProcess(q, "resource", rate_fun = "state_bad_rhs")
    bad_spec <- prepare_state_capabilities(bad)$resource
    expect_error(
        evaluate_resource_equation(bad, bad_spec, bad$initial_n,
                                   bad$initial_n_pp, list(), list(), 0),
        "exact resource state"
    )
    q$resource_dynamics <- "resource_semichemostat"
    expect_false(prepare_state_capabilities(q)$resource$equation_current)
})


test_that("Zero logistic capacity is a domain limit, not a structural mask", {
    p <- setResource(NS_params_small, resource_dynamics = "resource_logistic",
                     balance = FALSE)
    p$cc_pp[1] <- 0
    p$initial_n_pp[1] <- 0
    p$rr_pp[1] <- 1
    spec <- prepare_state_capabilities(p)$resource
    expect_true(spec$dynamic_mask[1])
    expect_match(resource_equation_domain(p, spec, spec$dynamic_mask),
                 "zero-capacity")
    n_pp <- p$initial_n_pp
    n_pp[1] <- 1e-6
    rates <- list(resource_mort = n_pp * 0 + 0.3)
    expect_error(evaluate_resource_equation(p, spec, p$initial_n, n_pp,
                                            list(), rates, 0),
                 "positive capacity")
    p$rr_pp[p$cc_pp <= 0] <- 0
    expect_null(resource_equation_domain(p, spec, spec$dynamic_mask))
    derivative <- evaluate_resource_equation(p, spec, p$initial_n, n_pp,
                                             list(), rates, 0)
    expect_equal(unname(derivative[1]), -0.3e-6)
})
