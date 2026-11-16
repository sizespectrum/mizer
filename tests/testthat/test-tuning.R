test_that("Built-in tuning separates fixed inputs and evolving consumers", {
    plan <- prepare_tuning(NS_params_small)
    expect_identical(plan$hold_states, "resource")
    expect_identical(plan$evolve_states, "consumers")
    expect_true(plan$hold_recruitment)
    expect_equal(plan$recruitment, tuning_rate_context(NS_params_small)$rdd)
    expect_named(plan$verification, c("consumers", "resource", "reproduction"))
})

test_that("Declared components retain unrelated built-in tuning support", {
    rlang::local_bindings(
        tuning_zero_rate = function(n_other, component,
                                    ...) {
            n_other[[component]] * 0
        },
        tuning_pool_rates = function(params, n_other, component = "pool",
                                     ...) {
            list(change = n_other[[component]] * 0)
        },
        .env = globalenv()
    )
    strategy <- tuning_strategy_defaults()
    strategy$evolve <- "component:pool"
    p <- setComponent(NS_params_small, "pool", 1,
        rate_fun = "tuning_zero_rate",
        stationary = TRUE, equation_rates_fun = "tuning_pool_rates",
        tuning = list(default = "evolve", strategies = list(evolve = strategy))
    )
    plan <- prepare_tuning(p)
    expect_setequal(plan$evolve_states, c("consumers", "component:pool"))
    expect_error(prepare_tuning(p, solver = "newton"), "cannot evolve")
    targets <- capture_tuning_targets(p, plan)
    expect_named(targets$equations$`component:pool`, "change")
})

test_that("Prescribed role alone does not claim temporal stationarity", {
    p <- setComponent(NS_params_small, "input", 1,
        dynamics_fun = "constant_other", role = "prescribed"
    )
    expect_error(prepare_tuning(p), "stationary.*component:input")
    p <- setComponent(p, "input", stationary = TRUE)
    plan <- prepare_tuning(p)
    expect_setequal(plan$hold_states, c("resource", "component:input"))
})

test_that("Effective rate replacements invalidate built-in tuning", {
    rlang::local_bindings(tuning_custom_rdi = function(params, ...) {
        rep(
            1,
            nrow(params$species_params)
        )
    }, .env = globalenv())
    p <- NS_params_small
    p$rates_funcs$RDI <- "tuning_custom_rdi"
    expect_null(builtin_reproduction_law(p))
    expect_error(prepare_tuning(p), "stationary")
    q <- NS_params_small
    q$rates_funcs$RDD <- "noRDD"
    expect_error(prepare_tuning(q), "stationary")
})

test_that("Strategy fields are checked rather than silently ignored", {
    expect_error(
        normalise_tuning_strategy(list(guess = TRUE)),
        "supported fields"
    )
    invalid_hold <- list(hold = list(states = character(), rates = "rdi"))
    expect_error(normalise_tuning_strategy(invalid_hold), "supported held rate")
    expect_error(
        normalise_tuning_strategy(list(action = "balance")),
        "requires balance_fun"
    )
    invalid_settings <- list(settings = list(code = function() 1))
    expect_error(normalise_tuning_strategy(invalid_settings), "serialisable")
    invalid_adjust <- list(adjust = list(x = list(path = "cc_pp")))
    expect_error(normalise_tuning_strategy(invalid_adjust), "no-action")
})

test_that("Equation-rate invariants compare individual cells", {
    expected <- c(1e12, rep(1, 1000))
    actual <- expected
    actual[[5]] <- 2
    expect_false(tuning_values_equal(actual, expected))
    expect_true(tuning_values_equal(c(Inf, 0), c(Inf, 0)))
    expect_false(tuning_values_equal(1, Inf))
})

test_that("Independent verification rejects state and undeclared changes", {
    p <- NS_params_small
    plan <- prepare_tuning(p, owners = "reproduction", solver = "balance")
    targets <- capture_tuning_targets(p, plan)
    bad <- p
    bad$initial_n[1, 1] <- bad$initial_n[1, 1] * 1.001
    expect_error(verify_tuning_candidate(bad, plan, targets), "target state")
    bad <- p
    bad$species_params$w_mat[1] <- bad$species_params$w_mat[1] * 1.001
    expect_error(verify_tuning_candidate(bad, plan, targets), "outside.*paths")
})

test_that("Joint ownership conflicts fail before balancing", {
    p <- NS_params_small
    strategy <- tuning_strategy_defaults()
    strategy$hold$states <- "resource"
    p <- setTuningStrategy(p, "coupled",
        owners = c("resource", "reproduction"),
        definition = strategy, stationary = TRUE,
        equation_rates_fun = "builtin_prescribed_equation_rates"
    )
    expect_error(prepare_tuning(p,
                     joint = "coupled",
                     strategies = c(resource = "capacity")
                 ), "conflict")
    p <- setTuningStrategy(p, "second",
        owners = c("resource", "reproduction"),
        definition = strategy, stationary = TRUE,
        equation_rates_fun = "builtin_prescribed_equation_rates"
    )
    expect_error(prepare_tuning(p, joint = c("coupled", "second")), "overlap")
})

test_that("Balancing verifies separate equation terms despite cancellation", {
    rlang::local_bindings(
        tuning_terms = function(params, ...) {
            list(
                source = params$other_params$pool$source,
                loss = params$other_params$pool$loss
            )
        },
        tuning_bad_balance = function(params, ...) {
            params$other_params$pool$source <- 3
            params$other_params$pool$loss <- 3
            params
        }, .env = globalenv()
    )
    strategy <- tuning_strategy_defaults()
    strategy$action <- "balance"
    strategy$hold$states <- "component:pool"
    strategy$balance_fun <- "tuning_bad_balance"
    strategy$adjust <- list(
        source = list(path = c("other_params", "pool", "source"), lower = 0),
        loss = list(path = c("other_params", "pool", "loss"), lower = 0)
    )
    p <- setComponent(NS_params_small, "pool", 1,
        dynamics_fun = "constant_other",
        component_params = list(source = 2, loss = 2), stationary = TRUE,
        equation_rates_fun = "tuning_terms",
        tuning = list(
            default = "balance",
            strategies = list(balance = strategy)
        )
    )
    plan <- prepare_tuning(p, owners = "component:pool", solver = "balance")
    targets <- capture_tuning_targets(p, plan)
    expect_error(
        finish_tuning(p, plan, targets),
        "effective equation rate.*component:pool"
    )
    expect_equal(p$other_params$pool, list(source = 2, loss = 2))
})

test_that("Resource capacity targets survive dependency refresh", {
    p <- NS_params_small
    p$initial_n <- p$initial_n * 1.2
    plan <- prepare_tuning(p, owners = "resource", solver = "balance")
    targets <- capture_tuning_targets(p, plan)
    balanced <- finish_tuning(p, plan, targets)
    expect_equal(tuning_state(balanced), tuning_state(p))
    actual <- builtin_resource_equation_rates(balanced,
        n_pp = balanced$initial_n_pp, rates = tuning_rate_context(balanced)
    )
    expect_equal(actual$source, actual$loss, tolerance = 1e-10)
})

test_that("Declared bounds are enforced on returned candidates", {
    p <- NS_params_small
    plan <- prepare_tuning(p, owners = "resource", solver = "balance")
    targets <- capture_tuning_targets(p, plan)
    bad <- p
    bad$cc_pp[[1L]] <- -1
    expect_error(verify_tuning_candidate(bad, plan, targets), "bounds")
})

test_that("One joint balancer owns coupled resource and component targets", {
    calls <- new.env(parent = emptyenv())
    calls$joint <- 0L
    rlang::local_bindings(
        unused_local_balancer = function(...) stop("local balancer was called"),
        joint_pool_rates = function(params, n_pp, rates, ...) {
            bin <- params$other_params$pool$bin
            list(resource = builtin_resource_equation_rates(params, n_pp,
                                                            rates),
                 `component:pool` = list(source = params$cc_pp[[bin]],
                                         loss = params$other_params$pool$loss))
        },
        joint_pool_targets = function(params, n_pp, rates, ...) {
            resource <- target_resource_strategy(params, n_pp, rates)
            bin <- params$other_params$pool$bin
            capacity <- resource$source[[bin]] / params$rr_pp[[bin]]
            list(resource = resource,
                 `component:pool` = list(source = capacity, loss = capacity))
        },
        balance_joint_pool = function(params, n_pp, rates, ...) {
            calls$joint <- calls$joint + 1L
            params <- balance_resource_strategy(params, rates, n_pp)
            bin <- params$other_params$pool$bin
            params$other_params$pool$loss <- params$cc_pp[[bin]]
            params
        }, .env = globalenv()
    )
    p <- NS_params_small
    p$initial_n <- p$initial_n * 1.2
    bin <- which(p$rr_pp > 0 & p$initial_n_pp > 0)[[1L]]
    local <- tuning_strategy_defaults()
    local$action <- "balance"
    local$hold$states <- "component:pool"
    local$balance_fun <- "unused_local_balancer"
    local$adjust <- list(capacity = list(path = "cc_pp", lower = 0))
    p <- setComponent(
        p, "pool", 1, dynamics_fun = "constant_other",
        component_params = list(bin = bin, loss = p$cc_pp[[bin]]),
        stationary = TRUE,
        tuning = list(default = "local", strategies = list(local = local))
    )
    resource <- p$capabilities$processes$resource$tuning
    resource$strategies$capacity$balance_fun <- "unused_local_balancer"
    p <- setProcess(p, "resource", tuning = resource)
    joint <- local
    joint$hold$states <- c("resource", "component:pool")
    joint$balance_fun <- "balance_joint_pool"
    joint$target_fun <- "joint_pool_targets"
    joint$adjust$loss <- list(path = c("other_params", "pool", "loss"),
                              lower = 0)
    p <- setTuningStrategy(p, "coupled",
                           owners = c("resource", "component:pool"),
                           definition = joint,
                           equation_rates_fun = "joint_pool_rates",
                           stationary = TRUE,
                           default = TRUE)
    standalone <- prepare_tuning(p, owners = "reproduction", solver = "balance")
    expect_named(standalone$verification_joint, "coupled")
    expect_no_error(balanceRDD(p))
    expect_identical(calls$joint, 0L)
    plan <- prepare_tuning(p)
    expect_false("resource" %in% names(plan$selected))
    expect_false("component:pool" %in% names(plan$selected))
    targets <- capture_tuning_targets(p, plan)
    balanced <- finish_tuning(p, plan, targets)
    expect_identical(calls$joint, 1L)
    expect_equal(tuning_state(balanced), tuning_state(p))
    expect_equal(balanced$other_params$pool$loss, balanced$cc_pp[[bin]])
    expect_equal(tuning_equation_rates(balanced, plan,
                                       tuning_rate_context(balanced), t = 0),
                 targets$equations)
})

test_that("Provider metadata and unrelated components do not stale support", {
    p <- NS_params_small
    comment(p$rates_funcs) <- "a model annotation"
    expect_no_error(prepare_tuning(p))
    p <- setComponent(p, "first", 1, dynamics_fun = "constant_other",
                      role = "prescribed", stationary = TRUE)
    p <- setComponent(p, "second", 1, dynamics_fun = "constant_other",
                      role = "prescribed", stationary = TRUE)
    expect_no_error(plan <- prepare_tuning(p))
    expect_setequal(plan$hold_states,
                    c("resource", "component:first", "component:second"))
})

test_that("Undeclared frozen arrays cannot hide behind dependency refresh", {
    p <- NS_params_small
    plan <- prepare_tuning(p, owners = "reproduction", solver = "balance")
    targets <- capture_tuning_targets(p, plan)
    bad <- p
    bad$mu_b[1, 1] <- bad$mu_b[1, 1] + 1
    comment(bad$mu_b) <- "set manually"
    expect_error(verify_tuning_candidate(bad, plan, targets), "outside.*paths")
    bad <- p
    comment(bad$mu_b) <- "set manually"
    expect_error(verify_tuning_candidate(bad, plan, targets), "outside.*paths")
})

test_that("Resource parameter adjustments refresh unfrozen dependencies", {
    p <- NS_params_small
    comment(p$cc_pp) <- NULL
    plan <- prepare_tuning(p, owners = "reproduction", solver = "balance")
    plan$selected$reproduction$strategy$adjust$kappa <- list(
        path = c("resource_params", "kappa"), lower = 0,
        upper = Inf, species = FALSE
    )
    changed <- p
    changed$resource_params$kappa <- 2 * p$resource_params$kappa
    refreshed <- refresh_tuning_dependencies(changed, plan, p)
    expect_equal(unname(refreshed$cc_pp), 2 * p$resource_params$kappa *
                     p$w_full ^ -p$resource_params$lambda *
                     (p$w_full <= p$resource_params$w_pp_cutoff))
    targets <- capture_tuning_targets(p, plan)
    expect_error(verify_tuning_candidate(changed, plan, targets),
                 "effective equation rate.*resource")

    comment(p$cc_pp) <- "set manually"
    changed <- p
    changed$resource_params$kappa <- 2 * p$resource_params$kappa
    refreshed <- refresh_tuning_dependencies(changed, plan, p)
    expect_equal(refreshed$cc_pp, p$cc_pp)
})

test_that("Author defaults survive unrelated provider updates", {
    p <- NS_params_small
    tuning <- p$capabilities$processes$reproduction$tuning
    tuning$default <- "preserve_erepro"
    p <- setProcess(p, "reproduction", tuning = tuning)
    p <- setResource(p, resource_dynamics = "resource_logistic")
    p <- setComponent(p, "input", 1, dynamics_fun = "constant_other",
                      role = "prescribed", stationary = TRUE)
    expect_identical(p$capabilities$processes$reproduction$tuning, tuning)
    expect_identical(prepare_tuning(p)$selected$reproduction$id,
                     "preserve_erepro")
    changed <- setReproduction(p, RDD = "RickerRDD")
    expect_identical(changed$capabilities$processes$reproduction$tuning, tuning)
    expect_error(prepare_tuning(changed), "stationary")
})

test_that("Invalid preservation requests and structural adjustments stop", {
    expect_error(prepare_tuning(NS_params_small,
                                preserve = c("erepro", "R_max")),
                 "one supported quantity")
    p <- NS_params_small
    tuning <- p$capabilities$processes$reproduction$tuning
    tuning$strategies$preserve_erepro$adjust$grid <- list(path = "dw")
    expect_error(setProcess(p, "reproduction", tuning = tuning),
                 "grid structure")
    tuning <- p$capabilities$processes$resource$tuning
    tuning$strategies$capacity$hold$rates <- "consumers/recruitment"
    p <- setProcess(p, "resource", tuning = tuning)
    expect_error(prepare_tuning(p), "requires reproduction ownership")
})

test_that("Partial balancing preserves untouched supplied metadata", {
    p <- NS_params_small
    selected <- as.character(p$species_params$species[[1L]])
    plan <- prepare_tuning(p, owners = "reproduction", solver = "balance",
                           species = selected,
                           strategies = c(reproduction = "set_parameters"))
    targets <- capture_tuning_targets(p, plan)
    bad <- p
    bad$given_species_params$erepro <- rep(NA_real_, nrow(p$species_params))
    bad$given_species_params$erepro[[2L]] <- 0.123
    expect_error(verify_tuning_candidate(bad, plan, targets),
                 "unselected given parameters")
    p$other_params$note <- list(amount = 1)
    plan <- prepare_tuning(p, owners = "reproduction", solver = "balance")
    targets <- capture_tuning_targets(p, plan)
    bad <- p
    names(bad$other_params$note) <- "renamed"
    expect_error(verify_tuning_candidate(bad, plan, targets), "outside.*paths")
})

test_that("Extension parameter changes refresh dependent cached rates", {
    rlang::local_bindings(
        setParams.tuning_cache = function(object, ...) {
            params <- NextMethod()
            scale <- params$other_params$pool$scale
            params$search_vol <- params$search_vol * scale
            params
        },
        tuning_cached_rates = function(params, ...) {
            list(search = params$search_vol)
        },
        tuning_stale_balance = function(params, ...) {
            params$other_params$pool$scale <- 2
            params
        }, .env = globalenv()
    )
    p <- NS_params_small
    p$extensions <- c(tuning_cache = NA_character_)
    strategy <- tuning_strategy_defaults()
    strategy$action <- "balance"
    strategy$hold$states <- "component:pool"
    strategy$balance_fun <- "tuning_stale_balance"
    strategy$adjust <- list(scale = list(
                                         path = c("other_params", "pool",
                                                  "scale")))
    p <- setComponent(
        p, "pool", 1, dynamics_fun = "constant_other",
        component_params = list(scale = 1), stationary = TRUE,
        equation_rates_fun = "tuning_cached_rates",
        tuning = list(default = "bad", strategies = list(bad = strategy))
    )
    plan <- prepare_tuning(p, owners = "component:pool", solver = "balance")
    targets <- capture_tuning_targets(p, plan)
    expect_error(finish_tuning(p, plan, targets), "effective equation rate")
})
