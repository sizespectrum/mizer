test_that("Built-in processes explicitly support each transformation", {
    p <- NS_params_small
    for (operation in transformation_operations()) {
        plan <- prepare_transformation(p, operation)
        expect_setequal(
            names(plan$selected),
            c("consumers", "resource", "reproduction")
        )
        expect_equal(apply_transformation_hooks(p, plan), p)
    }
    expect_error(prepare_transformation(p, "invented"), "Unknown")
    plan <- prepare_transformation(p, "size_grid")
    candidate <- p
    candidate$capabilities$processes$resource$transformations <- NULL
    expect_error(apply_transformation_hooks(candidate, plan), "lost.*providers")
})

test_that("Missing component support differs from explicit no-action", {
    p <- setComponent(NS_params_small, "pool", list(stock = 2, label = "warm"),
        dynamics_fun = "constant_other"
    )
    expect_error(prepare_transformation(p, "scale_abundance"), "component:pool")
    p <- setComponent(p, "pool", transformations = list(
        scale_abundance = list(action = "none")
    ))
    plan <- prepare_transformation(p, "scale_abundance", list(factor = 2))
    candidate <- p
    candidate$initial_n <- 2 * candidate$initial_n
    result <- apply_transformation_hooks(candidate, plan)
    expect_equal(result$initial_n_other$pool, p$initial_n_other$pool)
    expect_error(prepare_transformation(p, "size_grid"), "component:pool")
})

test_that("Authored hooks selectively transform component state", {
    rlang::local_bindings(
        transform_stock = function(params, before, component, details, ...) {
            params$initial_n_other[[component]]$stock <-
                before$initial_n_other[[component]]$stock * details$factor
            params
        }, .env = globalenv()
    )
    p <- setComponent(NS_params_small, "pool",
        list(stock = 2, temperature = 10, label = "warm"),
        dynamics_fun = "constant_other", transformations = list(
            scale_abundance = list(
                action = "update", transform_fun = "transform_stock"
            )
        )
    )
    plan <- prepare_transformation(p, "scale_abundance", list(factor = 3))
    result <- apply_transformation_hooks(p, plan)
    expect_equal(
        result$initial_n_other$pool,
        list(stock = 6, temperature = 10, label = "warm")
    )
    expect_equal(p$initial_n_other$pool$stock, 2)
    expect_true("transform_stock" %in% capability_callback_names(p))
})

test_that("Preflight resolves all hooks and rejects undeclared writes", {
    rlang::local_bindings(
        transform_bad = function(params, ...) {
            params$initial_n <- 2 * params$initial_n
            params
        }, .env = globalenv()
    )
    p <- setComponent(NS_params_small, "pool", 1,
        dynamics_fun = "constant_other", transformations = list(
            size_grid = list(action = "update", transform_fun = "missing_hook")
        )
    )
    expect_error(prepare_transformation(p, "size_grid"), "missing_hook")
    p <- setComponent(p, "pool", transformations = list(
        size_grid = list(action = "update", transform_fun = "transform_bad")
    ))
    plan <- prepare_transformation(p, "size_grid")
    expect_error(apply_transformation_hooks(p, plan), "outside.*paths")
    expect_equal(p$initial_n_other$pool, 1)
    p <- setComponent(p, "pool", transformations = list(size_grid = list(
        action = "update", transform_fun = "transform_bad", paths = list("w")
    )))
    expect_error(prepare_transformation(p, "size_grid"), "responsibility")
})

test_that("Joint hooks resolve shared ownership without local double updates", {
    rlang::local_bindings(
        unused_transform = function(...) stop("local hook was called"),
        transform_shared = function(params, before, details, ...) {
            params$other_params$shared <- before$other_params$shared *
                details$factor
            params$initial_n_other$a <- params$other_params$shared
            params$initial_n_other$b <- params$other_params$shared
            params
        }, .env = globalenv()
    )
    spec <- list(
        action = "update", transform_fun = "unused_transform",
        paths = list(c("other_params", "shared"))
    )
    p <- setComponent(NS_params_small, "a", 1,
        dynamics_fun = "constant_other",
        transformations = list(scale_abundance = spec)
    )
    p <- setComponent(p, "b", 1,
        dynamics_fun = "constant_other",
        transformations = list(scale_abundance = spec)
    )
    p$other_params$shared <- 2
    expect_error(prepare_transformation(p, "scale_abundance"), "overlap")
    spec$transform_fun <- "transform_shared"
    p <- setTransformation(p, "shared", c("component:a", "component:b"),
        "scale_abundance",
        definition = spec
    )
    plan <- prepare_transformation(p, "scale_abundance", list(factor = 4))
    expect_false(any(c("component:a", "component:b") %in% names(plan$selected)))
    transformed <- apply_transformation_hooks(p, plan)
    expect_equal(transformed$other_params$shared, 8)
    expect_equal(transformed$initial_n_other, list(a = 8, b = 8))
    expect_equal(p$initial_n_other, list(a = 1, b = 1))
    expect_true("transform_shared" %in% capability_callback_names(p))
    removed <- setTransformation(p, "shared", definition = NULL)
    expect_length(removed$capabilities$transformations, 0)
})

test_that("Process hooks can update owned fields without changing layouts", {
    rlang::local_bindings(
        resource_transform = function(params, details, ...) {
            params$cc_pp <- params$cc_pp * details$factor
            params
        }, resource_bad_layout = function(params, ...) {
            params$cc_pp <- unname(params$cc_pp)
            params
        }, .env = globalenv()
    )
    spec <- list(
        action = "update", transform_fun = "resource_transform",
        paths = list("cc_pp")
    )
    p <- setProcess(NS_params_small, "resource",
        transformations = list(scale_abundance = spec)
    )
    plan <- prepare_transformation(p, "scale_abundance", list(factor = 2))
    result <- apply_transformation_hooks(p, plan)
    expect_equal(result$cc_pp, p$cc_pp * 2)
    spec$transform_fun <- "resource_bad_layout"
    p <- setProcess(
        p, "resource", transformations = list(scale_abundance = spec)
    )
    plan <- prepare_transformation(p, "scale_abundance", list(factor = 2))
    expect_error(apply_transformation_hooks(p, plan), "layout")
})

test_that("Provider changes retain explicit declarations", {
    p <- NS_params_small
    declarations <- p$capabilities$processes$resource$transformations
    p <- setProcess(p, "resource", transformations = declarations)
    p <- setResource(p, resource_dynamics = "resource_logistic")
    expect_identical(
        p$capabilities$processes$resource$transformations,
        declarations
    )
    expect_error(
        prepare_transformation(p, "scale_abundance"), "stale.*resource"
    )
    automatic <- setResource(NS_params_small,
        resource_dynamics = "resource_logistic"
    )
    expect_no_error(prepare_transformation(automatic, "scale_abundance"))
})

test_that("Transformation declarations reject ambiguous or executable claims", {
    p <- NS_params_small
    bad <- list(size_grid = list(action = "none", transform_fun = "anything"))
    expect_error(setProcess(p, "resource", transformations = bad), "no-action")
    bad <- list(size_grid = list())
    expect_error(setProcess(p, "resource", transformations = bad), "explicit")
    bad <- list(size_grid = list(
        action = "update", transform_fun = function() 1
    ))
    expect_error(setProcess(p, "resource", transformations = bad), "named")
    bad <- list(size_grid = list(
        action = "update", transform_fun = "anything",
        settings = list(bad = new.env())
    ))
    expect_error(
        setProcess(p, "resource", transformations = bad), "serialisable"
    )
})

test_that("A failing hook cannot mutate reference state", {
    rlang::local_bindings(
        transform_fail = function(params, before, component, ...) {
            params$initial_n_other[[component]]$amount <- 10
            before$initial_n_other[[component]]$amount <- 20
            stop("authored failure")
        }, .env = globalenv()
    )
    state <- new.env(parent = emptyenv())
    state$amount <- 1
    p <- setComponent(NS_params_small, "pool", state,
        dynamics_fun = "constant_other", transformations = list(
            size_grid = list(
                action = "update", transform_fun = "transform_fail"
            )
        )
    )
    plan <- prepare_transformation(p, "size_grid")
    expect_error(apply_transformation_hooks(p, plan), "authored failure")
    expect_equal(p$initial_n_other$pool$amount, 1)
    expect_equal(plan$before$initial_n_other$pool$amount, 1)
})

test_that("Custom construction providers need support before setters run", {
    calls <- new.env(parent = emptyenv())
    calls$count <- 0L
    rlang::local_bindings(
        setExtMort.construction_rho = function(params, ...) {
            calls$count <- calls$count + 1L
            stopifnot(identical(
                dim(params$other_params$rho),
                dim(params$initial_n)
            ))
            NextMethod()
        }, transform_construction_rho = function(params, details, ...) {
            params$other_params$rho <- matrix(1,
                nrow = length(details$new_species),
                ncol = length(details$new_w),
                dimnames = dimnames(params$initial_n)
            )
            params
        }, .env = globalenv()
    )
    p <- recordExtension(NS_params_small, "construction_rho")
    p$other_params$rho <- p$initial_n * 0 + 1
    expect_error(
        adjustSizeGrid(p, max_w = max(p$w) * 2),
        "transformation support.*consumers"
    )
    expect_identical(calls$count, 0L)
    p <- setProcess(p, "consumers", transformations = list(size_grid = list(
        action = "update", transform_fun = "transform_construction_rho",
        paths = list(c("other_params", "rho"))
    )))
    changed <- adjustSizeGrid(p, max_w = max(p$w) * 2)
    expect_gt(calls$count, 0L)
    expect_identical(dim(changed$other_params$rho), dim(changed$initial_n))
    expect_identical(dim(p$other_params$rho), dim(p$initial_n))
})

test_that("Only authored updates are restored after core reconstruction", {
    rlang::local_bindings(
        authored_capacity = function(params, ...) {
            params$cc_pp <- params$cc_pp * 2
            params
        }, .env = globalenv()
    )
    p <- setProcess(NS_params_small, "resource", transformations = list(
        size_grid = list(
            action = "update", transform_fun = "authored_capacity",
            paths = list("cc_pp")
        )
    ))
    plan <- prepare_transformation(p, "size_grid")
    authored <- apply_transformation_hooks(p, plan)
    updates <- capture_transformation_updates(p, authored, plan)
    expect_length(updates, 1)
    expect_identical(updates[[1L]]$path, "cc_pp")
    reconstructed <- p
    reconstructed$initial_n <- reconstructed$initial_n * 3
    restored <- restore_transformation_updates(reconstructed, updates)
    expect_equal(restored$cc_pp, authored$cc_pp)
    expect_equal(restored$initial_n, reconstructed$initial_n)
    expect_length(capture_transformation_updates(p, p, plan), 0)
    spec <- list(
        action = "update", transform_fun = "authored_capacity",
        paths = list("ft_pred_kernel_e")
    )
    p <- setProcess(p, "consumers", transformations = list(size_grid = spec))
    expect_error(prepare_transformation(p, "size_grid"), "derived FFT caches")
})

test_that("Assessment patches do not adopt unrelated support", {
    rlang::local_bindings(
        custom_resource_patch = function(params, ...) params$initial_n_pp,
        resource_patch_assessment = function(...) data.frame(),
        .env = globalenv()
    )
    p <- setResource(NS_params_small,
                     resource_dynamics = "custom_resource_patch")
    old_transform <- p$capabilities$processes$resource$transformation_binding
    old_tuning <- p$capabilities$processes$resource$tuning_binding
    p <- setProcess(p, "resource",
                    convergence_fun = "resource_patch_assessment")
    expect_identical(p$capabilities$processes$resource$transformation_binding,
                     old_transform)
    expect_identical(
        p$capabilities$processes$resource$tuning_binding, old_tuning
    )
    expect_error(
        prepare_transformation(p, "scale_abundance"), "stale.*resource"
    )
    expect_error(prepare_tuning(p), "stationary.*resource")
})

test_that("Direct kernel hooks must supply valid model geometry", {
    rlang::local_bindings(
        invalid_kernel_transform = function(params, ...) {
            params$pred_kernel <- array(1, c(1, 1, 1))
            comment(params$pred_kernel) <- "set manually"
            params
        }, .env = globalenv()
    )
    p <- setProcess(NS_params_small, "consumers", transformations = list(
        scale_rates = list(
            action = "update", transform_fun = "invalid_kernel_transform",
            paths = list("pred_kernel")
        )
    ))
    plan <- prepare_transformation(p, "scale_rates", list(factor = 2))
    expect_error(apply_transformation_hooks(p, plan), "predator/prey grids")
})
