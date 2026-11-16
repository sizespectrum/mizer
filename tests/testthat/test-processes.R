process_test_assessment <- function(params, ...) {
    data.frame(measure = "residual", value = 0, tolerance = 0.05,
               units = "1/year")
}
process_test_resource <- function(params, n_pp, ...) n_pp
assign("process_test_assessment", process_test_assessment, envir = globalenv())
assign("process_test_resource", process_test_resource, envir = globalenv())

test_that("Built-in processes have bound assessment declarations", {
    p <- NS_params_small
    expect_identical(assessment_owners(p), c("consumers", "resource"))
    expect_named(getProcess(p), c("consumers", "resource", "reproduction"))
    for (owner in c("consumers", "resource")) {
        capability <- assessment_capability(p, owner)
        expect_true(capability$current)
        expect_identical(capability$role, "dynamic")
        expect_type(capability$convergence_fun, "character")
    }
    expect_error(getProcess(p, "recruitment"), "process")
})

test_that("Custom resources need explicit assessment adoption", {
    p <- setResource(NS_params_small,
                     resource_dynamics = "process_test_resource",
                     balance = FALSE)
    expect_null(getProcess(p, "resource")$convergence_fun)
    expect_identical(assessment_capability(p, "resource")$role, "dynamic")
    expect_match(assessment_capability(p, "resource")$reason, "No convergence")
    q <- setProcess(p, "resource", convergence_fun = "process_test_assessment")
    expect_identical(getProcess(q, "resource")$convergence_fun,
                     "process_test_assessment")
    expect_true(assessment_capability(q, "resource")$current)
    expect_null(getProcess(p, "resource")$convergence_fun)
    expect_identical(q$initial_n_pp, p$initial_n_pp)
    expect_identical(q$resource_dynamics, p$resource_dynamics)
})

test_that("Process definitions replace while direct arguments patch", {
    p <- setProcess(NS_params_small, "resource",
                    convergence_fun = "process_test_assessment")
    q <- setProcess(p, "resource", role = "prescribed")
    expect_identical(getProcess(q, "resource")$convergence_fun,
                     "process_test_assessment")
    expect_identical(getProcess(q, "resource")$role, "prescribed")
    replaced <- setProcess(q, "resource", definition = list(role = "dynamic"))
    expect_null(getProcess(replaced, "resource")$convergence_fun)
    adopted <- setProcess(replaced, "resource",
                          definition = list(role = "prescribed"),
                          role = "dynamic",
                          convergence_fun = "process_test_assessment")
    expect_identical(getProcess(adopted, "resource")$role, "dynamic")
    expect_identical(getProcess(adopted, "resource")$convergence_fun,
                     "process_test_assessment")
    cleared <- setProcess(adopted, "resource", convergence_fun = NULL)
    expect_null(getProcess(cleared, "resource")$convergence_fun)
    expect_error(setProcess(p, "resource", definition = NULL), "definition")
    expect_error(setProcess(p, "resource", definition = list(binding = list())),
                 "definition")
    expect_error(setProcess(p, "resource", role = "partial"), "role")
    expect_error(setProcess(p, "resource",
                            convergence_fun = "missing_assessment"),
                 "missing_assessment")
})

test_that("Stale prescribed claims do not exclude changed dynamics", {
    p <- setResource(NS_params_small, resource_dynamics = "resource_constant")
    expect_identical(assessment_capability(p, "resource")$role, "prescribed")
    q <- p
    q$resource_dynamics <- "process_test_resource"
    capability <- assessment_capability(q, "resource")
    expect_false(capability$current)
    expect_identical(capability$role, "dynamic")
    expect_match(capability$reason, "does not match")

    changed <- setResource(p, resource_dynamics = "process_test_resource",
                           balance = FALSE)
    expect_identical(getProcess(changed, "resource")$role, "dynamic")
    expect_null(getProcess(changed, "resource")$convergence_fun)
    expect_identical(getProcess(p, "resource")$role, "prescribed")
})

test_that("Rate parameters and rate choices do not replace fish equations", {
    p <- NS_params_small
    binding <- getProcess(p, "consumers")$binding
    p$rates_funcs$Mort <- "a_different_rate"
    expect_identical(process_provider_binding(p, "consumers"), binding)
    expect_true(assessment_capability(p, "consumers")$current)
    p$extensions <- c(test_process_extension = NA_character_)
    expect_false(assessment_capability(p, "consumers")$current)
})

test_that("Effective S3 dispatch changes invalidate process declarations", {
    p <- NS_params_small
    original <- process_provider_binding(p, "consumers")
    class(p) <- c("process_test_subclass", class(p))
    expect_false(identical(process_provider_binding(p, "consumers"), original))
    expect_false(assessment_capability(p, "consumers")$current)

    local_method <- function(params, ...) params
    # nolint next: object_name_linter. This is an S3 method name.
    assign("project.process_test_subclass", local_method, envir = globalenv())
    withr::defer(rm("project.process_test_subclass", envir = globalenv()))
    methods <- process_dispatch_binding(p)$methods
    expect_true(any(grepl("project.process_test_subclass@", methods,
                          fixed = TRUE)))
    expect_null(default_process_registry(p)$consumers$convergence_fun)
})

test_that("Calling setProcess explicitly adopts a changed provider", {
    p <- setProcess(NS_params_small, "resource",
                    convergence_fun = "process_test_assessment")
    p$resource_dynamics <- "process_test_resource"
    expect_false(assessment_capability(p, "resource")$current)
    adopted <- setProcess(p, "resource")
    expect_true(assessment_capability(adopted, "resource")$current)
    expect_identical(getProcess(adopted, "resource")$convergence_fun,
                     "process_test_assessment")
})

test_that("Restoring a temporary resource override retains its assessment", {
    for (provider in c("resource_semichemostat", "process_test_resource")) {
        p <- setResource(NS_params_small, resource_dynamics = provider,
                         balance = FALSE)
        p <- setProcess(p, "resource",
                        convergence_fun = "process_test_assessment")
        specification <- p$capabilities$processes$resource
        p$resource_dynamics <- "resource_constant"
        expect_false(assessment_capability(p, "resource")$current)
        restored <- setResource(p, resource_dynamics = provider, balance = FALSE)
        expect_identical(restored$capabilities$processes$resource, specification)
        expect_true(assessment_capability(restored, "resource")$current)
    }
})

test_that("Component roles stay distinct from temporary overrides", {
    p <- setComponent(NS_params_small, "pool", 1,
                      dynamics_fun = "constant_other")
    expect_identical(assessment_capability(p, "component:pool")$role, "dynamic")
    expect_null(assessment_capability(p, "component:pool")$convergence_fun)
    q <- setComponent(p, "pool", role = "prescribed")
    expect_identical(assessment_capability(q, "component:pool")$role,
                     "prescribed")
    changed <- setComponent(q, "pool", dynamics_fun = "stats::median")
    expect_identical(getComponent(changed, "pool")$role, "dynamic")
    temporary <- q
    temporary$other_dynamics$pool <- "stats::median"
    expect_false(assessment_capability(temporary, "component:pool")$current)
    expect_identical(assessment_capability(temporary, "component:pool")$role,
                     "dynamic")
    expect_error(assessment_capability(p, "pool"), "Unknown assessment owner")
})

test_that("Process registry validation rejects malformed bindings", {
    registry <- new_capability_registry()
    registry$processes <-
        list(resource = default_process_registry(NS_params_small)$resource)
    registry$processes$resource["binding"] <- list(NULL)
    expect_match(capability_registry_errors(registry), "provider binding")
    registry$processes$resource$binding <- new.env(parent = emptyenv())
    expect_match(capability_registry_errors(registry), "provider binding")
    registry$processes <- list(joint = list())
    expect_match(capability_registry_errors(registry),
                 "consumers, resource and/or reproduction")
})

test_that("Level registration has an independent provider binding", {
    rlang::local_bindings(
        independent_level = function(params, ...) {
            rep(0.5, nrow(params$species_params))
        }, .env = globalenv()
    )
    p <- NS_params_small
    p$rates_funcs$RDD <- "noRDD"
    tuning_binding <- p$capabilities$processes$reproduction$tuning_binding
    transforms <- p$capabilities$processes$reproduction$transformation_binding
    p <- setProcess(p, "reproduction", level_fun = "independent_level")
    definition <- p$capabilities$processes$reproduction
    expect_true(level_binding_current(p, definition))
    expect_false(tuning_binding_current(p, definition, "reproduction"))
    expect_identical(definition$tuning_binding, tuning_binding)
    # Unrelated explicit patches cannot re-adopt a stale explicit transform.
    expect_identical(definition$transformation_binding, transforms)
    expect_error(prepare_tuning(p), "stationary")
    adopted <- setProcess(p, "reproduction")
    expect_true(tuning_binding_current(
        adopted, adopted$capabilities$processes$reproduction, "reproduction"
    ))
})
