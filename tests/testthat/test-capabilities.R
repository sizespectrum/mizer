test_that("Stored capability schemas reject ambiguous or executable data", {
    registry <- new_capability_registry()
    expect_true(capability_registry_errors(registry))

    registry$schema_version <- 3L
    expect_match(capability_registry_errors(registry), "Unsupported")
    registry$schema_version <- 2L
    registry$components <- list(pool = list(dynamics_fun = function() 1))
    expect_match(capability_registry_errors(registry), "function name")
    registry$components <- list(pool = list(dynamics_fun = "unloaded_callback"))
    # Availability is checked at registration/use, outside the fingerprint gate.
    expect_true(capability_registry_errors(registry))
    registry$components$pool$parameters <- new.env(parent = emptyenv())
    expect_match(capability_registry_errors(registry), "Unknown")

    registry$components <- setNames(list(list(), list()), c("pool", "pool"))
    expect_match(capability_registry_errors(registry), "uniquely named")
    registry$components <- setNames(list(list()), NA_character_)
    expect_match(capability_registry_errors(registry), "uniquely named")
    registry$components <- list()
    registry$processes <- list(resource = list())
    expect_match(capability_registry_errors(registry), "not supported yet")
})

test_that("Definitions fail early on misspellings and unnamed callbacks", {
    assemble <- function(definition) {
        assemble_component_definition(definition = definition, replace = TRUE)
    }
    expect_error(assemble(list(dynamic_fun = "constant_other")),
                 "Unknown component definition fields")
    expect_error(assemble(list("constant_other")), "unique")
    expect_error(assemble(list(dynamics_fun = "constant_other",
                               dynamics_fun = "constant_other")), "unique")
    expect_error(assemble(list(dynamics_fun = "missing_component_callback")),
                 "missing_component_callback")
    expect_error(assemble(list(dynamics_fun = function(...) 1)),
                 "function name")
    expect_error(assemble(NULL), "plain named list")
    expect_error(assemble(list(dynamics_fun = NULL)), "requires `dynamics_fun`")
})

test_that("Legacy snapshots preserve only known component registrations", {
    p <- NS_params_small
    p$other_dynamics <- list(pool = "unloaded_component")
    p$other_encounter <- list(pool = "unloaded_encounter",
                              free_food = "unloaded_free_contribution")
    p$initial_n_other <- list(pool = list(amount = 1))
    p$other_params <- list(pool = list(capacity = 2))
    p$linecolour <- c(p$linecolour, pool = "orange")
    registry <- legacy_capability_registry(p)

    expect_named(registry$components, "pool")
    definition <- registry$components$pool
    expect_identical(definition$dynamics_fun, "unloaded_component")
    expect_identical(definition$encounter_fun, "unloaded_encounter")
    expect_identical(definition$colour, "orange")
    expect_null(definition$mort_fun)
    expect_null(definition$linetype)
    expect_null(definition$rate_fun)
    expect_null(definition$integration_order)
    expect_false(any(c("initial_value", "component_params") %in%
                         names(definition)))
    expect_true(capability_registry_errors(registry))
    expect_length(registry$processes, 0)
})

test_that("Equations and custom steps are separate capabilities", {
    equations <- assemble_component_definition(
        definition = list(rate_fun = "stats::median"), replace = TRUE
    )
    expect_null(equations$dynamics_fun)
    expect_identical(equations$rate_fun, "stats::median")
    expect_null(equations$integration_order)
    expect_true(component_definition_errors(equations))

    both <- assemble_component_definition(
        existing = equations,
        overrides = list(dynamics_fun = "constant_other",
                         integration_order = 2L)
    )
    expect_identical(both$dynamics_fun, "constant_other")
    expect_identical(both$rate_fun, "stats::median")
    expect_identical(both$integration_order, 2L)
    expect_error(assemble_component_definition(
        existing = both, overrides = list(dynamics_fun = NULL, rate_fun = NULL)
    ), "dynamics_fun.*rate_fun")
})

test_that("Declared integration order is valid only for a custom step", {
    invalid <- list(0, -1L, 0.5, Inf, NA_real_, NaN, "2", TRUE, c(1L, 2L),
                    2 + 0i, .Machine$integer.max + 1)
    for (order in invalid) {
        expect_match(component_definition_errors(
            list(dynamics_fun = "constant_other", integration_order = order)
        ), "positive whole number")
    }
    expect_match(component_definition_errors(
        list(rate_fun = "stats::median", integration_order = 2L)
    ), "requires `dynamics_fun`")
    expect_true(component_definition_errors(
        list(dynamics_fun = "constant_other", integration_order = 2L)
    ))
})

test_that("Changing a custom step clears its inherited order", {
    original <- assemble_component_definition(
        definition = list(dynamics_fun = "constant_other",
                          integration_order = 2L),
        replace = TRUE
    )
    changed <- assemble_component_definition(
        existing = original, overrides = list(dynamics_fun = "stats::median")
    )
    expect_null(changed$integration_order)
    same <- assemble_component_definition(
        existing = original, overrides = list(dynamics_fun = "constant_other")
    )
    expect_identical(same$integration_order, 2L)
    explicit <- assemble_component_definition(
        existing = original,
        overrides = list(dynamics_fun = "stats::median", integration_order = 1L)
    )
    expect_identical(explicit$integration_order, 1L)
    overridden <- assemble_component_definition(
        definition = original, replace = TRUE,
        overrides = list(dynamics_fun = "stats::median")
    )
    expect_null(overridden$integration_order)
    removed <- assemble_component_definition(
        existing = original,
        overrides = list(dynamics_fun = NULL, rate_fun = "stats::median")
    )
    expect_null(removed$integration_order)
    expect_null(removed$dynamics_fun)
})

test_that("Schema one upgrades keep stored choices without resolving code", {
    legacy <- list(
        schema_version = 1L,
        components = list(pool = list(dynamics_fun = "unloaded_step",
                                      mort_fun = "unloaded_mortality",
                                      colour = "orange", linetype = "dashed",
                                      encounter_fun = NULL)),
        processes = list()
    )
    upgraded <- upgrade_capability_registry(legacy)
    expect_identical(upgraded$schema_version, 2L)
    expect_identical(legacy$schema_version, 1L)
    expect_named(upgraded$components$pool,
                 c("dynamics_fun", "rate_fun", "encounter_fun", "mort_fun",
                   "colour", "linetype", "integration_order"))
    expect_identical(upgraded$components$pool$dynamics_fun, "unloaded_step")
    expect_identical(upgraded$components$pool$mort_fun, "unloaded_mortality")
    expect_identical(upgraded$components$pool$colour, "orange")
    expect_identical(upgraded$components$pool$linetype, "dashed")
    expect_null(upgraded$components$pool$rate_fun)
    expect_null(upgraded$components$pool$integration_order)
    expect_true(capability_registry_errors(upgraded))
    expect_identical(upgrade_capability_registry(upgraded), upgraded)
})

test_that("Schema upgrades retain malformed and unknown data for validation", {
    unknown <- new_capability_registry()
    unknown$schema_version <- 99L
    expect_identical(upgrade_capability_registry(unknown), unknown)
    expect_match(capability_registry_errors(unknown), "Unsupported")
    expect_null(upgrade_capability_registry(NULL))

    invalid <- new_capability_registry()
    invalid$schema_version <- 1L
    invalid$components <- "not a list"
    upgraded <- upgrade_capability_registry(invalid)
    expect_match(capability_registry_errors(upgraded), "uniquely named list")
    invalid$components <- list(pool = list(dynamics_fun = "unloaded_step",
                                           unknown_field = "do not discard"))
    upgraded <- upgrade_capability_registry(invalid)
    expect_identical(upgraded$components$pool$unknown_field, "do not discard")
    expect_match(capability_registry_errors(upgraded), "Unknown")
})
