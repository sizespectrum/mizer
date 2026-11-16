test_that("Stored capability schemas reject ambiguous or executable data", {
    registry <- new_capability_registry()
    expect_true(capability_registry_errors(registry))

    registry$schema_version <- 99L
    expect_match(capability_registry_errors(registry), "Unsupported")
    registry$schema_version <- 7L
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
    expect_match(capability_registry_errors(registry), "Invalid specification")
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
    expect_error(assemble(list(dynamics_fun = "constant_other",
                               level_fun = "stats::median")),
                 "Unknown component definition fields.*level_fun")
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
    expect_named(registry$processes, c("consumers", "resource", "reproduction"))
    expect_identical(definition$role, "dynamic")
    expect_null(definition$convergence_fun)
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

test_that("Dynamics replacements require new assessment declarations", {
    original <- assemble_component_definition(
        definition = list(dynamics_fun = "constant_other",
                          convergence_fun = "stats::median",
                          role = "prescribed"), replace = TRUE
    )
    changed <- assemble_component_definition(
        existing = original, overrides = list(dynamics_fun = "stats::median")
    )
    expect_null(changed$convergence_fun)
    expect_identical(changed$role, "dynamic")
    rate_changed <- assemble_component_definition(
        existing = original, overrides = list(rate_fun = "stats::median")
    )
    expect_null(rate_changed$convergence_fun)
    expect_identical(rate_changed$role, "dynamic")
    unchanged <- assemble_component_definition(
        existing = original, overrides = list(dynamics_fun = "constant_other")
    )
    expect_identical(unchanged, original)
    explicit <- assemble_component_definition(
        existing = original,
        overrides = list(dynamics_fun = "stats::median",
                         convergence_fun = "stats::median", role = "prescribed")
    )
    expect_identical(explicit$convergence_fun, "stats::median")
    expect_identical(explicit$role, "prescribed")
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
    expect_identical(upgraded$schema_version, 7L)
    expect_identical(legacy$schema_version, 1L)
    expect_named(upgraded$components$pool,
                 c("dynamics_fun", "rate_fun", "encounter_fun", "mort_fun",
                   "colour", "linetype", "integration_order",
                   "convergence_fun", "role", "stationary",
                   "equation_rates_fun",
                   "tuning", "tuning_binding", "transformations",
                   "transformation_binding", "biomass_fun", "biomass_filters",
                   "prescribed_mask_fun",
                   "prescribed_binding", "autonomous", "temporal_binding"))
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
    invalid$tuning <- NULL
    invalid$transformations <- NULL
    invalid$components <- "not a list"
    upgraded <- upgrade_capability_registry(invalid)
    expect_match(capability_registry_errors(upgraded), "uniquely named list")
    invalid$components <- list(pool = list(dynamics_fun = "unloaded_step",
                                           unknown_field = "do not discard"))
    upgraded <- upgrade_capability_registry(invalid)
    expect_identical(upgraded$components$pool$unknown_field, "do not discard")
    expect_match(capability_registry_errors(upgraded), "Unknown")
})

test_that("Schema four upgrades add transformation support without inference", {
    p <- setComponent(NS_params_small, "pool", 1,
                      dynamics_fun = "constant_other")
    old <- p$capabilities
    old$schema_version <- 4L
    old$transformations <- NULL
    fields <- c("transformations", "transformation_binding", "level_binding")
    old$components <- lapply(old$components, function(definition) {
        definition[setdiff(names(definition), fields)]
    })
    old$processes <- lapply(old$processes, function(definition) {
        definition[setdiff(names(definition), fields)]
    })
    upgraded <- upgrade_capability_registry(old, params = p)
    expect_identical(upgraded$schema_version, 7L)
    expect_identical(upgraded$tuning, old$tuning)
    expect_identical(upgraded$components$pool[names(old$components$pool)],
                     old$components$pool)
    expect_null(upgraded$components$pool$transformations)
    expect_named(upgraded$processes$consumers$transformations,
                 transformation_operations())
    expect_true(capability_registry_errors(upgraded))
    p$capabilities <- upgraded
    changed <- setReproduction(p, RDD = "RickerRDD")
    expect_identical(changed$capabilities$processes$reproduction$level_fun,
                     "reproduction_level_RickerRDD")
})

test_that("Measurement fields validate without executing callbacks", {
    definition <- component_definition_defaults()
    definition$dynamics_fun <- "unloaded_step"
    definition$biomass_fun <- "unloaded_measurement"
    definition$biomass_filters <- c("min_w", "use_cutoff")
    expect_true(component_definition_errors(definition))
    for (filters in list("unknown", c("min_w", "min_w"), NA_character_, 1)) {
        definition$biomass_filters <- filters
        expect_match(component_definition_errors(definition), "biomass_filters")
    }
    definition$biomass_filters <- character()
    definition$biomass_fun <- function() 1
    expect_match(component_definition_errors(definition), "function name")
})

test_that("Schema five upgrades preserve declared measurement support", {
    p <- setComponent(NS_params_small, "pool", 1,
                      dynamics_fun = "constant_other")
    old <- p$capabilities
    old$schema_version <- 5L
    old$components <- lapply(old$components, function(definition) {
        absent <- c("biomass_fun", "biomass_filters", "prescribed_mask_fun",
                    "prescribed_binding", "autonomous", "temporal_binding")
        definition[setdiff(names(definition), absent)]
    })
    old$processes <- lapply(old$processes, function(definition) {
        definition[setdiff(names(definition),
                           c("mortality_fun", "mortality_binding"))]
    })
    p$capabilities <- old
    upgraded <- upgrade_capability_registry(old, params = p)
    expect_identical(upgraded$schema_version, 7L)
    expect_identical(upgraded$tuning, old$tuning)
    expect_identical(upgraded$transformations, old$transformations)
    expect_identical(upgraded$components$pool[names(old$components$pool)],
                     old$components$pool)
    expect_null(upgraded$components$pool$biomass_fun)
    expect_identical(upgraded$components$pool$biomass_filters, character())
    expect_identical(upgraded$processes$consumers$mortality_fun,
                     "builtin_mortality_breakdown")
    expect_true(capability_registry_errors(upgraded))
    p$rates_funcs$Mort <- "unloaded_custom_mortality"
    unsupported <- upgrade_capability_registry(old, params = p)
    expect_null(unsupported$processes$consumers$mortality_fun)
    expect_identical(
        upgrade_capability_registry(upgraded, params = p), upgraded
    )
})

test_that("Schema six keeps roles and adds honest temporal support", {
    p <- setComponent(NS_params_small, "pool", 1,
                      dynamics_fun = "constant_other", role = "prescribed",
                      stationary = TRUE)
    old <- p$capabilities
    old$schema_version <- 6L
    fields <- c("prescribed_mask_fun", "prescribed_binding", "autonomous",
                "temporal_binding")
    old$components$pool[fields] <- NULL
    for (owner in names(old$processes)) {
        missing <- c(fields, "rate_fun", "equation_binding")
        old$processes[[owner]][missing] <- NULL
    }
    upgraded <- upgrade_capability_registry(old, p)
    expect_identical(upgraded$schema_version, 7L)
    expect_identical(upgraded$components$pool$role, "prescribed")
    expect_null(upgraded$components$pool$autonomous)
    expect_null(upgraded$components$pool$prescribed_mask_fun)
    expect_true(upgraded$processes$consumers$autonomous)
    expect_identical(upgraded$processes$resource$rate_fun,
                     "equation_resource_semichemostat")
    expect_true(capability_registry_errors(upgraded))
    expect_identical(upgrade_capability_registry(upgraded, p), upgraded)
    restored <- p
    restored$capabilities <- upgraded
    restored <- setResource(restored, resource_dynamics = "resource_logistic",
                            balance = FALSE)
    expect_true(prepare_state_capabilities(restored)$resource$mask_current)
    p$rates_funcs$RDD <- "unknown_forcing"
    unknown <- upgrade_capability_registry(old, p)
    expect_null(unknown$processes$consumers$autonomous)
    expect_null(unknown$processes$resource$autonomous)
})


test_that("State declarations reject malformed private bindings and scalars", {
    for (field in c("prescribed_binding", "temporal_binding")) {
        definition <- list(dynamics_fun = "constant_other")
        definition[field] <- list("not-a-binding")
        expect_match(component_definition_errors(definition), "Invalid stored")
    }
    for (value in list(NA, 1, "yes", c(TRUE, FALSE))) {
        definition <- list(dynamics_fun = "constant_other", autonomous = value)
        expect_match(component_definition_errors(definition), "autonomous")
    }
})


test_that("Schema migration leaves malformed process records for validation", {
    p <- NS_params_small
    old <- p$capabilities
    old$schema_version <- 6L
    old$processes$resource <- "invalid"
    upgraded <- upgrade_capability_registry(old, p)
    expect_match(capability_registry_errors(upgraded), "Invalid specification")
})
