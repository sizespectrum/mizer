test_that("Stored capability schemas reject ambiguous or executable data", {
    registry <- new_capability_registry()
    expect_true(capability_registry_errors(registry))

    registry$schema_version <- 2L
    expect_match(capability_registry_errors(registry), "Unsupported")
    registry$schema_version <- 1L
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
                 "was not found")
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
    expect_false(any(c("initial_value", "component_params") %in%
                         names(definition)))
    expect_true(capability_registry_errors(registry))
    expect_length(registry$processes, 0)
})
