# Stored capability specifications. Component equations and custom steps are
# separate capabilities. Assessment declarations state quantitative convergence
# support and distinguish dynamic state from prescribed inputs.

new_capability_registry <- function() {
    list(schema_version = 4L, components = list(), processes = list(),
         tuning = list())
}

# Upgrade data only: extension callbacks may not be available until the model's
# required packages have been loaded. Leave malformed fields for the structural
# validator, rather than silently dropping them while adding the new defaults.
upgrade_capability_registry <- function(registry, params = NULL) {
    if (!is.list(registry) || is.object(registry) ||
            !capability_names_valid(registry) ||
            !setequal(names(registry), c("schema_version", "components",
                                         "processes")) ||
            !is.numeric(registry$schema_version) ||
            length(registry$schema_version) != 1L ||
            is.na(registry$schema_version) ||
            !registry$schema_version %in% 1:3) return(registry)
    registry$schema_version <- 4L
    registry$tuning <- list()
    if (!is.list(registry$components) ||
            !capability_names_valid(registry$components)) return(registry)
    for (component in names(registry$components)) {
        old <- registry$components[[component]]
        if (!is.list(old) || !capability_names_valid(old)) next
        upgraded <- component_definition_defaults()
        upgraded[c("colour", "linetype")] <- list(NULL, NULL)
        upgraded[names(old)] <- old
        registry$components[[component]] <- upgraded
    }
    if (!is.null(params) && is.list(registry$processes)) {
        defaults <- default_process_registry(params)
        for (owner in names(defaults)) {
            old <- registry$processes[[owner]]
            if (is.null(old)) {
                registry$processes[[owner]] <- defaults[[owner]]
            } else {
                missing <- setdiff(names(defaults[[owner]]), names(old))
                old[missing] <- defaults[[owner]][missing]
                registry$processes[[owner]] <- old
            }
        }
    }
    registry
}

component_definition_defaults <- function() {
    list(dynamics_fun = NULL, rate_fun = NULL, encounter_fun = NULL,
         mort_fun = NULL, colour = "grey", linetype = "solid",
         integration_order = NULL, convergence_fun = NULL, role = "dynamic",
         stationary = NULL, equation_rates_fun = NULL, tuning = NULL,
         tuning_binding = NULL)
}

component_definition_fields <- function() {
    c(names(component_definition_defaults()), "level_fun")
}

capability_names_valid <- function(x) {
    if (!length(x)) return(TRUE)
    ids <- names(x)
    !is.null(ids) && !anyNA(ids) && all(nzchar(ids)) && !anyDuplicated(ids)
}

# Validate data without resolving functions: their availability can change
# independently of the model, and must not be cached by validParams().
component_definition_errors <- function(definition, required = TRUE) {
    if (!is.list(definition) || is.object(definition)) {
        return("A component definition must be a plain named list.")
    }
    fields <- names(definition)
    if (!capability_names_valid(definition)) {
        return("Component definition fields must have unique, non-empty names.")
    }
    unknown <- setdiff(fields, component_definition_fields())
    if (length(unknown)) {
        return(paste0("Unknown component definition fields: ",
                      paste(unknown, collapse = ", "), "."))
    }
    errors <- character()
    if (required && is.null(definition$dynamics_fun) &&
            is.null(definition$rate_fun)) {
        errors <- c(errors, "Component requires `dynamics_fun` or `rate_fun`.")
    }
    callbacks <- c("dynamics_fun", "rate_fun", "encounter_fun", "mort_fun",
                   "convergence_fun", "equation_rates_fun", "level_fun")
    for (field in intersect(fields, callbacks)) {
        value <- definition[[field]]
        valid_name <- is.character(value) && length(value) == 1L &&
            !is.na(value) && nzchar(value)
        if (!is.null(value) && !valid_name) {
            errors <- c(errors, paste0("`", field,
                                       "` must be a function name or NULL."))
        }
    }
    if ("role" %in% fields) {
        valid_role <- is.character(definition$role) &&
            length(definition$role) == 1L && !is.na(definition$role) &&
            definition$role %in% c("dynamic", "prescribed")
        if (!valid_role) {
            errors <- c(errors, "`role` must be 'dynamic' or 'prescribed'.")
        }
    }
    if (!is.null(definition$stationary) &&
            (!is.logical(definition$stationary) ||
                 length(definition$stationary) != 1L ||
                 is.na(definition$stationary))) {
        errors <- c(errors, "stationary must be TRUE, FALSE or NULL.")
    }
    tuning_errors <- tuning_block_errors(definition$tuning)
    if (!isTRUE(tuning_errors)) errors <- c(errors, tuning_errors)
    if (!tuning_data_valid(definition$tuning_binding)) {
        errors <- c(errors, "Invalid stored tuning binding.")
    }
    for (field in intersect(fields, c("colour", "linetype"))) {
        value <- definition[[field]]
        valid_style <- (is.character(value) || is.numeric(value)) &&
            length(value) == 1L && !is.na(value)
        if (!is.null(value) && !valid_style) {
            errors <- c(errors, paste0("`", field,
                                       "` must be a single value or NULL."))
        }
    }
    order <- definition$integration_order
    if (!is.null(order)) {
        valid_order <- is.numeric(order) && !is.complex(order) &&
            length(order) == 1L &&
            !is.na(order) && is.finite(order) && order > 0 &&
            order <= .Machine$integer.max && order == floor(order)
        if (!valid_order) {
            errors <- c(errors,
                        "`integration_order` must be a positive whole number.")
        }
        if (required && is.null(definition$dynamics_fun)) {
            errors <- c(errors, paste0("`integration_order` requires ",
                                       "`dynamics_fun` (a custom step)."))
        }
    }
    if (length(errors)) errors else TRUE
}

capability_registry_errors <- function(registry) {
    fields <- c("schema_version", "components", "processes", "tuning")
    if (!is.list(registry) || is.object(registry) ||
            length(registry) != length(fields) ||
            !setequal(names(registry), fields)) {
        return(paste0("`capabilities` must be a registry containing ",
                      "schema_version, components and processes."))
    }
    version <- registry$schema_version
    if (!is.numeric(version) || length(version) != 1L ||
            is.na(version) || version != 4L) {
        return("Unsupported capabilities schema_version; update mizer.")
    }
    components <- registry$components
    if (!is.list(components) || is.object(components) ||
            !capability_names_valid(components)) {
        return("`capabilities$components` must be a uniquely named list.")
    }
    process_errors <- process_registry_errors(registry$processes)
    errors <- if (isTRUE(process_errors)) character() else process_errors
    joint_errors <- joint_tuning_errors(registry$tuning)
    if (!isTRUE(joint_errors)) errors <- c(errors, joint_errors)
    for (component in names(components)) {
        result <- component_definition_errors(components[[component]])
        if (!isTRUE(result)) {
            errors <- c(errors, paste0("Component '", component, "': ", result))
        }
    }
    if (length(errors)) errors else TRUE
}

# Read slots directly: this also runs inside upgradeParams(), where public
# accessors could re-enter validation before the object has been upgraded.
legacy_component_definition <- function(params, component) {
    style <- function(values) {
        if (component %in% names(values)) values[[component]] else NULL
    }
    dynamics_fun <- params$other_dynamics[[component]]
    if (identical(dynamics_fun, "component_rate_dynamics")) {
        dynamics_fun <- NULL
    }
    list(dynamics_fun = dynamics_fun,
         rate_fun = NULL,
         encounter_fun = params$other_encounter[[component]],
         mort_fun = params$other_mort[[component]],
         colour = style(params$linecolour),
         linetype = style(params$linetype),
         integration_order = NULL, convergence_fun = NULL, role = "dynamic",
         stationary = NULL, equation_rates_fun = NULL, tuning = NULL,
         tuning_binding = NULL)
}

legacy_capability_registry <- function(params) {
    registry <- new_capability_registry()
    for (component in names(params$other_dynamics)) {
        registry$components[[component]] <-
            legacy_component_definition(params, component)
    }
    registry$processes <- default_process_registry(params)
    registry
}

existing_component_definition <- function(params, component) {
    definition <- params$capabilities$components[[component]]
    if (!is.null(definition)) {
        return(definition)
    }
    if (component %in% names(params$other_dynamics)) {
        return(legacy_component_definition(params, component))
    }
    NULL
}

assemble_component_definition <- function(existing = NULL, definition = NULL,
                                          overrides = list(), replace = FALSE) {
    check <- function(value, required = FALSE) {
        errors <- component_definition_errors(value, required = required)
        if (!isTRUE(errors)) stop(paste(errors, collapse = "\n"), call. = FALSE)
    }
    check(overrides)
    if (replace) check(definition)
    result <- component_definition_defaults()
    if (!replace && !is.null(existing)) {
        check(existing, required = TRUE)
        result[names(existing)] <- existing
    }
    # Single-bracket assignment preserves explicit NULL entries. modifyList()
    # drops them by default, conflating removal with omission.
    if (replace) result[names(definition)] <- definition
    previous_dynamics <- result$dynamics_fun
    previous_rate <- result$rate_fun
    result[names(overrides)] <- overrides
    if ("dynamics_fun" %in% names(overrides) &&
            !identical(previous_dynamics, result$dynamics_fun) &&
            !"integration_order" %in% names(overrides)) {
        result["integration_order"] <- list(NULL)
    }
    provider_changed <- !identical(previous_dynamics, result$dynamics_fun) ||
        !identical(previous_rate, result$rate_fun)
    if (provider_changed && !"role" %in% names(overrides)) {
        result$role <- "dynamic"
    }
    if (provider_changed && !"convergence_fun" %in% names(overrides)) {
        result["convergence_fun"] <- list(NULL)
    }
    for (field in c("stationary", "equation_rates_fun", "tuning",
                    "tuning_binding")) {
        if (provider_changed && !field %in% names(overrides)) result[field] <-
            list(NULL)
    }
    check(result, required = TRUE)
    for (field in c("dynamics_fun", "rate_fun", "encounter_fun", "mort_fun",
                    "convergence_fun", "equation_rates_fun", "level_fun")) {
        fun <- result[[field]]
        if (!is.null(fun)) {
            context <- paste0("component `", field, "`")
            resolve_model_callback(fun, context = context)
        }
    }
    result[names(component_definition_defaults())]
}

# Plot setters are explicit updates too. Keep their accepted values in the
# snapshot so that a later component patch cannot undo a style-only update.
update_component_styles <- function(params, values, field) {
    # Legacy S4 objects have no stored definitions. Their upgrade will capture
    # the updated styles; keep these lightweight setters usable before then.
    if (isS4(params)) return(params)
    components <- names(params$capabilities$components)
    for (component in intersect(names(values), components)) {
        params$capabilities$components[[component]][field] <-
            list(values[[component]])
    }
    params
}
