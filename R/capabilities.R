# Stored capability specifications. The first schema covers existing component
# callbacks; process capabilities and additional callbacks are added together
# with their runtime consumers in subsequent development previews.

new_capability_registry <- function() {
    list(schema_version = 1L, components = list(), processes = list())
}

component_definition_fields <- function() {
    c("dynamics_fun", "encounter_fun", "mort_fun", "colour", "linetype")
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
    if (required && is.null(definition$dynamics_fun)) {
        errors <- c(errors, "A component definition requires `dynamics_fun`.")
    }
    callbacks <- c("dynamics_fun", "encounter_fun", "mort_fun")
    for (field in intersect(fields, callbacks)) {
        value <- definition[[field]]
        valid_name <- is.character(value) && length(value) == 1L &&
            !is.na(value) && nzchar(value)
        if (!is.null(value) && !valid_name) {
            errors <- c(errors, paste0("`", field,
                                       "` must be a function name or NULL."))
        }
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
    if (length(errors)) errors else TRUE
}

capability_registry_errors <- function(registry) {
    fields <- c("schema_version", "components", "processes")
    if (!is.list(registry) || is.object(registry) ||
            length(registry) != length(fields) ||
            !setequal(names(registry), fields)) {
        return(paste0("`capabilities` must be a registry containing ",
                      "schema_version, components and processes."))
    }
    version <- registry$schema_version
    if (!is.numeric(version) || length(version) != 1L ||
            is.na(version) || version != 1L) {
        return("Unsupported capabilities schema_version; update mizer.")
    }
    components <- registry$components
    if (!is.list(components) || is.object(components) ||
            !capability_names_valid(components)) {
        return("`capabilities$components` must be a uniquely named list.")
    }
    if (!identical(registry$processes, list())) {
        return("Process capability specifications are not supported yet.")
    }
    errors <- character()
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
    list(dynamics_fun = params$other_dynamics[[component]],
         encounter_fun = params$other_encounter[[component]],
         mort_fun = params$other_mort[[component]],
         colour = style(params$linecolour),
         linetype = style(params$linetype))
}

legacy_capability_registry <- function(params) {
    registry <- new_capability_registry()
    for (component in names(params$other_dynamics)) {
        registry$components[[component]] <-
            legacy_component_definition(params, component)
    }
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
    result <- list(dynamics_fun = NULL, encounter_fun = NULL, mort_fun = NULL,
                   colour = "grey", linetype = "solid")
    if (!replace && !is.null(existing)) {
        check(existing, required = TRUE)
        result[names(existing)] <- existing
    }
    # Single-bracket assignment preserves explicit NULL entries. modifyList()
    # drops them by default, conflating removal with omission.
    if (replace) result[names(definition)] <- definition
    result[names(overrides)] <- overrides
    check(result, required = TRUE)
    for (field in c("dynamics_fun", "encounter_fun", "mort_fun")) {
        fun <- result[[field]]
        if (!is.null(fun) && !is.function(get0(fun, mode = "function"))) {
            stop("`", field, "` needs to be the name of a function: '",
                 fun, "' was not found.", call. = FALSE)
        }
    }
    result[component_definition_fields()]
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
