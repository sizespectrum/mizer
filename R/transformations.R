# Explicit transformation capabilities. Core operations construct candidates;
# extension hooks update declared owned data before dependent calculations.

transformation_operations <- function() {
    c(
        "add_species", "remove_species", "reorder_species", "rename_species",
        "size_grid", "scale_abundance", "scale_rates"
    )
}

transformation_defaults <- function() {
    list(
        action = "none", transform_fun = NULL, paths = list(),
        settings = list()
    )
}

transformation_spec_errors <- function(specification) {
    if (!is.list(specification) || is.object(specification) ||
            !capability_names_valid(specification) ||
            any(!names(specification) %in% names(transformation_defaults())) ||
            !"action" %in% names(specification)) {
        return("A transformation needs an explicit action and known fields.")
    }
    spec <- transformation_defaults()
    spec[names(specification)] <- specification
    if (!is.character(spec$action) || length(spec$action) != 1L ||
            is.na(spec$action) || !spec$action %in% c("none", "update")) {
        return("Transformation action must be 'none' or 'update'.")
    }
    if (spec$action == "none" &&
            (!is.null(spec$transform_fun) || length(spec$paths) ||
                 length(spec$settings))) {
        return("A no-action transformation cannot declare a hook or updates.")
    }
    if (spec$action == "update" &&
            (!is.character(spec$transform_fun) ||
                 length(spec$transform_fun) != 1L ||
                 is.na(spec$transform_fun) ||
                 !nzchar(spec$transform_fun))) {
        return("An update transformation requires a named transform_fun.")
    }
    if (!is.list(spec$paths) || is.object(spec$paths) ||
        !all(vapply(spec$paths, function(path) {
            is.character(path) && length(path) && !anyNA(path) &&
                all(nzchar(path))
        }, logical(1)))) {
        return("Transformation paths must be a list of character vectors.")
    }
    if (!is.list(spec$settings) || !tuning_data_valid(spec$settings)) {
        return("Transformation settings must contain serialisable data.")
    }
    TRUE
}

transformation_map_errors <- function(transformations) {
    if (is.null(transformations)) {
        return(TRUE)
    }
    if (!is.list(transformations) || is.object(transformations) ||
            !capability_names_valid(transformations) ||
            any(!names(transformations) %in% transformation_operations())) {
        return("Transformations must name supported operations.")
    }
    errors <- lapply(transformations, transformation_spec_errors)
    errors <- errors[!vapply(errors, isTRUE, logical(1))]
    if (length(errors)) unlist(errors, use.names = FALSE) else TRUE
}

normalise_transformation <- function(specification) {
    errors <- transformation_spec_errors(specification)
    if (!isTRUE(errors)) stop(errors, call. = FALSE)
    result <- transformation_defaults()
    result[names(specification)] <- specification
    result
}

# nolint next: object_length_linter. Shared capability helper.
transformation_provider_binding <- function(params) {
    binding <- tuning_provider_binding(params)
    # Reproduction balancing does not define structural or scaling behaviour.
    binding$balance_methods <- NULL
    setters <- list(
        consumers = c(
            "setParams", "setInteraction",
            setdiff(set_params_setters, "setReproduction"),
            "setColours", "setLinetypes"
        ),
        resource = "setResource", reproduction = "setReproduction"
    )
    binding$construction_methods <- lapply(setters, function(generics) {
        tuning_method_names(params, generics)
    })
    binding
}

transformation_binding_current <- function(params, definition, owner) {
    stored <- definition$transformation_binding
    if (is.null(stored)) {
        return(FALSE)
    }
    stored$automatic <- NULL
    current <- transformation_provider_binding(params)
    if (!startsWith(owner, "component:")) {
        stored$construction_methods <- stored$construction_methods[owner]
        current$construction_methods <- current$construction_methods[owner]
    }
    identical(
        relevant_provider_binding(stored, owner),
        relevant_provider_binding(current, owner)
    )
}

# The core owns its documented built-in structural and scaling updates.
# Custom processes must explicitly adopt that behaviour or provide a hook.
# nolint next: object_length_linter. Shared capability helper.
default_provider_transformations <- function(params, owner) {
    known <- builtin_tuning_rates(
        params, reproduction = owner == "reproduction"
    )
    binding <- transformation_provider_binding(params)
    known <- known && all(endsWith(
        binding$construction_methods[[owner]],
        "@mizer"
    ))
    if (owner == "resource") {
        known <- known && sub("^mizer::", "", params$resource_dynamics) %in%
            c("resource_constant", "resource_semichemostat",
              "resource_logistic")
    } else if (owner == "reproduction") {
        known <- known && !is.null(builtin_reproduction_law(params))
    }
    support <- if (known) {
        stats::setNames(
            lapply(transformation_operations(), function(operation) {
                transformation_defaults()
            }), transformation_operations()
        )
    } else {
        NULL
    }
    binding$automatic <- TRUE
    list(transformations = support, transformation_binding = binding)
}

# Built-in support is refreshed only when the old recipe was automatic.
# Explicit author definitions stay stored and require explicit re-adoption.
# nolint next: object_length_linter. Shared capability helper.
refresh_builtin_transformations <- function(params) {
    for (owner in c("consumers", "resource", "reproduction")) {
        old <- params$capabilities$processes[[owner]]
        if (is.null(old) ||
                transformation_binding_current(params, old, owner)) {
            next
        }
        automatic <- isTRUE(old$transformation_binding$automatic)
        if (!automatic) next
        defaults <- default_provider_transformations(params, owner)
        if (!is.null(defaults$transformations)) {
            old[names(defaults)] <- defaults
            params$capabilities$processes[[owner]] <- old
        }
    }
    params
}

joint_transformation_errors <- function(registry) {
    if (!is.list(registry) || is.object(registry) ||
            !capability_names_valid(registry)) {
        return("Joint transformations must be a uniquely named list.")
    }
    for (entry in registry) {
        if (!is.list(entry) || !capability_names_valid(entry) ||
            !setequal(names(entry), c(
                "owners", "operation",
                "definition", "binding"
            )) ||
            !is.character(entry$owners) || !length(entry$owners) ||
            anyNA(entry$owners) || any(!nzchar(entry$owners)) ||
            anyDuplicated(entry$owners) ||
            !is.character(entry$operation) || length(entry$operation) != 1L ||
            is.na(entry$operation) ||
            !entry$operation %in% transformation_operations() ||
            !tuning_data_valid(entry$binding)) {
            return("Invalid joint transformation declaration.")
        }
        errors <- transformation_spec_errors(entry$definition)
        if (!isTRUE(errors)) {
            return(errors)
        }
    }
    TRUE
}

transformation_default_paths <- function(owners) {
    components <- substring(owners[startsWith(owners, "component:")], 11L)
    unlist(lapply(components, function(component) {
        list(c("initial_n_other", component), c("other_params", component))
    }), recursive = FALSE)
}

transformation_path_allowed <- function(path, owners) {
    root <- path[[1L]]
    components <- substring(owners[startsWith(owners, "component:")], 11L)
    if (root == "initial_n_other") {
        return(length(path) >= 2L && path[[2L]] %in% components)
    }
    if (root == "other_params") {
        return(length(path) >= 2L)
    }
    if (root == "capabilities") {
        if (length(path) < 3L) return(FALSE)
        component <- path[[2L]] == "components" && path[[3L]] %in% components
        process <- path[[2L]] == "processes" && path[[3L]] %in% owners
        return(component || process)
    }
    if (root %in% c("species_params", "given_species_params")) {
        return(length(path) == 2L &&
                   any(c("consumers", "reproduction") %in% owners) &&
                   !path[[2L]] %in% c("species", "w_min", "w_max"))
    }
    resource <- c("initial_n_pp", "cc_pp", "rr_pp", "resource_params")
    if (root %in% resource) {
        return("resource" %in% owners)
    }
    consumer <- c(
        "initial_n", "maturity", "psi", "intake_max", "search_vol",
        "metab", "mu_b", "ext_encounter", "ext_diffusion",
        "ft_pred_kernel_e", "ft_pred_kernel_p", "ft_pred_kernel_d",
        "pred_kernel", "selectivity", "catchability", "interaction"
    )
    if (root %in% consumer) {
        return("consumers" %in% owners)
    }
    # Remaining core bookkeeping and discretisation are owned by the operation.
    !root %in% c(
        "metadata", "mizer_version", "extensions", "time_created",
        "time_modified", "w", "dw", "w_full", "dw_full", "w_min_idx",
        "gear_params", "initial_effort", "sc", "other_dynamics",
        "other_encounter", "other_mort", "rates_funcs",
        "resource_dynamics", "A", "linecolour", "linetype", "ft_mask",
        "use_predation_diffusion", "second_order_w", "defaults_edition"
    )
}

transformation_paths <- function(specification, owners) {
    if (specification$action == "none") {
        return(list())
    }
    paths <- c(transformation_default_paths(owners), specification$paths)
    derived <- c("ft_pred_kernel_e", "ft_pred_kernel_p", "ft_pred_kernel_d")
    fft <- vapply(paths, function(path) path[[1L]] %in% derived, logical(1))
    if (any(fft)) {
        stop("Declare pred_kernel or its parameters, not derived FFT caches.")
    }
    keys <- vapply(paths, paste, character(1), collapse = "/")
    paths <- paths[!duplicated(keys)]
    if (any(!vapply(paths, transformation_path_allowed, logical(1),
                owners = owners
            ))) {
        stop("A transformation path is outside its owners' responsibility.",
            call. = FALSE
        )
    }
    # A parent path already grants the more specific permission.
    keys <- vapply(paths, paste, character(1), collapse = "/")
    keep <- vapply(seq_along(keys), function(i) {
        !any(startsWith(keys[[i]], paste0(keys[-i], "/")))
    }, logical(1))
    paths[keep]
}

prepare_transformation <- function(before, operation, details = list()) {
    before <- validParams(before)
    if (!is.character(operation) || length(operation) != 1L ||
            is.na(operation) || !operation %in% transformation_operations()) {
        stop("Unknown model transformation operation.", call. = FALSE)
    }
    if (!is.list(details)) stop("Transformation details must be a list.")
    owners <- tuning_owners(before)
    binding <- transformation_provider_binding(before)
    selected <- list()
    covered <- character()
    for (id in names(before$capabilities$transformations)) {
        joint <- before$capabilities$transformations[[id]]
        if (joint$operation != operation) next
        if (any(!joint$owners %in% owners) ||
                length(intersect(covered, joint$owners))) {
            stop("Joint transformation owners overlap or are unavailable.")
        }
        if (!identical(joint$binding, binding)) {
            stop("Joint transformation support is stale; register it again.")
        }
        selected[[paste0("joint:", id)]] <- list(
            owners = joint$owners, specification = joint$definition
        )
        covered <- c(covered, joint$owners)
    }
    for (owner in setdiff(owners, covered)) {
        definition <- tuning_definition(before, owner)
        spec <- definition$transformations[[operation]]
        if (is.null(spec) ||
                !transformation_binding_current(before, definition, owner)) {
            stop("Missing or stale transformation support for '", owner,
                "' during ", operation, ".",
                call. = FALSE
            )
        }
        selected[[owner]] <- list(owners = owner, specification = spec)
    }
    claimed <- character()
    for (id in names(selected)) {
        node <- selected[[id]]
        node$specification <- normalise_transformation(node$specification)
        node$paths <- transformation_paths(node$specification, node$owners)
        paths <- vapply(node$paths, paste, character(1), collapse = "/")
        for (path in paths) {
            if (any(path == claimed | startsWith(path, paste0(claimed, "/")) |
                        startsWith(claimed, paste0(path, "/")))) {
                stop("Transformation paths overlap; declare joint ownership.")
            }
        }
        claimed <- c(claimed, paths)
        if (node$specification$action == "update") {
            node$fun <- resolve_model_callback(
                node$specification$transform_fun,
                paste0("Transformation ", id)
            )
        }
        selected[[id]] <- node
    }
    list(
        before = clone_params_runtime(before), operation = operation,
        details = unserialize(serialize(details, NULL)), selected = selected,
        binding = binding
    )
}

transformation_path_entry <- function(object, path) {
    for (part in path) {
        if (!is.list(object) || !part %in% names(object)) {
            return(list(exists = FALSE, value = NULL))
        }
        object <- object[[part]]
    }
    list(exists = TRUE, value = object)
}

replace_transformation_path <- function(object, path, entry) {
    key <- path[[1L]]
    if (length(path) == 1L) {
        if (entry$exists) {
            object[key] <- list(entry$value)
        } else {
            object[key] <- NULL
        }
    } else if (entry$exists || key %in% names(object)) {
        object[[key]] <- replace_transformation_path(
            object[[key]], path[-1L],
            entry
        )
    }
    object
}

capture_transformation_updates <- function(before_hooks, after_hooks, plan) {
    updates <- list()
    for (node in plan$selected) {
        for (path in node$paths) {
            before <- transformation_path_entry(before_hooks, path)
            after <- transformation_path_entry(after_hooks, path)
            if (!isTRUE(all.equal(before, after, tolerance = 0,
                                  check.attributes = TRUE))) {
                updates[[length(updates) + 1L]] <- list(path = path,
                                                        entry = after)
            }
        }
    }
    unserialize(serialize(updates, NULL))
}

restore_transformation_updates <- function(candidate, updates) {
    for (update in updates) {
        entry <- unserialize(serialize(update$entry, NULL))
        candidate <- replace_transformation_path(candidate, update$path, entry)
    }
    candidate
}

apply_transformation_hooks <- function(candidate, plan) {
    candidate <- clone_params_runtime(candidate)
    if (!identical(transformation_provider_binding(candidate), plan$binding) ||
            !isTRUE(all.equal(candidate$capabilities,
                plan$before$capabilities,
                tolerance = 0,
                check.attributes = TRUE
            ))) {
        stop(
            "The transformation candidate lost its stored providers or ",
            "capability specifications."
        )
    }
    for (id in names(plan$selected)) {
        node <- plan$selected[[id]]
        if (node$specification$action == "none") next
        expected <- clone_params_runtime(candidate)
        owner <- if (length(node$owners) == 1L) node$owners else NULL
        component <- if (!is.null(owner) && startsWith(owner, "component:")) {
            substring(owner, 11L)
        } else {
            NULL
        }
        process <- if (is.null(component)) owner else NULL
        args <- list(
            before = clone_params_runtime(plan$before),
            operation = plan$operation,
            details = unserialize(serialize(plan$details, NULL)),
            owners = node$owners, owner = owner, component = component,
            process = process,
            settings = unserialize(serialize(node$specification$settings, NULL))
        )
        candidate <- invoke_tuning_callback(
            node$fun, candidate, args,
            paste0("Transformation ", id)
        )
        if (!inherits(candidate, "MizerParams")) {
            stop("A transformation hook must return a MizerParams object.")
        }
        kernel <- candidate$pred_kernel
        kernel_changed <- !isTRUE(all.equal(kernel, expected$pred_kernel,
                                            tolerance = 0,
                                            check.attributes = TRUE))
        empty_kernel <- isTRUE(all.equal(kernel, array(), tolerance = 0,
                                         check.attributes = TRUE))
        if (kernel_changed && !empty_kernel) {
            kernel_dim <- c(nrow(candidate$species_params), length(candidate$w),
                            length(candidate$w_full))
            kernel_names <- list(
                sp = candidate$species_params$species,
                w_pred = as.character(signif(candidate$w, 3)),
                w_prey = as.character(signif(candidate$w_full, 3))
            )
            if (!is.numeric(kernel) || is.complex(kernel) ||
                    !identical(dim(kernel), kernel_dim) ||
                    any(!is.finite(kernel)) || any(kernel < 0) ||
                    (!is.null(dimnames(kernel)) &&
                         !identical(dimnames(kernel), kernel_names))) {
                stop("An authored pred_kernel must match the model's species ",
                     "and predator/prey grids and be finite and nonnegative.")
            }
        }
        core_arrays <- c(
            "initial_n", "initial_n_pp", "cc_pp", "rr_pp",
            "maturity", "psi", "intake_max", "search_vol",
            "metab", "mu_b", "ext_encounter", "ext_diffusion",
            "selectivity", "catchability", "interaction"
        )
        shape <- function(value) {
            list(
                dim = dim(value), dimnames = dimnames(value),
                names = names(value), length = length(value)
            )
        }
        if (any(!vapply(core_arrays, function(field) {
            identical(shape(candidate[[field]]), shape(expected[[field]]))
        }, logical(1)))) {
            stop("A transformation hook changed core state or rate layout.")
        }
        for (path in node$paths) {
            entry <- transformation_path_entry(candidate, path)
            expected <- replace_transformation_path(expected, path, entry)
        }
        expected$time_modified <- candidate$time_modified
        if (!isTRUE(all.equal(candidate, expected,
                tolerance = 0,
                check.attributes = TRUE
            ))) {
            stop("Transformation hook '", id,
                "' changed data outside its declared paths.",
                call. = FALSE
            )
        }
    }
    candidate
}

#' Register a joint model transformation
#'
#' Local transformation declarations belong to [setComponent()] or
#' [setProcess()]. Use a joint declaration when one hook must update data shared
#' by several owners. It replaces their local declarations for its operation;
#' overlapping owners or update paths are errors. Hooks run on the constructed
#' candidate before dependent rates are rebuilt or the complete model is
#' validated. Missing support is an error; `action = "none"` explicitly declares
#' that the core operation needs no additional update from this provider.
#'
#' @param params A MizerParams object.
#' @param name Name identifying the joint declaration.
#' @param owners Canonical owners: consumers, resource, reproduction or
#'   component:name.
#' @param operation One of add_species, remove_species, reorder_species,
#'   rename_species, size_grid, scale_abundance or scale_rates.
#' @param definition A list with action (none or update), transform_fun (a named
#'   callback for update), paths (a list of literal character-vector paths) and
#'   serialisable settings. NULL removes the named joint declaration.
#' @return The updated MizerParams object.
#' @export
setTransformation <- function(params, name, owners = NULL, operation = NULL,
                              definition) {
    params <- validParams(params)
    if (!is.character(name) || length(name) != 1L || is.na(name) ||
            !nzchar(name)) {
        stop("A transformation needs one nonempty name.")
    }
    if (is.null(definition)) {
        params$capabilities$transformations[[name]] <- NULL
        return(params)
    }
    entry <- list(
        owners = owners, operation = operation,
        definition = normalise_transformation(definition),
        binding = transformation_provider_binding(params)
    )
    errors <- joint_transformation_errors(list(entry = entry))
    if (!isTRUE(errors)) stop(errors, call. = FALSE)
    if (any(!owners %in% tuning_owners(params))) {
        stop("Unknown transformation owner.")
    }
    transformation_paths(entry$definition, owners)
    params$capabilities$transformations[[name]] <- entry
    params$time_modified <- lubridate::now()
    params
}
