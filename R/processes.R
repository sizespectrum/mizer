# Process assessment declarations are bound to the provider that was explicitly
# adopted. Comparing identifiers does not freeze callback implementations.

process_dispatch_binding <- function(params) {
    classes <- unique(c(class(params),
                        names(extensionRequirements(params$extensions))))
    methods <- unlist(lapply(c("project", "project_simple"), function(generic) {
        unlist(lapply(classes, function(class_name) {
            method <- getS3method(generic, class_name, optional = TRUE)
            if (is.null(method)) return(character())
            paste0(generic, ".", class_name, "@",
                   environmentName(environment(method)))
        }), use.names = FALSE)
    }), use.names = FALSE)
    list(classes = classes, methods = methods)
}

process_provider_binding <- function(params, process) {
    binding <- process_dispatch_binding(params)
    binding$provider <- if (process == "resource") {
        params$resource_dynamics
    } else {
        "mizer_consumer_equations"
    }
    binding
}

process_registry_errors <- function(processes) {
    if (!is.list(processes) || is.object(processes) ||
            !capability_names_valid(processes) ||
            any(!names(processes) %in% c("consumers", "resource"))) {
        return("`capabilities$processes` must name consumers and/or resource.")
    }
    errors <- character()
    for (process in names(processes)) {
        specification <- processes[[process]]
        if (!is.list(specification) || is.object(specification) ||
                !capability_names_valid(specification) ||
                !setequal(names(specification),
                          c("convergence_fun", "role", "binding"))) {
            errors <- c(errors, paste0("Invalid specification for process '",
                                       process, "'."))
            next
        }
        fragment <- specification[c("convergence_fun", "role")]
        result <- component_definition_errors(fragment, required = FALSE)
        if (!isTRUE(result)) errors <- c(errors, result)
        binding <- specification$binding
        binding_fields <- c("classes", "methods", "provider")
        if (!is.list(binding) || is.object(binding) ||
            !capability_names_valid(binding) ||
            !setequal(names(binding), binding_fields) ||
            !all(vapply(binding, function(x) {
                is.character(x) && !anyNA(x) && all(nzchar(x))
            }, logical(1))) || length(binding$provider) != 1L) {
            errors <- c(errors, paste0("Invalid provider binding for process '",
                                       process, "'."))
        }
    }
    if (length(errors)) errors else TRUE
}

# Only core projection methods can receive automatic built-in declarations.
# Extensions adopt equivalent support explicitly with setProcess().
default_process_registry <- function(params) {
    dispatch <- process_dispatch_binding(params)
    expected <- c("project.MizerParams@mizer",
                  "project_simple.MizerParams@mizer")
    builtin_projection <- identical(dispatch$methods, expected)
    consumers <- list(
        convergence_fun = if (builtin_projection) {
            "builtin_consumer_assessment"
        } else {
            NULL
        },
        role = "dynamic",
        binding = process_provider_binding(params, "consumers")
    )
    resource_name <- sub("^mizer::", "", params$resource_dynamics)
    builtin_resource <- builtin_projection && resource_name %in%
        c("resource_semichemostat", "resource_logistic", "resource_constant")
    prescribed_resource <- builtin_resource &&
        resource_name == "resource_constant"
    resource <- list(
        convergence_fun = if (builtin_resource && !prescribed_resource) {
            "builtin_resource_assessment"
        } else {
            NULL
        },
        role = if (prescribed_resource) {
            "prescribed"
        } else {
            "dynamic"
        },
        binding = process_provider_binding(params, "resource")
    )
    list(consumers = consumers, resource = resource)
}

refresh_resource_assessment <- function(params, previous) {
    current_binding <- process_provider_binding(params, "resource")
    stored_binding <- params$capabilities$processes$resource$binding
    if (!identical(previous, params$resource_dynamics) &&
            !identical(stored_binding, current_binding)) {
        params$capabilities$processes$resource <-
            default_process_registry(params)$resource
    }
    params
}

assessment_owners <- function(params) {
    components <- names(params$other_dynamics)
    c("consumers", "resource",
      if (length(components)) paste0("component:", components) else character())
}

assessment_capability <- function(params, owner) {
    if (!owner %in% assessment_owners(params)) {
        stop("Unknown assessment owner '", owner, "'.", call. = FALSE)
    }
    component <- startsWith(owner, "component:")
    if (component) {
        id <- substring(owner, nchar("component:") + 1L)
        definition <- params$capabilities$components[[id]]
        expected <- definition$dynamics_fun
        if (is.null(expected)) expected <- "component_rate_dynamics"
        current <- !is.null(definition) &&
            identical(params$other_dynamics[[id]], expected)
    } else {
        definition <- params$capabilities$processes[[owner]]
        current <- !is.null(definition) &&
            identical(definition$binding,
                      process_provider_binding(params, owner))
    }
    role <- if (current) definition$role else "dynamic"
    if (is.null(role)) role <- "dynamic"
    reason <- if (is.null(definition)) {
        "No assessment capability registered."
    } else if (!current) {
        paste0("The assessment declaration does not match ",
               "the current dynamics provider.")
    } else if (role == "dynamic" && is.null(definition$convergence_fun)) {
        "No convergence callback declared."
    } else {
        NULL
    }
    list(owner = owner, role = role,
         convergence_fun = definition$convergence_fun,
         current = current, reason = reason)
}

#' Declare assessment support for a built-in model process
#'
#' Register a convergence callback or a prescribed-input declaration for the
#' current consumers or resource process. The declaration is stored in the model
#' together with the identity of its current dynamics provider. Replacing that
#' provider or changing the projection dispatch chain invalidates the
#' declaration until support is explicitly adopted again. Changing parameters
#' or consumer rate functions leaves the identity of the fish equations intact.
#' Calling `setProcess()` explicitly adopts the current provider, including when
#' omitted arguments retain the previous callback and role.
#'
#' A convergence callback returns a data frame with the columns `measure`,
#' `value`, `tolerance` and `units`. It supplies measurements, not verdicts;
#' [getSteadyAssessment()] compares each absolute value with its tolerance.
#' Callback arguments include the model, state, complete rate list, derivative,
#' time and process name; accept `...` for arguments not needed by the callback.
#'
#' A new `definition` replaces the previous callback and role before explicit
#' arguments are applied. Without `definition`, omitted arguments are retained.
#' The default role is `"dynamic"`. A `"prescribed"` role declares the whole
#' state block as an input rather than a dynamic variable. Partial masks are not
#' supported by this interface. Supplying `convergence_fun = NULL` removes its
#' callback. State and equations are not changed by this registration.
#'
#' @param params A MizerParams object.
#' @param process Either `"consumers"` or `"resource"`.
#' @param definition A named list with `convergence_fun` and/or `role`.
#' @param convergence_fun Name of an assessment callback, optionally qualified
#'   as `package::exported_function`, or `NULL`.
#' @param role Either `"dynamic"` or `"prescribed"`.
#' @return `setProcess()` returns the updated model. `getProcess()` returns
#'   the callback, role, reusable `definition` and recorded provider `binding`.
#'   With no `process` argument, `getProcess()` returns all stored processes.
#' @seealso [setComponent()], [getSteadyAssessment()]
#' @export
#' @family extension tools
setProcess <- function(params, process, definition, convergence_fun, role) {
    params <- validParams(params)
    if (!is.character(process) || length(process) != 1L || is.na(process) ||
            !process %in% c("consumers", "resource")) {
        stop("`process` must be 'consumers' or 'resource'.", call. = FALSE)
    }
    fields <- c("convergence_fun", "role")
    assembled <- list(convergence_fun = NULL, role = "dynamic")
    existing <- params$capabilities$processes[[process]]
    if (missing(definition)) {
        if (!is.null(existing)) assembled[fields] <- existing[fields]
    } else {
        if (!is.list(definition) || is.object(definition) ||
                !capability_names_valid(definition) ||
                any(!names(definition) %in% fields)) {
            stop("A process definition contains only convergence_fun and role.",
                 call. = FALSE)
        }
        assembled[names(definition)] <- definition
    }
    if (!missing(convergence_fun)) assembled["convergence_fun"] <-
        list(convergence_fun)
    if (!missing(role)) assembled["role"] <- list(role)
    errors <- component_definition_errors(assembled, required = FALSE)
    if (!isTRUE(errors)) stop(paste(errors, collapse = "\n"), call. = FALSE)
    if (!is.null(assembled$convergence_fun)) {
        resolve_model_callback(assembled$convergence_fun,
                               paste0("Assessment for process '", process, "'"))
    }
    assembled$binding <- process_provider_binding(params, process)
    params$capabilities$processes[[process]] <- assembled
    params$time_modified <- lubridate::now()
    params
}

#' @rdname setProcess
#' @export
getProcess <- function(params, process) {
    params <- validParams(params)
    if (missing(process)) {
        return(lapply(names(params$capabilities$processes), function(id) {
            getProcess(params, id)
        }) |> stats::setNames(names(params$capabilities$processes)))
    }
    if (!is.character(process) || length(process) != 1L || is.na(process) ||
            !process %in% c("consumers", "resource")) {
        stop("`process` must be 'consumers' or 'resource'.", call. = FALSE)
    }
    specification <- params$capabilities$processes[[process]]
    if (is.null(specification)) return(NULL)
    list(convergence_fun = specification$convergence_fun,
         role = specification$role,
         definition = specification[c("convergence_fun", "role")],
         binding = specification$binding)
}
