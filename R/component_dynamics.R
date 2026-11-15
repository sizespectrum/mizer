# Resolve stored identifiers without retaining namespaces or function bodies in
# the saved model. Availability is checked when an operation is prepared.
resolve_model_callback <- function(name, context = "Model callback") {
    if (!is.character(name) || length(name) != 1L ||
            is.na(name) || !nzchar(name)) {
        stop(context, " must be a function name.", call. = FALSE)
    }
    parts <- strsplit(name, "::", fixed = TRUE)[[1L]]
    if (grepl("::", name, fixed = TRUE)) {
        if (length(parts) != 2L || any(!nzchar(parts)) ||
                grepl(":", parts[[2L]], fixed = TRUE)) {
            stop(context, " must use 'package::exported_function': '",
                 name, "'.", call. = FALSE)
        }
        if (!requireNamespace(parts[[1L]], quietly = TRUE)) {
            stop(context, " requires package '", parts[[1L]],
                 "' for '", name, "'.", call. = FALSE)
        }
        fun <- tryCatch(getExportedValue(parts[[1L]], parts[[2L]]),
                        error = function(e) NULL)
    } else {
        fun <- get0(name, mode = "function")
    }
    if (!is.function(fun)) {
        stop(context, " '", name, "' was not found as an available function.",
             call. = FALSE)
    }
    fun
}

# A temporary override, e.g. constant_other during tuning, takes priority over
# a stored equation. Declaring an equation never bypasses that active override.
active_component_rate_name <- function(params, component) {
    definition <- existing_component_definition(params, component)
    expected <- definition$dynamics_fun
    if (is.null(expected)) expected <- "component_rate_dynamics"
    if (identical(params$other_dynamics[[component]], expected)) {
        definition$rate_fun
    } else {
        NULL
    }
}

prepare_component_rates <- function(params) {
    ids <- names(params$other_dynamics)
    result <- setNames(vector("list", length(ids)), ids)
    for (component in ids) {
        name <- active_component_rate_name(params, component)
        if (!is.null(name)) {
            result[[component]] <- resolve_model_callback(
                name = name,
                context = paste0("Equation for component '", component, "'")
            )
        }
    }
    result
}

checked_component_rate <- function(fun, params, n, n_pp, n_other, rates, t,
                                   component, ...) {
    state <- n_other[[component]]
    if (!is.numeric(state) || any(!is.finite(state))) {
        stop("Equation-based component '", component,
             "' needs a finite numeric state.", call. = FALSE)
    }
    derivative <- fun(params, n = n, n_pp = n_pp, n_other = n_other,
                      rates = rates, t = t, component = component, ...)
    if (!is.numeric(derivative) || any(!is.finite(derivative))) {
        stop("Equation for component '", component,
             "' must return a finite numeric derivative.", call. = FALSE)
    }
    if (length(derivative) != length(state) ||
            !identical(dim(derivative), dim(state)) ||
            !identical(dimnames(derivative), dimnames(state)) ||
            !identical(names(derivative), names(state))) {
        stop("Equation for component '", component,
             "' must preserve the state's dimensions and names.", call. = FALSE)
    }
    derivative
}

make_component_rate_step <- function(rate_fun) {
    force(rate_fun)
    step <- function(params, n, n_pp, n_other, rates, t, dt, component, ...) {
        derivative <- checked_component_rate(
            rate_fun, params, n, n_pp, n_other, rates, t, component, ...
        )
        n_other[[component]] + dt * derivative
    }
    attr(step, "component_rate") <- rate_fun
    step
}

prepare_component_dynamics <- function(params) {
    ids <- names(params$other_dynamics)
    result <- setNames(vector("list", length(ids)), ids)
    for (component in ids) {
        name <- params$other_dynamics[[component]]
        if (identical(name, "component_rate_dynamics")) {
            rate_name <- active_component_rate_name(params, component)
            if (is.null(rate_name)) {
                stop("Component '", component,
                     "' has no active equation. Update it with setComponent().",
                     call. = FALSE)
            }
            rate_fun <- resolve_model_callback(
                name = rate_name,
                context = paste0("Equation for component '", component, "'")
            )
            result[[component]] <- make_component_rate_step(rate_fun)
        } else {
            result[[component]] <- resolve_model_callback(
                name = name,
                context = paste0("Step for component '", component, "'")
            )
        }
    }
    result
}

# The legacy registry still determines state storage. This named bridge also
# supports callers that resolve it directly; project() prepares its rate once
# and performs the Heun corrector at the full predicted model state.
component_rate_dynamics <- function(params, n, n_pp, n_other, rates, t, dt,
                                    component, ...) {
    step <- prepare_component_dynamics(params)[[component]]
    step(params, n = n, n_pp = n_pp, n_other = n_other, rates = rates,
         t = t, dt = dt, component = component, ...)
}

invoke_model_step <- function(fun, args, predicted_state = NULL,
                              endpoint_rates = NULL) {
    optional <- list(predicted_state = predicted_state,
                     endpoint_rates = endpoint_rates)
    accepted <- names(formals(fun))
    if (!"..." %in% accepted) {
        optional <- optional[intersect(names(optional), accepted)]
    }
    # Steps historically receive the model positionally, so their first formal
    # can be called p/object rather than params. Keep that calling convention.
    arg_names <- names(args)
    if (length(args) && !is.null(arg_names) &&
            identical(arg_names[[1L]], "params")) {
        arg_names[[1L]] <- ""
        names(args) <- arg_names
    }
    # These are supplied by the numerical stage, not a stale value in `...`.
    if (!is.null(arg_names)) args <- args[!arg_names %in% names(optional)]
    do.call(fun, c(args, optional))
}

component_integrators <- function(params, method) {
    method <- normalise_project_method(method)
    ids <- names(params$other_dynamics)
    if (is.null(ids)) ids <- character()
    result <- data.frame(component = ids, method = rep("custom", length(ids)),
                         order = rep(NA_integer_, length(ids)),
                         callback = vapply(params$other_dynamics, identity,
                                           character(1)),
                         stringsAsFactors = FALSE, row.names = NULL)
    for (i in seq_along(ids)) {
        definition <- existing_component_definition(params, ids[[i]])
        active <- params$other_dynamics[[ids[[i]]]]
        if (identical(active, "component_rate_dynamics")) {
            result$method[[i]] <- if (method == "euler") "euler" else "heun"
            result$order[[i]] <- if (method == "euler") 1L else 2L
            result$callback[[i]] <- definition$rate_fun
        } else if (identical(active, definition$dynamics_fun) &&
                       !is.null(definition$integration_order)) {
            result$order[[i]] <- as.integer(definition$integration_order)
        }
    }
    result
}
