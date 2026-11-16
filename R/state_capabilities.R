# Structural inputs, temporal dependence and resource equations are separate
# provider-bound capabilities. None is inferred from a zero observed rate.

state_binding_valid <- function(binding) {
    if (is.null(binding)) return(TRUE)
    is.list(binding) && !is.object(binding) &&
        capability_names_valid(binding) && tuning_data_valid(binding)
}

process_state_errors <- function(definition, process) {
    if (process == "reproduction" &&
            (!is.null(definition$prescribed_mask_fun) ||
                 !is.null(definition$autonomous) ||
                 !is.null(definition$rate_fun) ||
                 !is.null(definition$prescribed_binding) ||
                 !is.null(definition$temporal_binding))) {
        return(paste0("Reproduction has no state capabilities; declare ",
                      "temporal support on the consumers equations."))
    }
    if (process != "resource" &&
            (!is.null(definition$rate_fun) ||
                 !is.null(definition$equation_binding))) {
        return("A process rate_fun belongs to the resource equations.")
    }
    if (!state_binding_valid(definition$equation_binding)) {
        return("Invalid stored resource equation binding.")
    }
    TRUE
}

state_provider_binding <- function(params, owner) {
    if (!startsWith(owner, "component:")) {
        return(process_provider_binding(params, owner))
    }
    component <- substring(owner, 11L)
    list(dispatch = process_dispatch_binding(params),
         dynamics_fun = params$other_dynamics[[component]],
         rate_fun = params$capabilities$components[[component]]$rate_fun)
}

state_binding_current <- function(params, definition, owner, field) {
    binding <- definition[[field]]
    if (is.null(binding)) return(FALSE)
    binding$automatic <- NULL
    identical(binding, state_provider_binding(params, owner))
}

temporal_provider_binding <- function(params, owner) {
    binding <- relevant_provider_binding(tuning_provider_binding(params), owner)
    # Balancing changes no equation and therefore is not a temporal provider.
    binding$balance_methods <- NULL
    binding$other_encounter <- lapply(params$other_encounter, as.character)
    binding$other_mort <- lapply(params$other_mort, as.character)
    binding
}

temporal_binding_current <- function(params, definition, owner) {
    binding <- definition$temporal_binding
    if (is.null(binding)) return(FALSE)
    binding$automatic <- NULL
    identical(binding, temporal_provider_binding(params, owner))
}

builtin_autonomous_rates <- function(params) {
    law <- sub("^mizer::", "", params$rates_funcs$RDD)
    builtin_tuning_rates(params, reproduction = TRUE) && law %in%
        c("BevertonHoltRDD", "RickerRDD", "SheperdRDD", "constantRDD", "noRDD")
}

default_provider_state <- function(params, owner) {
    state_binding <- state_provider_binding(params, owner)
    temporal_binding <- temporal_provider_binding(params, owner)
    resource <- owner == "resource"
    law <- sub("^mizer::", "", params$resource_dynamics)
    known_resource <- builtin_consumer_projection(params) && law %in%
        c("resource_semichemostat", "resource_logistic", "resource_constant")
    known_time <- builtin_autonomous_rates(params) &&
        (owner == "consumers" || (resource && known_resource))
    list(prescribed_mask_fun = NULL,
         prescribed_binding = c(state_binding, list(automatic = TRUE)),
         autonomous = if (known_time) TRUE else NULL,
         temporal_binding = c(temporal_binding, list(automatic = TRUE)),
         rate_fun = if (resource && known_resource) {
             paste0("equation_", law)
         } else {
             NULL
         },
         equation_binding = if (resource) {
             c(state_binding, list(automatic = TRUE))
         } else {
             NULL
         })
}

refresh_builtin_state <- function(params) {
    for (owner in c("consumers", "resource")) {
        definition <- params$capabilities$processes[[owner]]
        if (is.null(definition)) next
        defaults <- default_provider_state(params, owner)
        for (family in c("prescribed", "temporal", "equation")) {
            field <- paste0(family, "_binding")
            if (!isTRUE(definition[[field]]$automatic)) next
            current <- if (family == "temporal") {
                temporal_binding_current(params, definition, owner)
            } else {
                state_binding_current(params, definition, owner, field)
            }
            if (current) next
            if (family == "temporal" && !is.null(defaults$autonomous)) {
                fields <- c("autonomous", field)
            } else if (family == "equation" && !is.null(defaults$rate_fun)) {
                fields <- c("rate_fun", field)
            } else if (family == "prescribed" &&
                           builtin_consumer_projection(params) &&
                           (owner == "consumers" ||
                                !is.null(defaults$rate_fun))) {
                fields <- c("prescribed_mask_fun", field)
            } else {
                next
            }
            definition[fields] <- defaults[fields]
        }
        params$capabilities$processes[[owner]] <- definition
    }
    params
}

state_mask_like <- function(state, value = FALSE) {
    mask <- rep(value, length(state))
    if (!is.null(dim(state))) {
        dim(mask) <- dim(state)
        dimnames(mask) <- dimnames(state)
    }
    names(mask) <- names(state)
    mask
}

state_mask_valid <- function(mask, state) {
    is.logical(mask) && !anyNA(mask) && length(mask) == length(state) &&
        identical(dim(mask), dim(state)) &&
        identical(dimnames(mask), dimnames(state)) &&
        identical(names(mask), names(state))
}

prepare_state_capabilities <- function(params, fixed = character()) {
    params <- validParams(params)
    owners <- assessment_owners(params)
    if (!is.character(fixed) || anyNA(fixed) || anyDuplicated(fixed) ||
            any(!fixed %in% owners)) {
        stop("`fixed` must contain unique canonical state-owner IDs.")
    }
    result <- setNames(vector("list", length(owners)), owners)
    for (owner in owners) {
        component <- startsWith(owner, "component:")
        id <- if (component) substring(owner, 11L) else owner
        definition <- if (component) params$capabilities$components[[id]] else
            params$capabilities$processes[[owner]]
        state <- if (component) params$initial_n_other[[id]] else if (
            owner == "consumers"
        ) params$initial_n else params$initial_n_pp
        numeric_state <- is.numeric(state) && !is.complex(state)
        structural <- if (numeric_state) state_mask_like(state, TRUE) else TRUE
        if (owner == "consumers" && builtin_consumer_projection(params)) {
            structural[] <- col(structural) >= params$w_min_idx &
                col(structural) <= support_top_idx(params)
        }
        prescribed <- if (numeric_state) state_mask_like(state) else FALSE
        current <- state_binding_current(params, definition, owner,
                                         "prescribed_binding")
        reason <- if (!current) {
            "The prescribed-input declaration does not match the provider."
        } else {
            NULL
        }
        role <- definition$role %||% "dynamic"
        held <- owner %in% fixed
        if (!held && current) {
            if (role == "prescribed") {
                prescribed[] <- TRUE
            } else if (!is.null(definition$prescribed_mask_fun)) {
                if (!numeric_state) {
                    stop("Partial prescribed masks require numeric state for ",
                         owner, ".")
                }
                prescribed <- invoke_tuning_callback(
                    definition$prescribed_mask_fun, params,
                    list(owner = owner, component = if (component) id else NULL,
                         process = if (!component) owner else NULL),
                    paste0("Prescribed-input mask for ", owner)
                )
                if (!state_mask_valid(prescribed, state)) {
                    stop("Prescribed-input mask for ", owner,
                         " must be logical without NA and exactly match ",
                         "the state dimensions and names.")
                }
            }
        }
        temporal_current <- temporal_binding_current(params, definition, owner)
        autonomous <- if (temporal_current) definition$autonomous else NULL
        temporal_reason <- if (!temporal_current) {
            "The temporal declaration does not match the current provider."
        } else if (is.null(autonomous)) {
            "No autonomous/nonautonomous declaration is registered."
        } else if (!autonomous) {
            "The provider is explicitly nonautonomous."
        } else {
            NULL
        }
        equation_current <- owner == "resource" && state_binding_current(
            params, definition, owner, "equation_binding"
        )
        rate_fun <- if (equation_current) definition$rate_fun else NULL
        equation_reason <- if (owner != "resource") NULL else if (
            !equation_current
        ) {
            "The resource equation declaration does not match its provider."
        } else if (is.null(rate_fun)) {
            "No resource equation callback is registered."
        } else {
            NULL
        }
        contributes <- component &&
            (!is.null(params$other_encounter[[id]]) ||
                 !is.null(params$other_mort[[id]]))
        dynamic <- structural & !prescribed
        if (held) dynamic[] <- FALSE
        result[[owner]] <- list(
            owner = owner, role = role, fixed = held,
            numeric_state = numeric_state, structural_mask = structural,
            prescribed_mask = prescribed, dynamic_mask = dynamic,
            mask_current = current, mask_reason = reason,
            autonomous = autonomous, temporal_current = temporal_current,
            temporal_required = !held || contributes,
            temporal_reason = temporal_reason,
            rate_fun = rate_fun, equation_current = equation_current,
            equation_reason = equation_reason
        )
    }
    result
}

evaluate_resource_equation <- function(params, specification, n, n_pp,
                                       n_other, rates, t) {
    if (!isTRUE(specification$equation_current) ||
            is.null(specification$rate_fun)) {
        stop(specification$equation_reason %||%
                 "No current resource equation callback is registered.")
    }
    derivative <- invoke_tuning_callback(
        specification$rate_fun, params,
        list(n = n, n_pp = n_pp, n_other = n_other, rates = rates, t = t),
        "Resource equation"
    )
    if (!is.numeric(derivative) || is.complex(derivative) ||
            any(!is.finite(derivative)) ||
            length(derivative) != length(n_pp) ||
            !identical(dim(derivative), dim(n_pp)) ||
            !identical(dimnames(derivative), dimnames(n_pp)) ||
            !identical(names(derivative), names(n_pp))) {
        stop("Resource equation must return finite numeric derivatives with ",
             "the exact resource state dimensions and names.")
    }
    derivative
}

# nolint next: object_length_linter. Descriptive stored callback identifier.
equation_resource_semichemostat <- function(params, n_pp, rates, ...) {
    params$rr_pp * (params$cc_pp - n_pp) - rates$resource_mort * n_pp
}

equation_resource_logistic <- function(params, n_pp, rates, ...) {
    result <- -rates$resource_mort * n_pp
    positive <- n_pp != 0 & params$rr_pp != 0
    if (any(positive & params$cc_pp <= 0)) {
        stop("Nonzero resource abundance requires positive capacity.")
    }
    density_factor <- 1 - n_pp[positive] / params$cc_pp[positive]
    result[positive] <- (params$rr_pp[positive] * density_factor -
                             rates$resource_mort[positive]) * n_pp[positive]
    result
}

equation_resource_constant <- function(params, n_pp, ...) n_pp * 0

# The logistic step can collapse positive density instantaneously at zero
# capacity. That map has no finite continuous equation on such a neighbourhood.
# Keep the coordinates and report the domain limitation rather than declaring
# them prescribed from observed zero density or zero capacity.
resource_equation_domain <- function(params, specification, dynamic_mask) {
    if (isTRUE(specification$equation_current) &&
            identical(specification$rate_fun, "equation_resource_logistic") &&
            any(dynamic_mask & params$cc_pp <= 0 & params$rr_pp != 0)) {
        return(paste0("The logistic resource equation has no finite ",
                      "neighbourhood at active zero-capacity bins with ",
                      "nonzero replenishment."))
    }
    NULL
}
