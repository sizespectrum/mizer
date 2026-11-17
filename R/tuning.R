# Author-declared tuning transactions. Stored declarations contain data only;
# resolved functions and captured targets exist only for one operation.

tuning_data_valid <- function(x) {
    if (is.null(x)) {
        return(TRUE)
    }
    if (is.list(x) && !is.object(x)) {
        return(capability_names_valid(x) && all(vapply(
            x, tuning_data_valid,
            logical(1)
        )))
    }
    is.atomic(x) && !is.object(x) && !is.complex(x) &&
        !is.raw(x) && !anyNA(x)
}

tuning_strategy_defaults <- function() {
    list(
        action = "none", hold = list(states = character(), rates = character()),
        evolve = character(), adjust = list(), preserve = list(),
        balance_fun = NULL, target_fun = NULL,
        solvers = c("project", "newton", "balance"),
        partial_species = TRUE, settings = list()
    )
}

tuning_parameter_path_allowed <- function(path) {
    protected <- c("initial_n", "initial_n_pp", "initial_n_other",
                   "capabilities", "rates_funcs", "other_dynamics",
                   "resource_dynamics",
                   "extensions", "mizer_version", "w", "w_full",
                   "dw", "dw_full",
                   "w_min_idx", "w_max_idx", "second_order_w", "sc", "ft_mask",
                   "use_predation_diffusion")
    !path[[1L]] %in% protected &&
        !(length(path) > 1L && path[[1L]] == "species_params" &&
              path[[2L]] %in% c("species", "w_min", "w_max"))
}

tuning_strategy_errors <- function(strategy) {
    fields <- names(tuning_strategy_defaults())
    if (!is.list(strategy) || is.object(strategy) ||
            !capability_names_valid(strategy) ||
            any(!names(strategy) %in% fields)) {
        return("A tuning strategy must be a named list of supported fields.")
    }
    x <- tuning_strategy_defaults()
    x[names(strategy)] <- strategy
    ids <- function(value) {
        is.character(value) && !anyNA(value) &&
            all(nzchar(value)) && !anyDuplicated(value)
    }
    if (!identical(x$action, "none") && !identical(x$action, "balance")) {
        return("Tuning `action` must be 'none' or 'balance'.")
    }
    if (!is.list(x$hold) || !capability_names_valid(x$hold) ||
            !setequal(names(x$hold), c("states", "rates")) ||
            !ids(x$hold$states) || !ids(x$hold$rates) || !ids(x$evolve) ||
            any(!x$hold$rates %in% "consumers/recruitment")) {
        return(paste0("Tuning hold/evolve declarations require unique IDs; ",
                      "the supported held rate is 'consumers/recruitment'."))
    }
    if (!ids(x$solvers) || any(!x$solvers %in% c(
        "project", "newton",
        "balance"
    )) ||
        !is.logical(x$partial_species) || length(x$partial_species) != 1L ||
        is.na(x$partial_species)) {
        return("Invalid tuning solver or partial_species declaration.")
    }
    for (field in c("balance_fun", "target_fun")) {
        value <- x[[field]]
        if (!is.null(value) && (!is.character(value) || length(value) != 1L ||
                                    is.na(value) || !nzchar(value))) {
            return(paste0("Tuning `", field, "` must name a function."))
        }
    }
    if (x$action == "none" && length(x$adjust)) {
        return("A no-action strategy cannot authorise parameter adjustments.")
    }
    if ((x$action == "balance" && is.null(x$balance_fun)) ||
            (x$action == "none" && !is.null(x$balance_fun))) {
        return("Balance requires balance_fun; no-action must omit it.")
    }
    for (field in c("adjust", "preserve")) {
        descriptors <- x[[field]]
        if (!is.list(descriptors) || is.object(descriptors) ||
                !capability_names_valid(descriptors)) {
            return(paste0(
                "Tuning `", field,
                "` must be a named descriptor list."
            ))
        }
        for (descriptor in descriptors) {
            allowed <- if (field == "adjust") {
                c("path", "lower", "upper", "species")
            } else {
                c("path", "fun")
            }
            if (!is.list(descriptor) || is.object(descriptor) ||
                    !capability_names_valid(descriptor) ||
                    any(!names(descriptor) %in% allowed)) {
                return(paste0("Invalid tuning ", field, " descriptor."))
            }
            path <- descriptor$path
            if (is.null(path)) {
                if (field == "adjust" || !is.character(descriptor$fun) ||
                        length(descriptor$fun) != 1L || is.na(descriptor$fun) ||
                        !nzchar(descriptor$fun)) {
                    return("A descriptor needs a path or value function.")
                }
            } else if (!is.character(path) || !length(path) || anyNA(path) ||
                           any(!nzchar(path)) || !is.null(descriptor$fun)) {
                return("A parameter path must be a nonempty character vector.")
            }
            if (field == "adjust") {
                if (!tuning_parameter_path_allowed(path)) {
                    return("Cannot adjust state, providers or grid structure.")
                }
                for (bound in c("lower", "upper")) {
                    value <- descriptor[[bound]]
                    if (!is.null(value) && (!is.numeric(value) ||
                                                is.complex(value) ||
                                                !length(value) ||
                                                anyNA(value))) {
                        return("Parameter bounds must be numeric without NAs.")
                    }
                }
                by_species <- descriptor$species
                if (!is.null(by_species) && (!is.logical(by_species) ||
                                                 length(by_species) != 1L ||
                                                 is.na(by_species))) {
                    return("The species selector must be TRUE or FALSE.")
                }
            }
        }
    }
    if (!is.list(x$settings) || !tuning_data_valid(x$settings)) {
        return("Tuning settings must contain serialisable named data only.")
    }
    TRUE
}

normalise_tuning_strategy <- function(strategy) {
    errors <- tuning_strategy_errors(strategy)
    if (!isTRUE(errors)) stop(errors, call. = FALSE)
    result <- tuning_strategy_defaults()
    result[names(strategy)] <- strategy
    result
}

tuning_block_errors <- function(tuning) {
    if (is.null(tuning)) {
        return(TRUE)
    }
    if (!is.list(tuning) || is.object(tuning) ||
            !capability_names_valid(tuning) ||
            !setequal(names(tuning), c("default", "strategies")) ||
            !is.list(tuning$strategies) ||
            !capability_names_valid(tuning$strategies) ||
            !is.character(tuning$default) || length(tuning$default) != 1L ||
            is.na(tuning$default) ||
            !tuning$default %in% names(tuning$strategies)) {
        return("Tuning requires named strategies and one declared default.")
    }
    errors <- lapply(tuning$strategies, tuning_strategy_errors)
    errors <- errors[!vapply(errors, isTRUE, logical(1))]
    if (length(errors)) unlist(errors, use.names = FALSE) else TRUE
}

tuning_method_names <- function(params, generics) {
    classes <- unique(c(
        class(params),
        names(extensionRequirements(params$extensions))
    ))
    unlist(lapply(generics, function(generic) {
        unlist(lapply(classes, function(id) {
            method <- getS3method(generic, id, optional = TRUE)
            if (is.null(method)) {
                return(character())
            }
            paste0(generic, ".", id, "@", environmentName(environment(method)))
        }), use.names = FALSE)
    }), use.names = FALSE)
}

tuning_provider_binding <- function(params) {
    list(
        dispatch = process_dispatch_binding(params),
        rate_functions = lapply(params$rates_funcs, as.character),
        balance_methods = tuning_method_names(params, "balanceRDD"),
        rate_methods = tuning_method_names(params, c(
            "projectRates", "projectEncounter", "projectFeedingLevel",
            "projectEReproAndGrowth", "projectERepro", "projectEGrowth",
            "projectDiffusion", "projectPredRate", "projectPredMort",
            "projectFMort", "projectMort", "projectRDI", "projectRDD",
            "projectResourceMort"
        )),
        resource_dynamics = as.character(params$resource_dynamics),
        other_dynamics = lapply(params$other_dynamics, as.character),
        component_equations = lapply(
            params$capabilities$components, `[[`,
            "rate_fun"
        ),
        other_encounter = lapply(params$other_encounter, as.character),
        other_mort = lapply(params$other_mort, as.character)
    )
}

relevant_provider_binding <- function(binding, owner) {
    components <- names(binding$other_dynamics)
    if (owner %in% c("consumers", "reproduction")) {
        binding$resource_dynamics <- NULL
    }
    own <- if (startsWith(owner, "component:")) {
        substring(owner, 11L)
    } else {
        character()
    }
    binding$other_dynamics <- binding$other_dynamics[
        intersect(own, components)
    ]
    binding$component_equations <-
        binding$component_equations[intersect(
            own,
            names(binding$component_equations)
        )]
    for (slot in c("other_encounter", "other_mort")) {
        keep <- setdiff(names(binding[[slot]]), setdiff(components, own))
        binding[[slot]] <- binding[[slot]][keep]
    }
    for (field in c("other_dynamics", "component_equations",
                    "other_encounter", "other_mort")) {
        if (!length(binding[[field]])) binding[[field]] <- list()
    }
    binding
}

tuning_binding_current <- function(params, definition, owner) {
    stored <- definition$tuning_binding
    if (is.null(stored)) return(FALSE)
    current <- tuning_provider_binding(params)
    identical(relevant_provider_binding(stored, owner),
              relevant_provider_binding(current, owner))
}

builtin_tuning_rates <- function(params, reproduction = FALSE) {
    expected <- list(
        Rates = "mizerRates", Encounter = "mizerEncounter",
        FeedingLevel = "mizerFeedingLevel",
        EReproAndGrowth = "mizerEReproAndGrowth",
        PredRate = "mizerPredRate", PredMort = "mizerPredMort",
        FMort = "mizerFMort",
        Mort = "mizerMort", ERepro = "mizerERepro", EGrowth = "mizerEGrowth",
        Diffusion = "mizerDiffusion", ResourceMort = "mizerResourceMort"
    )
    if (reproduction) expected$RDI <- "mizerRDI"
    named <- all(vapply(names(expected), function(id) {
        identical(params$rates_funcs[[id]], expected[[id]])
    }, logical(1)))
    methods <- tuning_provider_binding(params)$rate_methods
    if (!reproduction) methods <- methods[!grepl("^projectRD[ID]\\.", methods)]
    named && all(endsWith(methods, "@mizer")) &&
        !length(freeRateContributions(params, "other_mort")) &&
        !length(freeRateContributions(params, "other_encounter")) &&
        builtin_consumer_projection(params)
}

builtin_no_action_tuning <- function(owner, prescribed = FALSE) {
    strategy <- tuning_strategy_defaults()
    if (prescribed) strategy$hold$states <- owner else strategy$evolve <- owner
    list(default = "unchanged", strategies = list(unchanged = strategy))
}

default_provider_tuning <- function(params, owner) {
    result <- list(
        stationary = NULL, equation_rates_fun = NULL, tuning = NULL,
        tuning_binding = tuning_provider_binding(params), level_fun = NULL,
        level_binding = c(tuning_provider_binding(params),
                          list(automatic = TRUE))
    )
    known <- builtin_tuning_rates(params,
        reproduction = owner == "reproduction"
    )
    if (owner == "consumers" && known) {
        result$stationary <- TRUE
        result$equation_rates_fun <- "builtin_consumer_equation_rates"
        result$tuning <- builtin_no_action_tuning(owner)
    } else if (owner == "resource" && known) {
        name <- sub("^mizer::", "", params$resource_dynamics)
        if (name == "resource_constant") {
            result$stationary <- TRUE
            result$equation_rates_fun <- "builtin_prescribed_equation_rates"
            result$tuning <- builtin_no_action_tuning(owner, prescribed = TRUE)
        } else if (name %in% c("resource_semichemostat", "resource_logistic")) {
            result$stationary <- TRUE
            result$equation_rates_fun <- "builtin_resource_equation_rates"
            strategy <- tuning_strategy_defaults()
            strategy$action <- "balance"
            strategy$hold$states <- "resource"
            strategy$adjust <- list(capacity = list(
                path = "cc_pp", lower = 0,
                upper = Inf, species = FALSE
            ))
            strategy$balance_fun <- "balance_resource_strategy"
            strategy$target_fun <- "target_resource_strategy"
            result$tuning <- list(
                default = "capacity",
                strategies = list(capacity = strategy)
            )
        }
    } else if (owner == "reproduction" && known &&
        sub("^mizer::", "", params$rates_funcs$RDD) %in% c(
            "BevertonHoltRDD", "RickerRDD",
            "SheperdRDD", "constantRDD", "noRDD"
        )) {
        result$stationary <- TRUE
        result$equation_rates_fun <- "builtin_reproduction_equation_rates"
        result$tuning <- builtin_reproduction_tuning(params)
        result$level_fun <- paste0(
            "reproduction_level_",
            sub("^mizer::", "", params$rates_funcs$RDD)
        )
    }
    result
}

automatic_tuning_definition <- function(params, definition, owner) {
    binding <- definition$tuning_binding
    if (is.null(binding)) return(FALSE)
    previous <- params
    previous$rates_funcs <- binding$rate_functions
    previous$resource_dynamics <- binding$resource_dynamics
    previous$other_dynamics <- binding$other_dynamics
    previous$other_encounter <- binding$other_encounter
    previous$other_mort <- binding$other_mort
    class(previous) <- binding$dispatch$classes
    defaults <- default_provider_tuning(previous, owner)
    fields <- c("stationary", "equation_rates_fun", "tuning")
    !is.null(defaults$tuning) && identical(definition[fields], defaults[fields])
}

refresh_builtin_tuning <- function(params) {
    for (owner in c("consumers", "resource", "reproduction")) {
        old <- params$capabilities$processes[[owner]]
        if (is.null(old)) next
        if (tuning_binding_current(params, old, owner)) next
        if (!automatic_tuning_definition(params, old, owner)) next
        defaults <- default_provider_tuning(params, owner)
        # A stale explicit declaration remains stale unless this is a known
        # built-in route with a matching built-in adapter.
        if (!is.null(defaults$tuning)) {
            fields <- c("stationary", "equation_rates_fun", "tuning",
                        "tuning_binding")
            old[fields] <- defaults[fields]
            params$capabilities$processes[[owner]] <- old
        }
    }
    params <- refresh_builtin_levels(params)
    params <- refresh_builtin_mortality(params)
    params <- refresh_builtin_state(params)
    refresh_builtin_transformations(params)
}

invoke_tuning_callback <- function(fun, params, args = list(), context) {
    if (is.character(fun)) fun <- resolve_model_callback(fun, context)
    args <- c(list(params), args)
    formals <- names(formals(fun))
    if (!"..." %in% formals) {
        args <- args[c(
            TRUE,
            names(args)[-1L] %in% formals
        )]
    }
    tryCatch(do.call(fun, args), error = function(error) {
        stop(context, ": ", conditionMessage(error), call. = FALSE)
    })
}

# nolint next: object_length_linter. Stored descriptive callback identifier.
builtin_consumer_equation_rates <- function(params, rates, ...) {
    list(
        growth = rates$e_growth, diffusion = rates$diffusion,
        mortality = rates$mort, recruitment = rates$rdd
    )
}

# nolint next: object_length_linter. Stored descriptive callback identifier.
builtin_reproduction_equation_rates <- function(params, rates, ...) {
    list(recruitment = rates$rdd)
}

# nolint next: object_length_linter. Stored descriptive callback identifier.
builtin_prescribed_equation_rates <- function(params, ...) list()

# nolint next: object_length_linter. Stored descriptive callback identifier.
builtin_resource_equation_rates <- function(params, n_pp, rates, ...) {
    law <- sub("^mizer::", "", params$resource_dynamics)
    if (law == "resource_semichemostat") {
        list(
            source = params$rr_pp * params$cc_pp,
            loss = (params$rr_pp + rates$resource_mort) * n_pp
        )
    } else if (law == "resource_logistic") {
        positive <- n_pp > 0
        if (any(positive & params$cc_pp <= 0)) {
            stop("Positive resource abundance requires positive capacity.")
        }
        loss <- n_pp * 0
        density_loss <- params$rr_pp[positive] * n_pp[positive] /
            params$cc_pp[positive]
        loss[positive] <- (density_loss + rates$resource_mort[positive]) *
            n_pp[positive]
        list(growth = params$rr_pp * n_pp, loss = loss)
    } else {
        stop("No built-in effective resource equation adapter.")
    }
}

target_resource_strategy <- function(params, n_pp, rates, ...) {
    law <- sub("^mizer::", "", params$resource_dynamics)
    actual <- builtin_resource_equation_rates(params, n_pp, rates)
    if (law == "resource_semichemostat") {
        actual$source <- actual$loss
    } else {
        actual$loss <- actual$growth
    }
    actual
}

balance_resource_strategy <- function(params, rates, n_pp, ...) {
    law <- sub("^mizer::", "", params$resource_dynamics)
    if (law == "resource_semichemostat") {
        numerator <- (params$rr_pp + rates$resource_mort) * n_pp
        invalid <- params$rr_pp == 0 & numerator != 0
        if (any(invalid)) stop("Resource capacity cannot balance zero renewal.")
        capacity <- params$cc_pp
        valid <- params$rr_pp != 0
        capacity[valid] <- numerator[valid] / params$rr_pp[valid]
    } else if (law == "resource_logistic") {
        denominator <- params$rr_pp - rates$resource_mort
        if (any(n_pp > 0 & denominator <= 0)) {
            stop("Resource capacity cannot balance mortality at this state.")
        }
        capacity <- params$cc_pp
        positive <- n_pp > 0
        capacity[positive] <-
            params$rr_pp[positive] * n_pp[positive] / denominator[positive]
    } else {
        stop("No built-in resource balancing strategy.")
    }
    params$cc_pp[] <- capacity
    comment(params$cc_pp) <- NULL
    params
}

tuning_owners <- function(params) c(assessment_owners(params), "reproduction")

tuning_definition <- function(params, owner) {
    if (startsWith(owner, "component:")) {
        params$capabilities$components[[substring(owner, 11L)]]
    } else {
        params$capabilities$processes[[owner]]
    }
}

tuning_state <- function(params) {
    list(
        n = params$initial_n, n_pp = params$initial_n_pp,
        n_other = params$initial_n_other
    )
}

tuning_rate_context <- function(params, effort = params$initial_effort, t = 0) {
    functions <- projectRateFunctions(params)
    functions$Rates(params,
        n = params$initial_n, n_pp = params$initial_n_pp,
        n_other = params$initial_n_other, effort = effort, t = t,
        rates_fns = functions
    )
}

tuning_path_value <- function(params, path) {
    value <- params
    for (part in path) {
        if (!is.list(value) || !part %in% names(value)) {
            stop("Unknown tuning parameter path: ", paste(path, collapse = "/"),
                call. = FALSE
            )
        }
        value <- value[[part]]
    }
    value
}

tuning_path_replace <- function(params, path, value) {
    if (length(path) == 1L) {
        params[path] <- list(value)
    } else {
        params[[path[[1L]]]] <- tuning_path_replace(
            params[[path[[1L]]]], path[-1L], value
        )
    }
    params
}

tuning_context_args <- function(params, rates, t, effort) {
    c(tuning_state(params), list(rates = rates, t = t, effort = effort))
}

# Prepared callbacks never escape into stored model data.
prepare_tuning <- function(params, strategies = NULL, joint = NULL,
                           solver = "project", species = NULL,
                           effort = params$initial_effort, t = 0,
                           preserve = NULL, owners = NULL) {
    params <- clone_params_runtime(validParams(params))
    if (!solver %in% c(
        "project", "newton",
        "balance"
    )) {
        stop("Unsupported tuning solver.")
    }
    all_owners <- tuning_owners(params)
    requested <- if (is.null(owners)) all_owners else owners
    if (any(!requested %in% all_owners) ||
            anyDuplicated(requested)) {
        stop("Invalid tuning owners.")
    }
    if (!is.null(strategies) && (!is.character(strategies) ||
                                     anyNA(strategies) ||
                                     is.null(names(strategies)) ||
                                     anyDuplicated(names(strategies)) ||
                                     any(!names(strategies) %in% requested))) {
        stop("Strategies must name requested tuning owners.")
    }
    if (!is.null(preserve)) {
        if (!is.character(preserve) || length(preserve) != 1L ||
                is.na(preserve) || !nzchar(preserve)) {
            stop("Preserve must name one supported quantity.")
        }
        if (!"reproduction" %in% requested) {
            stop("Preserve requires a reproduction tuning responsibility.")
        }
        if ("reproduction" %in% names(strategies)) {
            stop("Select a reproduction strategy or preserve, not both.")
        }
        strategies["reproduction"] <- paste0("preserve_", preserve)
    }
    species_names <- as.character(params$species_params$species)
    if (is.null(species)) species <- species_names
    if (!is.character(species) || anyNA(species) || anyDuplicated(species) ||
            any(!species %in% species_names)) {
        stop("Unknown or repeated tuning species.")
    }
    partial <- length(species) != length(species_names)
    binding <- tuning_provider_binding(params)
    selected <- list()
    covered <- character()
    joint_registry <- params$capabilities$tuning
    defaults <- names(joint_registry)[vapply(joint_registry,
                                             function(x) isTRUE(x$default),
                                             logical(1))]
    related <- vapply(joint_registry[defaults], function(x) {
        length(intersect(x$owners, requested)) > 0L
    }, logical(1))
    verification_joint <- list()
    verification_covered <- character()
    for (id in defaults[!related]) {
        definition <- joint_registry[[id]]
        if (any(!definition$owners %in% all_owners) ||
                length(intersect(verification_covered, definition$owners))) {
            stop("Joint verification owners overlap or are unavailable.")
        }
        if (!identical(definition$binding, binding) ||
                !isTRUE(definition$stationary)) {
            stop("Joint equation-rate support is stale or not stationary.")
        }
        verification_joint[[id]] <- list(
            owners = definition$owners, joint = TRUE,
            equation_rates_fun = resolve_model_callback(
                definition$equation_rates_fun, "Joint verification rates"
            )
        )
        verification_covered <- c(verification_covered, definition$owners)
    }
    if (is.null(joint)) joint <- defaults[related]
    if (is.null(joint)) joint <- character()
    if (!is.character(joint) || anyNA(joint) || anyDuplicated(joint) ||
            any(!joint %in% names(joint_registry))) {
        stop("Unknown or repeated joint tuning strategy.")
    }
    for (id in joint) {
        definition <- joint_registry[[id]]
        if (any(!definition$owners %in% requested) ||
                length(intersect(covered, definition$owners)) ||
                length(intersect(names(strategies), definition$owners))) {
            stop("Joint owners overlap or conflict with local selections.")
        }
        if (!identical(definition$binding, binding) ||
                !isTRUE(definition$stationary)) {
            stop(
                "Joint tuning strategy '", id,
                "' lacks current stationary support."
            )
        }
        spec <- normalise_tuning_strategy(definition$strategy)
        selected[[paste0("joint:", id)]] <- list(
            id = id, owners = definition$owners,
            strategy = spec, equation_rates_fun = definition$equation_rates_fun,
            joint = TRUE
        )
        covered <- c(covered, definition$owners)
    }
    for (owner in setdiff(requested, covered)) {
        definition <- tuning_definition(params, owner)
        if (is.null(definition) || !tuning_binding_current(
            params, definition,
            owner
        ) ||
            !isTRUE(definition$stationary)) {
            stop("Missing current stationary tuning support for ", owner, ".",
                call. = FALSE
            )
        }
        block <- definition$tuning
        if (is.null(block) && identical(definition$role, "prescribed")) {
            block <- builtin_no_action_tuning(owner, prescribed = TRUE)
        }
        if (is.null(block)) {
            stop("Missing tuning strategy for ", owner, ".",
                call. = FALSE
            )
        }
        id <- if (owner %in% names(strategies)) strategies[[owner]] else NULL
        if (is.null(id)) id <- block$default
        if (!id %in% names(block$strategies)) {
            stop("Unsupported tuning strategy '", id, "' for ", owner, ".")
        }
        selected[[owner]] <- list(
            id = id, owners = owner,
            strategy = normalise_tuning_strategy(block$strategies[[id]]),
            equation_rates_fun = definition$equation_rates_fun, joint = FALSE
        )
    }
    # Even owners not being balanced need equation coverage: another balancer
    # can alter their inputs. No-action support is not inferred from silence.
    verification <- list()
    for (owner in setdiff(all_owners, verification_covered)) {
        joint_id <- names(selected)[vapply(selected, function(x) {
            x$joint &&
                owner %in% x$owners
        }, logical(1))]
        if (length(joint_id)) next
        definition <- tuning_definition(params, owner)
        if (is.null(definition) || !tuning_binding_current(
            params, definition,
            owner
        ) ||
            !isTRUE(definition$stationary)) {
            stop("Missing current stationary tuning support for ", owner, ".")
        }
        fun <- definition$equation_rates_fun
        if (is.null(fun) && identical(definition$role, "prescribed")) {
            fun <-
                "builtin_prescribed_equation_rates"
        }
        if (is.null(fun)) {
            stop(
                "Missing effective equation-rate support for ",
                owner, "."
            )
        }
        verification[[owner]] <- resolve_model_callback(
            fun,
            paste0("Equation rates for ", owner)
        )
    }
    held <- evolved <- paths <- character()
    hold_recruitment <- FALSE
    rates <- tuning_rate_context(params, effort, t)
    if (!is.numeric(rates$rdd) || is.complex(rates$rdd) ||
        length(rates$rdd) != length(species_names) ||
        any(!is.finite(rates$rdd)) || any(rates$rdd < 0) ||
        (!is.null(names(rates$rdd)) && !identical(
            names(rates$rdd),
            species_names
        ))) {
        stop("Effective recruitment requires one finite nonnegative value ",
             "per species in model order.")
    }
    args <- tuning_context_args(params, rates, t, effort)
    for (id in names(selected)) {
        node <- selected[[id]]
        spec <- node$strategy
        if (!solver %in% spec$solvers) {
            stop(
                "Tuning strategy '", id,
                "' does not support solver '", solver, "'."
            )
        }
        if (partial && !spec$partial_species) {
            stop(
                "Tuning strategy '", id,
                "' does not support partial species requests."
            )
        }
        states <- setdiff(node$owners, "reproduction")
        if (length(spec$hold$rates) && !"reproduction" %in% node$owners) {
            stop("Holding recruitment requires reproduction ownership.")
        }
        if (any(!c(spec$hold$states, spec$evolve) %in% states) ||
                anyDuplicated(c(spec$hold$states, spec$evolve)) ||
                !setequal(c(spec$hold$states, spec$evolve), states)) {
            stop(
                "Tuning strategy '", id,
                "' must partition owned states into held and evolving states."
            )
        }
        held <- c(held, spec$hold$states)
        evolved <- c(evolved, spec$evolve)
        hold_recruitment <- hold_recruitment || length(spec$hold$rates) > 0L
        for (field in c("balance_fun", "target_fun")) {
            if (!is.null(spec[[field]])) {
                node[[field]] <-
                    resolve_model_callback(spec[[field]], paste0(
                        "Tuning ", field,
                        " for ", id
                    ))
            }
        }
        if (node$joint) {
            node$equation_rates_fun <-
                resolve_model_callback(
                    node$equation_rates_fun,
                    paste0("Joint equation rates for ", id)
                )
        }
        node$preserved <- lapply(spec$preserve, function(descriptor) {
            if (!is.null(descriptor$path)) {
                tuning_path_value(params, descriptor$path)
            } else {
                invoke_tuning_callback(
                    descriptor$fun, params, args,
                    paste0("Preserved value for ", id)
                )
            }
        })
        for (descriptor in spec$adjust) {
            if (!tuning_parameter_path_allowed(descriptor$path)) {
                stop("Tuning cannot change state, providers or grid structure.")
            }
            path <- paste(descriptor$path, collapse = "/")
            overlap <- vapply(paths, function(old) {
                old == path ||
                    startsWith(old, paste0(path, "/")) || startsWith(
                    path,
                    paste0(old, "/")
                )
            }, logical(1))
            if (any(overlap)) {
                stop(
                    "Shared adjustable parameter '", path,
                    "' requires one joint strategy."
                )
            }
            paths <- c(paths, path)
            value <- tuning_path_value(params, descriptor$path)
            if (!is.numeric(value) || is.complex(value) ||
                    anyNA(value)) {
                stop("Adjustable parameters must be numeric without NAs.")
            }
            if (isTRUE(descriptor$species) &&
                    length(value) != length(species_names)) {
                stop("A species parameter must have one value per species.")
            }
            lower <- descriptor$lower
            upper <- descriptor$upper
            for (bound in list(lower, upper)) {
                if (!is.null(bound) &&
                    !length(bound) %in% c(
                        1L,
                        length(value)
                    )) {
                    stop("Parameter bounds have incompatible lengths.")
                }
            }
            if (!is.null(lower) && !is.null(upper) &&
                    any(lower > upper)) {
                stop("Parameter lower bounds exceed upper bounds.")
            }
        }
        selected[[id]] <- node
    }
    if (solver != "balance" &&
            "consumers" %in% held) {
        stop("Tuning solvers cannot hold consumers fixed.")
    }
    if (solver == "newton" && any(startsWith(evolved, "component:"))) {
        stop("Newton cannot evolve extension components.")
    }
    list(
        selected = selected, verification = verification,
        verification_joint = verification_joint, binding = binding,
        hold_states = held, evolve_states = evolved,
        hold_recruitment = hold_recruitment, recruitment = rates$rdd,
        rates = rates, species = species, all_species = species_names,
        effort = effort, t = t, solver = solver
    )
}

checked_equation_rates <- function(value, context) {
    if (!is.list(value) || is.object(value) || !capability_names_valid(value) ||
        !all(vapply(value, function(x) {
            is.numeric(x) && !is.complex(x) &&
                all(is.finite(x))
        }, logical(1)))) {
        stop(context,
            " must return a named list of finite numeric equation rates.",
            call. = FALSE
        )
    }
    value
}

tuning_equation_rates <- function(params, plan, rates, t) {
    args <- tuning_context_args(params, rates, t, plan$effort)
    result <- lapply(names(plan$verification), function(owner) {
        checked_equation_rates(invoke_tuning_callback(
            plan$verification[[owner]],
            params, c(args, list(
                owner = owner, process = owner,
                component = if (startsWith(owner, "component:")) {
                    substring(owner, 11L)
                } else {
                    NULL
                }
            )), paste0("Equation rates for ", owner)
        ), owner)
    })
    names(result) <- names(plan$verification)
    for (node in c(plan$selected, plan$verification_joint)) {
        if (node$joint) {
            values <- invoke_tuning_callback(
                node$equation_rates_fun, params,
                c(args, list(owners = node$owners)), "Joint equation rates"
            )
            if (!is.list(values) || !setequal(
                names(values),
                node$owners
            )) {
                stop("Joint equation rates must cover their declared owners.")
            }
            for (owner in node$owners) {
                result[[owner]] <-
                    checked_equation_rates(values[[owner]], owner)
            }
        }
    }
    result
}

capture_tuning_targets <- function(params, plan, rates = NULL, t = 0,
                                   rdd_target = NULL) {
    params <- clone_params_runtime(params)
    if (!identical(
        tuning_provider_binding(params),
        plan$binding
    )) {
        stop("Tuning providers changed before target capture.")
    }
    if (is.null(rates)) rates <- tuning_rate_context(params, plan$effort, t)
    if (is.null(rdd_target) && plan$hold_recruitment) {
        rdd_target <-
            plan$recruitment
    }
    if (!is.null(rdd_target)) {
        if (!is.numeric(rdd_target) || is.complex(rdd_target) ||
                length(rdd_target) != length(plan$all_species) ||
                any(!is.finite(rdd_target)) || any(rdd_target < 0)) {
            stop("Recruitment targets must be finite nonnegative values ",
                 "for every species.")
        }
        rates$rdd <- rdd_target
    }
    targets <- tuning_equation_rates(params, plan, rates, t)
    args <- tuning_context_args(params, rates, t, plan$effort)
    for (node in plan$selected) {
        if (!is.null(node$target_fun)) {
            value <- invoke_tuning_callback(
                node$target_fun, params,
                c(args, list(
                    owners = node$owners, preserved = node$preserved,
                    species = plan$species, strategy = node$id,
                    settings = node$strategy$settings
                )),
                "Tuning target construction"
            )
            if (!node$joint) value <- setNames(list(value), node$owners)
            if (!is.list(value) || !setequal(
                names(value),
                node$owners
            )) {
                stop("A target callback must cover its declared owners.")
            }
            for (owner in node$owners) {
                targets[[owner]] <-
                    checked_equation_rates(value[[owner]], owner)
            }
        }
    }
    # The standard consumer equations always receive an independent invariant,
    # even when an extension supplies its own broader equation-rate adapter.
    list(
        reference = clone_params_runtime(params),
        state = unserialize(serialize(tuning_state(params), NULL)),
        rates = rates,
        equations = targets, core = builtin_consumer_equation_rates(
            params,
            rates
        ),
        rdd_target = rates$rdd, t = t
    )
}

apply_tuning_balancers <- function(params, plan, targets) {
    params <- clone_params_runtime(params)
    for (node in plan$selected) {
        if (node$strategy$action == "none") next
        if (!node$joint && identical(node$owners, "reproduction")) {
            context <- list(
                plan = plan, targets = targets,
                rates = targets$rates
            )
            params <- balanceRDD_dispatch(params,
                rdd_target = targets$rdd_target,
                effort = plan$effort, t = targets$t, .balance_context = context
            )
            if (!inherits(params, "MizerParams")) {
                stop("A balancer must return a MizerParams candidate.")
            }
            next
        }
        args <- tuning_context_args(
            params, targets$rates, targets$t,
            plan$effort
        )
        settings <- node$strategy$settings
        reserved <- intersect(names(settings), c(
            names(args), "targets",
            "rdd_target", "strategy", "preserved", "species", "settings"
        ))
        if (length(reserved)) {
            stop("Tuning settings conflict with reserved callback context.")
        }
        args <- c(args, list(
            targets = targets$equations[node$owners],
            rdd_target = targets$rdd_target, strategy = node$id,
            preserved = node$preserved, species = plan$species,
            settings = settings
        ), settings)
        params <- invoke_tuning_callback(
            node$balance_fun, params, args,
            paste0("Balancer for ", paste(node$owners, collapse = ", "))
        )
        if (!inherits(params, "MizerParams")) {
            stop("A balancer must return a MizerParams candidate.")
        }
    }
    params
}

tuning_values_equal <- function(a, b, tolerance = 1e-8) {
    if (!identical(dim(a), dim(b)) || !identical(dimnames(a), dimnames(b)) ||
            !identical(names(a), names(b)) || length(a) != length(b)) {
        return(FALSE)
    }
    if (!is.numeric(a) || !is.numeric(b)) {
        return(identical(a, b))
    }
    if (anyNA(a) || anyNA(b)) {
        return(FALSE)
    }
    equal <- a == b
    remaining <- !equal
    if (any(!is.finite(a[remaining])) || any(!is.finite(b[remaining]))) {
        return(FALSE)
    }
    limit <- 1e-12 + tolerance * abs(b[remaining])
    all(abs(a[remaining] - b[remaining]) <= limit)
}

refresh_tuning_dependencies <- function(params, plan, reference) {
    species_changed <- character()
    resource_changed <- character()
    rebuild <- FALSE
    for (node in plan$selected) for (descriptor in node$strategy$adjust) {
        path <- descriptor$path
        if (tuning_values_equal(tuning_path_value(params, path),
                                tuning_path_value(reference, path))) next
        if (length(path) == 2L && path[[1L]] == "species_params") {
            species_changed <- c(species_changed, path[[2L]])
        } else if (path[[1L]] == "resource_params") {
            resource_changed <- c(resource_changed, path[[2L]])
        } else if (!path[[1L]] %in% c("cc_pp", "rr_pp")) {
            # Extension setters may cache quantities from these parameters.
            # Rebuild both candidates independently, never trust stale arrays.
            rebuild <- TRUE
        }
    }
    required <- recalculation_species_params(reference, params$species_params)
    if (rebuild || length(intersect(species_changed, required))) {
        params <- setParams(params, info_level = 0)
    }
    if (length(resource_changed)) {
        params <- with_info_level(
            setResource(params, balance = FALSE,
                resource_params_changed = unique(resource_changed)
            ), info_level = 0
        )
    }
    validParams(params, info_level = 0)
}

normalise_tuning_metadata <- function(x) {
    attributes <- attributes(x)
    # Comment presence freezes rate arrays. Only the annotation text is
    # immaterial; adding or removing a freeze requires an authorised path.
    if (!is.null(attributes$comment)) attributes$comment <- "present"
    if (is.list(x)) {
        # numeric_version and other list classes can return their own class
        # from [[; recurse on raw storage, then restore every other attribute.
        x <- lapply(unclass(x), normalise_tuning_metadata)
    }
    attributes(x) <- attributes
    x
}

verify_tuning_candidate <- function(
    candidate, plan, targets,
    tolerance = 1e-8
) {
    reference <- targets$reference
    if (!identical(
        tuning_provider_binding(candidate),
        plan$binding
    )) {
        stop("A balancer changed effective providers.")
    }
    for (state in names(targets$state)) {
        if (!isTRUE(all.equal(tuning_state(candidate)[[state]],
                targets$state[[state]],
                check.attributes = TRUE,
                tolerance = 0
            ))) {
            stop(
                "A balancer changed the target state '",
                state, "'."
            )
        }
    }
    # Rebuild known dependent arrays before verification, preserving explicit
    # frozen-array semantics. Resource capacity is an authorised direct target.
    candidate <- refresh_tuning_dependencies(candidate, plan, reference)
    expected <- reference
    for (node in plan$selected) {
        for (descriptor in node$strategy$adjust) {
            value <- tuning_path_value(candidate, descriptor$path)
            original <- tuning_path_value(reference, descriptor$path)
            if (!is.numeric(value) || is.complex(value) || anyNA(value) ||
                    !identical(dim(value), dim(original)) ||
                    length(value) != length(original)) {
                stop("An adjustable parameter changed shape or type.")
            }
            if ((!is.null(descriptor$lower) && any(value < descriptor$lower)) ||
                    (!is.null(descriptor$upper) &&
                         any(value > descriptor$upper))) {
                stop("A balanced parameter violates its declared bounds.")
            }
            if (isTRUE(descriptor$species)) {
                untouched <- !plan$all_species %in% plan$species
                if (!tuning_values_equal(
                    value[untouched],
                    original[untouched]
                )) {
                    stop("A balancer changed an unselected species parameter.")
                }
            }
            expected <- tuning_path_replace(expected, descriptor$path, value)
            # Table setters also record the authorised parameter as given.
            if (length(descriptor$path) == 2L &&
                    descriptor$path[[1L]] == "species_params") {
                column <- descriptor$path[[2L]]
                if (isTRUE(descriptor$species)) {
                    before <- reference$given_species_params[[column]]
                    after <- candidate$given_species_params[[column]]
                    if (is.null(before)) before <- rep(NA_real_, length(value))
                    if (is.null(after)) after <- rep(NA_real_, length(value))
                    if (!isTRUE(all.equal(before[untouched], after[untouched],
                                          tolerance = 0,
                                          check.attributes = FALSE))) {
                        stop("A balancer changed unselected given parameters.")
                    }
                }
                expected$given_species_params[[column]] <-
                    candidate$given_species_params[[column]]
            }
        }
    }
    # Derive the authorised model independently from the original and the
    # allowed parameters. This rejects hand-edited frozen rate arrays rather
    # than merely whitelisting every array that a setter might rebuild.
    expected <- refresh_tuning_dependencies(expected, plan, reference)
    expected$time_modified <- candidate$time_modified
    if (!identical(tuning_provider_binding(candidate), plan$binding)) {
        stop("Rebuilding the balanced model changed effective providers.")
    }
    for (state in names(targets$state)) {
        if (!isTRUE(all.equal(tuning_state(candidate)[[state]],
                targets$state[[state]],
                check.attributes = TRUE, tolerance = 0
            ))) {
            stop(
                "Rebuilding the balanced model changed the target state '",
                state, "'."
            )
        }
    }
    # Ignore annotation text while preserving freeze status and all structure.
    if (!isTRUE(all.equal(normalise_tuning_metadata(candidate),
                          normalise_tuning_metadata(expected),
                          tolerance = 0, check.attributes = TRUE))) {
        stop("A balancer changed parameters outside its declared ",
             "adjustable paths.", call. = FALSE)
    }
    rates <- tuning_rate_context(candidate, plan$effort, targets$t)
    args <- tuning_context_args(candidate, rates, targets$t, plan$effort)
    for (node in plan$selected) {
        for (id in names(node$strategy$preserve)) {
            descriptor <- node$strategy$preserve[[id]]
            value <- if (!is.null(descriptor$path)) {
                tuning_path_value(candidate, descriptor$path)
            } else {
                invoke_tuning_callback(
                    descriptor$fun, candidate, args,
                    "Preserved tuning value"
                )
            }
            if (!tuning_values_equal(
                value, node$preserved[[id]],
                tolerance
            )) {
                stop("Balancer did not preserve '", id, "'.")
            }
        }
    }
    actual <- tuning_equation_rates(candidate, plan, rates, targets$t)
    actual$core <- builtin_consumer_equation_rates(candidate, rates)
    expected_rates <- c(targets$equations, list(core = targets$core))
    if (!setequal(
        names(actual),
        names(expected_rates)
    )) {
        stop("Equation-rate ownership changed during balancing.")
    }
    for (owner in names(expected_rates)) {
        wanted <- expected_rates[[owner]]
        found <- actual[[owner]]
        if (!setequal(
            names(wanted),
            names(found)
        )) {
            stop("Equation-rate fields changed for ", owner, ".")
        }
        for (rate in names(wanted)) {
            if (!tuning_values_equal(
                found[[rate]], wanted[[rate]],
                tolerance
            )) {
                stop("Balancing changed effective equation rate '", owner,
                    "/", rate, "'.",
                    call. = FALSE
                )
            }
        }
    }
    candidate
}

finish_tuning <- function(params, plan, targets) {
    verify_tuning_candidate(
        apply_tuning_balancers(params, plan, targets),
        plan, targets
    )
}

joint_tuning_errors <- function(registry) {
    if (!is.list(registry) || is.object(registry) ||
            !capability_names_valid(registry)) {
        return("Joint tuning registry must be a named list.")
    }
    for (entry in registry) {
        if (!is.list(entry) || !capability_names_valid(entry) ||
            !setequal(names(entry), c(
                "owners", "stationary",
                "equation_rates_fun", "strategy", "binding", "default"
            )) ||
            !is.character(entry$owners) || !length(entry$owners) ||
            anyNA(entry$owners) || any(!nzchar(entry$owners)) ||
            anyDuplicated(entry$owners) ||
            !is.logical(entry$stationary) || length(entry$stationary) != 1L ||
            is.na(entry$stationary) ||
            !is.character(entry$equation_rates_fun) ||
            length(entry$equation_rates_fun) != 1L ||
            is.na(entry$equation_rates_fun) ||
            !nzchar(entry$equation_rates_fun) ||
            !is.logical(entry$default) || length(entry$default) != 1L ||
            is.na(entry$default) || !tuning_data_valid(entry$binding)) {
            return("Invalid joint tuning definition.")
        }
        errors <- tuning_strategy_errors(entry$strategy)
        if (!isTRUE(errors)) {
            return(errors)
        }
    }
    TRUE
}

#' Register a joint tuning strategy
#'
#' A joint strategy owns all listed responsibilities. Selecting it suppresses
#' their local defaults; overlapping owners or adjustable parameter paths are
#' errors. Its equation-rate callback returns one named rate list per owner.
#' The declaration is bound to the current effective dispatch, including named
#' rate replacements and extension methods.
#'
#' @param params A MizerParams object.
#' @param name Unique name for the joint strategy.
#' @param owners Canonical tuning owners: consumers, resource, reproduction,
#'   or component:name.
#' @param definition A tuning strategy specification.
#' @param equation_rates_fun Callback returning effective equation quantities.
#' @param stationary Explicit declaration that this strategy supports stationary
#'   tuning of its processes. A prescribed state alone does not establish this.
#' @param default Whether tuning selects this joint strategy by default.
#' @return The updated model.
#' @export
setTuningStrategy <- function(params, name, owners, definition,
                              equation_rates_fun, stationary, default = FALSE) {
    params <- validParams(params)
    if (!is.character(name) || length(name) != 1L || is.na(name) ||
            !nzchar(name)) {
        stop("A joint strategy needs one nonempty name.")
    }
    if (any(!owners %in% tuning_owners(params))) {
        stop("Unknown joint tuning owner.")
    }
    entry <- list(
        owners = owners, stationary = stationary,
        equation_rates_fun = equation_rates_fun,
        strategy = normalise_tuning_strategy(definition),
        binding = tuning_provider_binding(params), default = default
    )
    errors <- joint_tuning_errors(setNames(list(entry), name))
    if (!isTRUE(errors)) stop(errors, call. = FALSE)
    resolve_model_callback(
        equation_rates_fun,
        paste0("Joint equation rates for ", name)
    )
    params$capabilities$tuning[[name]] <- entry
    params$time_modified <- lubridate::now()
    params
}

capability_callback_names <- function(params) {
    strategy_callbacks <- function(strategy) {
        c(
            strategy$balance_fun, strategy$target_fun,
            unlist(lapply(strategy$preserve, `[[`, "fun"), use.names = FALSE)
        )
    }
    provider_callbacks <- function(definition) {
        c(
            unlist(definition[intersect(names(definition), c(
                "dynamics_fun", "rate_fun",
                "encounter_fun", "mort_fun", "convergence_fun",
                "equation_rates_fun", "level_fun", "biomass_fun",
                "mortality_fun", "prescribed_mask_fun"
            ))], use.names = FALSE),
            unlist(lapply(definition$transformations, `[[`, "transform_fun"),
                   use.names = FALSE),
            unlist(lapply(definition$tuning$strategies, strategy_callbacks),
                use.names = FALSE
            )
        )
    }
    unique(c(
        unlist(lapply(c(
            params$capabilities$components,
            params$capabilities$processes
        ), provider_callbacks), use.names = FALSE),
        unlist(lapply(params$capabilities$transformations, function(entry) {
            entry$definition$transform_fun
        }), use.names = FALSE),
        unlist(lapply(params$capabilities$tuning, function(entry) {
            c(entry$equation_rates_fun, strategy_callbacks(entry$strategy))
        }), use.names = FALSE)
    ))
}

builtin_reproduction_law <- function(params) {
    if (!builtin_tuning_rates(params, reproduction = TRUE)) {
        return(NULL)
    }
    law <- sub("^mizer::", "", params$rates_funcs$RDD)
    if (law %in% c(
        "BevertonHoltRDD", "RickerRDD", "SheperdRDD",
        "constantRDD", "noRDD"
    )) {
        law
    } else {
        NULL
    }
}
