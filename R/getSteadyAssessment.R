#' Assess the steady state of the whole model or a conditional subsystem
#'
#' `r lifecycle::badge("experimental")`
#' Returns a structured assessment of every dynamical owner in the requested
#' scope. The default includes the consumers, resource and all components.
#' Missing convergence support is reported explicitly; it never counts as a
#' successful check. A conditional assessment requires an explicit `fixed`
#' list and describes the subsystem with those owners held at their current
#' values. Owners declared prescribed are recorded separately. Structural
#' prescribed-coordinate masks exclude only those cells from the assessment;
#' the remaining cells of a partially prescribed owner still require coverage.
#'
#' The default `criterion = "state"` measures the integral of the absolute
#' rate of change of each consumer spectrum and the resource spectrum, weighted
#' by individual mass. Changes in different size classes therefore cannot
#' cancel. `criterion = "biomass"` instead takes the absolute value after
#' integrating and explicitly permits such cancellation. Both built-in checks
#' compare against `atol + tol * integral(w * abs(n))`. They use the model's
#' quadrature scheme and report values and tolerances in grams per year.
#'
#' Extensions supply named quantitative measures with their own units and
#' tolerances. Mizer compares `abs(value) <= tolerance` for every row and
#' preserves signed values in the report. `tol` and `atol` apply only to the
#' built-in checks. `tolerances` overrides measures from any provider.
#'
#' @param params A MizerParams object.
#' @param tol Relative tolerance for built-in checks, in 1/year.
#' @param effort Fishing effort at which to assess the state.
#' @param scope `"whole"` or `"conditional"`.
#' @param fixed Owners held fixed for a conditional assessment. Use
#'   `"consumers"`, `"resource"`, or `"component:name"` for a component named
#'   `name`. This must be empty for whole-model assessments and nonempty for
#'   conditional assessments.
#' @param criterion `"state"` for the absolute state residual, or `"biomass"`
#'   for net biomass drift in the built-in checks.
#' @param tolerances Optional named numeric vector overriding individual
#'   tolerances. Names are `owner/measure`, as shown in the `measures` table.
#' @param atol Absolute tolerance for built-in checks, in grams per year.
#' @param t Time at which to evaluate the model, in years.
#' @param dt Step length made available to extension convergence callbacks.
#' @param ... Additional arguments made available to convergence callbacks.
#' @return A `MizerSteadyAssessment` list. `verdict` is `TRUE` when every
#'   required owner is covered and passes, `FALSE` when any measured check
#'   fails, and `NA` when no check fails but required support is missing.
#'   The list records `scope`, `fixed`, `prescribed`, `required`, `covered`,
#'   `missing`, `time` and `criterion`. `prescribed_masks` records the exact
#'   prescribed coordinates for each owner. `unavailable` gives each missing
#'   owner, capability and reason. `measures` has columns `owner`, `measure`,
#'   `value`, `tolerance`, `units` and `passed`. A scope with no required
#'   dynamic owners is steady by construction and has no measures.
#' @seealso [isSteady()], [getSteadyResidual()]
#' @export
getSteadyAssessment <- function(params, tol = 0.05,
                                effort = params@initial_effort,
                                scope = c("whole", "conditional"),
                                fixed = character(),
                                criterion = c("state", "biomass"),
                                tolerances = NULL, atol = 1e-12,
                                t = 0, dt = 1e-4, ...) {
    UseMethod("getSteadyAssessment")
}

#' @rdname getSteadyAssessment
#' @usage NULL
#' @export
getSteadyAssessment.MizerParams <- function(params, tol = 0.05,
                                            effort = params@initial_effort,
                                            scope = c("whole", "conditional"),
                                            fixed = character(),
                                            criterion = c("state", "biomass"),
                                            tolerances = NULL, atol = 1e-12,
                                            t = 0, dt = 1e-4, ...) {
    prepared <- prepare_steady_assessment(
        params, tol = tol, scope = scope, fixed = fixed,
        criterion = criterion, tolerances = tolerances, atol = atol, dt = dt
    )
    steady_assessment(params, prepared = prepared, effort = effort, t = t, ...)
}

assessment_state_structure <- function(state) {
    list(type = typeof(state), length = length(state),
         dimensions = dim(state), dimnames = dimnames(state),
         names = names(state), class = class(state),
         children = if (is.list(state)) {
             lapply(state, assessment_state_structure)
         } else {
             NULL
         })
}

assessment_model_key <- function(params) {
    # Structural masks may read arbitrary model parameters, but never current
    # state values. Retain state structure while permitting runtime evolution.
    model <- params
    for (field in c("initial_n", "initial_n_pp", "initial_n_other")) {
        model[[field]] <- assessment_state_structure(params[[field]])
    }
    model$initial_effort <- NULL
    digest::digest(list(model = model,
                        dispatch = process_dispatch_binding(params)),
                   algo = "xxhash64")
}

# Resolve and validate once per operation, without evaluating any equations or
# assessment callbacks. The prepared functions belong to runtime only.
prepare_steady_assessment <- function(params, tol = 0.05,
                                      scope = c("whole", "conditional"),
                                      fixed = character(),
                                      criterion = c("state", "biomass"),
                                      tolerances = NULL, atol = 1e-12,
                                      dt = 1e-4,
                                      rate_constraints = character()) {
    params <- validParams(params)
    scope <- match.arg(scope)
    criterion <- match.arg(criterion)
    for (arg in c("tol", "atol", "dt")) {
        value <- get(arg)
        if (!is.numeric(value) || is.complex(value) || length(value) != 1L ||
                !is.finite(value) || value < 0 || (arg == "dt" && value == 0)) {
            stop("`", arg, "` must be a finite ",
                 if (arg == "dt") "positive" else "nonnegative",
                 " number.", call. = FALSE)
        }
    }
    owners <- assessment_owners(params)
    if (!is.character(fixed) || anyNA(fixed) || anyDuplicated(fixed) ||
            any(!fixed %in% owners)) {
        stop("`fixed` must contain unique canonical owner IDs: ",
             paste(owners, collapse = ", "), ".", call. = FALSE)
    }
    if (!is.character(rate_constraints) || anyNA(rate_constraints) ||
            anyDuplicated(rate_constraints) ||
            any(!rate_constraints %in% "consumers/recruitment")) {
        stop("Unsupported internal rate constraints.", call. = FALSE)
    }
    if ((scope == "whole" && (length(fixed) || length(rate_constraints))) ||
            (scope == "conditional" && !length(fixed) &&
                 !length(rate_constraints))) {
        stop("Use `scope = \"conditional\"` with an explicit nonempty ",
             "`fixed` list; whole-model assessments require an empty list.",
             call. = FALSE)
    }
    if (!is.null(tolerances) &&
            (!is.numeric(tolerances) || is.complex(tolerances) ||
                 any(!is.finite(tolerances)) || any(tolerances < 0) ||
                 is.null(names(tolerances)) || anyNA(names(tolerances)) ||
                 any(!nzchar(names(tolerances))) ||
                 anyDuplicated(names(tolerances)))) {
        stop("`tolerances` must be a named finite nonnegative numeric vector ",
             "with unique owner/measure names.", call. = FALSE)
    }
    specifications <- setNames(lapply(owners, function(owner) {
        assessment_capability(params, owner)
    }), owners)
    state_capabilities <- prepare_state_capabilities(params, fixed = fixed)
    prescribed <- owners[vapply(state_capabilities, function(specification) {
        !specification$fixed && isTRUE(specification$mask_current) &&
            !any(specification$dynamic_mask)
    }, logical(1))]
    required <- setdiff(owners, c(fixed, prescribed))
    callbacks <- list()
    derivatives <- list()
    unavailable <- data.frame(owner = character(), capability = character(),
                              reason = character(), stringsAsFactors = FALSE)
    for (owner in required) {
        specification <- specifications[[owner]]
        state <- state_capabilities[[owner]]
        if (!isTRUE(state$mask_current)) {
            unavailable <- rbind(unavailable, data.frame(
                owner = owner, capability = "prescribed_inputs",
                reason = state$mask_reason, stringsAsFactors = FALSE
            ))
            next
        }
        if (!isTRUE(specification$current) ||
                is.null(specification$convergence_fun)) {
            unavailable <- rbind(unavailable, data.frame(
                owner = owner, capability = "convergence",
                reason = specification$reason, stringsAsFactors = FALSE
            ))
        } else {
            callbacks[[owner]] <- resolve_model_callback(
                specification$convergence_fun,
                paste0("Convergence assessment for '", owner, "'")
            )
            if (startsWith(owner, "component:")) {
                component <- substring(owner, nchar("component:") + 1L)
                name <- active_component_rate_name(params, component)
                if (!is.null(name)) {
                    derivatives[[owner]] <- resolve_model_callback(
                        name, paste0("Equation for '", owner, "'")
                    )
                }
            }
        }
    }
    list(scope = scope, fixed = fixed, rate_constraints = rate_constraints,
         prescribed = prescribed,
         prescribed_masks = lapply(state_capabilities, `[[`, "prescribed_mask"),
         state_capabilities = state_capabilities,
         required = required,
         covered = if (length(callbacks)) names(callbacks) else character(),
         missing = unavailable$owner, unavailable = unavailable,
         callbacks = callbacks, derivatives = derivatives,
         specifications = specifications,
         tol = tol, criterion = criterion, tolerances = tolerances,
         atol = atol, dt = dt,
         builtin_projection = !is.null(
             default_process_registry(params)$consumers$convergence_fun
         ),
         model_key = assessment_model_key(params))
}

empty_assessment_measures <- function() {
    data.frame(owner = character(), measure = character(), value = numeric(),
               tolerance = numeric(), units = character(), passed = logical(),
               stringsAsFactors = FALSE)
}

# All arguments other than the positional model are optional context: strict
# callbacks may declare just those they use, without accepting `...`.
invoke_convergence_callback <- function(fun, args, owner) {
    accepted <- names(formals(fun))
    if (!"..." %in% accepted) {
        args <- args[c(TRUE, names(args)[-1L] %in% accepted)]
    }
    tryCatch(do.call(fun, args), error = function(error) {
        stop("Convergence assessment for '", owner, "' failed: ",
             conditionMessage(error), call. = FALSE)
    })
}

checked_assessment_measures <- function(value, owner) {
    fields <- c("measure", "value", "tolerance", "units")
    valid <- is.data.frame(value) && nrow(value) > 0L &&
        all(fields %in% names(value)) && !anyDuplicated(names(value))
    if (valid) {
        valid <- all(vapply(value[fields], function(column) {
            is.null(dim(column)) && length(column) == nrow(value)
        }, logical(1)))
    }
    if (valid) {
        valid <- is.character(value$measure) && !anyNA(value$measure) &&
            all(nzchar(value$measure)) && !anyDuplicated(value$measure) &&
            is.numeric(value$value) && !is.complex(value$value) &&
            all(is.finite(value$value)) &&
            is.numeric(value$tolerance) && !is.complex(value$tolerance) &&
            all(is.finite(value$tolerance)) && all(value$tolerance >= 0) &&
            is.character(value$units) && !anyNA(value$units) &&
            all(nzchar(value$units))
    }
    if (!valid) {
        stop("Convergence assessment for '", owner,
             "' must return a nonempty data frame with unique character ",
             "`measure`, finite numeric `value`, finite nonnegative ",
             "`tolerance`, and nonempty character `units` columns.",
             call. = FALSE)
    }
    data.frame(owner = owner, value[fields],
               passed = abs(value$value) <= value$tolerance,
               stringsAsFactors = FALSE, row.names = NULL)
}

# Explicit state arguments avoid copying the model at every settling check.
steady_assessment <- function(params, prepared = NULL,
                              n = params@initial_n,
                              n_pp = params@initial_n_pp,
                              n_other = params@initial_n_other,
                              effort = params@initial_effort, t = 0,
                              rates_fns = projectRateFunctions(params), ...) {
    params <- validParams(params)
    if (is.null(prepared)) prepared <- prepare_steady_assessment(params)
    if (!identical(prepared$model_key, assessment_model_key(params))) {
        stop("The model's parameters or providers changed after assessment ",
             "preparation; ",
             "prepare a new assessment.", call. = FALSE)
    }
    if (!is.numeric(t) || is.complex(t) || length(t) != 1L || !is.finite(t)) {
        stop("`t` must be a finite number.", call. = FALSE)
    }
    covered <- prepared$covered
    resource <- prepared$state_capabilities$resource
    evaluate_resource <- "resource" %in% covered &&
        isTRUE(resource$equation_current) && !is.null(resource$rate_fun)
    evaluated <- if (length(covered)) {
        steady_rates(
            params, effort = effort, dt = prepared$dt, n = n, n_pp = n_pp,
            n_other = n_other, t = t, evaluate_other = FALSE,
            evaluate_consumers = prepared$builtin_projection &&
                "consumers" %in% covered,
            evaluate_resource = FALSE, rates_fns = rates_fns
        )
    } else {
        list(params = params, n = n, n_pp = n_pp, n_other = n_other,
             dNdt = NULL, dn_pp_dt = NULL, other = list(),
             rates = NULL, time = t)
    }
    params <- evaluated$params
    if (evaluate_resource) {
        evaluated$dn_pp_dt <- evaluate_resource_equation(
            params, resource, n = n, n_pp = n_pp, n_other = n_other,
            rates = evaluated$rates, t = t
        )
    }
    measures <- empty_assessment_measures()
    for (owner in covered) {
        derivative <- NULL
        args <- list(params, n = n, n_pp = n_pp, n_other = n_other,
                     rates = evaluated$rates, t = t, dt = prepared$dt,
                     tol = prepared$tol, atol = prepared$atol,
                     criterion = prepared$criterion,
                     prescribed_mask = prepared$prescribed_masks[[owner]],
                     dynamic_mask = prepared$state_capabilities[[owner]]$
                         dynamic_mask)
        if (startsWith(owner, "component:")) {
            component <- substring(owner, nchar("component:") + 1L)
            args$component <- component
            rate_fun <- prepared$derivatives[[owner]]
            if (!is.null(rate_fun)) {
                derivative <- tryCatch(checked_component_rate(
                    rate_fun, params, n = n, n_pp = n_pp, n_other = n_other,
                    rates = evaluated$rates, t = t, component = component
                ), error = function(error) {
                    stop("Equation for '", owner, "' failed: ",
                         conditionMessage(error), call. = FALSE)
                })
            }
        } else {
            args$process <- owner
            derivative <- if (owner == "consumers") {
                evaluated$dNdt
            } else {
                evaluated$dn_pp_dt
            }
        }
        args["derivative"] <- list(derivative)
        extra <- list(...)
        if (length(extra)) {
            reserved <- intersect(names(extra), names(args))
            if (length(reserved)) {
                stop("Assessment context arguments cannot be overridden in ",
                     "`...`: ", paste(reserved, collapse = ", "), ".",
                     call. = FALSE)
            }
            args <- c(args, extra)
        }
        value <- invoke_convergence_callback(prepared$callbacks[[owner]],
                                             args, owner)
        measures <- rbind(measures, checked_assessment_measures(value, owner))
    }
    overrides <- prepared$tolerances
    ids <- paste(measures$owner, measures$measure, sep = "/")
    if (anyDuplicated(ids)) {
        stop("Assessment owner/measure IDs must be unique.", call. = FALSE)
    }
    if (length(overrides)) {
        unknown <- setdiff(names(overrides), ids)
        if (length(unknown)) {
            stop("Unknown assessment tolerance IDs: ",
                 paste(unknown, collapse = ", "), ".", call. = FALSE)
        }
        index <- match(names(overrides), ids)
        measures$tolerance[index] <- overrides
        measures$passed <- abs(measures$value) <= measures$tolerance
    }
    verdict <- if (any(!measures$passed)) {
        FALSE
    } else if (length(prepared$missing)) {
        NA
    } else {
        TRUE
    }
    report <- structure(list(
        verdict = verdict, scope = prepared$scope, fixed = prepared$fixed,
        prescribed = prepared$prescribed,
        prescribed_masks = prepared$prescribed_masks,
        required = prepared$required,
        covered = covered, missing = prepared$missing,
        unavailable = prepared$unavailable, measures = measures,
        time = t, criterion = prepared$criterion
    ), class = "MizerSteadyAssessment")
    if (length(prepared$rate_constraints)) {
        report$rate_constraints <- prepared$rate_constraints
    }
    attr(report, "drift_report") <- steady_drift_report(params,
                                                        evaluated = evaluated)
    report
}

builtin_consumer_assessment <- function(params, n, derivative, tol, atol,
                                        criterion, dynamic_mask = NULL, ...) {
    if (is.null(derivative)) {
        stop("No built-in consumer derivative is available for this provider.",
             call. = FALSE)
    }
    if (!is.null(dynamic_mask)) {
        n[!dynamic_mask] <- 0
        derivative[!dynamic_mask] <- 0
    }
    abundance <- as.numeric(sizeIntegral(params, weighting = params@w,
                                         n = abs(n)))
    value <- as.numeric(sizeIntegral(params, weighting = params@w,
                                     n = if (criterion == "state") {
                                         abs(derivative)
                                     } else {
                                         derivative
                                     }))
    data.frame(measure = as.character(params@species_params$species),
               value = abs(value), tolerance = atol + tol * abundance,
               units = "g/year", stringsAsFactors = FALSE)
}

builtin_resource_assessment <- function(params, n_pp, derivative, tol, atol,
                                        criterion, dynamic_mask = NULL, ...) {
    if (is.null(derivative)) {
        stop("No built-in resource derivative is available for this provider.",
             call. = FALSE)
    }
    if (!is.null(dynamic_mask)) {
        n_pp[!dynamic_mask] <- 0
        derivative[!dynamic_mask] <- 0
    }
    weight <- bin_average_weight(params@w_full, params) * params@dw_full
    value <- if (criterion == "state") {
        sum(abs(derivative) * weight)
    } else {
        abs(sum(derivative * weight))
    }
    data.frame(measure = "resource", value = value,
               tolerance = atol + tol * sum(abs(n_pp) * weight),
               units = "g/year", stringsAsFactors = FALSE)
}

assessment_ratio <- function(assessment) {
    measures <- assessment$measures
    if (!nrow(measures)) return(if (length(assessment$missing)) NA_real_ else 0)
    ratio <- abs(measures$value) / measures$tolerance
    ratio[measures$value == 0 & measures$tolerance == 0] <- 0
    max(ratio)
}
