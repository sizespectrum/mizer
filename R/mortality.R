# Mortality decomposition is diagnostic support, independent of the rate used
# by the dynamical equations. Bind it to the two effective outer assemblers.
mortality_provider_binding <- function(params) {
    providers <- lapply(params$rates_funcs[c("Rates", "Mort")], as.character)
    generics <- character()
    for (rate in names(providers)) {
        name <- providers[[rate]]
        if (identical(name, paste0("mizer", rate)) ||
                sub("^mizer::", "", name) == paste0("project", rate)) {
            generics <- c(generics, paste0("project", rate))
        }
    }
    list(
        rate_functions = providers,
        rate_methods = tuning_method_names(params, generics),
        contributions = lapply(params$other_mort, as.character)
    )
}

mortality_binding_current <- function(params, definition) {
    binding <- definition$mortality_binding
    if (is.null(binding)) {
        return(FALSE)
    }
    binding$automatic <- NULL
    identical(binding, mortality_provider_binding(params))
}

builtin_mortality_provider <- function(params) {
    binding <- mortality_provider_binding(params)
    known <- all(vapply(names(binding$rate_functions), function(rate) {
        sub("^mizer::", "", binding$rate_functions[[rate]]) %in%
            c(paste0("mizer", rate), paste0("project", rate))
    }, logical(1)))
    known && all(binding$rate_methods %in%
        c(
            "projectRates.MizerParams@mizer",
            "projectMort.MizerParams@mizer"
        ))
}

default_provider_mortality <- function(params) {
    list(
        mortality_fun = if (builtin_mortality_provider(params)) {
            "builtin_mortality_breakdown"
        } else {
            NULL
        },
        mortality_binding = c(
            mortality_provider_binding(params),
            list(automatic = TRUE)
        )
    )
}

refresh_builtin_mortality <- function(params) {
    definition <- params$capabilities$processes$consumers
    if (is.null(definition) || mortality_binding_current(params, definition) ||
            !isTRUE(definition$mortality_binding$automatic)) {
        return(params)
    }
    defaults <- default_provider_mortality(params)
    if (!is.null(defaults$mortality_fun)) {
        definition[names(defaults)] <- defaults
        params$capabilities$processes$consumers <- definition
    }
    params
}

#' Diagnose mortality by cause
#'
#' Calculate mortality contributions separately from the authoritative total
#' mortality used in the dynamical equations. The getter evaluates the complete
#' effective rate pipeline once at each requested state, then checks a declared
#' complete decomposition against its total mortality.
#'
#' Built-in mortality is decomposed into predation, fishing, external mortality
#' and each additive contribution registered with [setComponent()] or
#' [other_mort()]. Replacing the outer `Rates` or `Mort` provider requires an
#' explicitly adopted decomposition through
#' `setProcess(params, "consumers", mortality_fun = "myBreakdown")`.
#' Unsupported decomposition returns an incomplete report with a reason and the
#' authoritative total; it does not assign unexplained mortality to a cause.
#'
#' A `mortality_fun` receives the model, `n`, `n_pp`, `n_other`, `effort`, `t`,
#' the complete `rates` list, `total` and `process = "consumers"`. It returns a
#' list with `causes` (a named list of finite numeric species-by-size matrices)
#' and `complete` (one logical value). Optional `metadata` is a data frame with
#' `id` and any of `label`, `owner`, `colour`, `linetype`; `id` must name each
#' cause once. Optional `unavailable` is a data frame with `id`, `owner` and
#' `reason` columns. It is only allowed for incomplete decompositions.
#' Callback errors and malformed results are errors. A complete decomposition
#' must satisfy `abs(total - sum(causes)) <= atol + rtol * abs(total)` at every
#' species, size and saved time.
#'
#' @param object A MizerParams or MizerSim object.
#' @param rtol Relative tolerance for the decomposition check. Default `1e-8`.
#' @param atol Absolute tolerance in 1/year. Default `1e-12`.
#' @param ... For MizerParams, `n`, `n_pp`, `n_other` and `effort` default to
#'   the model's initial values and `t` defaults to zero. For MizerSim,
#'   `time_range` selects saved times (all by default), `t` is an alias for
#'   `time_range`, and `drop = TRUE` drops dimensions of length one, as in
#'   [getMort()]. Simulation states and effort are taken from each saved time.
#'   Other arguments are passed to rate and diagnostic callbacks.
#' @return A `MizerMortalityBreakdown` list with `total`, a named list of
#'   `causes`, plotting `metadata`, logical `complete`, an `unavailable` data
#'   frame including time and reasons, and `residual` (total minus the sum of
#'   supplied causes). `total`, `causes` and `residual` use the same shapes
#'   as [getMort()]. The residual is a diagnostic difference, not an attributed
#'   cause. Simulation causes and metadata must remain consistent through time.
#' @seealso [getMort()], [setProcess()]
#' @family summary functions
#' @export
#' @examples
#' mortality <- getMortByCause(NS_params)
#' mortality$complete
#' mortality$metadata
#' plot(mortality$total)
getMortByCause <- function(object, rtol = 1e-8, atol = 1e-12, ...) {
    UseMethod("getMortByCause")
}

#' @rdname getMortByCause
#' @usage NULL
#' @export
getMortByCause.MizerParams <- function(object, rtol = 1e-8, atol = 1e-12,
                                       n = object$initial_n,
                                       n_pp = object$initial_n_pp,
                                       n_other = object$initial_n_other,
                                       effort = object$initial_effort,
                                       t = 0, ...) {
    # Defaults are lazy: bind them to the upgraded object, which may have gained
    # initial state fields absent from a saved S4 model.
    object <- validParams(object)
    params <- object
    prepared <- prepare_mortality_breakdown(params, rtol, atol)
    result <- evaluate_mortality_breakdown(
        params, prepared, n, n_pp, n_other, effort, t, ...
    )
    wrap <- function(value, label) {
        ArraySpeciesBySize(value,
            value_name = label, units = "1/year",
            params = params, representation = "average"
        )
    }
    result$total <- wrap(result$total, "Total mortality")
    result$residual <- wrap(result$residual, "Mortality decomposition residual")
    for (id in names(result$causes)) {
        label <- result$metadata$label[match(id, result$metadata$id)]
        result$causes[[id]] <- wrap(result$causes[[id]], label)
    }
    result
}

#' @rdname getMortByCause
#' @usage NULL
#' @export
getMortByCause.MizerSim <- function(object, rtol = 1e-8, atol = 1e-12,
                                    ..., time_range, drop = TRUE, t) {
    if (isS4(object) || needs_upgrading(object)) {
        object <- validSim(object)
    }
    params <- validParams(object$params)
    prepared <- prepare_mortality_breakdown(params, rtol, atol)
    if (!is.logical(drop) || length(drop) != 1L || is.na(drop)) {
        stop("`drop` must be TRUE or FALSE.", call. = FALSE)
    }
    if (missing(time_range) && !missing(t)) time_range <- t
    selected <- which(get_sim_rate_time_elements(object, time_range))
    if (!length(selected)) stop("No saved times were selected.", call. = FALSE)
    results <- lapply(selected, function(index) {
        slice <- get_sim_rate_slice(object, index)
        do.call(
            evaluate_mortality_breakdown,
            c(list(params = params, prepared = prepared), slice, list(...))
        )
    })
    reference <- results[[1L]]
    ids <- names(reference$causes)
    for (result in results) {
        if (!identical(names(result$causes), ids) ||
                !identical(result$metadata, reference$metadata)) {
            stop("Mortality causes and their metadata must be stable across ",
                "the selected simulation times.",
                call. = FALSE
            )
        }
    }
    times <- dimnames(object$n)$time[selected]
    wrap <- function(values, label) {
        stack_mortality_arrays(values, params, times, drop, label)
    }
    reference$total <- wrap(lapply(results, `[[`, "total"), "Total mortality")
    reference$residual <- wrap(
        lapply(results, `[[`, "residual"),
        "Mortality decomposition residual"
    )
    for (id in ids) {
        label <- reference$metadata$label[match(id, reference$metadata$id)]
        reference$causes[[id]] <- wrap(lapply(results, function(x) {
            x$causes[[id]]
        }), label)
    }
    reference$complete <- all(vapply(results, `[[`, logical(1), "complete"))
    reference$unavailable <- do.call(
        rbind, lapply(results, `[[`, "unavailable")
    )
    rownames(reference$unavailable) <- NULL
    reference
}

prepare_mortality_breakdown <- function(params, rtol, atol) {
    for (tolerance in list(rtol = rtol, atol = atol)) {
        if (!is.numeric(tolerance) || is.complex(tolerance) ||
                length(tolerance) != 1L || !is.finite(tolerance) ||
                tolerance < 0) {
            stop("`rtol` and `atol` must be finite nonnegative numbers.",
                call. = FALSE
            )
        }
    }
    definition <- params$capabilities$processes$consumers
    reason <- NULL
    if (is.null(definition$mortality_fun)) {
        reason <- "No mortality decomposition callback declared."
    } else if (!mortality_binding_current(params, definition)) {
        reason <- "Mortality decomposition does not match the current provider."
    }
    callback <- if (is.null(reason)) {
        resolve_model_callback(
            definition$mortality_fun,
            "Mortality decomposition callback"
        )
    } else {
        NULL
    }
    list(
        rates_fns = projectRateFunctions(params), callback = callback,
        reason = reason, rtol = rtol, atol = atol
    )
}

evaluate_mortality_breakdown <- function(params, prepared, n, n_pp, n_other,
                                         effort, t, ...) {
    if (!is.numeric(t) || length(t) != 1L || !is.finite(t)) {
        stop("`t` must be one finite time.", call. = FALSE)
    }
    rates <- prepared$rates_fns$Rates(
        params,
        n = n, n_pp = n_pp, n_other = n_other, effort = effort, t = t,
        rates_fns = prepared$rates_fns, ...
    )
    if (!is.list(rates)) {
        stop("The effective Rates provider must return a list.", call. = FALSE)
    }
    total <- checked_mortality_array(
        rates$mort, params,
        "Effective total mortality"
    )
    unavailable <- data.frame(
        id = character(), owner = character(),
        reason = character(), time = numeric()
    )
    if (is.null(prepared$callback)) {
        result <- list(causes = list(), complete = FALSE)
        unavailable <- data.frame(
            id = "mortality", owner = "consumers",
            reason = prepared$reason, time = t
        )
    } else {
        result <- invoke_tuning_callback(
            prepared$callback, params,
            c(list(
                n = n, n_pp = n_pp, n_other = n_other, effort = effort,
                t = t, rates = rates, total = total,
                process = "consumers"
            ), list(...)),
            "Mortality decomposition for consumers"
        )
    }
    if (!is.list(result) || is.object(result) ||
        !capability_names_valid(result) ||
        !all(c("causes", "complete") %in% names(result)) ||
        any(!names(result) %in% c(
            "causes", "complete", "metadata",
            "unavailable"
        )) ||
        !is.logical(result$complete) || length(result$complete) != 1L ||
        is.na(result$complete) || !is.list(result$causes) ||
        is.object(result$causes) ||
        !capability_names_valid(result$causes)) {
        stop("Mortality decomposition must return named `causes` and one ",
            "logical `complete` value.",
            call. = FALSE
        )
    }
    for (id in names(result$causes)) {
        result$causes[[id]] <- checked_mortality_array(
            result$causes[[id]], params, paste0("Mortality cause '", id, "'")
        )
    }
    metadata <- mortality_metadata(result$metadata, names(result$causes))
    if (!is.null(result$unavailable)) {
        missing <- result$unavailable
        if (!is.data.frame(missing) ||
            !all(c("id", "owner", "reason") %in% names(missing)) ||
            any(!vapply(missing[c("id", "owner", "reason")], function(x) {
                is.character(x) && !anyNA(x) && all(nzchar(x))
            }, logical(1)))) {
            stop("Unavailable mortality causes need `id`, `owner` and ",
                "`reason` character columns.",
                call. = FALSE
            )
        }
        missing <- missing[c("id", "owner", "reason")]
        missing$time <- rep(t, nrow(missing))
        unavailable <- rbind(unavailable, missing)
    }
    if (result$complete && nrow(unavailable)) {
        stop("A complete mortality decomposition cannot have unavailable ",
            "causes.",
            call. = FALSE
        )
    }
    if (!result$complete && !nrow(unavailable)) {
        unavailable <- data.frame(
            id = "mortality", owner = "consumers",
            reason = "The callback declares an incomplete decomposition.",
            time = t
        )
    }
    supplied <- Reduce(`+`, result$causes, init = total * 0)
    residual <- total - supplied
    if (any(!is.finite(supplied)) || any(!is.finite(residual))) {
        stop("The sum of mortality causes and its residual must be finite.",
             call. = FALSE)
    }
    if (result$complete &&
            any(abs(residual) > prepared$atol + prepared$rtol * abs(total))) {
        stop("Complete mortality decomposition disagrees with the effective ",
            "total mortality at time ", t, ".",
            call. = FALSE
        )
    }
    structure(list(
        total = total, causes = result$causes, metadata = metadata,
        complete = result$complete, unavailable = unavailable,
        residual = residual
    ), class = "MizerMortalityBreakdown")
}

checked_mortality_array <- function(value, params, context) {
    expected <- dim(params$initial_n)
    if (!is.numeric(value) || is.complex(value) ||
            !identical(dim(value), expected) || any(!is.finite(value))) {
        stop(context, " must be a finite numeric species-by-size matrix.",
            call. = FALSE
        )
    }
    expected_names <- dimnames(params$initial_n)
    provided <- dimnames(value)
    for (dimension in seq_along(expected)) {
        matches <- identical(provided[[dimension]], expected_names[[dimension]])
        if (!is.null(provided[[dimension]]) && !matches) {
            stop(context, " has incompatible species or size names.",
                call. = FALSE
            )
        }
    }
    array(as.numeric(value), dim = expected, dimnames = expected_names)
}

mortality_metadata <- function(metadata, ids) {
    defaults <- data.frame(
        id = ids %||% character(),
        label = ids %||% character(),
        owner = rep("consumers", length(ids)),
        colour = rep("grey", length(ids)),
        linetype = rep("solid", length(ids))
    )
    if (is.null(metadata)) {
        return(defaults)
    }
    if (!is.data.frame(metadata) || !"id" %in% names(metadata) ||
        any(!names(metadata) %in% names(defaults)) ||
        anyDuplicated(metadata$id) || !setequal(metadata$id, ids) ||
        any(!vapply(metadata, function(x) {
            is.character(x) && !anyNA(x) && all(nzchar(x))
        }, logical(1)))) {
        stop("Mortality metadata must identify every cause exactly once ",
            "with character plotting fields.",
            call. = FALSE
        )
    }
    metadata <- metadata[match(ids, metadata$id), , drop = FALSE]
    rownames(metadata) <- NULL
    defaults[names(metadata)] <- metadata
    defaults
}

builtin_mortality_breakdown <- function(params, n, n_pp, n_other, rates, t,
                                        effort, total, process, ...) {
    causes <- list(
        predation = rates$pred_mort, fishing = rates$f_mort,
        external = params$mu_b
    )
    metadata <- data.frame(
        id = names(causes),
        label = c("Predation", "Fishing", "External"),
        owner = rep("consumers", 3),
        colour = c("#D55E00", "blue", "purple"),
        linetype = rep("solid", 3)
    )
    for (index in seq_len(nrow(metadata))) {
        label <- metadata$label[[index]]
        if (label %in% names(params$linecolour)) {
            metadata$colour[[index]] <- params$linecolour[[label]]
        }
        if (label %in% names(params$linetype)) {
            metadata$linetype[[index]] <- as.character(params$linetype[[label]])
        }
    }
    for (name in names(params$other_mort)) {
        component <- name %in% names(params$other_dynamics)
        id <- paste0(if (component) "component:" else "contribution:", name)
        fun <- resolve_model_callback(
            params$other_mort[[name]],
            paste0(
                "Mortality contribution '", name,
                "'"
            )
        )
        # Match the additive callback's inputs in projectMort(): the aggregate
        # diagnostic receives the full context but does not reinterpret terms.
        causes[[id]] <- fun(
            params = params, n = n, n_pp = n_pp,
            n_other = n_other, t = t, component = name, ...
        )
        metadata <- rbind(metadata, data.frame(
            id = id, label = name, owner = id,
            colour = if (name %in% names(params$linecolour)) {
                params$linecolour[[name]]
            } else {
                "grey"
            },
            linetype = if (name %in% names(params$linetype)) {
                as.character(params$linetype[[name]])
            } else {
                "solid"
            }
        ))
    }
    list(causes = causes, complete = TRUE, metadata = metadata)
}

stack_mortality_arrays <- function(values, params, times, drop, label) {
    dimensions <- dim(params$initial_n)
    array <- array(unlist(values, use.names = FALSE),
        dim = c(dimensions, length(values))
    )
    array <- aperm(array, c(3L, 1L, 2L))
    dimnames(array) <- c(list(time = times), dimnames(params$initial_n))
    result <- array[, , , drop = drop]
    if (length(dim(result)) == 3L) {
        result <- ArrayTimeBySpeciesBySize(
            result,
            value_name = label, units = "1/year", params = params,
            representation = "average"
        )
    } else if (is.matrix(result) && names(dimnames(result))[[1L]] == "sp") {
        result <- ArraySpeciesBySize(result,
            value_name = label,
            units = "1/year", params = params,
            representation = "average"
        )
    }
    result
}

#' Print a mortality decomposition report
#'
#' @param x A MizerMortalityBreakdown object from [getMortByCause()].
#' @param ... Unused.
#' @return The report, invisibly.
#' @keywords internal
#' @export
print.MizerMortalityBreakdown <- function(x, ...) {
    cat(
        "Mortality decomposition:",
        if (x$complete) "complete" else "incomplete", "[1/year]\n"
    )
    if (length(x$causes)) {
        print(x$metadata, row.names = FALSE)
    }
    if (nrow(x$unavailable)) {
        print(x$unavailable, row.names = FALSE)
    }
    invisible(x)
}
