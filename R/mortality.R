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
#' @param predation_by_species Whether to refine built-in predation into one
#'   cause per predator species. The refinement is used only when its sum
#'   agrees with the effective predation mortality at every species and size;
#'   otherwise the aggregate predation cause is retained. For the built-in
#'   callback and its public alias, any selected simulation time requiring
#'   aggregate predation makes every selected time use that grouping. Custom
#'   callbacks are responsible for stable grouping across times.
#'   Default `FALSE`.
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
getMortByCause <- function(object, rtol = 1e-8, atol = 1e-12, ...,
                           predation_by_species = FALSE) {
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
                                       t = 0, ...,
                                       predation_by_species = FALSE) {
    # Defaults are lazy: bind them to the upgraded object, which may have gained
    # initial state fields absent from a saved S4 model.
    object <- validParams(object)
    params <- object
    prepared <- prepare_mortality_breakdown(
        params, rtol, atol, predation_by_species
    )
    result <- evaluate_mortality_breakdown(
        params, prepared, n, n_pp, n_other, effort, t, ...
    )
    attr(result, "predation_aggregate") <- NULL
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
                                    ..., time_range, drop = TRUE, t,
                                    predation_by_species = FALSE) {
    if (isS4(object) || needs_upgrading(object)) {
        object <- validSim(object)
    }
    params <- validParams(object$params)
    prepared <- prepare_mortality_breakdown(
        params, rtol, atol, predation_by_species
    )
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
    times <- dimnames(object$n)$time[selected]
    results <- uniform_mortality_predation(results, prepared, as.numeric(times))
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

prepare_mortality_breakdown <- function(params, rtol, atol,
                                        predation_by_species = FALSE) {
    if (!is.logical(predation_by_species) ||
            length(predation_by_species) != 1L || is.na(predation_by_species)) {
        stop("`predation_by_species` must be TRUE or FALSE.", call. = FALSE)
    }
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
        reason = reason, rtol = rtol, atol = atol,
        predation_by_species = predation_by_species
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
                process = "consumers",
                predation_by_species = prepared$predation_by_species,
                rtol = prepared$rtol, atol = prepared$atol
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
    residual <- mortality_residual(total, result$causes, result$complete,
                                   prepared$rtol, prepared$atol, t)
    report <- structure(list(
        total = total, causes = result$causes, metadata = metadata,
        complete = result$complete, unavailable = unavailable,
        residual = residual
    ), class = "MizerMortalityBreakdown")
    if (identical(prepared$callback, builtin_mortality_breakdown)) {
        attr(report, "predation_aggregate") <-
            attr(result, "predation_aggregate")
    }
    report

}

mortality_residual <- function(total, causes, complete, rtol, atol, t) {
    supplied <- Reduce(`+`, causes, init = total * 0)
    residual <- total - supplied
    if (any(!is.finite(supplied)) || any(!is.finite(residual))) {
        stop("The sum of mortality causes and its residual must be finite.",
             call. = FALSE)
    }
    if (complete && any(abs(residual) > atol + rtol * abs(total))) {
        stop("Complete mortality decomposition disagrees with the effective ",
             "total mortality at time ", t, ".", call. = FALSE)
    }
    residual
}

# A built-in refinement must use the same grouping across selected times. Keep
# the original aggregate privately until that decision; never rerun Rates or
# reconstruct its exact value by summing approximate diagnostic contributions.
uniform_mortality_predation <- function(results, prepared, times) {
    fallback <- isTRUE(prepared$predation_by_species) &&
        identical(prepared$callback, builtin_mortality_breakdown) &&
        any(vapply(results, function(result) {
            "predation" %in% names(result$causes)
        }, logical(1)))
    lapply(seq_along(results), function(index) {
        result <- results[[index]]
        aggregate <- attr(result, "predation_aggregate")
        if (fallback && !is.null(aggregate)) {
            other <- !startsWith(names(result$causes), "predator:") &
                names(result$causes) != "predation"
            ids <- names(result$causes)[other]
            result$causes <- c(list(predation = aggregate$value),
                               result$causes[ids])
            rows <- match(ids, result$metadata$id)
            result$metadata <- rbind(aggregate$metadata,
                                     result$metadata[rows, ])
            rownames(result$metadata) <- NULL
            result$residual <- mortality_residual(
                result$total, result$causes, result$complete,
                prepared$rtol, prepared$atol, times[[index]]
            )
        }
        attr(result, "predation_aggregate") <- NULL
        result
    })
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
                                        effort, total, process,
                                        predation_by_species = FALSE,
                                        rtol = 1e-8, atol = 1e-12, ...) {
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
    aggregate <- NULL
    if (predation_by_species) {
        aggregate <- list(
            value = checked_mortality_array(rates$pred_mort, params,
                                            "Effective predation mortality"),
            metadata = style_mortality_metadata(params, metadata[1L, ])
        )
        predators <- mortality_predator_causes(params, rates, rtol, atol)
        if (!is.null(predators)) {
            causes <- c(predators, causes[c("fishing", "external")])
            metadata <- rbind(data.frame(
                id = names(predators),
                label = as.character(params$species_params$species),
                owner = rep("consumers", length(predators)),
                colour = rep("#D55E00", length(predators)),
                linetype = rep("solid", length(predators))
            ), metadata[metadata$id != "predation", ])
        }
    }
    metadata <- style_mortality_metadata(params, metadata)
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
    result <- list(causes = causes, complete = TRUE, metadata = metadata)
    attr(result, "predation_aggregate") <- aggregate
    result
}

style_mortality_metadata <- function(params, metadata) {
    for (index in seq_len(nrow(metadata))) {
        label <- metadata$label[[index]]
        if (label %in% names(params$linecolour)) {
            metadata$colour[[index]] <- params$linecolour[[label]]
        }
        if (label %in% names(params$linetype)) {
            metadata$linetype[[index]] <- as.character(params$linetype[[label]])
        }
    }
    metadata
}

# Refine a known contribution, never replace the authoritative aggregate with
# an assumed decomposition. All inputs come from the same full rate evaluation.
mortality_predator_causes <- function(params, rates, rtol, atol) {
    predation <- checked_mortality_array(rates$pred_mort, params,
                                         "Effective predation mortality")
    pred_rate <- rates$pred_rate
    expected <- c(nrow(params$initial_n), length(params$w_full))
    if (!is.numeric(pred_rate) || is.complex(pred_rate) ||
            !identical(dim(pred_rate), expected) ||
            any(!is.finite(pred_rate))) {
        stop("Predation rates must be a finite numeric predator-by-size ",
             "matrix on the full size grid.", call. = FALSE)
    }
    provided <- dimnames(pred_rate)
    if (!is.null(provided[[1L]]) &&
            !identical(provided[[1L]],
                       as.character(params$species_params$species))) {
        stop("Predation rates have incompatible predator names.", call. = FALSE)
    }
    index <- seq.int(length(params$w_full) - length(params$w) + 1L,
                     length(params$w_full))
    species <- as.character(params$species_params$species)
    causes <- setNames(lapply(seq_along(species), function(predator) {
        array(as.numeric(outer(params$interaction[predator, ],
                               pred_rate[predator, index])),
              dim = dim(params$initial_n),
              dimnames = dimnames(params$initial_n))
    }), paste0("predator:", species))
    supplied <- Reduce(`+`, causes, init = predation * 0)
    if (any(!is.finite(supplied))) {
        stop("Predator-specific mortality contributions must be finite.",
             call. = FALSE)
    }
    if (any(abs(supplied - predation) > atol + rtol * abs(predation))) {
        return(NULL)
    }
    causes
}

#' Build the standard mortality diagnostic from supplied rates
#'
#' An author helper for a registered `mortality_fun` callback. It partitions
#' the supplied effective rates into predation, fishing, external mortality and
#' registered additive contributions without evaluating the rate pipeline
#' again. An extension that changes the outer mortality assembler can call this
#' helper, append its additional terms, and return the resulting specification.
#'
#' This helper returns an unchecked callback specification. The complete sum is
#' checked against the effective total by [getMortByCause()] after the author's
#' callback returns, so terms added by that callback participate in the check.
#'
#' @param params A MizerParams object.
#' @param n,n_pp,n_other Current consumer, resource and component states.
#' @param rates The complete effective rate list supplied to the callback.
#' @param t Current model time.
#' @param effort Current fishing effort.
#' @param total Effective total mortality supplied to the callback.
#' @param process The diagnostic owner, `"consumers"`.
#' @inheritParams getMortByCause
#' @param ... Further callback context.
#' @return A raw callback specification with `causes`, `metadata` and
#'   `complete = TRUE`. Authors must append every additional mortality term
#'   before returning it as a complete decomposition.
#' @seealso [setProcess()], [getMortByCause()]
#' @export
mizerMortalityBreakdown <- builtin_mortality_breakdown

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
