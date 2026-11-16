#' Balance the active reproduction law while preserving dynamical rates
#'
#' Adjusts only parameters authorised by a stored, provider-bound tuning
#' strategy. By default the current effective recruitment is retained, even
#' when the initial state is not steady. `rdd_target` supplies explicit targets
#' for a tuning or calibration operation; unselected species retain their
#' original recruitment and parameters.
#'
#' The complete effective rate pipeline is evaluated before and after the
#' candidate change. Mizer checks unchanged state, declared parameter bounds
#' and preserved quantities, and every effective rate entering the equations.
#' Unsupported providers and infeasible controls are errors. There is no
#' numerical fallback for an unknown reproduction function.
#'
#' Built-in laws support explicit `erepro`, their density parameter (`R_max`,
#' `ricker_b` or `sheperd_b`), or `reproduction_level` in `...`. Supply at most
#' one of these controls. Shepherd also accepts `sheperd_c`. Ricker and
#' overcompensating Shepherd inversions use `limb = "ascending"` by default;
#' choose `limb = "descending"` explicitly for the other finite solution.
#' Built-in defaults preserve their law-specific reproduction level, except
#' unbounded Shepherd curves, which preserve `sheperd_b`. `noRDD` adjusts
#' `erepro`, and `constantRDD` adjusts `constant_reproduction`.
#'
#' Extension methods may transform targets and delegate with `NextMethod()`.
#' The outer transaction verifies the final result once, after the full method
#' chain. Methods must accept `...` to pass through runtime context. Providing
#' a method alone does not declare support: the effective provider also needs
#' an explicit bound tuning strategy.
#'
#' @param params A MizerParams object.
#' @param rdd_target Optional scalar, full vector, or named partial vector of
#'   recruitment targets. `NULL` preserves current effective recruitment.
#'   `NA` leaves a species unselected.
#' @param strategy Name of a registered reproduction strategy. By default the
#'   provider's declared default is used.
#' @param preserve Alternative to `strategy`, selecting the corresponding
#'   `preserve_...` strategy. Built-ins support `"erepro"`,
#'   `"reproduction_level"`, and their density parameter where applicable.
#' @param ... Law-specific controls described above, or controls accepted by
#'   an extension's balancing callback.
#' @param effort Fishing effort at the evaluation state.
#' @param t Time at which to evaluate the rates. Stationarity must be declared
#'   by the provider; it is never inferred by probing other times.
#' @param info_level Controls information and warnings; see
#'   [default_info_level()].
#' @return The verified MizerParams candidate. Errors leave the input unchanged.
#' @seealso [setBevertonHolt()], [reproduction_level()], [setProcess()]
#' @export
#' @family functions for setting parameters
balanceRDD <- function(params, rdd_target = NULL, strategy = NULL,
                       preserve = NULL, ...,
                       effort = params@initial_effort, t = 0,
                       info_level = default_info_level()) {
    with_info_level(info_level = info_level, {
        params <- validParams(params)
        controls <- list(...)
        if (length(controls) && (is.null(names(controls)) ||
                                     any(!nzchar(names(controls))))) {
            stop("All reproduction controls in ... must be named.",
                call. = FALSE
            )
        }
        reserved <- intersect(names(controls), c(
            "rdi", "rates", "n", "n_pp",
            "n_other",
            "targets", "preserved",
            "species", "settings",
            ".balance_context"
        ))
        if (length(reserved)) {
            stop("Balancing context cannot be supplied independently: ",
                paste(reserved, collapse = ", "), ".",
                call. = FALSE
            )
        }
        parameters <- c(
            "erepro", "R_max", "reproduction_level", "ricker_b",
            "sheperd_b"
        )
        supplied <- intersect(names(controls), parameters)
        supplied <- supplied[!vapply(controls[supplied], is.null, logical(1))]
        if (length(supplied) > 1L ||
                (length(supplied) && !is.null(preserve))) {
            stop("Choose one explicit reproduction control or one preserved ",
                "quantity, not both.",
                call. = FALSE
            )
        }
        selected <- as.character(params@species_params$species)
        requested <- reproduction_species_values(
            params, rdd_target,
            "rdd_target"
        )
        if (!is.null(requested)) {
            if (any(requested < 0, na.rm = TRUE)) {
                stop("`rdd_target` must be nonnegative or NA.", call. = FALSE)
            }
            selected <- names(requested)[!is.na(requested)]
        }
        for (field in intersect(names(controls), c(parameters, "sheperd_c"))) {
            controls[[field]] <- reproduction_species_values(
                params, controls[[field]], field,
                allow_infinite = field == "R_max"
            )
            if (field %in% supplied) {
                selected <- intersect(
                    selected,
                    names(controls[[field]])[!is.na(controls[[field]])]
                )
            }
        }
        if (!length(selected)) {
            return(params)
        }
        if (!is.null(builtin_reproduction_law(params)) && length(supplied)) {
            if (!is.null(strategy) && strategy != "set_parameters") {
                stop("Explicit parameter values require ",
                    "strategy = 'set_parameters'.",
                    call. = FALSE
                )
            }
            strategy <- "set_parameters"
        }
        plan <- prepare_tuning(
            params,
            owners = "reproduction", solver = "balance",
            strategies = if (is.null(strategy)) {
                NULL
            } else {
                c(reproduction = strategy)
            },
            preserve = preserve, species = selected, effort = effort, t = t
        )
        full_target <- stats::setNames(
            as.numeric(plan$rates$rdd),
            plan$all_species
        )
        if (!is.null(requested)) full_target[selected] <- requested[selected]
        # Keep the effective vector's original names for invariant comparison.
        names(full_target) <- names(plan$rates$rdd)
        targets <- capture_tuning_targets(params, plan,
            rates = plan$rates,
            t = t, rdd_target = full_target
        )
        dispatched_target <- stats::setNames(
            rep(NA_real_, length(full_target)),
            plan$all_species
        )
        dispatched_target[selected] <- as.numeric(full_target)[
            match(
                selected,
                plan$all_species
            )
        ]
        context <- list(
            plan = plan, targets = targets, rates = plan$rates,
            controls = controls
        )
        args <- c(
            list(params,
                rdd_target = dispatched_target,
                strategy = strategy,
                preserve = preserve, effort = effort, t = t
            ), controls,
            list(.balance_context = context)
        )
        candidate <- do.call(balanceRDD_dispatch, args)
        verify_tuning_candidate(candidate, plan, targets)
    })
}
# nolint next: object_name_linter.
balanceRDD_dispatch <- function(params, ...) UseMethod("balanceRDD")
#' @rdname balanceRDD
#' @usage NULL
#' @exportS3Method balanceRDD MizerParams
balanceRDD.MizerParams <- function(params, rdd_target = NULL, strategy = NULL,
                                   preserve = NULL, ...,
                                   effort = params@initial_effort, t = 0,
                                   .balance_context = NULL) {
    context <- .balance_context
    if (is.null(context)) {
        stop("Call balanceRDD() to prepare and verify ",
            "the balancing transaction.",
            call. = FALSE
        )
    }
    node <- context$plan$selected$reproduction
    if (is.null(node)) stop("No local reproduction strategy was selected.")
    if (node$strategy$action == "none") {
        return(params)
    }
    controls <- list(...)
    settings <- node$strategy$settings
    if (!is.null(preserve)) settings$preserve <- preserve
    settings[names(controls)] <- controls
    if (!is.null(builtin_reproduction_law(params))) {
        unsupported <- setdiff(names(settings), c(
            "preserve", "limb", "erepro",
            "R_max",
            "reproduction_level",
            "ricker_b", "sheperd_b",
            "sheperd_c"
        ))
        if (length(unsupported)) {
            stop("Unsupported reproduction controls: ",
                paste(unsupported,
                    collapse = ", "
                ), ".",
                call. = FALSE
            )
        }
    }
    args <- tuning_context_args(params, context$rates, t, effort)
    args <- c(args, list(
        rdd_target = rdd_target, strategy = node$id,
        preserved = node$preserved,
        species = context$plan$species,
        settings = settings
    ), settings)
    invoke_tuning_callback(
        node$balance_fun, params, args,
        "Reproduction balancing callback"
    )
}
# nolint next: object_length_linter.
require_reproduction_balance_support <- function(params, strategy = NULL,
                                                 preserve = NULL,
                                                 partial = FALSE) {
    plan <- prepare_tuning(
        params,
        owners = "reproduction", solver = "balance",
        strategies = if (is.null(strategy)) {
            NULL
        } else {
            c(reproduction = strategy)
        },
        preserve = preserve
    )
    if (partial && any(vapply(plan$selected, function(node) {
        !node$strategy$partial_species
    }, logical(1)))) {
        stop("The declared reproduction strategy does not support partial ",
            "species requests.",
            call. = FALSE
        )
    }
    invisible(plan)
}
# The built-in laws use their own inverse and level, never a numerical fit to
# an unknown callback. Their use is authorised by the stored provider binding.
reproduction_law_name <- function(params) {
    sub("^mizer::", "", params@rates_funcs$RDD)
}
reproduction_species_values <- function(params, values, argument,
                                        allow_infinite = FALSE) {
    if (is.null(values)) {
        return(NULL)
    }
    if (is.logical(values) && all(is.na(values))) values <- as.numeric(values)
    if (!is.numeric(values) || is.complex(values) ||
            any(is.nan(values)) ||
            (!allow_infinite && any(is.infinite(values)))) {
        stop("`", argument, "` must contain numeric values",
            if (!allow_infinite) " that are finite or NA" else " or NA",
            ".",
            call. = FALSE
        )
    }
    species <- as.character(params@species_params$species)
    if (is.null(names(values))) {
        if (length(values) == 1L) values <- rep(values, length(species))
        if (length(values) != length(species)) {
            stop("`", argument, "` needs one value, a vector of length ",
                length(species), ", or a named species vector.",
                call. = FALSE
            )
        }
        return(stats::setNames(values, species))
    }
    if (anyNA(names(values)) || any(!nzchar(names(values))) ||
            anyDuplicated(names(values)) || any(!names(values) %in% species)) {
        stop("`", argument, "` must have unique names of species in the model.",
            call. = FALSE
        )
    }
    result <- stats::setNames(rep(NA_real_, length(species)), species)
    result[names(values)] <- values
    result
}
set_reproduction_parameter <- function(params, field, selected, values,
                                       record = FALSE) {
    no_sp <- nrow(params@species_params)
    if (!field %in% names(params@species_params)) {
        params@species_params[[field]] <- rep(NA_real_, no_sp)
    }
    old <- params@species_params[[field]][selected]
    unchanged <- isTRUE(all.equal(as.numeric(values), as.numeric(old),
        tolerance = 1e-12
    ))
    if (unchanged && !record) {
        return(params)
    }
    if (!field %in% names(params@given_species_params)) {
        params@given_species_params[[field]] <- rep(NA_real_, no_sp)
    }
    params@species_params[[field]][selected] <- values
    params@given_species_params[[field]][selected] <- values
    if (field == "erepro") {
        unusual <- values > 1 & !unchanged
        if (any(unusual)) {
            signal_info("erepro", paste0(
                "The following species require an unrealistic value greater ",
                "than 1 for `erepro`: ",
                paste(params@species_params$species[selected][unusual],
                    collapse = ", "
                )
            ),
            level = 1, severity = "warning", unhandled = "show"
            )
        }
    }
    params
}
# The proportional relationship is part of the supported built-in RDI route.
# Its zero-efficiency case is evaluated from the same reproductive energy,
# rather than dividing zero by zero or perturbing an unknown RDI callback.
# nolint next: object_length_linter.
reproduction_efficiency_coefficient <- function(params, rates) {
    efficiency <- params@species_params$erepro
    result <- rates$rdi / efficiency
    zero <- efficiency == 0
    if (any(zero)) {
        unit <- params
        unit@species_params$erepro[] <- 1
        result[zero] <- mizerRDI(
            unit,
            n = params@initial_n, n_pp = params@initial_n_pp,
            n_other = params@initial_n_other,
            e_growth = rates$e_growth, mort = rates$mort,
            e_repro = rates$e_repro, diffusion = rates$diffusion
        )[zero]
    }
    if (any(!is.finite(result)) || any(result < 0)) {
        stop("The declared proportional reproduction route has an invalid ",
            "reproductive-efficiency coefficient.",
            call. = FALSE
        )
    }
    result
}
reproduction_balance_error <- function(species, reason) {
    stop("Cannot balance reproduction for '", species, "': ", reason,
        call. = FALSE
    )
}
ricker_inverse <- function(target, b, limb, species) {
    if (b == 0) {
        if (limb == "descending") {
            reproduction_balance_error(
                species, paste0(
                    "Ricker reproduction with ricker_b = 0 ",
                    "has no descending limb."
                )
            )
        }
        return(target)
    }
    if (target == 0) {
        if (limb == "descending") {
            reproduction_balance_error(
                species, paste0(
                    "zero recruitment has no finite ",
                    "descending-limb solution."
                )
            )
        }
        return(0)
    }
    log_q <- log(b) + log(target)
    if (log_q > -1 + 64 * .Machine$double.eps) {
        reproduction_balance_error(
            species, paste0(
                "the target exceeds the maximum of the ",
                "specified Ricker curve."
            )
        )
    }
    if (abs(log_q + 1) <= 64 * .Machine$double.eps) {
        return(1 / b)
    }
    if (limb == "ascending") {
        log_z <- stats::uniroot(function(u) u - exp(u) - log_q,
            c(log_q, 0),
            tol = 1e-12
        )$root
        return(target * exp(exp(log_z)))
    }
    upper <- log(-log_q + log(-log_q) + 2)
    log_z <- stats::uniroot(function(u) u - exp(u) - log_q,
        c(0, upper),
        tol = 1e-12
    )$root
    exp(log_z) / b
}
log_one_plus_exp <- function(value) {
    pmax(value, 0) + log1p(exp(-abs(value)))
}
shepherd_peak_factor <- function(shape) {
    exp((1 - 1 / shape) * log(shape - 1) - log(shape))
}
shepherd_inverse <- function(target, b, shape, limb, species) {
    if (b == 0) {
        if (limb == "descending") {
            reproduction_balance_error(
                species, paste0(
                    "Shepherd reproduction with sheperd_b = 0 ",
                    "has no descending limb."
                )
            )
        }
        return(target)
    }
    if (shape <= 1 && limb == "descending") {
        reproduction_balance_error(
            species, paste0(
                "Shepherd reproduction has no descending ",
                "limb when sheperd_c <= 1."
            )
        )
    }
    if (target == 0) {
        if (limb == "descending") {
            reproduction_balance_error(
                species, paste0(
                    "zero recruitment has no finite ",
                    "descending-limb solution."
                )
            )
        }
        return(0)
    }
    if (shape == 1) {
        level <- b * target
        if (level >= 1) {
            reproduction_balance_error(
                species, paste0(
                    "the target must be below the specified ",
                    "Shepherd asymptote."
                )
            )
        }
        return(target / (1 - level))
    }
    log_b <- log(b)
    log_target <- log(target)
    equation <- function(log_x) {
        log_x - log_one_plus_exp(shape * (log_b + log_x)) - log_target
    }
    if (shape > 1) {
        peak <- -log_b - log(shape - 1) / shape
        maximum <- log(shepherd_peak_factor(shape)) - log_b
        if (log_target > maximum + 64 * .Machine$double.eps) {
            reproduction_balance_error(
                species, paste0(
                    "the target exceeds the maximum of the ",
                    "specified Shepherd curve."
                )
            )
        }
        if (abs(log_target - maximum) <= 64 * .Machine$double.eps) {
            return(exp(peak))
        }
        interval <- if (limb == "ascending") {
            c(log_target, peak)
        } else {
            c(peak, max(
                peak + 1,
                (-log_target - shape * log_b) / (shape - 1) + 1
            ))
        }
    } else {
        upper <- max(
            -log_b, log_target + log(2),
            (log_target + shape * log_b + log(2)) / (1 - shape)
        )
        interval <- c(log_target, upper)
    }
    exp(stats::uniroot(equation, interval, tol = 1e-12)$root)
}
# Each law's level is recruitment divided by its own maximal recruitment.
# A zero constant has a zero maximum too, so its normalised level is undefined.
# nolint next: object_length_linter, object_name_linter.
reproduction_level_BevertonHoltRDD <- function(params, rates, ...) {
    stats::setNames(
        rates$rdd / params@species_params$R_max,
        params@species_params$species
    )
}
# nolint next: object_name_linter.
reproduction_level_RickerRDD <- function(params, rates, ...) {
    stats::setNames(
        exp(1) * params@species_params$ricker_b * rates$rdd,
        params@species_params$species
    )
}
# nolint next: object_name_linter.
reproduction_level_SheperdRDD <- function(params, rates, ...) {
    shape <- params@species_params$sheperd_c
    b <- params@species_params$sheperd_b
    value <- numeric(length(rates$rdd))
    value[shape == 1] <- (b * rates$rdd)[shape == 1]
    over <- shape > 1
    value[over] <- b[over] * rates$rdd[over] / shepherd_peak_factor(shape[over])
    stats::setNames(value, params@species_params$species)
}
# nolint next: object_name_linter.
reproduction_level_noRDD <- function(params, rates, ...) {
    stats::setNames(rep(0, length(rates$rdd)), params@species_params$species)
}
# nolint next: object_name_linter.
reproduction_level_constantRDD <- function(params, rates, ...) {
    stats::setNames(
        ifelse(rates$rdd == 0, NA_real_, 1),
        params@species_params$species
    )
}
# Common validation and parameter recording surround the law-specific algebra.
balance_reproduction_law <- function(params, rates, rdd_target, law,
                                     preserve = "erepro", preserved = list(),
                                     # nolint next: object_name_linter.
                                     erepro = NULL, R_max = NULL,
                                     reproduction_level = NULL,
                                     ricker_b = NULL,
                                     sheperd_b = NULL, sheperd_c = NULL,
                                     limb = c("ascending", "descending"),
                                     species = NULL, ...) {
    limb <- match.arg(limb)
    sp <- params@species_params
    targets <- reproduction_species_values(params, rdd_target, "rdd_target")
    if (!is.null(species)) targets[!names(targets) %in% species] <- NA_real_
    if (any(targets < 0, na.rm = TRUE)) {
        stop("`rdd_target` must be nonnegative or NA.", call. = FALSE)
    }
    controls <- list(
        erepro = erepro, R_max = R_max,
        reproduction_level = reproduction_level,
        ricker_b = ricker_b, sheperd_b = sheperd_b
    )
    controls <- controls[!vapply(controls, is.null, logical(1))]
    allowed <- switch(law,
        BevertonHoltRDD = c(
            "erepro",
            "R_max", "reproduction_level"
        ),
        RickerRDD = c("erepro", "ricker_b", "reproduction_level"),
        SheperdRDD = c("erepro", "sheperd_b", "reproduction_level"),
        noRDD = "reproduction_level",
        constantRDD = "reproduction_level"
    )
    if (length(controls) > 1L || any(!names(controls) %in% allowed)) {
        stop("Choose one supported reproduction control for ", law, ": ",
            paste(allowed, collapse = ", "), ".",
            call. = FALSE
        )
    }
    selector <- if (length(controls)) names(controls)[[1L]] else preserve
    if (is.null(selector)) selector <- "erepro"
    values <- if (length(controls)) {
        reproduction_species_values(params, controls[[1L]], selector,
            allow_infinite = selector == "R_max"
        )
    } else if (selector %in% names(preserved)) {
        reproduction_species_values(params, preserved[[selector]], selector,
            allow_infinite = selector == "R_max"
        )
    } else if (selector == "reproduction_level") {
        fun <- get(paste0("reproduction_level_", law), mode = "function")
        fun(params, rates)
    } else {
        sp[[selector]]
    }
    if (law == "noRDD" && selector != "reproduction_level") selector <- "none"
    if (law == "constantRDD" && selector != "reproduction_level") {
        selector <- "none"
    }
    if (!selector %in% c(allowed, "none")) {
        stop("Unsupported preserved quantity '", selector, "' for ", law, ".",
            call. = FALSE
        )
    }
    if (selector == "none" && !law %in% c("noRDD", "constantRDD")) {
        stop("The set_parameters strategy requires an explicit law control.",
            call. = FALSE
        )
    }
    if (selector != "none" && is.null(values)) {
        stop("Missing values for reproduction control '", selector, "'.",
            call. = FALSE
        )
    }
    selected <- which(!is.na(targets))
    if (length(controls)) selected <- intersect(selected, which(!is.na(values)))
    if (!length(selected)) {
        return(params)
    }
    if (!is.null(sheperd_c) && law != "SheperdRDD") {
        stop("`sheperd_c` only applies to Shepherd reproduction.",
            call. = FALSE
        )
    }
    shape <- sp$sheperd_c
    if (law == "SheperdRDD") {
        if (!is.null(sheperd_c)) {
            supplied_shape <- reproduction_species_values(
                params, sheperd_c,
                "sheperd_c"
            )
            changed <- which(!is.na(supplied_shape))
            shape[changed] <- supplied_shape[changed]
        }
        if (is.null(shape) || any(!is.finite(shape[selected])) ||
                any(shape[selected] <= 0)) {
            stop("`sheperd_c` must be finite and positive.", call. = FALSE)
        }
    }
    coefficient <- if (law == "constantRDD") {
        NULL
    } else {
        reproduction_efficiency_coefficient(params, rates)
    }
    for (i in selected) {
        species <- as.character(sp$species[[i]])
        target <- targets[[i]]
        value <- if (selector == "none") NULL else values[[i]]
        if (!is.null(value) && (is.na(value) || value < 0)) {
            reproduction_balance_error(
                species,
                paste0(
                    "the specified ",
                    selector, " must be nonnegative."
                )
            )
        }
        if (law == "constantRDD") {
            if (selector == "reproduction_level" &&
                    (target == 0 || value != 1)) {
                reproduction_balance_error(
                    species, paste0(
                        "positive constant reproduction has level ",
                        "1; zero has no level."
                    )
                )
            }
            params <- set_reproduction_parameter(
                params,
                "constant_reproduction",
                i, target
            )
            next
        }
        if (law == "noRDD") {
            if (selector == "reproduction_level" && value != 0) {
                reproduction_balance_error(species, "noRDD has level 0.")
            }
            required_rdi <- target
            density <- NULL
        } else {
            density_field <- switch(law,
                BevertonHoltRDD = "R_max",
                RickerRDD = "ricker_b",
                SheperdRDD = "sheperd_b"
            )
            density <- sp[[density_field]][[i]]
            if (selector == "erepro") {
                required_rdi <- coefficient[[i]] * value
                if (required_rdi < target) {
                    reproduction_balance_error(
                        species, paste0(
                            "the specified erepro cannot supply the ",
                            "target recruitment."
                        )
                    )
                }
                if (required_rdi == 0 && target == 0) {
                    params <- set_reproduction_parameter(
                        params, "erepro", i,
                        value,
                        record = length(controls) > 0L
                    )
                    if (law == "SheperdRDD" && !is.null(sheperd_c)) {
                        params <- set_reproduction_parameter(
                            params, "sheperd_c", i, shape[[i]],
                            record = TRUE
                        )
                    }
                    next
                }
                if (target == 0) {
                    reproduction_balance_error(
                        species, paste0(
                            "positive egg production needs a finite ",
                            "positive target."
                        )
                    )
                }
                density <- switch(law,
                    BevertonHoltRDD = if (required_rdi == target) {
                        Inf
                    } else {
                        target / (1 - target / required_rdi)
                    },
                    RickerRDD = {
                        (log(required_rdi) - log(target)) / required_rdi
                    },
                    SheperdRDD = {
                        ratio <- log(required_rdi) - log(target)
                        excess <- if (ratio > log(2)) {
                            ratio + log1p(-exp(-ratio))
                        } else {
                            log(expm1(ratio))
                        }
                        exp(excess / shape[[i]] - log(required_rdi))
                    }
                )
            } else {
                if (selector == "reproduction_level") {
                    if (target == 0 && value > 0) {
                        reproduction_balance_error(
                            species, paste0(
                                "a positive reproduction level requires ",
                                "positive recruitment."
                            )
                        )
                    }
                    if (target == 0 && value == 0) {
                        if (coefficient[[i]] > 0) {
                            params <- set_reproduction_parameter(
                                params, "erepro", i, 0
                            )
                        }
                        if (law == "SheperdRDD" && !is.null(sheperd_c)) {
                            params <- set_reproduction_parameter(
                                params,
                                "sheperd_c", i, shape[[i]],
                                record = TRUE
                            )
                        }
                        next
                    }
                    if (law == "BevertonHoltRDD") {
                        if (value >= 1) {
                            reproduction_balance_error(
                                species, paste0(
                                    "Beverton-Holt reproduction_level must be ",
                                    "below 1."
                                )
                            )
                        }
                        density <- if (value == 0) Inf else target / value
                    } else if (law == "RickerRDD") {
                        if (value > 1) {
                            reproduction_balance_error(
                                species, paste0(
                                    "Ricker reproduction_level must not ",
                                    "exceed 1."
                                )
                            )
                        }
                        density <- value / (exp(1) * target)
                    } else {
                        if (shape[[i]] < 1) {
                            reproduction_balance_error(
                                species, paste0(
                                    "Shepherd c < 1 has no informative ",
                                    "reproduction-level control."
                                )
                            )
                        }
                        if (value > 1 || (shape[[i]] == 1 && value == 1)) {
                            reproduction_balance_error(
                                species, paste0(
                                    "the specified Shepherd ",
                                    "reproduction_level is unattainable."
                                )
                            )
                        }
                        factor <- if (shape[[i]] == 1) {
                            1
                        } else {
                            shepherd_peak_factor(shape[[i]])
                        }
                        density <- factor * value / target
                    }
                } else {
                    density <- value
                }
                if (law == "BevertonHoltRDD") {
                    if (density <= target || density <= 0) {
                        reproduction_balance_error(
                            species, paste0(
                                "R_max must be strictly greater than the ",
                                "target recruitment."
                            )
                        )
                    }
                    required_rdi <- target / (1 - target / density)
                } else {
                    if (!is.finite(density) || density < 0) {
                        reproduction_balance_error(
                            species, paste0(
                                "the density-dependence parameter must be ",
                                "finite and nonnegative."
                            )
                        )
                    }
                    required_rdi <- if (law == "RickerRDD") {
                        ricker_inverse(target, density, limb, species)
                    } else {
                        shepherd_inverse(
                            target,
                            density, shape[[i]], limb, species
                        )
                    }
                }
            }
        }
        if (!is.null(density) && law != "BevertonHoltRDD" &&
                (!is.finite(density) || density < 0)) {
            reproduction_balance_error(
                species, paste0(
                    "the required density-dependence ",
                    "parameter is not finite and nonnegative."
                )
            )
        }
        if (!is.finite(required_rdi) || required_rdi < 0) {
            reproduction_balance_error(
                species, paste0(
                    "the requested curve has no finite ",
                    "nonnegative egg-production solution."
                )
            )
        }
        if (coefficient[[i]] == 0) {
            if (required_rdi > 0) {
                reproduction_balance_error(
                    species, paste0(
                        "there is no reproductive energy ",
                        "available to supply the target."
                    )
                )
            }
            efficiency <- sp$erepro[[i]]
        } else {
            efficiency <- required_rdi / coefficient[[i]]
        }
        if (!is.finite(efficiency)) {
            reproduction_balance_error(
                species, paste0(
                    "the required reproductive efficiency is ",
                    "not finite."
                )
            )
        }
        if (selector == "erepro") efficiency <- value
        params <- set_reproduction_parameter(
            params, "erepro", i, efficiency,
            record = length(controls) > 0L && selector == "erepro"
        )
        if (!is.null(density)) {
            params <- set_reproduction_parameter(
                params, density_field, i, density,
                record = length(controls) > 0L && selector == density_field
            )
        }
        if (law == "SheperdRDD" && !is.null(sheperd_c)) {
            params <- set_reproduction_parameter(params, "sheperd_c", i,
                shape[[i]],
                record = TRUE
            )
        }
    }
    params@time_modified <- lubridate::now()
    params
}
# nolint next: object_name_linter.
balance_BevertonHoltRDD <- function(params, rates, rdd_target, ...) {
    balance_reproduction_law(params, rates, rdd_target, "BevertonHoltRDD", ...)
}
# nolint next: object_name_linter.
balance_RickerRDD <- function(params, rates, rdd_target, ...) {
    balance_reproduction_law(params, rates, rdd_target, "RickerRDD", ...)
}
# nolint next: object_name_linter.
balance_SheperdRDD <- function(params, rates, rdd_target, ...) {
    balance_reproduction_law(params, rates, rdd_target, "SheperdRDD", ...)
}
# nolint next: object_name_linter.
balance_noRDD <- function(params, rates, rdd_target, ...) {
    balance_reproduction_law(params, rates, rdd_target, "noRDD", ...)
}
# nolint next: object_name_linter.
balance_constantRDD <- function(params, rates, rdd_target, ...) {
    balance_reproduction_law(params, rates, rdd_target, "constantRDD", ...)
}
builtin_reproduction_tuning <- function(params) {
    law <- reproduction_law_name(params)
    fields <- switch(law,
        BevertonHoltRDD = c("erepro", "R_max"),
        RickerRDD = c("erepro", "ricker_b"),
        SheperdRDD = c("erepro", "sheperd_b", "sheperd_c"),
        noRDD = "erepro",
        constantRDD = "constant_reproduction"
    )
    if (is.null(fields)) {
        return(NULL)
    }
    descriptor <- function(field) {
        list(
            path = c("species_params", field), lower = 0, upper = Inf,
            species = TRUE
        )
    }
    make_strategy <- function(preserve = NULL, change = fields) {
        if (law == "SheperdRDD") change <- union(change, "sheperd_c")
        kept <- list()
        if (!is.null(preserve)) {
            kept[[preserve]] <- if (preserve == "reproduction_level") {
                list(fun = paste0("reproduction_level_", law))
            } else {
                list(path = c("species_params", preserve))
            }
        }
        list(
            action = "balance",
            hold = list(states = character(), rates = "consumers/recruitment"),
            evolve = character(),
            adjust = stats::setNames(lapply(change, descriptor), change),
            preserve = kept, balance_fun = paste0("balance_", law),
            target_fun = NULL, solvers = c("project", "newton", "balance"),
            partial_species = TRUE,
            settings = list(
                preserve = if (is.null(preserve)) "none" else preserve,
                limb = "ascending"
            )
        )
    }
    strategies <- list(set_parameters = make_strategy())
    if (law %in% c("noRDD", "constantRDD")) {
        return(list(default = "set_parameters", strategies = strategies))
    }
    density <- fields[[2L]]
    strategies$preserve_erepro <- make_strategy("erepro", density)
    strategies$preserve_reproduction_level <-
        make_strategy("reproduction_level", c("erepro", density))
    strategies[[paste0("preserve_", density)]] <-
        make_strategy(density, "erepro")
    default <- if (law == "SheperdRDD" &&
                       any(params@species_params$sheperd_c < 1)) {
        "preserve_sheperd_b"
    } else {
        "preserve_reproduction_level"
    }
    list(default = default, strategies = strategies)
}

#' Reuse the built-in reproduction capability definition
#'
#' Returns the data-only balancing and reproduction-level capabilities for the
#' built-in law named in `params$rates_funcs$RDD`. This helper does not register
#' support or change the model. An extension author can explicitly adopt the
#' definition with [setProcess()] after verifying that its effective reproduction
#' route still implements that law and that density-independent reproduction
#' remains proportional to `erepro`.
#'
#' The recognised laws are Beverton-Holt, Ricker, Shepherd, constant recruitment
#' and density-independent recruitment. Custom mortality or encounter providers
#' may leave these reproduction relationships unchanged. A custom RDI, RDD or
#' complete-rate wrapper needs its own analysis: merely naming a built-in RDD
#' does not establish that the effective route preserves those relationships.
#' Unknown law names are rejected. Transformation support must be declared
#' separately by the extension author.
#'
#' @param params A MizerParams object.
#' @return A named list containing stationary, equation_rates_fun, tuning and
#'   level_fun, suitable as the definition argument for the reproduction
#'   process. The strategies reuse mizer's law-specific implementations.
#' @seealso [setProcess()], [balanceRDD()]
#' @family extension tools
#' @export
mizerReproductionDefinition <- function(params) {
    params <- validParams(params)
    law <- reproduction_law_name(params)
    if (length(law) != 1L || is.na(law) ||
            !law %in% c("BevertonHoltRDD", "RickerRDD", "SheperdRDD",
                         "constantRDD", "noRDD")) {
        stop("No built-in reproduction definition for the stored RDD law.",
             call. = FALSE)
    }
    list(stationary = TRUE,
         equation_rates_fun = "builtin_reproduction_equation_rates",
         tuning = builtin_reproduction_tuning(params),
         level_fun = paste0("reproduction_level_", law))
}
