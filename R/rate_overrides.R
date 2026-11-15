# Species names, rather than row positions, survive reordering and subsetting.
# An absent attribute on a commented array is the legacy whole-array override.
rate_override_species <- function(params, var) {
    x <- slot(params, var)
    if (is.null(comment(x))) return(character())
    species <- params@species_params$species
    selected <- attr(x, "mizer_override_species", exact = TRUE)
    if (is.null(selected)) species else intersect(species, selected)
}

rate_species_dim <- function(var) {
    if (var %in% c("selectivity", "catchability")) 2L else 1L
}

check_override_species <- function(params, species) {
    if (is.null(species)) return(params@species_params$species)
    if (!is.character(species) || !length(species) || anyNA(species) ||
            anyDuplicated(species) ||
            !all(species %in% params@species_params$species)) {
        stop("`species` must contain unique species names from the model.")
    }
    species
}

# Index all sizes (and gears) of the selected species, without dropping axes.
rate_species_index <- function(x, var, species) {
    index <- lapply(dim(x), seq_len)
    axis <- rate_species_dim(var)
    index[[axis]] <- match(species, dimnames(x)[[axis]])
    index
}

rate_slice <- function(x, var, species) {
    do.call(`[`, c(list(x), rate_species_index(x, var, species),
                   list(drop = FALSE)))
}

replace_rate_slice <- function(x, var, species, value) {
    do.call(`[<-`, c(list(x), rate_species_index(x, var, species),
                     list(value = value)))
}

mark_rate_override <- function(x, species, all_species) {
    attr(x, "mizer_override_species") <-
        if (length(species) && !setequal(species, all_species)) species
    if (!length(species)) comment(x) <- NULL
    else if (is.null(comment(x))) comment(x) <- "set manually"
    x
}

reset_rate_override <- function(params, var, species) {
    selected <- check_override_species(params, species)
    remaining <- setdiff(rate_override_species(params, var), selected)
    slot(params, var) <- mark_rate_override(slot(params, var), remaining,
                                            params@species_params$species)
    params
}

# Expand a species-only value to the existing array shape. Full arrays remain
# accepted, but only the explicitly selected species are taken from them.
prepare_rate_value <- function(params, var, value, species) {
    check_override_species(params, species)
    if (is.null(value)) return(NULL)
    if (is.null(species)) {
        attr(value, "mizer_override_species") <- NULL
        return(value)
    }
    original <- switch(var,
                       pred_kernel = pred_kernel(params),
                       psi = repro_prop(params),
                       slot(params, var))
    original <- unclass_rate(original)
    expected <- rate_slice(original, var, species)
    if (is.null(dim(value))) {
        if (length(species) != 1L || length(value) != length(expected)) {
            stop("Supply an array for the selected species, ",
                 "or a full rate array.")
        }
    } else if (identical(dim(value), dim(original))) {
        dn <- dimnames(value)[[rate_species_dim(var)]]
        original_names <- dimnames(original)[[rate_species_dim(var)]]
        if (!is.null(dn) && !identical(dn, original_names)) {
            stop("The species dimnames must match the model's species order.")
        }
        value <- rate_slice(value, var, species)
    } else {
        if (!identical(dim(value), dim(expected))) {
            stop("The array dimensions must match the selected species ",
                 "and size grid.")
        }
        dn <- dimnames(value)[[rate_species_dim(var)]]
        if (!is.null(dn) && !identical(dn, species)) {
            stop("The species dimnames must match `species` ",
                 "in the supplied order.")
        }
    }
    if (!is.null(dim(value))) {
        axis <- rate_species_dim(var)
        for (d in setdiff(seq_along(dim(expected)), axis)) {
            dn <- dimnames(value)[[d]]
            if (!is.null(dn) && !identical(dn, dimnames(expected)[[d]])) {
                stop("The gear and size dimnames must match ",
                     "the existing rate array.")
            }
        }
    }
    if (!is.numeric(value) || any(!is.finite(value))) {
        stop("Manual rate values must be finite numbers.")
    }
    original <- replace_rate_slice(original, var, species, value)
    mark_rate_override(original,
                       union(rate_override_species(params, var), species),
                       params@species_params$species)
}

# Preserve only the protected slices of a newly calculated array.
preserve_rate_overrides <- function(params, var, value) {
    selected <- rate_override_species(params, var)
    if (!length(selected)) return(value)
    if (!identical(dim(value), dim(slot(params, var)))) {
        stop("Reset the manual override of `", var,
             "` before changing its array dimensions.")
    }
    dimnames(value) <- dimnames(slot(params, var))
    old <- rate_slice(slot(params, var), var, selected)
    proposed <- rate_slice(value, var, selected)
    if (different(old, proposed)) {
        entry <- frozen_rate_params()[[var]]
        quantity <- entry$quantity
        reset_call <- entry$reset_call
        if (length(selected) < nrow(params@species_params)) {
            quantity <- paste0(quantity, " for ", toString(selected))
            selection <- paste(encodeString(selected, quote = '"'),
                               collapse = ", ")
            replacement <- paste0("reset = TRUE, species = c(", selection, ")")
            reset_call <- sub("reset = TRUE", replacement, reset_call,
                              fixed = TRUE)
        }
        signal_not_recalculated(var, quantity, reset_call,
                                entry$derived_from %||% "species parameters")
    }
    replace_rate_slice(value, var, selected, old)
}

#' Inspect manual rate overrides
#'
#' Reports the species whose rate arrays have been set manually. The listed
#' parameters are retained inputs for recalculation after a reset: they no
#' longer determine this rate for this species, but may still contribute to
#' other rates or derivations. Resource overrides have `NA` in `species`.
#'
#' @param params A MizerParams object.
#' @return A data frame with one row per overridden rate and species, columns
#'   `rate`, `species`, `parameters` (a list of parameter names), and `reset`
#'   (a call restoring parameter-driven calculation for that species).
#' @export
#' @family functions for setting parameters
#' @examples
#' params <- NS_params
#' metab(params, species = "Cod") <- 0.5 * w(params)^0.8
#' rate_overrides(params)
#' params <- setMetabolicRate(params, species = "Cod", reset = TRUE)
rate_overrides <- function(params) {
    rows <- list()
    for (var in names(frozen_rate_params())) {
        entry <- frozen_rate_params()[[var]]
        if (is.null(comment(slot(params, var)))) next
        resource <- identical(entry$derived_from, "resource parameters")
        species <- if (resource) NA_character_ else
            rate_override_species(params, var)
        columns <- if (resource) names(params@resource_params) else
            if (identical(entry$derived_from, "gear parameters"))
                names(params@gear_params) else names(params@species_params)
        for (sp in species) {
            reset <- entry$reset_call
            if (!is.na(sp)) {
                selection <- paste0("reset = TRUE, species = ",
                                    encodeString(sp, quote = '"'))
                reset <- sub("reset = TRUE", selection, reset, fixed = TRUE)
            }
            rate <- switch(var, mu_b = "ext_mort", psi = "repro_prop",
                           cc_pp = "resource_capacity", rr_pp = "resource_rate",
                           var)
            parameters <- intersect(rate_input_parameters(params, var, sp),
                                    columns)
            row <- data.frame(rate = rate, species = sp,
                              parameters = I(list(parameters)), reset = reset)
            rows[[length(rows) + 1L]] <- row
        }
    }
    if (length(rows)) return(do.call(rbind, rows))
    data.frame(rate = character(), species = character(),
               parameters = I(list()), reset = character())
}

# Kernel and selectivity functions can have custom parameter columns.
rate_input_parameters <- function(params, var, species) {
    entry <- frozen_rate_params()[[var]]
    if (var == "pred_kernel") {
        type <- params@species_params[species, "pred_kernel_type", drop = TRUE]
        fun <- paste0(type, "_pred_kernel")
        args <- tryCatch(names(formals(fun)), error = function(e) character())
        return(c("pred_kernel_type", setdiff(args, c("ppmr", "..."))))
    }
    if (var == "selectivity") {
        gp <- params@gear_params
        functions <- unique(gp$sel_func[gp$species == species])
        args <- unlist(lapply(functions, function(fun) {
            tryCatch(names(formals(fun)), error = function(e) character())
        }))
        args <- setdiff(unique(args), c("w", "species_params", "..."))
        return(c("sel_func", args))
    }
    entry$params
}

print_rate_overrides <- function(overrides) {
    if (!nrow(overrides)) return(invisible(NULL))
    cat("Manual rate overrides (parameters retained for reset):\n")
    for (i in seq_len(nrow(overrides))) {
        species <- overrides$species[i]
        label <- if (!is.na(species)) paste0(" for ", species) else ""
        cat("  ", overrides$rate[i], label,
            ": ", toString(overrides$parameters[[i]]), "\n", sep = "")
    }
    cat("These parameters may still contribute to other rates.\n",
        "Use rate_overrides(params)$reset for calls restoring calculation.\n",
        sep = "")
    invisible(NULL)
}
