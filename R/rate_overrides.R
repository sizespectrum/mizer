# The rate arrays that can be overridden for selected species only. The
# predation kernel is not among them: it can only be set as a whole, because a
# kernel overridden for some species would need explicit sums over the kernel
# for those species at every time step instead of the Fourier transform. The
# resource arrays have no species dimension.
species_override_rates <- c("intake_max", "search_vol", "metab", "mu_b",
                            "ext_encounter", "ext_diffusion", "maturity",
                            "psi", "selectivity", "catchability")

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

# Copy the override marker, that is the freeze comment together with the
# overridden species, from one array to another.
copy_override_marker <- function(to, from) {
    comment(to) <- comment(from)
    attr(to, "mizer_override_species") <-
        attr(from, "mizer_override_species", exact = TRUE)
    to
}

reset_rate_override <- function(params, var, species) {
    selected <- check_override_species(params, species)
    remaining <- setdiff(rate_override_species(params, var), selected)
    slot(params, var) <- mark_rate_override(slot(params, var), remaining,
                                            params@species_params$species)
    params
}

# The call that makes the rate `var` follow its parameters again for
# `species`. Without `species`, or for a rate that can only be overridden as a
# whole, it is the call that resets the whole array.
override_reset_call <- function(params, var, species = NULL) {
    call <- frozen_rate_params()[[var]]$reset_call
    if (is.null(species) || !var %in% species_override_rates ||
            setequal(species, params@species_params$species)) {
        return(call)
    }
    selection <- encodeString(species, quote = '"')
    if (length(selection) > 1) {
        selection <- paste0("c(", paste(selection, collapse = ", "), ")")
    }
    sub("reset = TRUE", paste0("reset = TRUE, species = ", selection), call,
        fixed = TRUE)
}

# Non-species dimnames of a supplied value must match those of the rate array.
check_other_dimnames <- function(dn, expected) {
    for (d in seq_along(expected)) {
        if (!is.null(dn[[d]]) && !identical(dn[[d]], expected[[d]])) {
            stop("The gear and size dimnames must match ",
                 "the existing rate array.")
        }
    }
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
    original <- unclass_rate(if (var == "psi") repro_prop(params) else
        slot(params, var))
    expected <- rate_slice(original, var, species)
    axis <- rate_species_dim(var)
    if (identical(dim(value), dim(original))) {
        dn <- dimnames(value)[[axis]]
        if (!is.null(dn) && !identical(dn, dimnames(original)[[axis]])) {
            stop("The species dimnames must match the model's species order.")
        }
        check_other_dimnames(dimnames(value)[-axis], dimnames(original)[-axis])
        value <- rate_slice(value, var, species)
    } else if (length(species) == 1L &&
               length(value) == length(expected) &&
               length(dim(value)) < length(dim(expected))) {
        # The slice of a single species without its species dimension: a
        # vector, or for the selectivity a gear x size matrix.
        if (!is.null(dim(value))) {
            if (!identical(dim(value), dim(expected)[-axis])) {
                stop("The array dimensions must match the selected species ",
                     "and size grid.")
            }
            check_other_dimnames(dimnames(value), dimnames(expected)[-axis])
        }
        value <- array(as.vector(value), dim = dim(expected),
                       dimnames = dimnames(expected))
    } else {
        if (is.null(dim(value))) {
            stop("Supply an array for the selected species, ",
                 "or a full rate array.")
        }
        if (!identical(dim(value), dim(expected))) {
            stop("The array dimensions must match the selected species ",
                 "and size grid.")
        }
        dn <- dimnames(value)[[axis]]
        if (!is.null(dn) && !identical(dn, species)) {
            stop("The species dimnames must match `species` ",
                 "in the supplied order.")
        }
        check_other_dimnames(dimnames(value)[-axis], dimnames(expected)[-axis])
    }
    # Infinite values are allowed as they are for whole arrays: an infinite
    # maximum intake rate, for example, means that all food is consumed. The
    # setter's own checks apply to them in the same way.
    if (!is.numeric(value) || anyNA(value)) {
        stop("Manual rate values must be numbers, not NA.")
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
    if (different(old, rate_slice(value, var, selected))) {
        entry <- frozen_rate_params()[[var]]
        quantity <- entry$quantity
        if (length(selected) < nrow(params@species_params)) {
            quantity <- paste0(quantity, " for ", toString(selected))
        }
        signal_not_recalculated(var, quantity,
                                override_reset_call(params, var, selected),
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
#'   (a call restoring parameter-driven calculation for that species). The
#'   predation kernel can only be set as a whole, so its `reset` call
#'   restores the calculation for all species.
#' @export
#' @family functions for setting parameters
#' @examples
#' params <- NS_params
#' metab(params, species = "Cod") <- 0.5 * w(params)^0.8
#' rate_overrides(params)
#' params <- setMetabolicRate(params, species = "Cod", reset = TRUE)
rate_overrides <- function(params) {
    entries <- frozen_rate_params()
    rows <- list()
    for (var in names(entries)) {
        entry <- entries[[var]]
        if (is.null(comment(slot(params, var)))) next
        resource <- identical(entry$derived_from, "resource parameters")
        species <- if (resource) NA_character_ else
            rate_override_species(params, var)
        columns <- if (resource) names(params@resource_params) else
            if (identical(entry$derived_from, "gear parameters"))
                names(params@gear_params) else names(params@species_params)
        rate <- switch(var, mu_b = "ext_mort", psi = "repro_prop",
                       cc_pp = "resource_capacity", rr_pp = "resource_rate",
                       var)
        for (sp in species) {
            reset <- override_reset_call(params, var, if (!is.na(sp)) sp)
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
