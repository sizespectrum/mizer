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
    original <- if (var == "psi") repro_prop(params) else slot(params, var)
    original <- unclass_rate(original)
    expected <- rate_slice(original, var, species)
    axis <- rate_species_dim(var)
    if (identical(dim(value), dim(original))) {
        gn <- dimnames(value)[[1]]
        if (axis == 2L && !is.null(gn) && !anyDuplicated(gn) &&
                setequal(gn, dimnames(original)[[1]])) {
            # The same gears in another order are matched by name
            value <- align_gears(value, dimnames(original)[[1]])
        }
        dn <- dimnames(value)[[axis]]
        if (!is.null(dn) && !identical(dn, dimnames(original)[[axis]])) {
            stop("The species dimnames must match the model's species order.")
        }
        check_other_dimnames(dimnames(value)[-axis], dimnames(original)[-axis])
        # Dimnames are optional. Where they are missing, the entries are in
        # the order of the rate array, which supplies the names.
        dn <- dimnames(value)
        if (is.null(dn)) dn <- vector("list", length(dim(value)))
        missing_dn <- vapply(dn, is.null, logical(1))
        dn[missing_dn] <- dimnames(original)[missing_dn]
        dimnames(value) <- dn
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

# Index the given gears and species (all species if `NULL`) of a fishing
# array by name.
fishing_index <- function(x, gears, species = NULL) {
    index <- lapply(dim(x), seq_len)
    index[[1]] <- match(gears, dimnames(x)[[1]])
    if (!is.null(species)) index[[2]] <- match(species, dimnames(x)[[2]])
    index
}

fishing_slice <- function(x, gears, species = NULL) {
    do.call(`[`, c(list(x), fishing_index(x, gears, species),
                   list(drop = FALSE)))
}

replace_fishing_slice <- function(x, gears, species, value) {
    do.call(`[<-`, c(list(x), fishing_index(x, gears, species),
                     list(value = value)))
}

# Bring a fishing array calculated for the gears in the gear parameters onto
# the gear axis `gears`, matching by name. A gear without parameters gets
# zeros.
align_gears <- function(x, gears) {
    if (identical(dimnames(x)[[1]], gears)) return(x)
    dn <- dimnames(x)
    dn[[1]] <- gears
    out <- array(0, dim = c(length(gears), dim(x)[-1]), dimnames = dn)
    common <- intersect(gears, dimnames(x)[[1]])
    replace_fishing_slice(out, common, NULL, fishing_slice(x, common))
}

# The gear axis that both fishing arrays share. A stored fishing array that
# has been set by hand as a whole keeps the gears it was set with, as it did
# before arrays could be set for selected species, and the other array must
# follow it. `vars` names the stored arrays that are not being replaced and so
# can impose their gears. Without such an array the gears are `gear_names`,
# those in the gear parameters. Returns the gears, together with the name of
# the array that imposed them, if any.
fishing_gears <- function(params, gear_names,
                          vars = c("selectivity", "catchability")) {
    for (var in vars) {
        if (length(rate_override_species(params, var)) ==
                nrow(params@species_params)) {
            return(list(gears = dimnames(slot(params, var))[[1]],
                        set_by = var))
        }
    }
    list(gears = gear_names, set_by = NULL)
}

# Check that a fishing array supplied by the user has the gears `gears` and
# put them into that order, matching by name. An array without gear names
# must already be in that order. When a frozen array imposes the gears, an
# array with other gears is rejected here, because the two arrays could not
# fit together.
fit_supplied_gears <- function(x, var, gears, set_by) {
    gn <- dimnames(x)[[1]]
    fits <- if (is.null(gn)) {
        dim(x)[[1]] == length(gears)
    } else {
        !anyDuplicated(gn) && setequal(gn, gears)
    }
    if (!fits && !is.null(set_by) && set_by != var) {
        stop("The gears of the ", var, " array must be ", toString(gears),
             ", because the ", set_by, " has been set manually as a whole ",
             "with these gears. Supply both arrays together, or call ",
             "`setFishing(params, reset = TRUE)` to calculate them from the ",
             "gear parameters again.")
    }
    if (fits && !is.null(gn) && !identical(gn, gears)) {
        x <- copy_override_marker(align_gears(x, gears), x)
    }
    x
}

signal_kept_override <- function(params, var, selected) {
    entry <- frozen_rate_params()[[var]]
    quantity <- entry$quantity
    if (length(selected) < nrow(params@species_params)) {
        quantity <- paste0(quantity, " for ", toString(selected))
    }
    signal_not_recalculated(var, quantity,
                            override_reset_call(params, var, selected),
                            entry$derived_from %||% "species parameters")
}

# Preserve only the protected slices of a newly calculated array. An array
# set by hand as a whole is kept as it is. The fishing arrays are given the
# gear axis `gears` (see `fishing_gears()`), and their protected slices are
# matched by gear name, because the gears in the gear parameters may have been
# reordered, added or removed since the array was set. The protected values
# for a gear that is no longer there are dropped with it.
preserve_rate_overrides <- function(params, var, value, gears = NULL) {
    selected <- rate_override_species(params, var)
    old <- slot(params, var)
    if (length(selected) &&
            length(selected) == nrow(params@species_params)) {
        if (different(value, old)) signal_kept_override(params, var, selected)
        return(old)
    }
    fishing <- var %in% c("selectivity", "catchability")
    if (fishing) {
        value <- align_gears(value, gears %||% dimnames(value)[[1]])
    }
    if (!length(selected)) return(value)
    if (fishing) {
        common <- intersect(dimnames(value)[[1]], dimnames(old)[[1]])
        dropped <- setdiff(dimnames(old)[[1]], common)
        kept <- fishing_slice(old, common, selected)
        lost <- fishing_slice(old, dropped, selected)
        if (different(kept, fishing_slice(value, common, selected)) ||
                any(lost != 0)) {
            signal_kept_override(params, var, selected)
        }
        return(replace_fishing_slice(value, common, selected, kept))
    }
    if (!identical(dim(value), dim(old))) {
        stop("Reset the manual override of `", var,
             "` before changing its array dimensions.")
    }
    # The calculated species x size arrays have the layout of the stored ones
    dimnames(value) <- dimnames(old)
    kept <- rate_slice(old, var, selected)
    if (different(kept, rate_slice(value, var, selected))) {
        signal_kept_override(params, var, selected)
    }
    replace_rate_slice(value, var, selected, kept)
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
