# Measurement of biomass independently of the component's dynamics.

biomass_requested_filters <- function(dots, use_cutoff) {
    allowed <- c("min_w", "max_w", "min_l", "max_l")
    if (length(dots) && (!capability_names_valid(dots) ||
                             any(!names(dots) %in% allowed))) {
        stop("Biomass filters must be named `min_w`, `max_w`, `min_l` or ",
             "`max_l` arguments.")
    }
    if (!is.logical(use_cutoff) || length(use_cutoff) != 1L ||
            is.na(use_cutoff)) {
        stop("`use_cutoff` must be TRUE or FALSE.")
    }
    # NULL size limits have their ordinary meaning of no requested limit.
    dots <- dots[!vapply(dots, is.null, logical(1))]
    if (!length(dots)) dots <- list()
    if (!is.null(dots$min_l)) dots$min_w <- NULL
    if (!is.null(dots$max_l)) dots$max_w <- NULL
    if (use_cutoff) dots$use_cutoff <- TRUE
    dots
}

prepare_biomass <- function(params, include, filters) {
    if (!is.character(include) || !length(include) || anyNA(include) ||
            any(!nzchar(include))) {
        stop("`include` must be a non-empty character vector of biomass ",
             "selectors.")
    }
    components <- names(params$initial_n_other)
    owners <- paste0("component:", components)
    if (!length(components)) owners <- character()
    selectors <- list()
    add_selector <- function(selectors, owner, id = NULL) {
        c(selectors, list(list(owner = owner, id = id)))
    }
    for (selector in unique(include)) {
        if (selector %in% c("species", "resource")) {
            selectors <- add_selector(selectors, selector)
        } else if (selector == "components") {
            for (owner in owners) {
                selectors <- add_selector(selectors, owner)
            }
        } else if (selector %in% owners) {
            selectors <- add_selector(selectors, selector)
        } else {
            matches <- owners[startsWith(selector, paste0(owners, "/"))]
            if (!length(matches)) {
                stop("Unknown biomass selector '", selector, "'.")
            }
            # An exact component name takes precedence over a series prefix.
            owner <- matches[[which.max(nchar(matches))]]
            if (nchar(selector) == nchar(owner) + 1L) {
                stop("A component biomass series identifier must not be empty.")
            }
            selectors <- add_selector(selectors, owner, selector)
        }
    }
    selected_owners <- unique(vapply(selectors, `[[`, "", "owner"))
    plans <- setNames(vector("list", length(selected_owners)), selected_owners)
    for (owner in selected_owners) {
        if (owner %in% c("species", "resource")) next
        component <- substring(owner, nchar("component:") + 1L)
        spec <- params$capabilities$components[[component]]
        reason <- NULL
        callback <- NULL
        unsupported <- setdiff(names(filters), spec$biomass_filters)
        if (is.null(spec$biomass_fun)) {
            reason <- paste0("No biomass callback is registered for ", owner,
                             ".")
        } else if (length(unsupported)) {
            reason <- paste0(owner, " does not support biomass filters: ",
                             paste(unsupported, collapse = ", "), ".")
        } else {
            callback <- tryCatch(resolve_model_callback(spec$biomass_fun),
                                 error = function(e) {
                                     stop("Cannot resolve biomass callback ",
                                          "for ", owner, ": ",
                                          conditionMessage(e), call. = FALSE)
                                 })
        }
        plans[[owner]] <- list(component = component, spec = spec,
                               callback = callback, reason = reason)
    }
    list(selectors = selectors, owners = selected_owners, plans = plans,
         filters = filters)
}

biomass_row <- function(id, value, label = id, owner = id, kind = "component",
                        colour = "grey", linetype = "solid", source = id,
                        available = TRUE, reason = NA_character_) {
    data.frame(id = id, label = label, owner = owner, kind = kind,
               colour = colour, linetype = linetype, source = source,
               value = value, available = available, reason = reason,
               stringsAsFactors = FALSE)
}

empty_biomass_rows <- function() {
    biomass_row(character(), numeric(), label = character(),
                owner = character(), kind = character(), colour = character(),
                linetype = character(), source = character(),
                available = logical(), reason = character())
}

unavailable_biomass <- function(owner, label, reason, kind = "component",
                                colour = "grey", linetype = "solid") {
    biomass_row(owner, NA_real_, label, owner, kind, colour, linetype,
                available = FALSE, reason = reason)
}

resource_biomass <- function(params, n_pp, filters) {
    style <- function(field, default) {
        value <- unname(params[[field]]["Resource"])
        if (length(value) && !is.na(value)) as.character(value) else default
    }
    colour <- style("linecolour", "green")
    linetype <- style("linetype", "solid")
    unavailable <- function(reason) {
        unavailable_biomass("resource", "Resource", reason, kind = "resource",
                            colour = colour, linetype = linetype)
    }
    if (isTRUE(filters$use_cutoff)) {
        return(unavailable("The resource has no species biomass cutoff."))
    }
    limits <- filters[c("min_w", "max_w", "min_l", "max_l")]
    limits <- limits[!vapply(limits, is.null, logical(1))]
    scalar <- vapply(limits, function(x) {
        is.numeric(x) && length(x) == 1L
    }, logical(1))
    if (any(!scalar)) {
        return(unavailable("Resource biomass requires scalar size filters."))
    }
    if (anyNA(unlist(limits, use.names = FALSE))) {
        stop("Resource biomass size limits must not be NA or NaN.")
    }
    min_w <- filters$min_w %||% min(params$w_full)
    max_w <- filters$max_w %||% max(params$w_full)
    if (!is.null(filters$min_l) || !is.null(filters$max_l)) {
        lengths <- unlist(filters[c("min_l", "max_l")], use.names = FALSE)
        if (any(lengths < 0)) {
            stop("Resource biomass length limits must be non-negative.")
        }
        mapping <- resource_length_params(params)
        valid_mapping <- vapply(mapping, function(value) {
            is.numeric(value) && length(value) == 1L && is.finite(value) &&
                value > 0
        }, logical(1))
        if (!all(valid_mapping)) {
            stop("Resource biomass requires positive finite weight-length ",
                 "parameters.")
        }
        if (!is.null(filters$min_l)) {
            min_w <- mapping$a * filters$min_l ^ mapping$b
        }
        if (!is.null(filters$max_l)) {
            max_w <- mapping$a * filters$max_l ^ mapping$b
        }
    }
    if (anyNA(c(min_w, max_w))) {
        stop("Resource biomass size limits could not be converted to weights.")
    }
    if (!is.numeric(n_pp) || length(n_pp) != length(params$w_full) ||
            any(!is.finite(n_pp)) || any(n_pp < 0)) {
        stop("Resource biomass requires a finite non-negative resource state ",
             "on the full size grid.")
    }
    mask <- params$w_full >= min_w & params$w_full <= max_w
    weighting <- bin_average_weight(params$w_full * mask, params)
    value <- as.numeric(n_pp %*% (weighting * params$dw_full))
    if (!is.finite(value) || value < 0) {
        stop("Resource biomass must be finite and non-negative.")
    }
    biomass_row("resource", value, "Resource", "resource", "resource",
                colour = colour, linetype = linetype)
}

component_biomass_rows <- function(value, owner, spec) {
    fail <- function(...) stop("Biomass callback for ", owner, ": ", ...)
    if (is.numeric(value) && is.null(dim(value))) {
        if (!capability_names_valid(value)) {
            fail("numeric results must have unique, non-empty series names.")
        }
        value <- data.frame(series = names(value), value = as.numeric(value),
                            stringsAsFactors = FALSE)
    }
    if (!is.data.frame(value) || !all(c("series", "value") %in% names(value)) ||
            !nrow(value)) {
        fail("return a named numeric vector or a data frame with `series` ",
             "and `value` columns and at least one series.")
    }
    if (!is.character(value$series) || anyNA(value$series) ||
            any(!nzchar(value$series)) || anyDuplicated(value$series)) {
        fail("series identifiers must be unique, non-empty strings.")
    }
    if (!is.numeric(value$value) || !is.null(dim(value$value))) {
        fail("`value` must be a numeric vector.")
    }
    optional <- c("label", "colour", "linetype", "source", "reason")
    for (field in intersect(optional, names(value))) {
        if (!is.character(value[[field]])) {
            fail(paste0("`", field, "` must be character."))
        }
    }
    available <- value$available %||% rep(TRUE, nrow(value))
    if (!is.logical(available) || anyNA(available)) {
        fail("`available` must contain TRUE or FALSE.")
    }
    if (any(available & (!is.finite(value$value) | value$value < 0))) {
        fail("available biomass values must be finite and non-negative.")
    }
    reason <- value$reason %||% rep(NA_character_, nrow(value))
    if (any(!available & (!is.na(value$value) | is.na(reason) |
                              !nzchar(reason)))) {
        fail("unavailable biomass needs an NA value and a non-empty reason.")
    }
    id <- paste0(owner, "/", value$series)
    fill <- function(field, default) {
        column <- value[[field]] %||% rep(NA_character_, nrow(value))
        missing <- is.na(column) | !nzchar(column)
        column[missing] <- rep(default, length.out = nrow(value))[missing]
        column
    }
    reason[available] <- NA_character_
    biomass_row(id, value$value, fill("label", value$series), owner,
                colour = fill("colour", spec$colour %||% "grey"),
                linetype = fill("linetype", spec$linetype %||% "solid"),
                source = fill("source", id), available = available,
                reason = reason)
}

evaluate_biomass <- function(params, plan, n, n_pp, n_other, t,
                             species_biomass = NULL) {
    rows <- setNames(vector("list", length(plan$owners)), plan$owners)
    for (owner in plan$owners) {
        if (owner == "species") {
            if (is.null(species_biomass)) {
                filters <- plan$filters
                cutoff <- isTRUE(filters$use_cutoff)
                filters$use_cutoff <- NULL
                if (cutoff && "biomass_cutoff" %in%
                        names(params$species_params)) {
                    filters <- list(min_w = cutoff_min_w(params))
                }
                arguments <- c(list(object = params, weighting = params$w,
                                    n = n), filters)
                species_biomass <- do.call(sizeIntegral, arguments)
            }
            if (any(!is.finite(species_biomass)) || any(species_biomass < 0)) {
                stop("Species biomass must be finite and non-negative.")
            }
            species <- names(species_biomass)
            source <- paste0("species:", species)
            colour <- unname(params$linecolour[species])
            linetype <- unname(params$linetype[species])
            rows[[owner]] <- biomass_row(species, as.numeric(species_biomass),
                                         species,
                                         owner = source,
                                         kind = "species", colour = colour,
                                         linetype = linetype,
                                         source = source)
        } else if (owner == "resource") {
            rows[[owner]] <- resource_biomass(params, n_pp, plan$filters)
        } else {
            component_plan <- plan$plans[[owner]]
            if (!is.null(component_plan$reason)) {
                spec <- component_plan$spec
                row <- unavailable_biomass(owner, component_plan$component,
                                           component_plan$reason,
                                           colour = spec$colour %||% "grey",
                                           linetype = spec$linetype %||%
                                               "solid")
                rows[[owner]] <- row
            } else {
                args <- list(n = n, n_pp = n_pp, n_other = n_other, t = t,
                             component = component_plan$component,
                             filters = plan$filters)
                callback <- component_plan$callback
                context <- paste("Biomass measurement for", owner)
                value <- invoke_tuning_callback(callback, params, args, context)
                rows[[owner]] <- component_biomass_rows(value, owner,
                                                        component_plan$spec)
            }
        }
    }
    selected <- empty_biomass_rows()
    for (selector in plan$selectors) {
        part <- rows[[selector$owner]]
        if (!is.null(selector$id)) {
            if (selector$id %in% part$id) {
                part <- part[part$id == selector$id, , drop = FALSE]
            } else if (nrow(part) == 1L && !part$available &&
                           part$id == selector$owner) {
                part$id <- part$source <- selector$id
            } else {
                stop("Unknown biomass series '", selector$id, "'.")
            }
        }
        overlap <- match(part$id, selected$id)
        if (any(!is.na(overlap) &
                    part$owner != selected$owner[overlap])) {
            stop("Biomass series identifiers conflict across model owners.")
        }
        selected <- rbind(selected, part[!part$id %in% selected$id, ,
                                         drop = FALSE])
    }
    if (anyDuplicated(selected$source)) {
        stop("Selected biomass series describe the same physical source: ",
             paste(unique(selected$source[duplicated(selected$source)]),
                   collapse = ", "), ".")
    }
    rownames(selected) <- NULL
    selected
}

assemble_biomass <- function(rows, times, params, simulation = FALSE) {
    fields <- c("id", "label", "owner", "kind", "colour", "linetype", "source")
    metadata <- rows[[1L]][fields]
    ids <- metadata$id
    values <- matrix(NA_real_, nrow = length(rows), ncol = length(ids),
                     dimnames = list(time = as.character(times), sp = ids))
    coverage <- vector("list", length(rows))
    for (i in seq_along(rows)) {
        current <- rows[[i]]
        if (!setequal(current$id, ids)) {
            stop("Component biomass series identifiers changed between saved ",
                 "times. Return the same series, using unavailable rows when ",
                 "a measurement is missing.")
        }
        current <- current[match(ids, current$id), , drop = FALSE]
        rownames(current) <- NULL
        if (!identical(current[fields], metadata)) {
            stop("Biomass series metadata changed between saved times.")
        }
        values[i, ] <- current$value
        coverage[[i]] <- data.frame(time = rep(times[[i]], length(ids)),
                                    id = ids, available = current$available,
                                    reason = current$reason,
                                    stringsAsFactors = FALSE)
    }
    coverage <- do.call(rbind, coverage)
    if (simulation) {
        result <- ArrayTimeBySpecies(values, value_name = "Biomass",
                                     units = "g", params = params)
    } else {
        result <- setNames(as.numeric(values[1L, ]), ids)
        result <- structure(result, class = c("MizerBiomass", "numeric"),
                            params = params, value_name = "Biomass",
                            units = "g", type = "value")
    }
    attr(result, "series") <- metadata
    attr(result, "coverage") <- coverage
    attr(result, "unavailable") <- coverage[!coverage$available, , drop = FALSE]
    result
}
