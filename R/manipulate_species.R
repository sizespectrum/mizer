#' Add new species
#'
#' Takes a \linkS4class{MizerParams} object and adds additional species with
#' given parameters to the ecosystem. It sets the initial values for these new
#' species to their steady-state solution in the given initial state of the
#' existing ecosystem. This will be close to the true steady state if the
#' abundances of the new species are sufficiently low. Hence the abundances of
#' the new species are set so that they are at most 1/100th of the resource
#' power law. Their reproductive efficiencies are set so as to keep them at that
#' low level.
#'
#' @param params A mizer params object for the original system.
#' @param species_params Data frame with the species parameters of the new
#'   species we want to add to the system.
#' @param interaction Interaction matrix. A square matrix giving either the
#'   interaction coefficients between all species or only those between the new
#'   species. In the latter case all interaction between an old and a new
#'   species are set to 1. If this argument is missing, all interactions
#'   involving a new species are set to 1.
#' @param gear_params Data frame with the gear parameters for the new
#'   species. If not provided then the new species will not be fished.
#' @param initial_effort A named vector with the effort for any new fishing gear
#'   introduced in `gear_params`. Not needed if the added species are only
#'   fished by already existing gear. Should not include effort values
#'   for existing gear. New gear for which no effort is set via this
#'   vector will have an initial effort of 0.
#' @param steady If `TRUE` (default), runs [steadySingleSpecies()] to
#'   initialise the new species at their single-species steady state and
#'   retuning their reproductive efficiencies. Set to `FALSE` when the caller
#'   (e.g. an extension package using `NextMethod()`) needs to make further
#'   changes to the params object before that steady-state calculation can be
#'   run successfully.
#' @param info_level Controls the amount of information messages that are shown
#'   when the function sets default values for parameters. Higher levels lead
#'   to more messages. Set to 0 to suppress all such messages.
#' @param ... Currently unused.
#'
#' @return An object of type \linkS4class{MizerParams}
#'
#' @details The resulting MizerParams object will use the same size grid where
#'   possible, but if one of the new species needs a larger range of w (either
#'   because a new species has an egg size smaller than those of existing
#'   species or a maximum size larger than those of existing species) then the
#'   grid will be expanded and all arrays will be enlarged accordingly.
#'
#'   If any of the rate arrays of the existing species had been set by the user
#'   to values other than those calculated as default from the species
#'   parameters, then these will be preserved. Only the rates for the new
#'   species will be calculated from their species parameters.
#'
#'   After adding the new species, the background species are not retuned and
#'   the system is not run to steady state. This could be done with
#'   [tuneSteadyState()].
#'   With Beverton-Holt reproduction, the new species receive a reproduction
#'   level of 1/4. Other reproduction laws retain their declared balancing
#'   strategy. Only the new species' recruitment is targeted; a coupled strategy
#'   that cannot honour a partial request must reject the operation.
#'
#' @examples
#' params <- newTraitParams()
#' species_params <- data.frame(
#'     species = "Mullet",
#'     w_max = 173,
#'     w_mat = 15,
#'     beta = 283,
#'     sigma = 1.8,
#'     h = 30,
#'     a = 0.0085,
#'     b = 3.11
#' )
#' params <- addSpecies(params, species_params)
#' plotSpectra(params)
#' @seealso [removeSpecies()], [renameSpecies()]
#' @export
#' @rdname addSpecies
addSpecies <- function(params, species_params,
                       gear_params = data.frame(), initial_effort,
                       interaction, steady = TRUE,
                       info_level = default_info_level(), ...) {
    UseMethod("addSpecies")
}

getPreservedParamsClass <- function(original_params, params) {
    original_class <- class(original_params)
    target_class <- class(params)
    if (identical(target_class, "MizerParams") &&
        !identical(original_class, "MizerParams")) {
        return(original_class)
    }
    target_class
}

restoreParamsClass <- function(params, target_class) {
    if (isS4(params)) {
        if (target_class[[1]] != "MizerParams") {
            params <- as(params, target_class[[1]])
        }
    } else {
        if (is.character(target_class)) {
            class(params) <- target_class
        }
    }
    params
}

copyPreservedParamsSlots <- function(params, old_params) {
    extra <- setdiff(names(old_params), names(params))
    params[extra] <- old_params[extra]
    params@other_dynamics <- old_params@other_dynamics
    params@initial_n_other <- old_params@initial_n_other
    params@other_encounter <- old_params@other_encounter
    params@other_mort <- old_params@other_mort
    params@other_params <- old_params@other_params
    params$capabilities <- old_params$capabilities
    params@rates_funcs <- old_params@rates_funcs
    params@use_predation_diffusion <- old_params@use_predation_diffusion
    params@second_order_w <- old_params@second_order_w
    params@defaults_edition <- defaults_edition(old_params)

    params@metadata <- old_params@metadata
    params@time_created <- old_params@time_created
    params@mizer_version <- old_params@mizer_version
    params@extensions <- old_params@extensions
    preserve_upgrade_guard(params, old_params)
}

copyParamsComments <- function(params, old_params) {
    comment(params) <- comment(old_params)
    for (slot in slotNames(params)) {
        comment(slot(params, slot)) <- comment(slot(old_params, slot))
    }
    params
}

# Give hooks explicit species and bin maps rather than inferring the meaning
# of arbitrary extension arrays from their dimensions.
structural_details <- function(before, after = before,
                               new_species = after$species_params$species,
                               rename = character(), ...) {
    old_species <- before$species_params$species
    renamed <- old_species
    renamed[match(names(rename), old_species)] <- unname(rename)
    origin <- before$w[[1L]]
    step <- log(before$w[[2L]] / origin)
    grid_map <- function(from, to) {
        match(round(log(from / origin) / step, 8),
              round(log(to / origin) / step, 8))
    }
    c(list(
        old_species = old_species, new_species = new_species,
        species_map = setNames(match(renamed, new_species), old_species),
        new_to_old = setNames(match(new_species, renamed), new_species),
        rename = rename,
        old_w = before$w, new_w = after$w,
        old_w_full = before$w_full, new_w_full = after$w_full,
        w_map = grid_map(before$w, after$w),
        new_w_to_old = grid_map(after$w, before$w),
        w_full_map = grid_map(before$w_full, after$w_full),
        new_w_full_to_old = grid_map(after$w_full, before$w_full)
    ), list(...))
}

# Core-owned axes are enumerated explicitly. Extension-owned objects are left
# for their declared transformation hooks, even when their dimensions match.
remap_core_species <- function(params, indices, drop_unused = FALSE) {
    p <- params
    old_species <- params$species_params$species
    species <- old_species[indices]
    for (field in c("linecolour", "linetype")) {
        extra <- setdiff(names(params[[field]]), old_species)
        p[[field]] <- params[[field]][c(species, extra)]
    }
    for (field in c("psi", "maturity", "initial_n", "intake_max",
                    "search_vol", "metab", "mu_b", "ext_encounter",
                    "ext_diffusion", "ft_mask")) {
        p[[field]] <- params[[field]][indices, , drop = FALSE]
    }
    for (field in c("ft_pred_kernel_e", "ft_pred_kernel_p",
                    "ft_pred_kernel_d")) {
        if (length(dim(params[[field]])) == 2L) {
            p[[field]] <- params[[field]][indices, , drop = FALSE]
        }
    }
    if (length(dim(params$pred_kernel)) == 3L) {
        p$pred_kernel <- params$pred_kernel[indices, , , drop = FALSE]
    }
    p$species_params <- params$species_params[indices, , drop = FALSE]
    p$given_species_params <-
        params$given_species_params[indices, , drop = FALSE]
    p$interaction <- params$interaction[indices, indices, drop = FALSE]
    p$selectivity <- params$selectivity[, indices, , drop = FALSE]
    p$catchability <- params$catchability[, indices, drop = FALSE]
    p$w_min_idx <- params$w_min_idx[indices]
    p$A <- params$A[indices]
    p$gear_params <- params$gear_params[
        params$gear_params$species %in% species, , drop = FALSE
    ]
    p$gear_params <- validGearParams(p$gear_params, p$species_params,
                                     edition = defaults_edition(p))
    if (drop_unused) {
        for (field in c("species_params", "given_species_params")) {
            keep <- colSums(is.na(p[[field]])) < nrow(p[[field]])
            p[[field]] <- p[[field]][, keep, drop = FALSE]
        }
    }
    copyParamsComments(p, params)
}

# Preserve only the explicitly known core axes at overlapping grid points.
restore_core_overlap <- function(candidate, before, details,
                                 preserve_species = details$old_species) {
    species <- which(details$old_species %in% preserve_species &
                         !is.na(details$species_map))
    target_species <- unname(details$species_map[species])
    bins <- which(!is.na(details$w_map))
    target_bins <- details$w_map[bins]
    resource_bins <- which(!is.na(details$w_full_map))
    target_resource <- details$w_full_map[resource_bins]
    for (field in c("initial_n", "psi", "maturity", "mu_b", "ext_encounter",
                    "ext_diffusion", "intake_max", "search_vol", "metab")) {
        candidate[[field]][target_species, target_bins] <-
            before[[field]][species, bins, drop = FALSE]
    }
    candidate$A[target_species] <- before$A[species]
    candidate$sc[target_bins] <- before$sc[bins]
    for (field in c("initial_n_pp", "cc_pp", "rr_pp")) {
        candidate[[field]][target_resource] <- before[[field]][resource_bins]
    }
    gears <- intersect(dimnames(before$selectivity)$gear,
                       dimnames(candidate$selectivity)$gear)
    candidate$selectivity[gears, target_species, target_bins] <-
        before$selectivity[gears, species, bins, drop = FALSE]
    candidate$catchability[gears, target_species] <-
        before$catchability[gears, species, drop = FALSE]
    candidate
}

# A manually supplied dense kernel has no implied extension beyond its grid.
# Retain its overlap and use zero in newly introduced bins until a hook or an
# explicit setter supplies another value.
map_frozen_core_arrays <- function(candidate, before, details) {
    if (length(dim(before$pred_kernel)) == 3L &&
            !is.null(comment(before$pred_kernel))) {
        kernel <- array(0, dim = c(dim(candidate$initial_n),
                                   length(candidate$w_full)),
                        dimnames = list(sp = candidate$species_params$species,
                                        w_pred = signif(candidate$w, 3),
                                        w_prey = signif(candidate$w_full, 3)))
        species <- which(!is.na(details$species_map))
        bins <- which(!is.na(details$w_map))
        prey <- which(!is.na(details$w_full_map))
        kernel[details$species_map[species], details$w_map[bins],
               details$w_full_map[prey]] <-
            before$pred_kernel[species, bins, prey, drop = FALSE]
        comment(kernel) <- comment(before$pred_kernel)
        candidate$pred_kernel <- kernel
    }
    for (field in c("selectivity", "catchability")) {
        comment(candidate[[field]]) <- comment(before[[field]])
    }
    candidate
}

#' @export
addSpecies.MizerParams <- function(params, species_params, gear_params = data.frame(),
                                   initial_effort = NULL, interaction = NULL,
                                   steady = TRUE,
                                   info_level = default_info_level(), ...) {
    # Collect the information signals raised while the species are added and
    # report them together at the end. Without this the two validation calls
    # below each report a misspelled column separately.
    with_info_level(info_level = info_level, {
    # check validity of parameters ----
    original_params <- params
    params <- validParams(params)
        before <- params
    target_class <- getPreservedParamsClass(original_params, params)
    given_species_params <- validGivenSpeciesParams(species_params)
    species_params <- validSpeciesParams(species_params)
    # The new species get the defaults of the model they are added to
    edition <- defaults_edition(params)
    gear_params <- validGearParams(gear_params, species_params,
                                   edition = edition)
    if (any(species_params$species %in% params@species_params$species)) {
        stop("You can not add species that are already there.")
    }
    if (steady) require_reproduction_balance_support(params, partial = TRUE)
    if (!is.null(comment(params@pred_kernel))) {
        stop("addSpecies() can not add species to a MizerParams object that ",
             "has its predation kernel protected by a comment.")
    }
    if (!is.null(comment(params@selectivity))) {
        stop("addSpecies() can not add species to a MizerParams object that ",
             "has its selectivity array protected by a comment.")
    }
    if (!is.null(comment(params@catchability))) {
        stop("addSpecies() can not add species to a MizerParams object that ",
             "has its catchability array protected by a comment.")
    }

    # set interaction ----
    no_old_sp <- nrow(params@species_params)
    old_sp <- 1:no_old_sp
    no_new_sp <- nrow(species_params)
    new_sp <- 1:no_new_sp + no_old_sp
    no_sp <- no_old_sp + no_new_sp
    if (missing(interaction)) {
        # keep existing interactions between old species and
        # set interactions involving new species to 1
        inter <- matrix(1, nrow = no_sp, ncol = no_sp)
        inter[old_sp, old_sp] <- params@interaction
    } else if (all(dim(interaction) == c(no_new_sp, no_new_sp))) {
        # keep existing interactions between old species,
        # set interactions involving an old and a new species to 1
        # and use supplied matrix for interaction among new species
        inter <- matrix(1, nrow = no_sp, ncol = no_sp)
        inter[old_sp, old_sp] <- params@interaction
        inter[new_sp, new_sp] <- interaction
    } else if (all(dim(interaction) != c(no_sp, no_sp))) {
        stop("Interaction matrix has invalid dimensions.")
    } else {
        inter <- interaction
    }

    # combine species params ----

    # Move linecolour and linetype into species_params
    params@species_params$linetype <-
        params@linetype[params@species_params$species]
    params@species_params$linecolour <-
        params@linecolour[params@species_params$species]

    # Make sure that all columns exist in both data frames
    missing <- setdiff(names(params@given_species_params), names(given_species_params))
    given_species_params[missing] <- NA
    missing <- setdiff(names(given_species_params), names(params@given_species_params))
    params@given_species_params[missing] <- NA

    missing <- setdiff(names(params@species_params), names(species_params))
    species_params[missing] <- NA
    missing <- setdiff(names(species_params), names(params@species_params))
    params@species_params[missing] <- NA

    # add the new species (with parameters described by species_params),
    # to make a larger species_params dataframe.
    combi_species_params <- rbind(params@species_params, species_params,
                                  stringsAsFactors = FALSE)
    combi_given_species_params <- rbind(params@given_species_params, given_species_params,
                                  stringsAsFactors = FALSE)

    # combine gear params ----
    if (!all(gear_params$species %in% species_params$species)) {
        stop("gear_params should only set gear parameters for new species.")
    }
    # Make sure that all columns exist in both data frames
    if (nrow(gear_params) > 0) {
        missing <- setdiff(names(params@gear_params), names(gear_params))
        gear_params[missing] <- NA
    }
    if (nrow(params@gear_params) > 0) {
        missing <- setdiff(names(gear_params), names(params@gear_params))
        params@gear_params[missing] <- NA
    }
    combi_gear_params <- rbind(params@gear_params, gear_params,
                               stringsAsFactors = FALSE)

    # expand grid ----
    # in case the new species need a bigger range of w
    # We need to make sure that the new grid that newMultispeciesParams()
    # will create contains the old grid as a subgrid.
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)
    max_w <- max(params@w)
    min_w <- min(params@w)
    new_max_w <- max_w
    new_min_w <- min_w
    new_no_w <- no_w
    extra_no_w <- 0  # extra bins added for smaller egg size
    if (max(species_params$w_max) > max(params@w) + .Machine$double.eps) {
        new_max_w <- max(species_params$w_max)
        dx <- log10(max_w / min_w) / (no_w - 1)
        new_no_w <- ceiling(log10(new_max_w / min_w) / dx) + 1
        new_max_w <- min_w * 10^(dx * (new_no_w - 1))
    }
    if (min(species_params$w_min) < min(params@w) - .Machine$double.eps) {
        new_min_w <- min(species_params$w_min)
        if (new_min_w < min(params@w_full)) {
            stop("The smallest egg size is too small.")
        }
        # We need to set the smallest egg size to a size on the existing grid
        # so that the new grid will be compatible
        sel_min <- combi_species_params$w_min == new_min_w
        new_min_w <- max(params@w_full[params@w_full <= new_min_w])
        combi_species_params$w_min[sel_min] <- new_min_w

        extra_no_w <- sum(params@w_full >= new_min_w) - no_w
        new_no_w <- new_no_w + extra_no_w
    }

    # new params object ----
    # use dataframe and global settings from params to make a new MizerParams
    # object. It is built with the defaults edition and the numerical scheme of
    # params, so that the arrays it computes afresh, such as the predation
    # kernels, use the quadrature of params.
        p <- with_defaults_edition(edition, new_multispecies_structure(
        combi_species_params,
        second_order_w = params@second_order_w,
        max_w = new_max_w,
        # for min_w_pp we choose something that will then be rounded down
        # to the existing smallest size when emptyParams() creates the new grid
        min_w_pp = (params@w_full[[1]] + params@w_full[[2]]) / 2,
        no_w = new_no_w,
        gear_params = combi_gear_params,
        kappa = params@resource_params$kappa,
        n = params@resource_params[["n"]],
        lambda = params@resource_params$lambda,
            w_pp_cutoff = params@resource_params$w_pp_cutoff
    ))
    p@given_species_params <- combi_given_species_params
        p@interaction[] <- inter

    # Set effort ----
    new_gear <- setdiff(unique(gear_params$gear),
                        unique(params@gear_params$gear))
    p@initial_effort[names(params@initial_effort)] <- params@initial_effort
    if (!missing(initial_effort)) {
        if (is.null(names(initial_effort))) {
            stop("The `initial_effort` must be a named list or vector, with one named entry for each new gear introduced in the `gear_params` argument. You should not provide an `initial_effort` argument if you did not introduce any new gear.")
        }
        if (!all(names(initial_effort) %in% new_gear)) {
            stop("The names of the `initial_effort` do not match the names of the new gears. You should not include effort values for existing gear. You should not provide an `initial_effort` argument if you did not introduce any new gear.")
        }
        p@initial_effort[names(initial_effort)] <- initial_effort
    }

    # Keep resource spectrum ----
    p@initial_n_pp[1:no_w_full] <- params@initial_n_pp
    p@cc_pp[1:no_w_full] <- params@cc_pp
    p@rr_pp[1:no_w_full] <- params@rr_pp
    p@resource_dynamics <- params@resource_dynamics
    p@resource_params <- params@resource_params

    # Copy old data ----
    # selector for old w bins inside new w
    old_w <- (extra_no_w + 1):(extra_no_w + no_w)
    p@A[old_sp] <- params@A
    p@psi[old_sp, old_w] <- params@psi
    p@maturity[old_sp, old_w] <- params@maturity
    p@sc[old_w] <- params@sc
    p@mu_b[old_sp, old_w] <- params@mu_b
    p@ext_encounter[old_sp, old_w] <- params@ext_encounter
    p@ext_diffusion[old_sp, old_w] <- params@ext_diffusion
    p@intake_max[old_sp, old_w] <- params@intake_max
    p@search_vol[old_sp, old_w] <- params@search_vol
    p@metab[old_sp, old_w] <- params@metab

    p <- copyPreservedParamsSlots(p, params)
        p@initial_n[old_sp, old_w] <- params@initial_n
        p <- restoreParamsClass(p, target_class)
        details <- structural_details(
            before, p, added = species_params$species, steady = steady
        )
        plan <- prepare_transformation(before, "add_species", details)
        before_hooks <- p
        p <- apply_transformation_hooks(p, plan)
        updates <- capture_transformation_updates(before_hooks, p, plan)
        p <- finish_multispecies_params(
            p, set_params_args = list(interaction = inter,
                                      initial_effort = p@initial_effort),
            resource_args = list(
                resource_rate = 10, resource_capacity = p@resource_params$kappa,
                resource_dynamics = p@resource_dynamics,
                lambda = p@resource_params$lambda, n = p@resource_params[["n"]],
                w_pp_cutoff = p@resource_params$w_pp_cutoff, balance = FALSE
            ), match_observations = FALSE, updates = updates,
            restore = function(candidate) {
                restore_core_overlap(candidate, params, details)
            }
        )

        p@resource_params <- params@resource_params
        p <- copyParamsComments(p, params)
        p <- restore_transformation_updates(p, updates)

    # The following does not affect the new species but preserves
    # any changes the user might have made in the original params object
    p <- setColours(p, params@linecolour)
    p <- setLinetypes(p, params@linetype)

    # we assume same background death for all species
    # p@mu_b[new_sp, ] <- rep(params@mu_b[1, ], each = no_new_sp)

    # initial solution ----

    if (steady) {
        # Turn off self-interaction among the new species, so we can determine the
        # growth rates, and death rates induced upon them by the pre-existing species
        p@interaction[new_sp, new_sp] <- 0

        # Compute solution for new species
        p <- steadySingleSpecies(p, species = new_sp)

        # set low abundance ----
        for (i in new_sp) {
            # Normalise solution so that it is never more than 1/100th of the
            # Sheldon spectrum.
            # We look at the maximum of abundance times w^lambda
            # because that is always an increasing function at small size.
            idx <- which.max(p@initial_n[i, ] * p@w^p@resource_params$lambda)
            if (p@initial_n[i, idx] == 0) {
                stop("The single-species steady state of ",
                     p@species_params$species[[i]],
                     " is zero at all sizes, so it can not be rescaled to a ",
                     "low abundance. Check the species parameters, in ",
                     "particular that the species can grow from w_min to ",
                     "w_max.")
            }
            p@initial_n[i, ] <- p@initial_n[i, ] *
                p@resource_params$kappa * p@w[idx]^(-p@resource_params$lambda) /
                p@initial_n[i, idx] / 100
            # TODO: remove next line after release of mizer 3.0
            p@A[i] <- sum(p@initial_n[i, ] * p@w * p@dw * p@maturity[i, ])
            p@species_params$is_background[i] <- FALSE
        }

        if (any(is.infinite(p@initial_n))) {
            stop("Candidate steady state holds infinities.")
        }
        if (any(is.na(p@initial_n) | is.nan(p@initial_n))) {
            stop("Candidate steady state holds non-numeric values.")
        }

        # Turn self interaction back on
        p@interaction[new_sp, new_sp] <- inter[new_sp, new_sp]

        target <- getRequiredRDD(p)[new_sp]
        if (identical(builtin_reproduction_law(p), "BevertonHoltRDD")) {
            p <- balanceRDD(p, rdd_target = target, reproduction_level = 1 / 4)
        } else {
            p <- balanceRDD(p, rdd_target = target)
        }
    }

        p <- validParams(p)
    p <- restoreParamsClass(p, target_class)
    })

    return(p)
}


#' Remove species
#'
#' This function removes all entries from the MizerParams object that
#' refer to the selected species. It does not recalculate the steady state for
#' the remaining species or retune their reproductive efficiency.
#'
#' Removing species changes the rates of the species that remain: the
#' predators of the removed species lose some of their food, so their encounter
#' rate drops, and the prey of the removed species lose some of their
#' predators, so their predation mortality drops. A model that was at steady
#' state is therefore no longer at steady state after the removal. If you want
#' the remaining species not to notice the removal, set `compensate = TRUE`, as
#' described in the section "Compensating for the removed species" below.
#'
#' If a gear was targeting only the removed species, then this function will
#' NOT remove that gear. If you want to also remove that gear then you can do
#' that by calling [setFishing()].
#'
#' @section Compensating for the removed species:
#' With `compensate = TRUE`, the external rates of the remaining species take
#' over what the removed species contributed to their rates:
#' * the food that the remaining species found in the removed species is added
#'   to their external encounter rate, see [setExtEncounter()],
#' * the predation mortality that the removed species imposed on the remaining
#'   species is added to their external mortality rate, see [setExtMort()],
#' * the diffusion that the removed prey species caused is added to the
#'   external diffusion rate, see [setExtDiffusion()]. This only matters if
#'   [use_predation_diffusion()] is switched on.
#'
#' Compensation is evaluated through the complete effective rate pipeline at
#' the initial state. After trimming species and running all declared hooks,
#' mizer checks growth, diffusion, total mortality and final recruitment of
#' every remaining species. A failed check stops the operation. These quantities
#' are preserved at this state; their responses to later changes may differ.
#'
#' For the built-in semichemostat and logistic resource, reduced predation is
#' balanced analytically by changing capacity, using the effective resource
#' mortality. This puts the resource at steady state and keeps its renewal rate.
#' A capacity set manually is replaced. A constant resource needs no adjustment.
#' Custom resource dynamics and other components require explicit removal
#' support; their hooks receive `details$compensate` and are responsible for
#' their own compensation. This is not a whole-model steady-state guarantee.
#'
#' The core adds each external-rate shortfall once. Custom callbacks that scale
#' or otherwise alter that response need an authored transformation hook with
#' permission to update the affected arrays. Mizer does not estimate a numerical
#' correction for an unknown callback.
#'
#' Each external rate that is changed is marked as set manually. A later
#' change to the species parameters that the rate is usually calculated from
#' (such as `E_ext`, `z0` or `D_ext`) therefore does not overwrite the
#' compensation; mizer tells you when it ignores such a change. External rates
#' to which the removed species contributed nothing are left unchanged and
#' unmarked.
#'
#' @param params A mizer params object for the original system.
#' @param species The species to be removed. A vector of species names, or a
#'   numeric vector of species indices, or a logical vector indicating for
#'   each species whether it is to be removed (TRUE) or not.
#' @param compensate If `TRUE`, the external rates of the remaining species
#'   are adjusted and checked to preserve growth, diffusion, mortality and
#'   recruitment at the initial state. Supported resource dynamics are balanced.
#'   See the section
#'   "Compensating for the removed species". Default is `FALSE`.
#' @param info_level Controls the amount of information messages that are
#'   shown. Higher levels lead to more messages. Set to 0 to suppress all such
#'   messages.
#' @param ... Currently unused.
#'
#' @return An object of type \linkS4class{MizerParams}
#' @seealso [addSpecies()], [renameSpecies()]
#' @export
#' @rdname removeSpecies
#' @examples
#' params <- NS_params
#' species_params(params)$species
#' params <- removeSpecies(params, c("Cod", "Haddock"))
#' species_params(params)$species
#'
#' # Removing Cod and Haddock moves the model off its steady state, unless the
#' # remaining species are compensated for their loss
#' rowSums(getSteadyResidual(removeSpecies(NS_params, c("Cod", "Haddock"))))
#' rowSums(getSteadyResidual(removeSpecies(NS_params, c("Cod", "Haddock"),
#'                                         compensate = TRUE)))
removeSpecies <- function(params, species, compensate = FALSE,
                          info_level = default_info_level(), ...) {
    UseMethod("removeSpecies")
}

#' @export
removeSpecies.MizerParams <- function(params, species, compensate = FALSE,
                                      info_level = default_info_level(), ...) {
    assert_that(is.flag(compensate), noNA(compensate))
    with_info_level(info_level = info_level, {
    params <- validParams(params)
    species <- valid_species_arg(params, species,
                                 return.logical = TRUE)
        keep <- which(!species)
        if (!length(keep)) {
            stop("At least one species must remain in the model.")
        }
        details <- structural_details(
            params, new_species = params$species_params$species[keep],
            removed = params$species_params$species[species],
            compensate = compensate
        )
        plan <- prepare_transformation(params, "remove_species", details)
    if (compensate) {
            compensation <- compensate_removal(params, species)
            params <- compensation$params
    }
        p <- remap_core_species(params, keep, drop_unused = TRUE)
        p <- apply_transformation_hooks(p, plan)

        p <- validParams(p)
        if (compensate) {
            verify_removal_compensation(p, compensation, keep)
        }

    p@time_modified <- lubridate::now()
    })

    return(p)
}

# Compute core compensation on the intact model so every extension still sees
# its original species dimensions. Author hooks then adapt the reduced model.
compensate_removal <- function(params, remove) {
    target <- tuning_rate_context(params)
    params@initial_n[remove, ] <- 0
    effective_rate <- function(field) {
        force(field)
        function(params) tuning_rate_context(params)[[field]]
    }
    params <- compensate_rate(params, "ext_encounter", remove,
                              target$encounter, effective_rate("encounter"))
    params <- compensate_rate(params, "mu_b", remove,
                              target$mort, effective_rate("mort"))
    params <- compensate_rate(params, "ext_diffusion", remove,
                              target$diffusion, effective_rate("diffusion"))
    rates <- tuning_rate_context(params)
    law <- sub("^mizer::", "", params$resource_dynamics)
    balance_resource <- law %in%
        c("resource_semichemostat", "resource_logistic") &&
        !tuning_values_equal(target$resource_mort, rates$resource_mort)
    if (balance_resource) {
        params <- balance_resource_strategy(params, rates = rates,
                                            n_pp = params$initial_n_pp)
    }
    list(params = params, rates = target, resource_balanced = balance_resource)
}

# Add each shortfall once. Unknown callback responses are the responsibility
# of their explicit transformation hooks, never inferred by numerical probes.
compensate_rate <- function(params, slot, remove, target, get_rate) {
    compensation <- unclass_rate(target) - unclass_rate(get_rate(params))
    compensation[remove, ] <- 0
    if (all(compensation == 0)) return(params)
    add_compensation(params, slot, compensation)
}

# Check the complete effective pipeline again after all hooks and trimming.
# Resource compensation retains the established built-in convention of putting
# it at steady state; this is not a whole-model steady-state guarantee.
verify_removal_compensation <- function(params, compensation, keep) {
    actual <- tuning_rate_context(params)
    for (field in c("e_growth", "diffusion", "mort", "rdd")) {
        target <- compensation$rates[[field]]
        target <- if (field == "rdd") {
            target[keep]
        } else {
            target[keep, , drop = FALSE]
        }
        if (!tuning_values_equal(actual[[field]], target)) {
            stop("Species-removal compensation did not preserve effective ",
                 field, " for the remaining species. Supply a removal hook ",
                 "that supports `compensate = TRUE`.", call. = FALSE)
        }
    }
    if (compensation$resource_balanced) {
        terms <- builtin_resource_equation_rates(
            params, n_pp = params$initial_n_pp, rates = actual
        )
        if (!tuning_values_equal(terms[[1L]], terms[[2L]])) {
            stop("Species-removal compensation did not balance the resource. ",
                 "Supply a removal hook that supports `compensate = TRUE`.",
                 call. = FALSE)
        }
    }
    invisible(TRUE)
}

#' Add the compensation for removed species to an external rate
#'
#' Adds `compensation` to the external rate array `slot` and marks that array
#' as set manually, so that a later recalculation from the species parameters
#' does not overwrite it.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param slot One of `"ext_encounter"`, `"mu_b"` or `"ext_diffusion"`.
#' @param compensation An array (species x size) to add to the slot.
#' @return The `MizerParams` object with the slot updated.
#' @noRd
add_compensation <- function(params, slot, compensation) {
    value <- slot(params, slot) + compensation
    comment(value) <- comment(slot(params, slot)) %||%
        "increased by removeSpecies() to compensate for the removed species"
    switch(slot,
           ext_encounter = setExtEncounter(params, ext_encounter = value),
           mu_b = setExtMort(params, ext_mort = value),
           ext_diffusion = setExtDiffusion(params, ext_diffusion = value))
}


#' Rename species
#'
#' Changes the names of species in a MizerParams object. This involves for
#' example changing the species dimension names of rate arrays appropriately.
#'
#' @param params A mizer params object
#' @param replace A named character vector, with new names as values, and old
#'   names as names.
#' @param ... Currently unused.
#'
#' @return An object of type \linkS4class{MizerParams}
#' @seealso [renameGear()]
#' @export
#' @rdname renameSpecies
#' @examples
#' replace <- c(Cod = "Kabeljau", Haddock = "Schellfisch")
#' params <- renameSpecies(NS_params, replace)
#' species_params(params)$species
renameSpecies <- function(params, replace, ...) {
    UseMethod("renameSpecies")
}

#' @export
renameSpecies.MizerParams <- function(params, replace, ...) {
    params <- validParams(params)
    before <- params
    replace[] <- as.character(replace)
    to_replace <- names(replace)
    species <- params@species_params$species
    wrong <- setdiff(names(replace), species)
    if (length(wrong) > 0) {
        stop(paste(wrong, collapse = ", "),
             " do not exist.")
    }
    names(species) <- species
    species[to_replace] <- replace
    names(species) <- NULL
    if (anyNA(species) || any(!nzchar(species)) || anyDuplicated(species)) {
        stop("Species names must be nonempty and unique.")
    }
    details <- structural_details(before, new_species = species,
                                  rename = replace)
    plan <- prepare_transformation(before, "rename_species", details)
    rownames(params@species_params) <- species
    params@species_params$species <- species
    rownames(params@given_species_params) <- species
    params@given_species_params$species <- species
    params@gear_params$species <- params@gear_params$species
    for (i in seq_len(nrow(params@gear_params))) {
        if (params@gear_params$species[[i]] %in% names(replace)) {
            params@gear_params$species[[i]] <-
                replace[[params@gear_params$species[[i]]]]
        }
    }
    params@gear_params <- validGearParams(params@gear_params,
                                          params@species_params,
                                          edition = defaults_edition(params))
    # rename line colours
    linenames <- names(params@linecolour)
    names(linenames) <- linenames
    linenames[to_replace] <- replace
    names(linenames) <- NULL
    names(params@linecolour) <- linenames
    # rename line types
    linenames <- names(params@linetype)
    names(linenames) <- linenames
    linenames[to_replace] <- replace
    names(linenames) <- NULL
    names(params@linetype) <- linenames

    names(params@w_min_idx) <- species
    dimnames(params@maturity)$sp <- species
    dimnames(params@psi)$sp <- species
    dimnames(params@initial_n)$sp <- species
    dimnames(params@intake_max)$sp <- species
    dimnames(params@search_vol)$sp <- species
    dimnames(params@metab)$sp <- species
    if (length(dim(params@ft_pred_kernel_e)) == 2) {
        dimnames(params@ft_pred_kernel_e)$sp <- species
        dimnames(params@ft_pred_kernel_p)$sp <- species
        dimnames(params@ft_pred_kernel_d)$sp <- species
    } else {
        dimnames(params@pred_kernel)$sp <- species
    }
    dimnames(params@mu_b)$sp <- species
    dimnames(params@ext_encounter)$sp <- species
    dimnames(params@ext_diffusion)$sp <- species
    dimnames(params@interaction)$predator <- species
    dimnames(params@interaction)$prey <- species
    dimnames(params@selectivity)$sp <- species
    dimnames(params@catchability)$sp <- species
    if (!is.null(dimnames(params@ft_mask)[[1]])) {
        dimnames(params@ft_mask)[[1]] <- species
    }

    params <- apply_transformation_hooks(params, plan)
    params <- validParams(params)

    params@time_modified <- lubridate::now()
    return(params)
}
#' Reorder species
#'
#' Reorders all core species axes and invokes registered extension hooks.
#' Every current species must occur exactly once in `species`.
#'
#' @param params A MizerParams object.
#' @param species Species names in the desired order, or a permutation of their
#'   current integer indices.
#' @param ... Unused.
#' @return The updated MizerParams object.
#' @seealso [addSpecies()], [removeSpecies()], [renameSpecies()]
#' @export
reorderSpecies <- function(params, species, ...) {
    UseMethod("reorderSpecies")
}

#' @rdname reorderSpecies
#' @usage NULL
#' @export
reorderSpecies.MizerParams <- function(params, species, ...) {
    params <- validParams(params)
    current <- params$species_params$species
    if (is.numeric(species) && !is.complex(species) &&
            all(is.finite(species)) && all(species == floor(species)) &&
            setequal(species, seq_along(current))) {
        species <- current[species]
    }
    if (!is.character(species) || anyNA(species) ||
            length(species) != length(current) || anyDuplicated(species) ||
            !setequal(species, current)) {
        stop("`species` must contain every species exactly once.")
    }
    details <- structural_details(params, new_species = species)
    plan <- prepare_transformation(params, "reorder_species", details)
    candidate <- remap_core_species(params, match(species, current))
    candidate <- apply_transformation_hooks(candidate, plan)
    candidate <- validParams(candidate)
    candidate$time_modified <- lubridate::now()
    candidate
}

#' Adjust the size grid
#'
#' This function adjusts the size grid in a [MizerParams] object to the desired
#' minimum and maximum size. It can both expand and truncate the grid.
#' If the grid is truncated, any data outside the new grid is discarded.
#' A warning is issued if there is non-negligible biomass in the discarded
#' size bins. Manually supplied dense predation kernels and fishing arrays
#' retain their existing overlap; newly introduced bins are zero unless an
#' explicit transformation hook supplies other values.
#'
#' @param params A [MizerParams] object.
#' @param new_min_w The new minimum size in the grid. Defaults to the minimum
#'   egg size of all species.
#' @param new_max_w The new maximum size in the grid. Defaults to the maximum
#'   asymptotic size of all species.
#' @param new_min_w_pp The new minimum size of the resource spectrum. Defaults
#'   to the current minimum of `w_full`.
#' @param preserve_species A vector of species names for which all rate arrays
#'   should be copied over to the new params object rather than being
#'   re-calculated from the species parameters. If missing, all species are
#'   preserved.
#' @param tol A numeric value specifying the tolerance for truncation losses.
#'   The following checks are made separately for each species and a warning is
#'   raised, listing the affected species, if the lost fraction exceeds this
#'   value for any of them: the fraction of the species' biomass lost, the
#'   fraction of the diet of the smallest individuals of the species lost to
#'   resource truncation, and the fraction of the diet of the largest
#'   individuals of the species lost to resource truncation. Defaults to
#'   `1e-6`.
#' @param ... Additional arguments.
#'
#' @return A new [MizerParams] object with the updated size grid.
#' @export
#' @rdname adjustSizeGrid
adjustSizeGrid <- function(params, ...) {
    UseMethod("adjustSizeGrid")
}

#' @export
#' @rdname adjustSizeGrid
adjustSizeGrid.MizerParams <- function(params,
                                       new_min_w = min(params@species_params$w_min),
                                       new_max_w = max(params@species_params$w_max),
                                       new_min_w_pp = min(params@w_full),
                                       preserve_species = params@species_params$species,
                                       tol = 1e-6,
                                       ...) {
    target_class <- class(params)
    params <- validParams(params)
    before <- params
    sp_sel <- valid_species_arg(params, preserve_species, return.logical = TRUE)
    min_w <- min(params@w)
    max_w <- max(params@w)
    # Short-circuit if no adjustment is needed
    if (abs(new_min_w - min_w) < .Machine$double.eps &&
        abs(new_max_w - max_w) < .Machine$double.eps &&
        abs(new_min_w_pp - min(params@w_full)) < .Machine$double.eps &&
        all(preserve_species == params@species_params$species)) {
        params@time_modified <- lubridate::now()
        return(params)
    }

    dx <- log10(params@w[2] / params@w[1])

    # Align/snap new_min_w
    if (new_min_w < min(params@w_full) - .Machine$double.eps) {
        stop("The smallest egg size is too small.")
    }
    k_min <- floor(log10(new_min_w / min_w) / dx + 1e-9)
    k_min_limit <- round(log10(min(params@w_full) / min_w) / dx)
    if (k_min < k_min_limit) {
        k_min <- k_min_limit
    }
    new_min_w <- min_w * 10^(k_min * dx)

    # Align/snap new_max_w
    k_max <- ceiling(log10(new_max_w / min_w) / dx - 1e-9)
    new_max_w <- min_w * 10^(k_max * dx)

    # Align/snap new_min_w_pp
    if (new_min_w_pp < min(params@w_full) - .Machine$double.eps) {
        stop("The smallest resource size is too small.")
    }
    k_min_pp <- floor(log10(new_min_w_pp / min_w) / dx + 1e-9)
    k_min_pp_limit <- round(log10(min(params@w_full) / min_w) / dx)
    if (k_min_pp < k_min_pp_limit) {
        k_min_pp <- k_min_pp_limit
    }
    new_min_w_pp <- min_w * 10^(k_min_pp * dx)

    # Safety checks
    if (new_min_w >= new_max_w) {
        stop("new_min_w must be smaller than new_max_w.")
    }
    if (new_min_w_pp >= new_min_w) {
        stop("new_min_w_pp must be smaller than new_min_w.")
    }
    too_small_species <- params@species_params$species[sp_sel & params@species_params$w_max < new_min_w]
    if (length(too_small_species) > 0) {
        stop("The following species have their upper size-grid boundary w_max smaller than the new minimum size: ",
             paste(too_small_species, collapse = ", "))
    }
    too_large_species <- params@species_params$species[sp_sel & params@species_params$w_min > new_max_w]
    if (length(too_large_species) > 0) {
        stop("The following species have their minimum size w_min larger than the new maximum size: ",
             paste(too_large_species, collapse = ", "))
    }

    # Calculate new number of bins
    new_no_w <- round(log10(new_max_w / new_min_w) / dx) + 1

    # Setup constructor parameters
    sp_params_for_grid <- params@species_params
    sp_params_for_grid$linetype <- params@linetype[params@species_params$species]
    sp_params_for_grid$linecolour <- params@linecolour[params@species_params$species]
    sp_params_for_grid$w_min <- pmax(sp_params_for_grid$w_min, new_min_w)
    sp_params_for_grid$w_max <- pmin(sp_params_for_grid$w_max, new_max_w)
    sp_params_for_grid$w_mat <- sp_params_for_grid$w_min + (sp_params_for_grid$w_max - sp_params_for_grid$w_min) * 0.5
    sp_params_for_grid$w_mat25 <- sp_params_for_grid$w_min + (sp_params_for_grid$w_mat - sp_params_for_grid$w_min) * 0.5
    sp_params_for_grid$w_repro_max <- sp_params_for_grid$w_max

    idx_min <- which.min(sp_params_for_grid$w_min)
    sp_params_for_grid$w_min[idx_min] <- new_min_w

    # Built with the defaults edition and the numerical scheme of params
    edition <- defaults_edition(params)
    p <- with_defaults_edition(edition, new_multispecies_structure(
        sp_params_for_grid,
        second_order_w = params@second_order_w,
        max_w = new_max_w,
        min_w = new_min_w,
        min_w_pp = new_min_w_pp * 10^(0.5 * dx),
        no_w = new_no_w,
        gear_params = params@gear_params,
        kappa = params@resource_params$kappa,
        n = params@resource_params[["n"]],
        lambda = params@resource_params$lambda,
        w_pp_cutoff = params@resource_params$w_pp_cutoff
    ))

    # Update species_params and given_species_params to align with the new grid boundaries
    new_sp_params <- params@species_params
    new_sp_params$w_min <- pmax(new_sp_params$w_min, new_min_w)
    new_sp_params$w_max <- pmin(new_sp_params$w_max, new_max_w)
    new_sp_params$w_mat <- pmin(new_sp_params$w_mat, new_sp_params$w_max * 0.9)
    new_sp_params$w_mat25 <- pmin(new_sp_params$w_mat25, new_sp_params$w_mat * 0.9)
    if ("w_repro_max" %in% names(new_sp_params)) {
        new_sp_params$w_repro_max <- pmin(new_sp_params$w_repro_max, new_sp_params$w_max)
    }
    p@species_params <- new_sp_params
    p@w_min_idx <- get_w_min_idx(new_sp_params, p@w)

    new_given_sp_params <- params@given_species_params
    if ("w_min" %in% names(new_given_sp_params)) {
        new_given_sp_params$w_min <- pmax(new_given_sp_params$w_min, new_min_w)
    }
    if ("w_max" %in% names(new_given_sp_params)) {
        new_given_sp_params$w_max <- pmin(new_given_sp_params$w_max, new_max_w)
    }
    if ("w_mat" %in% names(new_given_sp_params)) {
        new_given_sp_params$w_mat <- pmin(new_given_sp_params$w_mat, new_sp_params$w_max * 0.9)
    }
    if ("w_mat25" %in% names(new_given_sp_params)) {
        new_given_sp_params$w_mat25 <- pmin(new_given_sp_params$w_mat25, new_sp_params$w_mat * 0.9)
    }
    if ("w_repro_max" %in% names(new_given_sp_params)) {
        new_given_sp_params$w_repro_max <- pmin(new_given_sp_params$w_repro_max, new_sp_params$w_max)
    }
    p@given_species_params <- new_given_sp_params

    # Now calculate grid indices
    old_grid_k <- round(log10(params@w / min_w) / dx)
    new_grid_k <- round(log10(p@w / min_w) / dx)
    match_new <- match(old_grid_k, new_grid_k)
    old_idx <- which(!is.na(match_new))
    new_idx <- match_new[old_idx]

    old_full_k <- round(log10(params@w_full / min_w) / dx)
    new_full_k <- round(log10(p@w_full / min_w) / dx)
    match_full_new <- match(old_full_k, new_full_k)
    old_full_idx <- which(!is.na(match_full_new))
    new_full_idx <- match_full_new[old_full_idx]

    p@interaction[] <- params@interaction
    p@initial_effort <- params@initial_effort
    # Copy over overlapping data
    p@initial_n[sp_sel, new_idx] <- params@initial_n[sp_sel, old_idx]
    p@A[sp_sel] <- params@A[sp_sel]
    p@psi[sp_sel, new_idx] <- params@psi[sp_sel, old_idx]
    p@maturity[sp_sel, new_idx] <- params@maturity[sp_sel, old_idx]
    p@sc[new_idx] <- params@sc[old_idx]
    p@mu_b[sp_sel, new_idx] <- params@mu_b[sp_sel, old_idx]
    p@ext_encounter[sp_sel, new_idx] <- params@ext_encounter[sp_sel, old_idx]
    p@ext_diffusion[sp_sel, new_idx] <- params@ext_diffusion[sp_sel, old_idx]
    p@intake_max[sp_sel, new_idx] <- params@intake_max[sp_sel, old_idx]
    p@search_vol[sp_sel, new_idx] <- params@search_vol[sp_sel, old_idx]
    p@metab[sp_sel, new_idx] <- params@metab[sp_sel, old_idx]

    p@initial_n_pp[new_full_idx] <- params@initial_n_pp[old_full_idx]
    p@cc_pp[new_full_idx] <- params@cc_pp[old_full_idx]
    p@rr_pp[new_full_idx] <- params@rr_pp[old_full_idx]
    p@resource_dynamics <- params@resource_dynamics
    p@resource_params <- params@resource_params

    p <- copyPreservedParamsSlots(p, params)
    p <- restoreParamsClass(p, target_class)
    details <- structural_details(
        before, p, preserve_species = params$species_params$species[sp_sel]
    )
    plan <- prepare_transformation(before, "size_grid", details)
    p <- restore_core_overlap(p, params, details,
                              params$species_params$species[sp_sel])
    p <- map_frozen_core_arrays(p, params, details)
    before_hooks <- p
    p <- apply_transformation_hooks(p, plan)
    updates <- capture_transformation_updates(before_hooks, p, plan)
    p <- finish_multispecies_params(
        p, set_params_args = list(interaction = params@interaction,
                                  initial_effort = params@initial_effort),
        resource_args = list(
            resource_rate = 10, resource_capacity = p@resource_params$kappa,
            resource_dynamics = p@resource_dynamics,
            lambda = p@resource_params$lambda, n = p@resource_params[["n"]],
            w_pp_cutoff = p@resource_params$w_pp_cutoff, balance = FALSE
        ), match_observations = FALSE, updates = updates,
        restore = function(candidate) {
            restore_core_overlap(candidate, params, details,
                                 params$species_params$species[sp_sel])
        }
    )

    # Check for biomass/diet loss warnings
    truncated_idx <- setdiff(seq_along(params@w), old_idx)

    # Low-end resource truncation (below new_min_w_pp)
    truncated_full_low_idx <- which(params@w_full < min(p@w_full) - .Machine$double.eps)
    # High-end resource truncation (above new_max_w)
    truncated_full_high_idx <- which(params@w_full > max(p@w_full) + .Machine$double.eps)

    if (length(truncated_idx) > 0) {
        lost_fracs <- sapply(seq_along(params@species_params$species), function(sp_idx) {
            tot <- sum(params@initial_n[sp_idx, ] * params@w * params@dw)
            if (tot == 0) return(0)
            sum(params@initial_n[sp_idx, truncated_idx] * params@w[truncated_idx] * params@dw[truncated_idx]) / tot
        })
        names(lost_fracs) <- params@species_params$species
        warn_sp <- lost_fracs[lost_fracs > tol]
        if (length(warn_sp) > 0) {
            warning("Non-negligible species biomass was lost due to grid truncation: ",
                    paste(names(warn_sp), sprintf("(%g%%)", signif(warn_sp * 100, 3)), collapse = ", "))
        }
    }

    # Resource low-end truncation: check diet loss of smallest fish
    if (length(truncated_full_low_idx) > 0) {
        pk <- pred_kernel(params)
        encounter_old <- getEncounter(params)

        lost_diet_fracs <- sapply(seq_along(params@species_params$species), function(sp_idx) {
            w_egg_idx <- params@w_min_idx[sp_idx]
            tot_diet <- encounter_old[sp_idx, w_egg_idx]
            if (tot_diet <= 0) return(0)

            lost_enc <- params@search_vol[sp_idx, w_egg_idx] *
                params@species_params$interaction_resource[sp_idx] *
                sum(pk[sp_idx, w_egg_idx, truncated_full_low_idx] *
                    params@w_full[truncated_full_low_idx] *
                    params@dw_full[truncated_full_low_idx] *
                    params@initial_n_pp[truncated_full_low_idx])

            return(lost_enc / tot_diet)
        })
        names(lost_diet_fracs) <- params@species_params$species
        warn_diet <- lost_diet_fracs[lost_diet_fracs > tol]
        if (length(warn_diet) > 0) {
            warning("Non-negligible diet of smallest fish was lost due to resource truncation: ",
                    paste(names(warn_diet), sprintf("(%g%%)", signif(warn_diet * 100, 3)), collapse = ", "))
        }
    }

    # Resource high-end truncation: check diet loss of largest fish
    if (length(truncated_full_high_idx) > 0) {
        pk <- pred_kernel(params)
        encounter_old <- getEncounter(params)

        lost_enc_fracs <- sapply(seq_along(params@species_params$species), function(sp_idx) {
            # Largest size bin the species occupies
            w_top_idx <- max(which(params@w <= params@species_params$w_max[sp_idx]))
            tot_enc <- encounter_old[sp_idx, w_top_idx]
            if (tot_enc <= 0) return(0)

            lost_enc <- params@search_vol[sp_idx, w_top_idx] *
                params@species_params$interaction_resource[sp_idx] *
                sum(pk[sp_idx, w_top_idx, truncated_full_high_idx] *
                    params@w_full[truncated_full_high_idx] *
                    params@dw_full[truncated_full_high_idx] *
                    params@initial_n_pp[truncated_full_high_idx])

            return(lost_enc / tot_enc)
        })
        names(lost_enc_fracs) <- params@species_params$species
        warn_enc <- lost_enc_fracs[lost_enc_fracs > tol]
        if (length(warn_enc) > 0) {
            warning("Non-negligible diet of largest fish was lost due to resource truncation: ",
                    paste(names(warn_enc), sprintf("(%g%%)", signif(warn_enc * 100, 3)), collapse = ", "))
        }
    }

    p@resource_params <- params@resource_params

    # Restore core plotting metadata after rate construction.
    p <- setColours(p, params@linecolour)
    p <- setLinetypes(p, params@linetype)

    # Preserve comments
    p <- copyParamsComments(p, params)
    p <- restore_transformation_updates(p, updates)

    p <- validParams(p)

    # Restore plot metadata exactly
    p@linecolour <- params@linecolour
    p@linetype <- params@linetype

    p <- restoreParamsClass(p, target_class)

    return(p)
}

#' Expand the size grid
#'
#' `r lifecycle::badge("deprecated")`
#' This function expands the size grid in a [MizerParams] object to the desired
#' min and max size, preserving all existing species. The function is deprecated
#' because you can achieve the same more flexibly with `adjustSizeGrid()`.
#'
#' @param params A [MizerParams] object.
#' @param new_min_w The new minimum size in the grid. Defaults to the current minimum.
#' @param new_max_w The new maximum size in the grid. Defaults to the current maximum.
#' @param preserve_species A vector of species names for which all rate arrays
#'   should be copied over to the new params object rather than being
#'   re-calculated from the species parameters. If missing, all species are
#'   preserved.
#' @param ... Additional arguments (currently unused).
#'
#' @return A new [MizerParams] object with the updated size grid.
#' @export
#' @rdname expandSizeGrid
#' @concept deprecated
expandSizeGrid <- function(params, ...) {
    lifecycle::deprecate_warn("3.1.1", "expandSizeGrid()", "adjustSizeGrid()")
    UseMethod("expandSizeGrid")
}

#' @export
#' @rdname expandSizeGrid
expandSizeGrid.MizerParams <- function(params,
                                       new_min_w = min(params@w),
                                       new_max_w = max(params@w),
                                       preserve_species = params@species_params$species,
                                       ...) {
    min_w <- min(params@w)
    max_w <- max(params@w)
    if (new_min_w > min_w || new_max_w < max_w) {
        stop("`expandSizeGrid()` can only expand, not shrink the grid.")
    }
    adjustSizeGrid(params, new_min_w = new_min_w, new_max_w = new_max_w,
                   preserve_species = preserve_species, ...)
}

#' Rename gears
#'
#' Changes the names of gears in a MizerParams object. This involves for
#' example changing the gear dimension names of selectivity and catchability
#' arrays appropriately.
#'
#' @param params A mizer params object
#' @param replace A named character vector, with new names as values, and old
#'   names as names.
#' @param ... Currently unused.
#'
#' @return An object of type \linkS4class{MizerParams}
#' @seealso [renameSpecies()]
#' @export
#' @rdname renameGear
#' @examples
#' replace <- c(Industrial = "Trawl", Otter = "Beam_Trawl")
#' params <- renameGear(NS_params, replace)
#' gear_params(params)$gear
renameGear <- function(params, replace, ...) {
    UseMethod("renameGear")
}

#' @export
renameGear.MizerParams <- function(params, replace, ...) {
    params <- validParams(params)
    replace[] <- as.character(replace)
    to_replace <- names(replace)
    gears <- dimnames(params@selectivity)$gear
    wrong <- setdiff(names(replace), gears)
    if (length(wrong) > 0) {
        stop(paste(wrong, collapse = ", "),
             " do not exist.")
    }
    names(gears) <- gears
    gears[to_replace] <- replace
    names(gears) <- NULL

    # Update gear_params data frame
    for (i in seq_len(nrow(params@gear_params))) {
        if (params@gear_params$gear[[i]] %in% names(replace)) {
            params@gear_params$gear[[i]] <-
                replace[[params@gear_params$gear[[i]]]]
        }
    }
    params@gear_params <- validGearParams(params@gear_params,
                                          params@species_params,
                                          edition = defaults_edition(params))

    # Update dimension names in arrays
    dimnames(params@selectivity)$gear <- gears
    dimnames(params@catchability)$gear <- gears

    # Update initial_effort names
    gearnames <- names(params@initial_effort)
    names(gearnames) <- gearnames
    gearnames[to_replace] <- replace
    names(gearnames) <- NULL
    names(params@initial_effort) <- gearnames

    validObject(params)

    params@time_modified <- lubridate::now()
    return(params)
}
