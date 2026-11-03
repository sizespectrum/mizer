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
#'   The reproduction of the new species is balanced with [balanceRDD()], so
#'   that the model keeps its density-dependence function. With Beverton-Holt
#'   reproduction the new species will have a reproduction level of 1/4, this
#'   can then be changed with [reproduction_level<-()].
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
    params@other_dynamics <- old_params@other_dynamics
    params@initial_n_other <- old_params@initial_n_other
    params@other_encounter <- old_params@other_encounter
    params@other_mort <- old_params@other_mort
    params@other_params <- old_params@other_params
    params@rates_funcs <- old_params@rates_funcs
    params@use_predation_diffusion <- old_params@use_predation_diffusion
    params@second_order_w <- old_params@second_order_w
    params@defaults_edition <- defaults_edition(old_params)

    params@metadata <- old_params@metadata
    params@time_created <- old_params@time_created
    params@mizer_version <- old_params@mizer_version
    params@extensions <- old_params@extensions
    params
}

copyParamsComments <- function(params, old_params) {
    comment(params) <- comment(old_params)
    for (slot in slotNames(params)) {
        comment(slot(params, slot)) <- comment(slot(old_params, slot))
    }
    params
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
    p <- with_defaults_edition(edition, newMultispeciesParams(
        combi_species_params,
        interaction = inter,
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
        w_pp_cutoff = params@resource_params$w_pp_cutoff,
        info_level = info_level
    ))
    p@given_species_params <- combi_given_species_params

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

    # Preserve comments ----
    p <- copyParamsComments(p, params)

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

    # The following does not affect the new species but preserves
    # any changes the user might have made in the original params object
    p <- setColours(p, params@linecolour)
    p <- setLinetypes(p, params@linetype)

    # we assume same background death for all species
    # p@mu_b[new_sp, ] <- rep(params@mu_b[1, ], each = no_new_sp)

    # initial solution ----
    p@initial_n[old_sp, old_w] <- params@initial_n

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

        # Balance the reproduction of the new species. With Beverton-Holt
        # they get a reproduction level of 1/4. Other density dependence is
        # balanced in its own way, and the old species are left alone. The
        # extension class is needed for this, because an extension may change
        # the reproduction rate and take part in balancing it.
        p <- restoreParamsClass(p, target_class)
        rdd_target <- getRequiredRDD(p)
        rdd_target[-new_sp] <- NA
        if (uses_beverton_holt_balancer(p)) {
            p <- balanceRDD(p, rdd_target = rdd_target,
                            reproduction_level = species_vector(p, 1 / 4))
        } else {
            p <- balanceRDD(p, rdd_target = rdd_target)
        }
    }

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
#' All of this is calculated at the initial state stored in `params`. At that
#' state the remaining species then have exactly the rates they had before the
#' removal. In effect the removed species are still present, held at their
#' initial abundances, but they no longer respond to changes in the rest of the
#' model. As the model moves away from its initial state, the compensation no
#' longer matches what the removed species would have done.
#'
#' The removed species also no longer eat the resource, so the resource
#' capacity is rebalanced with [setResource()]: it is set so that at the
#' initial state the resource replenishes at exactly the rate at which it is now
#' consumed. The resource rate is kept. This replaces a resource capacity that
#' you had set manually. Rebalancing puts the resource at steady state, so it
#' keeps its rate of change only if that was zero. Thus if the model was at
#' steady state, it stays at steady state. If the resource was not at steady
#' state, it is put at steady state instead of keeping the rate at which it
#' was changing.
#'
#' The compensation is added to the external rates on the assumption that
#' they are added unchanged to the rates, as they are in mizer's own rate
#' functions. If you have replaced a rate function with [setRateFunction()],
#' mizer checks that the rates are restored. Where they are not, it rescales
#' the compensation by how much the rate responds to it, which restores a rate
#' that your function calculates by multiplying mizer's rate by a factor. If a
#' rate is still not restored, you are warned.
#'
#' Each external rate that is changed is marked as set manually. A later
#' change to the species parameters that the rate is usually calculated from
#' (such as `E_ext`, `z0` or `D_ext`) therefore does not overwrite the
#' compensation; mizer tells you when it ignores such a change. External rates
#' to which the removed species contributed nothing are left unchanged and
#' unmarked.
#'
#' The compensation covers the species and the resource. Other components
#' added with [setComponent()] that have dynamics of their own are not
#' compensated, and you are warned about them.
#'
#' @param params A mizer params object for the original system.
#' @param species The species to be removed. A vector of species names, or a
#'   numeric vector of species indices, or a logical vector indicating for
#'   each species whether it is to be removed (TRUE) or not.
#' @param compensate If `TRUE`, the external rates of the remaining species
#'   are adjusted so that the remaining species have the same rates as before
#'   the removal, and the resource is rebalanced. See the section
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
    if (compensate) {
        params <- compensate_removal(params, species)
    }
    keep <- !species
    p <- params

    # Select only the parts corresponding the species we keep
    p@linecolour <-
        params@linecolour[!(names(params@linecolour) %in%
                                params@species_params$species[species])]
    p@linetype <-
        params@linetype[!(names(params@linetype) %in%
                              params@species_params$species[species])]
    p@psi <- params@psi[keep, , drop = FALSE]
    p@maturity <- params@maturity[keep, , drop = FALSE]
    p@initial_n <- params@initial_n[keep, , drop = FALSE]
    p@intake_max <- params@intake_max[keep, , drop = FALSE]
    p@search_vol <- params@search_vol[keep, , drop = FALSE]
    p@metab <- params@metab[keep, , drop = FALSE]
    if (length(dim(params@ft_pred_kernel_e)) == 2) {
        p@ft_pred_kernel_e <- params@ft_pred_kernel_e[keep, , drop = FALSE]
    }
    if (length(dim(params@ft_pred_kernel_p)) == 2) {
        p@ft_pred_kernel_p <- params@ft_pred_kernel_p[keep, , drop = FALSE]
    }
    if (length(dim(params@ft_pred_kernel_d)) == 2) {
        p@ft_pred_kernel_d <- params@ft_pred_kernel_d[keep, , drop = FALSE]
    }
    if (length(dim(params@pred_kernel)) == 3) {
        p@pred_kernel <- params@pred_kernel[keep, , , drop = FALSE]
    }
    p@ft_mask <- params@ft_mask[keep, , drop = FALSE]
    p@mu_b <- params@mu_b[keep, , drop = FALSE]
    p@ext_encounter <- params@ext_encounter[keep, , drop = FALSE]
    p@ext_diffusion <- params@ext_diffusion[keep, , drop = FALSE]
    p@species_params <- p@species_params[keep, , drop = FALSE]
    p@given_species_params <- p@given_species_params[keep, , drop = FALSE]
    p@interaction <- params@interaction[keep, keep, drop = FALSE]
    p@selectivity <- params@selectivity[, keep, , drop = FALSE]
    p@catchability <- params@catchability[, keep, drop = FALSE]
    p@w_min_idx <- params@w_min_idx[keep]
    p@A <- params@A[keep]
    p@gear_params <- p@gear_params[p@gear_params$species %in%
                                       p@species_params$species, ]
    p@gear_params <- validGearParams(p@gear_params, p@species_params,
                                     edition = defaults_edition(p))

    # Drop species parameters with all values NA
    keep <- colSums(is.na(p@species_params)) < nrow(p@species_params)
    p@species_params <- p@species_params[, keep]
    keep <- colSums(is.na(p@given_species_params)) < nrow(p@given_species_params)
    p@given_species_params <- p@given_species_params[, keep]

    # Preserve comments
    for (slot in (slotNames(p))) {
        comment(slot(p, slot)) <- comment(slot(params, slot))
    }

    validObject(p)

    p@time_modified <- lubridate::now()
    })

    return(p)
}

#' Let the external rates take over the contributions of removed species
#'
#' Called by [removeSpecies()] with `compensate = TRUE`, before it drops the
#' species. Adds to the external encounter rate, external mortality and
#' external diffusion of every species what the species to be removed
#' contributed to its encounter rate, mortality and diffusion at the initial
#' state, and rebalances the resource capacity, so that at that state the
#' remaining species and the resource have the rates they had before.
#'
#' The contributions are calculated on the full model, with the abundances of
#' the species to be removed set to zero, rather than on the model after the
#' removal. The two give the same rates, but an extension that keeps data of
#' its own for each species trims that data in its `removeSpecies()` method
#' only after `NextMethod()` returns, so its rate functions are consistent only
#' on the full model. The removed species keep their zero abundances, which
#' [removeSpecies()] then drops together with the rest of their rows.
#'
#' The encounter rate is compensated first, because the feeding level depends
#' on it and the mortality, the diffusion and the resource mortality depend on
#' the feeding level. With the encounter rate restored, the remaining predators
#' feed exactly as before and the only change left in those rates is what the
#' removed species contributed.
#'
#' The resource is rebalanced rather than given back its rate of change, so a
#' resource that was not at steady state is put at steady state. Keeping its
#' rate of change would need a formula for each resource dynamics (for the
#' semichemostat the capacity shifts by a fixed amount, for the logistic its
#' inverse does), whereas balancing goes through the `balance_<dynamics>()`
#' function that each resource dynamics provides anyway.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param remove A logical vector saying for each species whether it is to be
#'   removed.
#' @return The `MizerParams` object with the species still present, but with
#'   zero abundance, and with the compensated external rates and resource
#'   capacity.
#' @noRd
compensate_removal <- function(params, remove) {
    old_encounter <- getEncounter(params)
    old_mort <- getMort(params)
    old_diffusion <- getDiffusion(params)
    old_resource_mort <- getResourceMort(params)

    params@initial_n[remove, ] <- 0

    params <- compensate_rate(params, "ext_encounter", remove, old_encounter,
                              getEncounter)
    params <- compensate_rate(params, "mu_b", remove, old_mort, getMort)
    params <- compensate_rate(params, "ext_diffusion", remove, old_diffusion,
                              getDiffusion)

    # The removed species no longer eat the resource. `resource_constant()`
    # needs no balancing because it keeps the resource where it is anyway.
    if (any(old_resource_mort != getResourceMort(params)) &&
        params@resource_dynamics != "resource_constant") {
        if (is.function(get0(paste0("balance_", params@resource_dynamics)))) {
            # A frozen capacity would block the balancing, as in `steady()`
            comment(params@cc_pp) <- NULL
            params <- setResource(params)
        } else {
            signal_info("cc_pp", paste0(
                "There is no `balance_", params@resource_dynamics, "()` ",
                "function, so the resource capacity could not be rebalanced ",
                "after the removed species stopped eating the resource. The ",
                "resource is therefore not compensated for their removal."),
                level = 1, severity = "warning", unhandled = "show")
        }
    }

    dynamic <- names(params@other_dynamics)[
        !vapply(params@other_dynamics,
                function(f) identical(f, "constant_other"), logical(1))]
    if (length(dynamic) > 0) {
        signal_info("other_components", paste0(
            "The compensation for the removed species covers the remaining ",
            "species and the resource only. The component",
            if (length(dynamic) > 1) "s " else " ",
            paste0("`", dynamic, "`", collapse = ", "),
            if (length(dynamic) > 1) " have" else " has",
            " dynamics of ", if (length(dynamic) > 1) "their" else "its",
            " own and may still be affected by the removal."),
            level = 1, severity = "warning", unhandled = "show")
    }

    params
}

#' Restore a rate of the remaining species through its external rate
#'
#' Adds to the external rate array `slot` the amount by which the rate that
#' `get_rate()` calculates has fallen short of `target`, and checks that this
#' restores the rate. An array to which the removed species contributed
#' nothing is left alone, so that it is not marked as set manually needlessly.
#' That comparison is exact: a species that contributes nothing changes
#' neither the inputs of a rate nor, therefore, its numerical round-off.
#'
#' Adding the shortfall restores the rate whenever the external rate is added
#' to it unchanged, as in all of mizer's own rate functions. A custom rate
#' function need not do that; one that returns `2 * mizerEncounter(...)`
#' doubles the compensation as well. So where the rate is not restored, the
#' compensation is divided by the response of the rate to it. This is exact
#' for any rate that at each size is an affine function of the external rate,
#' which covers a rate function that scales mizer's rate. If the rate is still
#' not restored, the user is warned.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param slot One of `"ext_encounter"`, `"mu_b"` or `"ext_diffusion"`.
#' @param remove A logical vector saying for each species whether it is to be
#'   removed.
#' @param target The rate (species x size) before the removal.
#' @param get_rate The function calculating the rate from `params`.
#' @param tol The tolerance for a restored rate, relative to the largest value
#'   of the rate of that species.
#' @return The `MizerParams` object with the slot updated.
#' @noRd
compensate_rate <- function(params, slot, remove, target, get_rate,
                            tol = 1e-8) {
    target <- unclass_rate(target)
    reduced <- unclass_rate(get_rate(params))
    compensation <- target - reduced
    compensation[remove, ] <- 0
    if (all(compensation == 0)) {
        return(params)
    }
    # TRUE where the rate of a remaining species is not restored
    not_restored <- function(rate) {
        off <- abs(rate - target) > tol * apply(abs(target), 1, max)
        off[remove, ] <- FALSE
        off
    }
    compensated <- add_compensation(params, slot, compensation)
    rate <- unclass_rate(get_rate(compensated))
    off <- not_restored(rate)
    if (any(off)) {
        # Rescale everywhere, not only where the rate is off by more than the
        # tolerance, which is relative to the largest value of the rate and so
        # would leave small values off. The two evaluations of the rate differ
        # only in the external rate, so the response is measured to the
        # round-off in the rate, about eps * rate. Where the compensation is
        # too small to measure it to 1e-4, it is left as it is, which leaves an
        # error that is negligible compared to the rate.
        response <- (rate - reduced) / compensation
        fix <- abs(compensation) > 1e4 * .Machine$double.eps * abs(target) &
            is.finite(response) & response != 0
        compensation[fix] <- compensation[fix] / response[fix]
        compensated <- add_compensation(params, slot, compensation)
        rate <- unclass_rate(get_rate(compensated))
        off <- not_restored(rate)
    }
    if (any(off)) {
        quantity <- switch(slot, ext_encounter = "encounter rate",
                           mu_b = "mortality rate",
                           ext_diffusion = "diffusion rate")
        species <- params@species_params$species[rowSums(off) > 0]
        signal_info(slot, paste0(
            "The compensation for the removed species did not restore the ",
            quantity, " of ", paste(species, collapse = ", "), ". The ",
            quantity, " is calculated by a function that does not respond ",
            "to the external ", quantity, " in a way the compensation can ",
            "match, so the remaining species do not keep their rates."),
            level = 1, severity = "warning", unhandled = "show")
    }
    compensated
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

    validObject(params)

    params@time_modified <- lubridate::now()
    return(params)
}
#' Adjust the size grid
#'
#' This function adjusts the size grid in a [MizerParams] object to the desired
#' minimum and maximum size. It can both expand and truncate the grid.
#' If the grid is truncated, any data outside the new grid is discarded.
#' A warning is issued if there is non-negligible biomass in the discarded
#' size bins.
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
    p <- with_defaults_edition(defaults_edition(params), newMultispeciesParams(
        sp_params_for_grid,
        interaction = params@interaction,
        second_order_w = params@second_order_w,
        max_w = new_max_w,
        min_w = new_min_w,
        min_w_pp = new_min_w_pp * 10^(0.5 * dx),
        no_w = new_no_w,
        gear_params = params@gear_params,
        initial_effort = params@initial_effort,
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

    # Preserve other slots and metadata
    p@given_species_params <- new_given_sp_params
    p <- copyPreservedParamsSlots(p, params)
    p <- setColours(p, params@linecolour)
    p <- setLinetypes(p, params@linetype)

    # Preserve comments
    p <- copyParamsComments(p, params)

    p <- validParams(p)

    # Restore plot metadata exactly
    p@species_params <- new_sp_params
    p@given_species_params <- new_given_sp_params
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
