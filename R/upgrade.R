#' Determine whether a MizerParams or MizerSim object needs to be upgraded
#'
#' Looks at the mizer version that was used to last update the object and
#' returns TRUE if changes since that version require an upgrade of the object.
#' You would not usually have to call this function. Upgrades are initiated
#' automatically by `validParams` and `validSim` when necessary.
#'
#' @param object A MizerParams or MizerSim object
#' @return TRUE or FALSE
#' @concept helper
needs_upgrading <- function(object) {
    if (is(object, "MizerParams")) {
        params <- object
    } else if (is(object, "MizerSim")) {
        params <- object@params
    } else {
        stop("The object you supplied is neither a MizerParams nor a MizerSim object.")
    }
    mizer_needs_upgrading(params) || extension_needs_upgrading(params)
}

#' Whether the core mizer slots of an object need upgrading
#'
#' @param params A MizerParams object.
#' @return TRUE or FALSE.
#' @keywords internal
mizer_needs_upgrading <- function(params) {
    !.hasSlot(params, "mizer_version") ||
        params@mizer_version < "3.0.0.9003"
}

#' Whether any extension recorded on an object needs upgrading
#'
#' Returns TRUE if, for any extension recorded in the object's `@extensions`
#' slot whose package is installed, the recorded version stamp is missing (`NA`)
#' or older than the installed version of that package. A missing stamp counts
#' as needing an upgrade so that objects created before extension-version
#' tracking are brought up to date (and stamped) on first use; this is safe
#' because [runExtensionUpgrades()] only calls an upgrade method if one is registered and
#' such methods are written to be idempotent.
#'
#' @param params A MizerParams object.
#' @return TRUE or FALSE.
#' @keywords internal
extension_needs_upgrading <- function(params) {
    if (!.hasSlot(params, "extensions")) return(FALSE)
    reqs <- extensionRequirements(params@extensions)
    if (length(reqs) == 0) return(FALSE)
    vers <- extensionVersions(params@extensions)
    for (name in names(reqs)) {
        installed <- tryCatch(utils::packageVersion(name),
                              error = function(e) NULL)
        if (is.null(installed)) next  # extension package not installed
        stamp <- vers[[name]]
        if (is.na(stamp) || package_version(stamp) < installed) {
            return(TRUE)
        }
    }
    FALSE
}

#' Upgrade the core slots of a MizerParams object
#'
#' This is the `MizerParams` method of the [upgrade()] generic and performs the
#' *core mizer* upgrade only. It is called from the orchestrator
#' [runExtensionUpgrades()], which also invokes any registered extension upgrade
#' methods.
#' You should never need to call it directly; use `validParams()` (or
#' `readParams()`) instead.
#'
#' @param object An old MizerParams object to be upgraded
#' @param ... Unused.
#'
#' @return The upgraded MizerParams object
#' @concept helper
#' @keywords internal
#' @exportS3Method utils::upgrade
#' @seealso [validParams()]
upgrade.MizerParams <- function(object, ...) {
    params <- object
    if (!mizer_needs_upgrading(params)) return(params)
    original_params <- params

    # Preserve time_created
    if (.hasSlot(params, "time_created")) {
        time_created <- params@time_created
    } else {
        time_created <- lubridate::now()
    }

    # We'll use the version to decide which upgrades are needed. Copy it to
    # a variable because params@mizer_version might get changed during the
    # upgrade
    if (!.hasSlot(params, "mizer_version")) {
        version <- "2.2.1"
    } else {
        version <- params@mizer_version
    }

    # Before version 2.3 ----
    if (version < "2.4.1.9001") {

        # In these old versions `w_inf` was the maximum-size parameter and the
        # weight grid was built up to it. `w_inf` is now only the von Bertalanffy
        # asymptotic size and `w_max` (defaulting to `1.5 * w_inf`) is the
        # computational boundary. To preserve the original grid we copy `w_inf`
        # into `w_max` before the object is rebuilt below, so that the boundary
        # stays at the old value rather than being extended to `1.5 * w_inf`.
        if (!("w_max" %in% names(params@species_params)) &&
            "w_inf" %in% names(params@species_params)) {
            params@species_params$w_max <- params@species_params$w_inf
        }

        if ("interaction_p" %in% names(params@species_params)) {
            params@species_params$interaction_resource <-
                params@species_params$interaction_p
            params@species_params$interaction_p <- NULL
            message("The 'interaction_p' column has been renamed to 'interaction_resource'.")
        }

        if (!.hasSlot(params, "gear_params")) {
            gear_params <- validGearParams(data.frame(), params@species_params)
        } else {
            gear_params <- params@gear_params
        }

        if (!.hasSlot(params, "ft_pred_kernel_e")) {
            stop("Objects from versions 0.3 and earlier can not be upgraded.")
        }

        params@species_params$w_min_idx <- NULL

        if (.hasSlot(params, "srr")) {
            if (is.function(params@srr)) {
                if ("constant_recruitment" %in% names(params@species_params)) {
                    RDD <- "constantRDD"
                    params@species_params$constant_reproduction <-
                        params@species_params$constant_recruitment
                    params@species_params$constant_recruitment <- NULL
                } else {
                    RDD <- "BevertonHoltRDD"
                    message('The density-dependent reproduction rate function has been set to "BevertonHoltRDD".')
                }
            } else {
                RDD <- switch(params@srr,
                              srrBevertonHolt = "BevertonHoltRDD",
                              srrNone = "noRDD",
                              srrConstant = "constantRDD",
                              srrRicker = "RickerRDD",
                              srrSheperd = "SheperdRDD")
            }
        } else if (.hasSlot(params, "rates_funcs")) {
            RDD <- params@rates_funcs[["RDD"]]
        } else {
            RDD <- "BevertonHoltRDD"
        }

        if (.hasSlot(params, "initial_effort") &&
            length(params@initial_effort) == length(unique(gear_params$gear))) {
            initial_effort <- params@initial_effort
        } else {
            gears <- as.character(unique(gear_params$gear))
            initial_effort <- rep(0, length(gears))
            names(initial_effort) <- gears
            message("Initial effort has been set to 0.")
        }

        if (.hasSlot(params, "metab")) {
            metab <- params@metab
        } else {
            metab <- params@std_metab + params@activity
        }

        if (.hasSlot(params, "pred_kernel") &&
            length(dim(params@pred_kernel)) == 3) {
            pred_kernel <- params@pred_kernel
        } else {
            pred_kernel <- NULL
        }

        if (.hasSlot(params, "maturity")) {
            maturity <- params@maturity
            repro_prop <- params@psi / params@maturity
            repro_prop[params@maturity == 0] <- 0
            comment(repro_prop) <- comment(params@psi)
        } else {
            maturity <- NULL
            repro_prop <- NULL
        }

        if (.hasSlot(params, "mu_b")) {
            mu_b <- params@mu_b
        } else {
            mu_b <- NULL
        }

        if (.hasSlot(params, "ext_encounter")) {
            ext_encounter <- params@ext_encounter
        } else {
            ext_encounter <- NULL
        }

        if ("r_max" %in% names(params@species_params)) {
            params@species_params$R_max <- params@species_params$r_max
            params@species_params$r_max <- NULL
            message("The 'r_max' column has been renamed to 'R_max'.")
        }

        if (.hasSlot(params, "p")) {
            params@species_params[["p"]] <- params@p
        } else if (!("p" %in% names(params@species_params))) {
            # No p in params object, so extract from metabolism
            p <- log(metab[, 2] / metab[, 1]) /
                log(params@w[[2]] / params@w[[1]])
            p[is.nan(p)] <- NA
            params@species_params[["p"]] <- p
        }
        if (.hasSlot(params, "q")) {
            params@species_params[["q"]] <- params@q
        } else if (!("q" %in% names(params@species_params))) {
            # No q in params object, so extract from search volume
            q <- log(params@search_vol[, 2] / params@search_vol[, 1]) /
                log(params@w[[2]] / params@w[[1]])
            params@species_params[["q"]] <- q
        }
        if (.hasSlot(params, "n")) {
            params@species_params[["n"]] <- params@n
        } else if (!("n" %in% names(params@species_params))) {
            # No n in params object, so extract from intake_max
            n <- log(params@intake_max[, 2] / params@intake_max[, 1]) /
                log(params@w[[2]] / params@w[[1]])
            n[is.nan(n)] <- NA
            params@species_params[["n"]] <- n
        }
        if (.hasSlot(params, "f0")) {
            params@species_params[["f0"]] <- params@f0
        }

        params@species_params$species <- as.character(params@species_params$species)

        pnew <- newMultispeciesParams(
            params@species_params,
            interaction = params@interaction,
            no_w = length(params@w),
            min_w = params@w[1],
            max_w = params@w[length(params@w)],
            min_w_pp = params@w_full[1] + 1e-16, # To make
            # sure that we don't add an extra bracket.
            resource_rate = params@rr_pp,
            resource_capacity = params@cc_pp,
            pred_kernel = pred_kernel,
            search_vol = params@search_vol,
            intake_max = params@intake_max,
            metab = metab,
            ext_mort = mu_b,
            ext_encounter = ext_encounter,
            maturity = maturity,
            repro_prop = repro_prop,
            RDD = RDD,
            gear_params = gear_params,
            initial_effort = initial_effort)

        pnew@psi <- params@psi

        if (.hasSlot(params, "linecolour")) {
            names(params@linecolour)[names(params@linecolour) == "Plankton"] <- "Resource"
            names(params@linetype)[names(params@linetype) == "Plankton"] <- "Resource"
            pnew@linecolour <- params@linecolour
            pnew@linetype <- params@linetype
        }
        if (.hasSlot(params, "initial_n")) {
            pnew@initial_n <- params@initial_n
            pnew@initial_n_pp <- params@initial_n_pp
        }
        if (.hasSlot(params, "diffusion")) {
            pnew@ext_diffusion[] <- slot(params, "diffusion")
        } else if (.hasSlot(params, "ext_diffusion")) {
            pnew@ext_diffusion <- params@ext_diffusion
        }
        if (.hasSlot(params, "initial_n_other")) {
            pnew@initial_n_other <- params@initial_n_other
        }

        if (.hasSlot(params, "sc")) {
            pnew@sc <- params@sc
        }
        if (.hasSlot(params, "other_dynamics")) {
            pnew@other_dynamics <- params@other_dynamics
            pnew@other_params <- params@other_params
        }
        if (.hasSlot(params, "other_encounter")) {
            pnew@other_encounter <- params@other_encounter
        }
        if (.hasSlot(params, "other_pred_mort")) {
            pnew@other_mort <- params@other_pred_mort
        } else if (.hasSlot(params, "other_mort")) {
            pnew@other_mort <- params@other_mort
        }
        if (.hasSlot(params, "plankton_dynamics")) {
            if (is.function(params@plankton_dynamics)) {
                pnew@resource_dynamics <- "resource_semichemostat"
                message('The resource dynamics function has been set to "resource_semichemostat".')
            } else {
                if (params@plankton_dynamics == "plankton_semichemostat") {
                    pnew@resource_dynamics <- "resource_semichemostat"
                } else if (params@plankton_dynamics == "plankton_constant") {
                    pnew@resource_dynamics <- "resource_constant"
                } else {
                    pnew@resource_dynamics <- params@plankton_dynamics
                }
            }
        } else if (.hasSlot(params, "resource_dynamics")) {
            pnew@resource_dynamics <- params@resource_dynamics
        }
        if (.hasSlot(params, "plankton_params")) {
            pnew@resource_params <- params@plankton_params
        }
        if (.hasSlot(params, "resource_params")) {
            pnew@resource_params <- params@resource_params
        } else if (.hasSlot(params, "lambda")) {
            # r_pp was not stored in params object, so has to be reconstructed
            maxidx <- max(which(pnew@cc_pp > 0))  # The largest resource index
            r_pp <- params@rr_pp[[maxidx]] / params@w_full[[maxidx]] ^ (params@n - 1)
            pnew@resource_params[["r_pp"]] <- r_pp
            pnew@resource_params[["lambda"]] <- params@lambda
            pnew@resource_params[["kappa"]] <- params@kappa
            pnew@resource_params[["n"]] <- params@n
            pnew@resource_params[["w_pp_cutoff"]] <- max(pnew@w_full[pnew@cc_pp > 0]) +
                1e-10 # to make sure we don't remove a bucket by mistake
        } else {
            # No resource parameters saved in params object, so need to
            # reconstruct from cc_pp and rr_pp
            minidx <- min(which(pnew@cc_pp > 0))  # The smallest resource index
            maxidx <- max(which(pnew@cc_pp > 0))  # The largest resource index
            lambda <- -log(params@cc_pp[[maxidx]] / params@cc_pp[[minidx]]) /
                log(params@w_full[[maxidx]] / params@w_full[[minidx]])
            kappa <- params@cc_pp[[maxidx]] * params@w_full[[maxidx]] ^ lambda
            n <- 1 + log(params@rr_pp[[maxidx]] / params@rr_pp[[minidx]]) /
                log(params@w_full[[maxidx]] / params@w_full[[minidx]])
            r_pp <- params@rr_pp[[maxidx]] / params@w_full[[maxidx]] ^ (n - 1)
            pnew@resource_params[["r_pp"]] <- r_pp
            pnew@resource_params[["lambda"]] <- lambda
            pnew@resource_params[["kappa"]] <- kappa
            pnew@resource_params[["n"]] <- n
            pnew@resource_params[["w_pp_cutoff"]] <- pnew@w_full[[maxidx]]
        }

        if (.hasSlot(params, "A")) {
            pnew@A <- params@A
            # Derive is_background from the old @A slot if not already present
            if (!("is_background" %in% names(pnew@species_params))) {
                pnew@species_params$is_background <- is.na(params@A)
            }
        }

        if (.hasSlot(params, "metadata")) {
            pnew@metadata <- params@metadata
            pnew@time_created <- params@time_created
            pnew@extensions <- params@extensions
        }

        if (.hasSlot(params, "rates_funcs")) {
            pnew@rates_funcs <- params@rates_funcs
        }
        if (.hasSlot(params, "given_species_params")) {
            pnew@given_species_params <- params@given_species_params
        }
        if (.hasSlot(params, "selectivity")) {
            pnew@selectivity <- params@selectivity
        }
        if (.hasSlot(params, "catchability")) {
            pnew@catchability <- params@catchability
        }
        if (.hasSlot(params, "use_predation_diffusion")) {
            pnew@use_predation_diffusion <- params@use_predation_diffusion
        }

        # renaming catch_observed to yield_observed in mizer 2.3
        if ("catch_observed" %in% names(pnew@species_params)) {
            pnew@species_params$yield_observed <-
                pnew@species_params$catch_observed
            pnew@species_params$catch_observed <- NULL
        }

        # Copy over all comments
        comment(pnew) <- comment(params)
        for (slot in slotNames(pnew)) {
            if (.hasSlot(params, slot)) {
                comment(slot(pnew, slot)) <- comment(slot(params, slot))
            }
        }

        params <- pnew
    }
    if (!.hasSlot(params, "ext_diffusion")) {
        mat1 <- array(0, dim = c(nrow(params@species_params), length(params@w)),
                      dimnames = list(sp = params@species_params$species,
                                      w = signif(params@w, 3)))
        if (.hasSlot(params, "diffusion")) {
            params@ext_diffusion <- slot(params, "diffusion")
        } else {
            params@ext_diffusion <- mat1
        }
    }

    # Before version 2.4 ----
    if (version < "2.3.1.9001") {
        # copy w_inf to w_max
        par_names <- names(params@species_params)
        if (!("w_max" %in% par_names)) {
            params@species_params$w_max <- params@species_params$w_inf
        }
        # Add linecolour and linetype for fishing and external mortality
        if (!"Fishing" %in% names(getColours(params))) {
            params <- setColours(params, c("Fishing" = "red"))
        }
        if (!"Fishing" %in% names(getLinetypes(params))) {
            params <- setLinetypes(params, c("Fishing" = "solid"))
        }
        if (!"External" %in% names(getColours(params))) {
            params <- setColours(params, c("External" = "grey"))
        }
        if (!"External" %in% names(getLinetypes(params))) {
            params <- setLinetypes(params, c("External" = "solid"))
        }
    }

    # For version since 2.4.1.9002 ----
    if (!.hasSlot(params, "given_species_params")) {
        params@given_species_params <- params@species_params
    }
    if (!("w_mat25" %in% names(params@species_params))) {
        params <- set_species_param_default(params, "w_mat25",
            params@species_params$w_mat / (3 ^ (1 / 10)))
    }

    # Add w_inf if missing (added in 3.0.0.9003). `w_inf`, the von Bertalanffy
    # asymptotic size, is now the required maximum-size parameter and is used to
    # set defaults for `w_repro_max`, `w_mat` and `z0`. For objects that predate
    # this change we derive it from `w_repro_max` (which played the role of the
    # growth-stopping size in earlier versions) or otherwise from `w_max`. The
    # stored `w_max` is left untouched so that the size grid is preserved.
    if (!("w_inf" %in% names(params@species_params))) {
        if ("w_repro_max" %in% names(params@species_params)) {
            params@species_params$w_inf <- params@species_params$w_repro_max
        } else {
            params@species_params$w_inf <- params@species_params$w_max
        }
    }

    # Derive is_background from the old @A slot if not already present
    if (!("is_background" %in% names(params@species_params))) {
        params@species_params$is_background <- is.na(params@A)
    }

    # Add Diffusion rate function if missing (added in 2.5.4.9122)
    if (is.null(params@rates_funcs[["Diffusion"]])) {
        params@rates_funcs[["Diffusion"]] <- "mizerDiffusion"
    }

    # Add use_predation_diffusion slot if missing (added in 2.5.4.9xxx)
    # Default to FALSE to preserve behaviour of previous mizer versions.
    if (!.hasSlot(params, "use_predation_diffusion")) {
        params@use_predation_diffusion <- FALSE
    }

    # Add/upgrade second_order_w slot (added in 3.0.x, changed to list in 3.0.x)
    # Default to first-order settings to preserve behaviour of previous mizer versions.
    if (!.hasSlot(params, "second_order_w")) {
        params@second_order_w <- list(flux = "upwind", bin_average = FALSE)
    }

    # Add ft_pred_kernel_d slot if missing (added in 3.0.0.9002). It mirrors
    # ft_pred_kernel_e except under second_order_w bin-averaging; objects from
    # earlier versions are first order, so the diffusion kernel equals the
    # encounter kernel. (A later setParams() rebuilds it correctly if needed.)
    if (!.hasSlot(params, "ft_pred_kernel_d") ||
        !identical(dim(params@ft_pred_kernel_d),
                   dim(params@ft_pred_kernel_e))) {
        params@ft_pred_kernel_d <- params@ft_pred_kernel_e
    }

    params@mizer_version <- packageVersion("mizer")
    params <- validParams(params, info_level = 0)
    params@time_modified <- lubridate::now()
    params@time_created <- time_created
    reinstateParamsClass(params, original_params)
}

#' Back-compatible wrapper for the core MizerParams upgrade
#'
#' Retained because it was an (unexported) internal entry point. New code should
#' rely on `validParams()` / `readParams()`, which orchestrate the core upgrade
#' together with any extension upgrades.
#'
#' @param params An old MizerParams object.
#' @return The upgraded MizerParams object.
#' @keywords internal
upgradeParams <- function(params) {
    upgrade.MizerParams(params)
}

#' Run the registered extension upgrade methods on an object
#'
#' For each extension recorded in the object's `@extensions` slot (processed
#' innermost-first) whose installed package version is newer than the recorded
#' stamp (or whose stamp is missing), calls the extension's `upgrade` method if
#' one is registered, then records the installed version as the new stamp. The
#' core mizer upgrade is *not* run here; see [upgrade.MizerParams()].
#'
#' Extension upgrade methods are looked up with
#' `getS3method("upgrade", name, optional = TRUE)`, must perform only their own
#' migration, must be idempotent, and must not call `NextMethod()`.
#'
#' @param params A MizerParams object.
#' @return The object with extension migrations applied and stamps refreshed.
#' @keywords internal
runExtensionUpgrades <- function(params) {
    reqs <- extensionRequirements(params@extensions)
    if (length(reqs) == 0) return(params)
    vers <- extensionVersions(params@extensions)
    # innermost extension first (the chain is stored outermost-first)
    for (name in rev(names(reqs))) {
        installed <- tryCatch(utils::packageVersion(name),
                              error = function(e) NULL)
        if (is.null(installed)) next  # extension package not installed
        stamp <- vers[[name]]
        if (!(is.na(stamp) || package_version(stamp) < installed)) next
        method <- getS3method("upgrade", name, optional = TRUE)
        if (!is.null(method)) {
            params <- method(params)
        }
        params <- recordExtension(params, name,
                                  version = as.character(installed))
    }
    params
}

reinstateParamsClass <- function(params, original_params) {
    original_class <- class(original_params)[[1]]
    if (identical(original_class, "MizerParams")) {
        return(params)
    }

    params <- as(params, original_class)
    extra_slots <- setdiff(slotNames(original_params), slotNames("MizerParams"))
    for (slot in extra_slots) {
        if (.hasSlot(params, slot)) {
            slot(params, slot) <- slot(original_params, slot)
        }
    }
    validObject(params)
    params
}

#' Upgrade a MizerSim object from earlier versions
#'
#' This is the `MizerSim` method of the [upgrade()] generic. It rebuilds the
#' simulation around an upgraded params object; the params upgrade (core mizer
#' and any extensions) is performed by the `validParams()` call below. You
#' should never need to call it directly; use `validSim()` (or `readSim()`).
#'
#' @param object An old MizerSim object to be upgraded
#' @param ... Unused.
#'
#' @return The upgraded MizerSim object
#' @concept helper
#' @keywords internal
#' @exportS3Method utils::upgrade
upgrade.MizerSim <- function(object, ...) {
    sim <- object
    assert_that(is(sim, "MizerSim"))
    if (!needs_upgrading(sim)) return(sim)

    t_dimnames <- dimnames(sim@n_pp)[[1]]
    no_gears <- dim(sim@effort)[[2]]
    new_sim <- MizerSim(validParams(sim@params, info_level = 0),
                        t_dimnames = as.numeric(t_dimnames))

    new_sim@n <- sim@n
    new_sim@n_pp <- sim@n_pp
    if (dim(sim@effort)[[1]] < dim(sim@n)[[1]]) {
        # This happened in version 0.4 because there was no effort for initial
        # time
        new_sim@effort[] <- rbind(rep(0, no_gears), sim@effort)
    } else {
        new_sim@effort[] <- sim@effort
    }
    if (.hasSlot(sim, "n_other")) {
        new_sim@n_other <- sim@n_other
    }
    if (.hasSlot(sim, "sim_params")) {
        new_sim@sim_params <- sim@sim_params
    }
    comment(new_sim) <- comment(sim)
    comment(new_sim@effort) <- comment(sim@effort)
    validObject(new_sim)
    new_sim
}

#' Back-compatible wrapper for the MizerSim upgrade
#'
#' Retained because it was an (unexported) internal entry point. New code should
#' rely on `validSim()` / `readSim()`.
#'
#' @param sim An old MizerSim object.
#' @return The upgraded MizerSim object.
#' @keywords internal
upgradeSim <- function(sim) {
    upgrade.MizerSim(sim)
}
