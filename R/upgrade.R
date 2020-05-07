#' Upgrade MizerParams object from earlier mizer versions
#' 
#' Occasionally during the development of new features for mizer, the
#' \linkS4class{MizerParams} object gains extra slots. MizerParams objects
#' created in older versions of mizer are then no longer valid in the new
#' version because of the missing slots. You need to upgrade them with
#' ```
#'params <- upgradeParams(params)
#' ```
#' where `params` should be replaced by the name of your MizerParams object.
#' This function adds the missing slots and fills them with default values. Any
#' object from version 0.4 onwards can be upgraded. Any old
#' \linkS4class{MizerSim} objects should be similarly updated with
#' [upgradeSim()]. This function uses [newMultispeciesParams()] to create a new
#' MizerParams object using the parameters extracted from the old MizerParams
#' object.
#' 
#' @section Backwards compatibility:
#' The internal numerics in mizer have changed over time, so there may be small
#' discrepancies between the results obtained with the upgraded object
#' in the new version and the original object in the old version. If it
#' is important for you to reproduce the exact results then you should install
#' the version of mizer with which you obtained the results. You can do this
#' with
#' ```
#' remotes::install_github("sizespectrum/mizer", ref = "v0.2")
#' ```
#' where you should replace "v0.2" with the version number you require. You can
#' see the list of available releases at 
#' <https://github.com/sizespectrum/mizer/tags>.
#' 
#' If you only have a serialised version of the old object, for example
#' created via [saveRDS()], and you get an error when trying to read it in
#' with [readRDS()] then unfortunately you will need to install the old version
#' of mizer first to read the params object into your workspace, then switch
#' to the current version and then call [upgradeParams()]. You can then save
#' the new version again with [saveRDS()].
#' 
#' @param params An old MizerParams object to be upgraded
#' 
#' @return The upgraded MizerParams object
#' @export
upgradeParams <- function(params) {
    
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
            if (!is.null(params@species_params$constant_recruitment)) {
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
    
    if (.hasSlot(params, "initial_effort")) {
        initial_effort <- params@initial_effort
    } else {
        initial_effort <- NULL
    }
    
    if (.hasSlot(params, "metab")) {
        metab <- params@metab
    } else {
        metab <- params@std_metab + params@activity
    }
    
    if (.hasSlot(params, "pred_kernel") && 
        length(dim(params@pred_kernel)) == 3) {
        pred_kernel <- params@pred_kernel
    } else pred_kernel <- NULL
    
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
    
    if ("r_max" %in% names(params@species_params)) {
        params@species_params$R_max <- params@species_params$r_max
        params@species_params$r_max <- NULL
        message("The 'r_max' column has been renamed to 'R_max'.")
    }
    
    if (.hasSlot(params, "p")) {
        params@species_params[["p"]] <- params@p
    } else if (is.null(params@species_params[["p"]])) {
        # No p in params object, so extract from metabolism
        p <- log(metab[, 2] / metab[, 1]) / 
            log(params@w[[2]] / params@w[[1]])
        p[is.nan(p)] <- NA
        params@species_params[["p"]] <- p
    }
    if (.hasSlot(params, "q")) {
        params@species_params[["q"]] <- params@q
    } else if (is.null(params@species_params[["q"]])) {
        # No q in params object, so extract from search volume
        q <- log(params@search_vol[, 2] / params@search_vol[, 1]) / 
            log(params@w[[2]] / params@w[[1]])
        params@species_params[["q"]] <- q
    }
    if (.hasSlot(params, "n")) {
        params@species_params[["n"]] <- params@n
    } else if (is.null(params@species_params[["n"]])) {
        # No n in params object, so extract from intake_max
        n <- log(params@intake_max[, 2] / params@intake_max[, 1]) / 
            log(params@w[[2]] / params@w[[1]])
        n[is.nan(n)] <- NA
        params@species_params[["n"]] <- n
    }
    if (.hasSlot(params, "f0")) {
        params@species_params[["f0"]] <- params@f0
    }
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
        z0 = mu_b,
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
        pnew@resource_params[["w_pp_cutoff"]] <- max(pnew@w_full[pnew@cc_pp > 0])
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
    }
    
    # Copy over all comments
    comment(pnew) <- comment(params)
    for (slot in slotNames(pnew)) {
        if (.hasSlot(params, slot)) {
            comment(slot(pnew, slot)) <- comment(slot(params, slot))
        }
    }
    
    return(pnew)
}

#' Upgrade MizerSim object from earlier mizer versions
#' 
#' Occasionally, during the development of new features for mizer, the
#' \linkS4class{MizerSim} class or the \linkS4class{MizerParams} class gains
#' extra slots. MizerSim objects created in older versions of mizer are then no
#' longer valid in the new version because of the missing slots. You need to
#' upgrade them with
#' ```
#' sim <- upgradeSim(sim)
#' ```
#' where `sim` should be replaced by the name of your MizerSim object.
#' 
#' This function adds the missing slots and fills them with default values. It
#' calls [upgradeParams()] to upgrade the MizerParams object inside the MizerSim
#' object. Any object from version 0.4 onwards can be upgraded.
#' 
#' @inheritSection upgradeParams Backwards compatibility
#' 
#' @param sim An old MizerSim object to be upgraded
#' 
#' @return The upgraded MizerSim object
#' @export
upgradeSim <- function(sim) {
    t_dimnames <- dimnames(sim@n_pp)[[1]]
    no_gears <- dim(sim@effort)[[2]]
    new_sim <- MizerSim(upgradeParams(sim@params),
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
    comment(new_sim) <- comment(sim)
    comment(new_sim@effort) <- comment(sim@effort)
    new_sim
}
