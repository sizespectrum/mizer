#' Species parameters
#' 
#' These functions allow you to get or set the species parameters stored in
#' a MizerParams object.
#' 
#' The `species_params` data frame holds species-specific parameters. The data
#' frame has one row for each species and one column for each species parameter.
#' There are a lot of species parameters and we will list them all below, but
#' most of them have sensible default values. The only required columns are
#' `species` for the species name and `w_inf` for its asymptotic size. However
#' if you have information about the values of other parameters then you should
#' include them in the `species_params` data frame.
#' 
#' There are some species parameters that are used to set up the
#' size-dependent parameters that are used in the mizer model:
#' 
#' * `gamma` and `q` are used to set the search volume, see [setSearchVolume()].
#' * `h` and `n` are used to set the maximum intake rate, see [setMaxIntakeRate()].
#' * `k`, `ks` and `p` are used to set activity and basic metabolic rate, 
#'   see [setMetabolicRate()].
#' * `z0` is used to set the external mortality rate, see [setExtMort()].
#' * `w_mat`, `w_mat25`, `w_inf` and `m` are used to set the allocation to
#'   reproduction, see [setReproduction()].
#' * `pred_kernel_type` specifies the shape of the predation kernel. The default
#'   is a "lognormal", for other options see the "Setting predation kernel"
#'   section in the help for [setPredKernel()].
#' * `beta` and `sigma` are parameters of the lognormal predation kernel, see
#'   [lognormal_pred_kernel()]. There will be other parameters if you are 
#'   using other predation kernel functions.
#'   
#' When you change one of the above species parameters in an already existing
#' MizerParams object using [species_params<-()], the new value will be
#' used to update the corresponding size-dependent rates automatically, unless
#' you have set those size-dependent rates manually, in which case the
#' corresponding species parameters will be ignored.
#' 
#' There are some species parameters that are used directly in the model
#' rather than being used for setting up size-dependent parameters:
#' 
#' * `alpha` is the assimilation efficiency, the proportion of the consumed
#'   biomass that can be used for growth, metabolism and reproduction, see
#'   the help for [getEReproAndGrowth()].
#' * `w_min` is the egg size.
#' * `interaction_resource` sets the interaction strength with the resource,
#'   see "Predation encounter" section in the help for [getEncounter()].
#' * `erepro` is the reproductive efficiency, the proportion of the energy
#'   invested into reproduction that is converted to egg biomass, see
#'   [getRDI()].
#' * `Rmax` is the parameter in the Beverton-Holt density dependence added to
#'   the reproduction, see [setBevertonHolt()]. There will be other such
#'   parameters if you use other density dependence functions, see the
#'   "Density dependence" section in the help for [setReproduction()].
#'
#' Two parameters are used only by functions that need to convert between
#' weight and length:
#' 
#' * `a` and `b` are the parameters in the allometric weight-length
#'   relationship \eqn{w = a l ^ b}.
#'   
#' Not all of species parameters have to be specified by the user. If they are
#' missing, [newMultispeciesParams()] will give them default values, sometimes
#' by using other species parameters. The parameters that are only used to
#' calculate default values for other parameters are:
#' 
#' * `k_vb` and `t0` are the von Bertalanffy growth parameters and are used
#'   together with the length-weight relationship exponent `b` and the egg
#'   size `w_min` to
#'   get a default value for the coefficient of the maximum intake rate `h`, 
#'   see [get_h_default()].
#' * `f0` is the feeding level and is used to get a default value for the
#'   coefficient of the search volume `gamma`, see [get_gamma_default()].
#' * `fc` is the critical feeding level below which the species can not 
#'   maintain itself. This is used to get a default value for the coefficient
#'   of the metabolic rate `ks`, see [get_ks_default()].
#'   
#' Note that these parameters are ignored if the parameters for which they allow
#' defaults to be calculated have instead been set explicitly. Also, these
#' parameters will only be used when setting up a new model with 
#' [newMultispeciesParams()]. Changing them later will have no effect 
#' because the default for the other parameters will not be recalculated.
#' 
#' There are other species parameters that are used in tuning the model to
#' observations:
#' 
#' * `biomass_observed` and `biomass_cutoff` allow you to specify for each
#'   species the total observed biomass above some cutoff size. This is
#'   used by [calibrateBiomass()] and [matchBiomasses()].
#' * `yield_observed` allows you to specify for each
#'   species the total annual fisheries yield. This is
#'   used by [calibrateYield()] and [matchYields()].
#' 
#' Finally there are two species parameters that control the way the species are
#' represented in plots:
#'
#' * `linecolour` specifies the colour and can be any valid R colour value.
#' * `linetype` specifies the line type ("solid", "dashed", "dotted", "dotdash",
#'    "longdash", "twodash" or "blank") 
#' 
#' Other species-specific information that is related to how the species is
#' fished is specified in a gear parameter data frame, see [gear_params()].
#' However in the case where each species is caught by only a single gear, 
#' this information can also optionally be provided as columns in the
#' `species_params` data frame and [newMultispeciesParams()] will transfer
#' them to the `gear_params` data frame. However changing these parameters later
#' in the species parameter data frame will have no effect.
#' 
#' You are allowed to include additional columns in the `species_params`
#' data frame. They will simply be ignored by mizer but will be stored in the
#' MizerParams object, in case your own code makes use of them.
#' 
#' @param params A MizerParams object
#' @export
#' @seealso [validSpeciesParams()]
#' @family functions for setting parameters
species_params <- function(params) {
    params@species_params
}

#' @rdname species_params
#' @param value A data frame with the species parameters
#' @export
`species_params<-` <- function(params, value) {
    value <- validSpeciesParams(value)
    params@species_params <- value
    suppressMessages(setParams(params))
}


#' Set a species parameter to a default value
#'
#' If the species parameter does not yet exist in the species parameter data
#' frame, then create it and fill it with the default. Otherwise use the default
#' only to fill in any NAs. Optionally gives a message if the parameter
#' did not already exist.
#' @param object Either a MizerParams object or a species parameter data frame
#' @param parname A string with the name of the species parameter to set
#' @param default A single default value or a vector with one default value for
#'   each species
#' @param message A string with a message to be issued when the parameter did
#'   not already exist
#' @return The `object` with an updated column in the species params data frame.
#' @export
#' @concept helper
set_species_param_default <- function(object, parname, default,
                                      message = NULL) {
    if (is(object, "MizerParams")) {
        species_params <- object@species_params
    } else {
        species_params <- object
    }
    assert_that(is.data.frame(species_params))
    assert_that(is.string(parname))
    no_sp <- nrow(species_params)
    if (length(default) == 1) {
        default <- rep(default, no_sp)
    }
    assert_that(length(default) == no_sp)
    if (!(parname %in% colnames(species_params))) {
        if (!missing(message)) {
            signal(message,
                    class = "info_about_default", var = parname, level = 3)
        }
        species_params[parname] <- default
    } else {
        # We do not like factors
        if (is.factor(species_params[[parname]])) {
            species_params[[parname]] <- as.character(species_params[[parname]])
        }
        missing <- is.na(species_params[[parname]])
        if (any(missing)) {
            species_params[missing, parname] <- default[missing]
        }
    }
    if  (is(object, "MizerParams")) {
        object@species_params <- species_params
        return(object)
    } else {
        return(species_params)
    }
}



#' Get default value for h
#' 
#' Sets `h` so that the species reaches maturity 
#' size at the age predicted by the von Bertalanffy growth curve parameters
#' `k_vb` and (optionally `t0`) taken from the species parameter
#' data frame. Also needs the exponent `b` from the length-weight
#' relationship \eqn{w = a l^b}. If this is not present in the species
#' parameter data frame it is set to \code{b = 3}.
#' @param params A MizerParams object
#' @return A vector with the values of h for all species
#' @export
#' @keywords internal
#' @concept helper
#' @family functions calculating defaults
get_h_default <- function(params) {
    assert_that(is(params, "MizerParams"))
    species_params <- params@species_params %>%
        set_species_param_default("f0", 0.6)
    if (!("h" %in% colnames(species_params))) {
        species_params$h <- rep(NA, nrow(species_params))
    }
    missing <- is.na(species_params$h)
    if (any(missing)) {
        assert_that(is.numeric(species_params$f0),
                    noNA(species_params$alpha),
                    "alpha" %in% names(species_params))
        signal("No h provided for some species, so using f0 and k_vb to calculate it.",
               class = "info_about_default", var = "h", level = 3)
        if (!("k_vb" %in% colnames(species_params))) {
            stop("Except I can't because there is no k_vb column in the species data frame")
        }
        if (any(is.na(species_params$k_vb[missing]))) {
            stop("Can not calculate defaults for h because some k_vb values are NA.")
        }
        if (!isTRUE(all.equal(species_params$n[missing], species_params$p[missing],
                              check.attributes = FALSE))) {
            signal("Because you have n != p, the default value for `h` is not very good.",
                   class = "info_about_default", var = "h", level = 1)
        }
        species_params <- species_params %>% 
            set_species_param_default("b", 3) %>% 
            set_species_param_default("t0", 0) %>% 
            set_species_param_default("fc", 0.2)
        w_mat <- species_params$w_mat
        w_inf <- species_params$w_inf
        w_min <- species_params$w_min
        b <- species_params$b
        k_vb <- species_params$k_vb
        n <- species_params$n
        age_mat <- -log(1 - (w_mat/w_inf)^(1/b)) / k_vb + species_params$t0
        h <- (w_mat^(1 - n) - w_min^(1 - n)) / age_mat / (1 - n) / 
            params@species_params$alpha / (species_params$f0 - species_params$fc)
        
        if (any(is.na(h[missing])) || any(h[missing] <= 0)) {
            stop("Could not calculate h.")
        }
        species_params$h[missing] <- h[missing]
    }
    return(species_params$h)
}


#' Get default value for gamma
#' 
#' Fills in any missing values for gamma so that fish feeding on a resource
#' spectrum described by the power law \eqn{\kappa w^{-\lambda}} achieve a
#' feeding level \eqn{f_0}. Only for internal use.
#' 
#' @param params A MizerParams object
#' @return A vector with the values of gamma for all species
#' @export
#' @concept helper
#' @family functions calculating defaults
get_gamma_default <- function(params) {
    assert_that(is(params, "MizerParams"))
    species_params <- params@species_params %>%
        set_species_param_default("f0", 0.6)
    if (!("gamma" %in% colnames(species_params))) {
        species_params$gamma <- rep(NA, nrow(species_params))
    }
    missing <- is.na(species_params$gamma)
    if (any(missing)) {
        assert_that(is.number(params@resource_params$lambda),
                    is.number(params@resource_params$kappa),
                    is.numeric(species_params$f0))
        signal("Using f0, h, lambda, kappa and the predation kernel to calculate gamma.",
                class = "info_about_default", var = "h", level = 1)
        if (!"h" %in% names(params@species_params) || 
            any(is.na(species_params$h))) {
            species_params$h <- get_h_default(params)
        }
        # Calculate available energy by setting search_volume
        # coefficient to 1
        params@species_params$gamma <- 1
        params <- setSearchVolume(params)
        # and setting a power-law prey spectrum
        params@initial_n[] <- 0
        if (defaults_edition() < 2) {
            # See issue #238
            params@species_params$interaction_resource <- 1
        }
        params@initial_n_pp[] <- params@resource_params$kappa * 
            params@w_full^(-params@resource_params$lambda)
        avail_energy <- getEncounter(params)[, length(params@w)] /
            params@w[length(params@w)] ^ 
            (2 + params@species_params$q - params@resource_params$lambda)
        # Now set gamma so that this available energy leads to f0
        gamma_default <- (species_params$h / avail_energy) * 
            (species_params$f0 / (1 - species_params$f0))
        # Only overwrite missing gammas with calculated values
        if (any(is.na(gamma_default[missing]))) {
            stop("Could not calculate gamma.")
        }
        species_params$gamma[missing] <- gamma_default[missing]
    }
    return(species_params$gamma)
}

#' Get default value for f0
#' 
#' Fills in any missing values for f0 so that if the prey abundance was
#' described by the power law \eqn{\kappa w^{-\lambda}} then the encounter rate
#' coming from the given `gamma` parameter would lead to the feeding level
#' \eqn{f_0}. This is thus doing the inverse of [get_gamma_default()].
#' Only for internal use.
#' 
#' For species for which no value for `gamma` is specified in the species
#' parameter data frame, the `f0` values is kept as provided in the species
#' parameter data frame or it is set to 0.6 if it is not provided.
#' 
#' @param params A MizerParams object
#' @return A vector with the values of f0 for all species
#' @export
#' @concept helper
#' @family functions calculating defaults
get_f0_default <- function(params) {
    assert_that(is(params, "MizerParams"))
    species_params <- params@species_params %>%
        set_species_param_default("f0", 0.6)
    if (!("gamma" %in% colnames(species_params))) {
        species_params$gamma <- rep(NA, nrow(species_params))
    }
    given <- !is.na(species_params$gamma)
    if (any(given)) {
        assert_that(is.number(params@resource_params$lambda),
                    is.number(params@resource_params$kappa),
                    is.numeric(species_params$gamma))
        if (!"h" %in% names(params@species_params) || 
            any(is.na(species_params$h))) {
            species_params$h <- get_h_default(params)
        }
        # Calculate available energy by setting a power-law prey spectrum
        params@initial_n[] <- 0
        params@species_params$interaction_resource <- 1
        params@initial_n_pp[] <- params@resource_params$kappa * 
            params@w_full^(-params@resource_params$lambda)
        avail_energy <- getEncounter(params)[, length(params@w)] /
            params@w[length(params@w)] ^ 
            (2 + params@species_params$q - params@resource_params$lambda)
        # Now set f0 so that this available energy leads to f0
        f0_default <- 1 / (species_params$h / avail_energy + 1)
        if (any(is.na(f0_default[given]))) {
            stop("Could not calculate f0.")
        }
        # Only overwrite f0 for species where gamma was given
        species_params$f0[given] <- f0_default[given]
    }
    return(species_params$f0)
}

#' Get default value for `ks`
#' 
#' Fills in any missing values for `ks` so that the critical feeding level needed
#' to sustain the species is as specified in the `fc` column in the species
#' parameter data frame. If that column is not provided the default critical
#' feeding level \eqn{f_c = 0.2} is used.
#' 
#' @param params A MizerParams object
#' @return A vector with the values of ks for all species
#' @export
#' @concept helper
#' @family functions calculating defaults
get_ks_default <- function(params) {
    assert_that(is(params, "MizerParams"))
    if (!"h" %in% names(params@species_params) ||
        any(is.na(params@species_params$h))) {
        params@species_params$h <- get_h_default(params)
    }
    params <- set_species_param_default(params, "fc", 0.2)
    sp <- params@species_params
    ks_default <- sp$fc * sp$alpha * sp$h * sp$w_mat^(sp$n - sp$p)
    
    message <- ("No ks column so calculating from critical feeding level.")
    params <- set_species_param_default(params, "ks", ks_default, message)
    if (any(is.na(params@species_params$ks) | 
            is.infinite(params@species_params$ks))) {
        stop("Could not calculate default values for the missing species ",
             "parameter ks. Got: ", params@species_params$ks)
    }
    return(params@species_params$ks)
}

#' Validate species parameter data frame
#' 
#' Check validity of species parameters and set defaults for missing but
#' required parameters
#' 
#' @param species_params The user-supplied species parameter data frame
#' @return A valid species parameter data frame
#' 
#' This function throws an error if 
#' * the `species` column does not exist or contains duplicates
#' * the `w_inf` column does not exist or contains NAs or is not numeric
#' 
#' It sets default values if any of the following are missing or NA
#' * `w_mat` is set to `w_inf/4`
#' * `w_min` is set to `0.001`
#' * `alpha` is set to `0.6`
#' * `interaction_resource` is set to `1`
#' 
#' Any `w_mat` that is given that is not smaller than `w_inf` is set to
#' `w_inf / 4`.
#' 
#' Any `w_mat25` that is given that is not smaller than `w_mat` is set to
#' `w_mat * 3^(-0.1)`.
#' 
#' The row names of the returned data frame will be the species names.
#' If `species_params` was provided as a tibble it is converted back to an
#' ordinary data frame.
#' 
#' @concept helper
#' @export
validSpeciesParams <- function(species_params) {
    assert_that(is.data.frame(species_params))
    # Convert a tibble back to an ordinary data frame
    sp <- as.data.frame(species_params,
                        stringsAsFactors = FALSE) # for old versions of R
    
    # check species ----
    if (!("species" %in% colnames(sp))) {
        stop("The species params dataframe needs a column 'species' with the species names")
    }
    species_names <- as.character(sp$species)
    sp$species <- species_names
    row.names(sp) <- species_names
    no_sp <- nrow(sp)
    if (length(unique(species_names)) != no_sp) {
        stop("The species parameter data frame has multiple rows for the same species")
    }
    
    ## For backwards compatibility, allow r_max instead of R_max
    if (!("R_max" %in% names(sp)) &&
        "r_max" %in% names(sp)) {
        names(sp)[names(sp) == "r_max"] <- "R_max"
    }
    
    # check w_inf ----
    if (!("w_inf" %in% colnames(sp))) {
        sp$w_inf <- rep(NA, no_sp)
    }
    missing <- is.na(sp$w_inf)
    if (any(missing)) {
        stop("You need to specify maximum sizes for all species.")
    }
    if (!is.numeric(sp$w_inf)) {
        stop("`w_inf` contains non-numeric values.")
    }
    
    # Defaults ----
    sp <- sp %>% 
        set_species_param_default("w_mat", sp$w_inf / 4) %>% 
        set_species_param_default("w_min", 0.001) %>% 
        set_species_param_default("alpha", 0.6) %>% 
        set_species_param_default("interaction_resource", 1)
    
    # check w_mat ----
    wrong <- sp$w_mat >= sp$w_inf
    if (any(wrong)) {
        message("For the species ", 
                paste(sp$species[wrong], collapse = ", "),
                " the value for `w_mat` is not smaller than that of `w_inf`.",
                " I have corrected that by setting it to about 25% of `w_mat.")
        sp$w_mat[wrong] <- sp$w_inf[wrong] / 4
    }
    
    # check w_mat25 ----
    # For w_mat25 it is o.k. if it is NA, but if given it must be 
    #  smaller than w_mat
    wrong <- !is.na(sp$w_mat25) & sp$w_mat25 >= sp$w_mat
    if (any(wrong)) {
        message("For the species ", 
                paste(sp$species[wrong], collapse = ", "),
                " the value for `w_mat25` is not smaller than that of `w_mat`.",
                " I have corrected that by setting it to about 90% of `w_mat.")
        sp$w_mat25[wrong] <- sp$w_mat[wrong]/(3^(1/10))
    }
    
    sp
}
