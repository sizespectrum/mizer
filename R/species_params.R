#' Species parameters
#' 
#' These functions allow you to get or set the species parameters stored in
#' a MizerParams object.
#' 
#' The `species_params` data frame holds species-specific parameters. The data
#' frame has one row for each species and one column for each species parameter.
#' There are a lot of species parameters and we will list them all below, but
#' most of them have sensible default values. The only required columns are
#' `species` for the species name and `w_max` for its maximum size. However
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
#' * `w_mat`, `w_mat25`, `w_max` and `m` are used to set the allocation to
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
#' If you have supplied these, then you can replace weight parameters
#' like `w_max`, `w_mat`, `w_mat25` and `w_min` by their corresponding length
#' parameters `l_max`, `l_mat`, `l_mat25` and `l_min`.
#'   
#' Not all species parameters have to be specified by the user. If they are
#' missing, [newMultispeciesParams()] will give them default values, sometimes
#' by using other species parameters. The parameters that are only used to
#' calculate default values for other parameters are:
#' 
#' * `f0` is the feeding level and is used to get a default value for the
#'   coefficient of the search volume `gamma`, see [get_gamma_default()].
#' * `fc` is the critical feeding level below which the species can not 
#'   maintain itself. This is used to get a default value for the coefficient
#'   `ks` of the metabolic rate, see [get_ks_default()].
#' * `age_mat` is the age at maturity and is used to get a default value for
#'   the coefficient `h` of the maximum intake rate, see [get_h_default()].
#'   
#' Note that these parameters are ignored if the parameters for which they allow
#' defaults to be calculated have instead been set explicitly. Also, these
#' parameters will only be used when setting up a new model with 
#' [newMultispeciesParams()]. Changing them later will have no effect 
#' because the default for the other parameters will not be recalculated.
#' 
#' In the past mizer also used the von Bertalanffy parameters `k_vb`, `w_inf`
#' and `t0` to determine a default for `h`. This is unreliable and is therefore
#' now deprecated.
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
#' Sets `h` so that the species reaches maturity size `w_mat` at the maturity
#' age `age_mat` if it feeds at feeding level `f0`.
#'
#' If `age_mat` is missing in the species parameter data frame, then it is
#' calculated from the von Bertalanffy growth curve parameters `k_vb` and
#' (optionally `t0`) taken from the species parameter data frame. This is not
#' reliable and a warning is issued.
#' @param params A MizerParams object or a species parameter data frame
#' @return A vector with the values of h for all species
#' @export
#' @keywords internal
#' @concept helper
#' @family functions calculating defaults
get_h_default <- function(params) {
    if (is(params, "MizerParams")) {
        species_params <- params@species_params
    } else {
        species_params <- validSpeciesParams(params)
    }
    species_params <- set_species_param_default(species_params, "f0", 0.6)
    if (!("h" %in% colnames(species_params))) {
        species_params[["h"]] <- rep(NA, nrow(species_params))
    }
    missing <- is.na(species_params[["h"]])
    if (any(missing)) {
        # The following should be assured by `validSpeciesParams()`
        assert_that(is.numeric(species_params$f0),
                    noNA(species_params$alpha),
                    "alpha" %in% names(species_params))
        signal("No h provided for some species, so using age at maturity to calculate it.",
                      class = "info_about_default", var = "h", level = 3)
        if (!isTRUE(all.equal(species_params$n[missing], species_params$p[missing],
                              check.attributes = FALSE))) {
            signal("Because you have n != p, the default value for `h` is not very good.",
                   class = "info_about_default", var = "h", level = 1)
        }
        species_params <- species_params %>% 
            set_species_param_default("fc", 0.2) %>% 
            set_species_param_default("n", 3 / 4) %>% 
            set_species_param_default(
                "age_mat", age_mat_vB(species_params),
                strwrap("Because the age at maturity is not known, I need to 
                        fall back to using von Bertalanffy parameters, where 
                        available, and this is not reliable.")
            )
        w_mat <- species_params$w_mat
        w_min <- species_params$w_min
        age_mat <- species_params$age_mat
        n <- species_params[["n"]]
        h <- (w_mat^(1 - n) - w_min^(1 - n)) / age_mat / (1 - n) / 
            species_params$alpha / (species_params$f0 - species_params$fc)
        
        if (any(is.na(h[missing])) || any(h[missing] <= 0)) {
            stop("Could not calculate h.")
        }
        species_params[missing, "h"] <- h[missing]
    }
    return(species_params[["h"]])
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
                class = "info_about_default", var = "gamma", level = 3)
        if (!"h" %in% names(params@species_params) || 
            any(is.na(species_params[["h"]]))) {
            species_params[["h"]] <- get_h_default(params)
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
            (2 + params@species_params[["q"]] - params@resource_params$lambda)
        # Now set gamma so that this available energy leads to f0
        gamma_default <- (species_params[["h"]] / avail_energy) * 
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
            any(is.na(species_params[["h"]]))) {
            species_params[["h"]] <- get_h_default(params)
        }
        # Calculate available energy by setting a power-law prey spectrum
        params@initial_n[] <- 0
        params@species_params$interaction_resource <- 1
        params@initial_n_pp[] <- params@resource_params$kappa * 
            params@w_full^(-params@resource_params$lambda)
        avail_energy <- getEncounter(params)[, length(params@w)] /
            params@w[length(params@w)] ^ 
            (2 + params@species_params[["q"]] - params@resource_params$lambda)
        # Now set f0 so that this available energy leads to f0
        f0_default <- 1 / (species_params[["h"]] / avail_energy + 1)
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
        any(is.na(params@species_params[["h"]]))) {
        params@species_params[["h"]] <- get_h_default(params)
    }
    params <- set_species_param_default(params, "fc", 0.2)
    sp <- params@species_params
    ks_default <- sp$fc * sp$alpha * sp[["h"]] * sp$w_mat^(sp[["n"]] - sp[["p"]])
    
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
#' required parameters.
#' 
#' @param species_params The user-supplied species parameter data frame
#' @return A valid species parameter data frame
#' 
#' This function throws an error if 
#' * the `species` column does not exist or contains duplicates
#' * the maximum size is not specified for all species
#' 
#' If a weight-based parameter is missing but the corresponding length-based
#' parameter is given, as well as the `a` and `b` parameters for length-weight
#' conversion, then the weight-based parameters are added. If both length and
#' weight are given, then weight is used and a warning is issued if the two are
#' inconsistent.
#' 
#' If a `w_inf` column is given but no `w_max` then the value from `w_inf` is
#' used. This is for backwards compatibility.
#' 
#' The function sets default values if any of the following are missing or NA
#' * `w_mat` is set to `w_max/4`
#' * `w_min` is set to `0.001`
#' * `alpha` is set to `0.6`
#' * `interaction_resource` is set to `1`
#' 
#' Some inconsistencies in the size parameters are resolved as follows:
#' * Any `w_mat` that is not smaller than `w_max` is set to `w_max / 4`.
#' * Any `w_mat25` that is not smaller than `w_mat` is set to NA.
#' * Any `w_min` that is not smaller than `w_mat` is set to `0.001` or 
#'   `w_mat /10`, whichever is smaller.
#' 
#' The row names of the returned data frame will be the species names.
#' If `species_params` was provided as a tibble it is converted back to an
#' ordinary data frame.
#' 
#' The function tests for some typical misspellings of parameter names, like
#' wrong capitalisation or missing underscores and issues a warning if it 
#' detects such a name.
#' 
#' Note that the species parameters returned by this function are not guaranteed
#' to produce a viable model. More checks of the parameters are performed by the
#' individual rate-setting functions (see [setParams()] for the list of these
#' functions).
#' @seealso species_params
#' @concept helper
#' @export
validSpeciesParams <- function(species_params) {
    assert_that(is.data.frame(species_params))
    # Convert a tibble back to an ordinary data frame
    sp <- as.data.frame(species_params,
                        stringsAsFactors = FALSE) # for old versions of R
    
    # Check for misspellings ----
    misspellings <- c("wmin", "wmax", "wmat", "wmat25", "w_mat_25", "Rmax",
                      "Species", "Gamma", "Beta", "Sigma", "Alpha",
                      "W_min", "W_max", "W_mat", "e_repro", "Age_mat")
    query <- intersect(misspellings, names(sp))
    if (length(query) > 0) {
        warning("Some column names in your species parameter data ",
                "frame are very close to standard parameter names: ",
                paste(query, collapse = ", "),
                ". Did you perhaps mis-spell the names?")
    }
    
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
    
    # Convert lengths to weights ----
    if (all(c("a", "b") %in% names(sp))) {
        sp <- sp %>%
            set_species_param_from_length("w_mat", "l_mat") %>%
            set_species_param_from_length("w_mat25", "l_mat25") %>%
            set_species_param_from_length("w_max", "l_max") %>%
            set_species_param_from_length("w_min", "l_min")
    }
    
    # check w_max ----
    if (!("w_max" %in% names(sp))) {
        # If old name `w_inf` is used, then copy over to `w_max`
        if ("w_inf" %in% names(sp)) {
            sp$w_max <- sp$w_inf
        } else {
            sp$w_max <- rep(NA, no_sp)
        }
    }
    missing <- is.na(sp$w_max)
    if (any(missing)) {
        stop("You need to specify maximum sizes for all species.")
    }
    if (!is.numeric(sp$w_max)) {
        stop("`w_max` contains non-numeric values.")
    }
    
    # Defaults ----
    sp <- sp %>% 
        set_species_param_default("w_mat", sp$w_max / 4) %>% 
        set_species_param_default("w_min", 0.001) %>% 
        set_species_param_default("alpha", 0.6) %>% 
        set_species_param_default("interaction_resource", 1)
    
    # check w_mat ----
    wrong <- sp$w_mat >= sp$w_max
    if (any(wrong)) {
        message("For the species ", 
                paste(sp$species[wrong], collapse = ", "),
                " the value for `w_mat` is not smaller than that of `w_max`.",
                " I have corrected that by setting it to about 25% of `w_mat.")
        sp$w_mat[wrong] <- sp$w_max[wrong] / 4
    }
    
    # check w_mat25 ----
    # For w_mat25 it is o.k. if it is NA, but if given it must be 
    #  smaller than w_mat
    wrong <- !is.na(sp$w_mat25) & sp$w_mat25 >= sp$w_mat
    if (any(wrong)) {
        message("For the species ", 
                paste(sp$species[wrong], collapse = ", "),
                " the value for `w_mat25` is not smaller than that of `w_mat`.",
                " I have corrected that by setting it to NA.")
        sp$w_mat25[wrong] <- NA
    }
    
    # check w_min ----
    wrong <- sp$w_min >= sp$w_mat
    if (any(wrong)) {
        sp$w_min[wrong] <- pmin(0.001, sp$w_mat[wrong] / 10)
        message("For the species ", 
                paste(sp$species[wrong], collapse = ", "),
                " the value for `w_min` is not smaller than that of `w_mat`.",
                " I have reduced the values.")
    }
    
    sp
}

# Set weight-based parameter from length-based parameter
set_species_param_from_length <- function(sp, pw, pl) {
    if (pl %in% names(sp)) {
        vw <- l2w(sp[[pl]], sp)
        sp <- set_species_param_default(sp, pw, vw)
        incons <- !is.na(sp[[pw]]) & !is.na(sp[[pl]]) &
            sp[[pw]] != vw
        if (any(incons)) {
            warning("For the following species your value for ", pw,
                    " and your value for ", pl, " are not consistent: ",
                    paste(sp$species[incons], collapse = ", "))
        }
    }
    sp
}
