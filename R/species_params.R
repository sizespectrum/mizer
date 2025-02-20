#' Species parameters
#' 
#' These functions allow you to get or set the species-specific parameters
#' stored in a MizerParams object.
#'
#' 
#' There are a lot of species parameters and we will list them all below, but
#' most of them have sensible default values. The only required columns are
#' `species` for the species name and `w_max` for its maximum size. However
#' if you have information about the values of other parameters then you should
#' provide them.
#' 
#' Mizer distinguishes between the species parameters that you have given
#' explicitly and the species parameters that have been calculated by mizer or
#' set to default values. You can retrieve the given species parameters with
#' `given_species_params()` and the calculated ones with
#' `calculated_species_params()`. You get all species_params with
#' `species_params()`.
#' 
#' If you change given species parameters with `given_species_params<-()` this
#' will trigger a re-calculation of the calculated species parameters, where
#' necessary. However if you change species parameters with `species_params<-()`
#' no recalculation will take place and furthermore your values could be
#' overwritten by a future recalculation triggered by a call to
#' `given_species_params<-()` . So in most use cases you will only want to use
#' `given_species_params<-()`.
#' 
#' There are some species parameters that are used to set up the
#' size-dependent parameters that are used in the mizer model:
#' 
#' * `gamma` and `q` are used to set the search volume, see [setSearchVolume()].
#' * `h` and `n` are used to set the maximum intake rate, see [setMaxIntakeRate()].
#' * `k`, `ks` and `p` are used to set activity and basic metabolic rate, 
#'   see [setMetabolicRate()].
#' * `z0` is used to set the external mortality rate, see [setExtMort()].
#' * `w_mat`, `w_mat25`, `w_repro_max` and `m` are used to set the allocation to
#'   reproduction, see [setReproduction()].
#' * `pred_kernel_type` specifies the shape of the predation kernel. The default
#'   is a "lognormal", for other options see the "Setting predation kernel"
#'   section in the help for [setPredKernel()].
#' * `beta` and `sigma` are parameters of the lognormal predation kernel, see
#'   [lognormal_pred_kernel()]. There will be other parameters if you are 
#'   using other predation kernel functions.
#'   
#' When you change one of the above species parameters using 
#' `given_species_params<-()` or `species_params<-()`, the new value will be
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
#' If you have supplied the `a` and `b` parameters, then you can replace weight
#' parameters like `w_max`, `w_mat`, `w_mat25`, `w_repro_max` and `w_min` by
#' their corresponding length parameters `l_max`, `l_mat`, `l_mat25`,
#' `l_repro_max` and `l_min`.
#'   
#' The parameters that are only used to calculate default values for other
#' parameters are:
#' 
#' * `f0` is the feeding level and is used to get a default value for the
#'   coefficient of the search volume `gamma`, see [get_gamma_default()].
#' * `fc` is the critical feeding level below which the species can not 
#'   maintain itself. This is used to get a default value for the coefficient
#'   `ks` of the metabolic rate, see [get_ks_default()].
#' * `age_mat` is the age at maturity and is used to get a default value for
#'   the coefficient `h` of the maximum intake rate, see [get_h_default()].
#'   
#' Note that setting these parameters with `species_params<-()` will have no
#' effect. You need to set them with `given_species_params<-()` in order to
#' trigger a re-calculation of the other species parameters.
#' 
#' In the past, mizer also used the von Bertalanffy parameters `k_vb`, `w_inf`
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
#' this information can also optionally be provided as species parameters and
#' [newMultispeciesParams()] will transfer them to the `gear_params` data frame.
#' However changing these parameters later in the species parameter data frames
#' will have no effect.
#' 
#' You are allowed to include additional columns in the species parameter
#' data frames. They will simply be ignored by mizer but will be stored in the
#' MizerParams object, in case your own code makes use of them.
#' 
#' @param params A MizerParams object
#' @return Data frame of species parameters
#' @export
#' @seealso [validSpeciesParams()], [setParams()]
#' @family functions for setting parameters
species_params <- function(params) {
    assert_that(is(params, "MizerParams"))
    params@species_params
}

#' @rdname species_params
#' @param value A data frame with the species parameters
#' @export
`species_params<-` <- function(params, value) {
    assert_that(is(params, "MizerParams"))
    value <- validSpeciesParams(value)
    if (!all(value$species == params@species_params$species)) {
        stop("The species names in the new species parameter data frame do not match the species names in the model.")
    }
    params@species_params <- value
    suppressMessages(setParams(params))
}


#' @rdname species_params
#' @export
given_species_params <- function(params) {
    assert_that(is(params, "MizerParams"))
    params@given_species_params
}

#' @rdname species_params
#' @export
`given_species_params<-` <- function(params, value) {
    assert_that(is(params, "MizerParams"))
    value <- validGivenSpeciesParams(value)
    if (!all(value$species == params@species_params$species)) {
        stop("The species names in the new species parameter data frame do not match the species names in the model.")
    }
    old_value <- params@given_species_params
    
    # Create data frame which contains only the values that have changed 
    common_columns <- intersect(names(value), names(params@given_species_params))
    new_columns <- setdiff(names(value), names(params@given_species_params))
    changes <- value[common_columns]
    changes[changes == params@given_species_params[common_columns]] <- NA
    # Remove columns that only contain NAs
    changes <- changes %>% select(where(~ !all(is.na(.))))
    # Add new columns
    changes <- cbind(changes, value[new_columns])
    
    # Give warnings when values are changed that will have no impact
    if ("gamma" %in% names(params@given_species_params) &
        "f0" %in% names(changes) &
        any(!is.na(params@given_species_params$gamma[!is.na(changes$f0)]))) {
        warning("You have specified some values for `f0` that are going to be ignored because values for `gamma` have already been given.")
    }
    if ("ks" %in% names(params@given_species_params) &
        "fc" %in% names(changes) &
        any(!is.na(params@given_species_params$ks[!is.na(changes$fc)]))) {
        warning("You have specified some values for `fc` that are going to be ignored because values for `ks` have already been given.")
    }
    if ("h" %in% names(params@given_species_params) &
        "age_mat" %in% names(changes) &
        any(!is.na(params@given_species_params$h[!is.na(changes$age_mat)]))) {
        warning("You have specified some values for `age_mat` that are going to be ignored because values for `h` have already been given.")
    }
    
    # Warn when user tries to change gear parameters
    if (any(c("catchability", "selectivity", "l50", "l25", "sel_func") %in% 
            names(changes))) {
        warning("To make changes to gears you should use `gear_params()<-`, not `species_params()`.")
    }
    if ("yield_observed" %in% names(changes)) {
        warning("To change the observed yield you should use `gear_params()<-`, not `species_params()`.")
    }
    
    params@given_species_params <- value
    params@species_params <- validSpeciesParams(value)
    suppressMessages(setParams(params))
}

#' @rdname species_params
#' @export
calculated_species_params <- function(params) {
    assert_that(is(params, "MizerParams"))
    # Identifying common columns
    common_cols <- intersect(names(params@species_params), 
                             names(params@given_species_params))
    # Copy df1 to new_df
    calculated <- params@species_params
    # remove the entries that are also in given_species_params
    for (col in common_cols) {
        calculated[[col]] <- replace_with_na(calculated[[col]], 
                                             params@given_species_params[[col]])
    }
    # Removing columns that only contain NAs
    calculated <- calculated %>%
        select(where(~ !all(is.na(.))))
    
    return(calculated)
}

# Function to replace overlapping entries with NA
replace_with_na <- function(x, y) {
    ifelse(is.na(y), x, NA)
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
#' 
#' If no growth information is given at all for a species, the default is set
#' to `h = 30`.
#' 
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
    assert_that("n" %in% names(species_params))
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
        
        species_params[missing, "h"] <- h[missing]
        
        # If no acceptable default could be calculated, set h=30
        missing <- is.na(species_params[["h"]]) | species_params[["h"]] <= 0
        if (any(missing)) {
            signal("For species where no growth information is available the parameter h has been set to h = 30.",
                   class = "info_about_default", var = "h", level = 3)
            species_params[missing, "h"] <- 30
        }
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
    assert_that(is(params, "MizerParams"),
                "n" %in% names(params@species_params),
                "p" %in% names(params@species_params))
    if (!"h" %in% names(params@species_params) ||
        any(is.na(params@species_params[["h"]]))) {
        params@species_params[["h"]] <- get_h_default(params)
    }
    params <- set_species_param_default(params, "fc", 0.2)
    sp <- params@species_params
    ks_default <- sp$fc * sp$alpha * sp[["h"]] * sp$w_mat^(sp[["n"]] - sp[["p"]])
    
    message <- ("No ks column so calculating from critical feeding level.")
    sp <- set_species_param_default(sp, "ks", ks_default, message)
    if (any(is.na(sp$ks) |  is.infinite(sp$ks))) {
        stop("Could not calculate default values for the missing species ",
             "parameter ks. Got: ", sp$ks)
    }
    return(sp$ks)
}
