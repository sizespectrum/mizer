#' Species parameters
#' 
#' This is the right place to document the use of species parameters in mizer.
#' 
#' @rdname setParams
#' @export
species_params <- function(params) {
    params@species_params
}

#' @rdname setParams
#' @param value A data frame with the species parameters
#' @export
`species_params<-` <- function(params, value) {
    value <- validSpeciesParams(value)
    params@species_params <- value
    setParams(params)
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
#' @keywords internal
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
            message(message)
        }
        species_params <- data.frame(species_params, default,
                                     stringsAsFactors = FALSE)
        colnames(species_params)[[ncol(species_params)]] <- parname
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
                    !is.null(species_params$alpha))
        message("Note: No h provided for some species, so using f0 and k_vb to calculate it.")
        if (!("k_vb" %in% colnames(species_params))) {
            stop("Except I can't because there is no k_vb column in the species data frame")
        }
        if (any(is.na(species_params$k_vb[missing]))) {
            stop("Can not calculate defaults for h because some k_vb values are NA.")
        }
        if (any(species_params$n[missing] != species_params$p[missing])) {
            message("Note: Because you have n != p, the default value is not very good.")
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
#' Fills in any missing values for gamma so that if the prey abundance was
#' described by the power law \eqn{\kappa w^{-\lambda}} then the encounter rate
#' would lead to the feeding level \eqn{f_0}. Only for internal use.
#' 
#' @param params A MizerParams object
#' @return A vector with the values of gamma for all species
#' @export
#' @keywords internal
#' @concept helper
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
        message("Note: Using f0, h, lambda, kappa and the predation kernel to calculate gamma.")
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
        params@species_params$interaction_p <- 1
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

#' Get default value for ks
#' 
#' Fills in any missing values for ks so that the critical feeding level needed
#' to sustain the species is as specified in the `fc` column in the species
#' parameter data frame. If that column is not provided the default critical
#' feeding level \eqn{f_c = 0.2} is used.
#' 
#' @param params A MizerParams object
#' @return A vector with the values of ks for all species
#' @export
#' @keywords internal
#' @concept helper
get_ks_default <- function(params) {
    assert_that(is(params, "MizerParams"))
    if (!"h" %in% names(params@species_params) ||
        any(is.na(params@species_params$h))) {
        params@species_params$h <- get_h_default(params)
    }
    params <- set_species_param_default(params, "fc", 0.2)
    sp <- params@species_params
    ks_default <- sp$fc * sp$alpha * sp$h * sp$w_mat^(sp$n - sp$p)
    
    message <- ("Note: No ks column so calculating from critical feeding level.")
    params <- set_species_param_default(params, "ks", ks_default, message)
    if (any(is.na(params@species_params$ks) | 
            is.infinite(params@species_params$ks))) {
        stop("Could not calculate default values for the missing species ",
             "parameter ks. Got: ", params@species_params$ks)
    }
    return(params@species_params$ks)
}


#' Check validity of species parameters and set defaults for missing but
#' required parameters
#' 
#' @param species_params The user-supplied species parameter data frame
#' @return A valid species parameter data frame
#' @concept("helper")
validSpeciesParams <- function(species_params) {
    assert_that(is.data.frame(species_params))
    
    if (!("species" %in% colnames(species_params))) {
        stop("The species params dataframe needs a column 'species' with the species names")
    }
    species_names <- as.character(species_params$species)
    row.names(species_params) <- species_names
    no_sp <- nrow(species_params)
    if (length(unique(species_names)) != no_sp) {
        stop("The species parameter data frame has multiple rows for the same species")
    }
    
    if (!("w_inf" %in% colnames(species_params))) {
        species_params$w_inf <- rep(NA, no_sp)
    }
    missing <- is.na(species_params$w_inf)
    if (any(missing)) {
        stop("You need to specify maximum sizes for all species.")
    }
    
    species_params <- species_params %>% 
        set_species_param_default("w_mat", species_params$w_inf / 4) %>% 
        set_species_param_default("w_min", 0.001) %>% 
        set_species_param_default("alpha", 0.6) %>% 
        set_species_param_default("interaction_p", 1)
    
    if (any(species_params$w_mat25 >= species_params$w_mat)) {
        idx <- which(species_params$w_mat25 >= species_params$w_mat)
        message("For the species ", species_params$species[idx],
                " the value for `w_mat25` is not smaller than that of `w_mat`.",
                " I have corrected that by setting it to about 90% of `w_mat.")
        species_params$w_mat25[idx] <- species_params$w_mat[idx]/(3^(1/10))
    }
    
    species_params
}
