#' Species parameters
#'
#' These functions allow you to get or set the species-specific parameters
#' stored in a MizerParams object.
#'
#'
#' There are a lot of species parameters and we will list them all below, but
#' most of them have sensible default values. The only required columns are
#' `species` for the species name and `w_inf` for its von Bertalanffy
#' asymptotic size. However if you have information about the values of other
#' parameters then you should provide them.
#'
#' Three species parameters describe maximum sizes and play distinct roles:
#'
#' * `w_inf` is the von Bertalanffy asymptotic size of an average individual.
#'   It is the required maximum-size parameter and is used to set default values
#'   for `w_max`, `w_repro_max` and `w_mat`.
#' * `w_repro_max` is the size at which a typical mature individual invests all
#'   of its available energy into reproduction, see [setReproduction()]. It is
#'   not a hard ceiling on size and defaults to `w_inf`.
#' * `w_max` is purely a computational boundary: it sets the upper end of the
#'   size grid and the range of plots. It defaults to `1.5 * w_inf`. For
#'   backwards compatibility, if `w_inf` is not supplied it is taken from
#'   `w_repro_max` or `w_max` instead.
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
#' * `z0`, `z_ext` and `d` are used to set the external mortality rate, see
#'   [setExtMort()].
#' * `E_ext` and `n` are used to set the external encounter rate, see
#'   [setExtEncounter()].
#' * `D_ext` and `n` are used to set the external diffusion rate, see
#'   [setExtDiffusion()].
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
#' * `R_max` is the parameter in the Beverton-Holt density dependence added to
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
#' parameters like `w_inf`, `w_max`, `w_mat`, `w_mat25`, `w_repro_max` and
#' `w_min` by their corresponding length parameters `l_inf`, `l_max`, `l_mat`,
#' `l_mat25`, `l_repro_max` and `l_min`.
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
#' * If `age_mat` is not supplied, mizer used the von Bertalanffy parameters
#'   `k_vb`, `w_inf` and `t0` as well as the weight-length exponent `b` to
#'   determine it. This is unreliable and is therefore not recommended.
#'
#' Changing these parameters with `species_params<-()` updates the stored
#' species parameter table and triggers a recalculation via [setParams()].
#' However they only affect model behaviour if the corresponding downstream
#' parameters are recalculated rather than kept at explicitly supplied values.
#' In typical workflows these quantities should therefore be changed via
#' `given_species_params<-()`.
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
#' @param object A MizerParams object, a MizerSim object or a data frame
#' @return `species_params()`: Data frame containing all species parameters
#'   currently stored in the model.
#'
#'   `species_params<-()`: Updates the full species parameter table after
#'   validating it with [validSpeciesParams()] and then recalculating the model
#'   parameters with [setParams()].
#'
#'   `given_species_params()`: Data frame containing the species parameter
#'   values that were supplied explicitly by the user.
#'
#'   `given_species_params<-()`: Updates the explicitly supplied species
#'   parameters after validating them with [validGivenSpeciesParams()] and then
#'   recalculating the full species parameter table and dependent model
#'   quantities.
#'
#'   `calculated_species_params()`: Data frame containing only those species
#'   parameter entries that are not explicit user input. Columns that would
#'   consist entirely of `NA` values are dropped.
#' @export
#' @seealso [validSpeciesParams()], [setParams()]
#' @family functions for setting parameters
species_params <- function(object) {
    UseMethod("species_params")
}

#' @rdname species_params
#' @usage NULL
#' @export
species_params.MizerParams <- function(object) {
    object@species_params
}

#' @rdname species_params
#' @usage NULL
#' @export
species_params.MizerSim <- function(object) {
    object@params@species_params
}

#' @rdname species_params
#' @usage NULL
#' @export
species_params.data.frame <- function(object) {
    class(object) <- c("species_params", setdiff(class(object), c("given_species_params", "species_params")))
    check_and_convert_species_params(object)
}

#' @rdname species_params
#' @usage NULL
#' @export
species_params.species_params <- function(object) {
    class(object) <- c("species_params", setdiff(class(object), c("given_species_params", "species_params")))
    check_and_convert_species_params(object)
}

#' @rdname species_params
#' @export
`species_params<-` <- function(object, value) {
    UseMethod("species_params<-")
}

#' @rdname species_params
#' @usage NULL
#' @export
`species_params<-.MizerParams` <- function(object, value) {
    value <- validSpeciesParams(value)
    if (!all(value$species == object@species_params$species)) {
        stop("The species names in the new species parameter data frame do not match the species names in the model.")
    }
    object@species_params <- value
    suppressMessages(setParams(object))
}

#' Test if an object is a species_params object
#'
#' @param x An object to test.
#' @return `TRUE` if `x` is a `species_params` object, `FALSE` otherwise.
#' @export
is.species_params <- function(x) {
    inherits(x, "species_params")
}

check_and_convert_species_params <- function(x) {
    # Check for misspellings
    misspellings <- c("wmin", "wmax", "wmat", "wmat25", "w_mat_25", "Rmax",
                      "Species", "Gamma", "Beta", "Sigma", "Alpha",
                      "W_min", "W_max", "W_mat", "e_repro", "Age_mat",
                      "w_max_mat")
    query <- intersect(misspellings, names(x))
    if (length(query) > 0) {
        warning("Some column names in your species parameter data ",
                "frame are very close to standard parameter names: ",
                paste(query, collapse = ", "),
                ". Did you perhaps mis-spell the names?")
    }

    # Auto convert length to weight if allometric parameters exist
    if (all(c("a", "b") %in% names(x))) {
        mappings <- list(
            c("w_mat", "l_mat"),
            c("w_mat25", "l_mat25"),
            c("w_repro_max", "l_repro_max"),
            c("w_inf", "l_inf"),
            c("w_max", "l_max"),
            c("w_min", "l_min")
        )
        for (m in mappings) {
            pw <- m[1]
            pl <- m[2]
            if (pl %in% names(x)) {
                # Convert the values
                vw <- l2w(x[[pl]], x)
                # If weight is missing or different, update it without
                # triggering recursive validation
                if (!(pw %in% names(x)) || any(is.na(x[[pw]])) || any(abs(x[[pw]] - vw) > 1e-10, na.rm = TRUE)) {
                    saved_class <- class(x)
                    class(x) <- "data.frame"
                    x[[pw]] <- vw
                    class(x) <- saved_class
                }
            }
        }
    }

    # Check w_mat < w_inf consistency
    if (all(c("w_mat", "w_inf") %in% names(x))) {
        wrong <- !is.na(x$w_mat) & !is.na(x$w_inf) & x$w_mat >= x$w_inf
        if (any(wrong)) {
            warning("For the species ",
                    paste(x$species[wrong], collapse = ", "),
                    " the value for `w_mat` is not smaller than that of `w_inf`.")
        }
    }

    x
}

#' @export
`$.species_params` <- function(x, name) {
    out <- NextMethod()
    if (!is.null(out) && !is.data.frame(out) && name != "species") {
        names(out) <- rownames(x)
    }
    out
}


#' @export
`[.species_params` <- function(x, i, j, ..., drop = FALSE) {
    out <- NextMethod("[")
    if (is.data.frame(out)) {
        class(out) <- class(x)
    }
    out
}

#' @export
`[<-.species_params` <- function(x, i, j, ..., value) {
    out <- NextMethod("[<-")
    class(out) <- class(x)
    check_and_convert_species_params(out)
}

#' @export
`[[<-.species_params` <- function(x, i, j, ..., value) {
    out <- NextMethod("[[<-")
    class(out) <- class(x)
    check_and_convert_species_params(out)
}

#' @export
`$<-.species_params` <- function(x, name, value) {
    out <- NextMethod("$<-")
    class(out) <- class(x)
    check_and_convert_species_params(out)
}

#' @export
print.species_params <- function(x, ...) {
    cat("An object of class \"species_params\" containing parameters for", nrow(x), "species:\n")
    core_cols <- c("species", "w_inf", "w_mat", "w_min", "alpha", "erepro")
    cols_to_show <- intersect(core_cols, names(x))
    extra_cols <- setdiff(names(x), core_cols)

    print(as.data.frame(x)[, cols_to_show, drop = FALSE], row.names = FALSE, ...)

    if (length(extra_cols) > 0) {
        cat("With", length(extra_cols), "other parameters:", paste(extra_cols, collapse = ", "), "\n")
    }
    invisible(x)
}

#' @export
print.given_species_params <- function(x, ...) {
    cat("An object of class \"given_species_params\" containing given parameters for", nrow(x), "species:\n")
    core_cols <- c("species", "w_inf", "w_mat", "w_min", "alpha", "erepro")
    cols_to_show <- intersect(core_cols, names(x))
    extra_cols <- setdiff(names(x), core_cols)

    print(as.data.frame(x)[, cols_to_show, drop = FALSE], row.names = FALSE, ...)

    if (length(extra_cols) > 0) {
        cat("With", length(extra_cols), "other parameters:", paste(extra_cols, collapse = ", "), "\n")
    }
    invisible(x)
}

#' @export
summary.species_params <- function(object, ...) {
    cat("Summary of species_params:\n")
    cat("Number of species:", nrow(object), "\n")
    num_cols <- names(object)[vapply(object, is.numeric, logical(1))]
    cat("Parameter ranges:\n")
    for (col in num_cols) {
        vals <- object[[col]]
        vals <- vals[is.finite(vals)]
        if (length(vals) > 0) {
            cat("  ", col, ": min = ", min(vals), ", max = ", max(vals), "\n", sep = "")
        }
    }
    invisible(object)
}


#' @rdname species_params
#' @export
given_species_params <- function(object) {
    UseMethod("given_species_params")
}

#' @rdname species_params
#' @usage NULL
#' @export
given_species_params.MizerParams <- function(object) {
    object@given_species_params
}

#' @rdname species_params
#' @usage NULL
#' @export
given_species_params.MizerSim <- function(object) {
    object@params@given_species_params
}

#' @rdname species_params
#' @usage NULL
#' @export
given_species_params.data.frame <- function(object) {
    class(object) <- c("given_species_params", "species_params", setdiff(class(object), c("given_species_params", "species_params")))
    check_and_convert_species_params(object)
}

#' @rdname species_params
#' @usage NULL
#' @export
given_species_params.given_species_params <- function(object) {
    check_and_convert_species_params(object)
}

#' Test if an object is a given_species_params object
#'
#' @param x An object to test.
#' @return `TRUE` if `x` is a `given_species_params` object, `FALSE` otherwise.
#' @export
is.given_species_params <- function(x) {
    inherits(x, "given_species_params")
}

#' @rdname species_params
#' @export
`given_species_params<-` <- function(object, value) {
    UseMethod("given_species_params<-")
}

#' @rdname species_params
#' @usage NULL
#' @export
`given_species_params<-.MizerParams` <- function(object, value) {
    params <- object
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
#' did not already exist. The signal has class `info_about_default`.
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
get_h_default <- function(object) {
    UseMethod("get_h_default")
}

#' @rdname get_h_default
#' @usage NULL
#' @export
get_h_default.MizerParams <- function(object) {
    get_h_default(object@species_params)
}

#' @rdname get_h_default
#' @usage NULL
#' @export
get_h_default.species_params <- function(object) {
    species_params <- object
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

#' @rdname get_h_default
#' @usage NULL
#' @export
get_h_default.data.frame <- function(object) {
    get_h_default(validSpeciesParams(object))
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
        params@initial_n_pp[] <- resource_power_law(
            params, params@resource_params$kappa,
            params@resource_params$lambda)
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
        params@initial_n_pp[] <- resource_power_law(
            params, params@resource_params$kappa,
            params@resource_params$lambda)
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
