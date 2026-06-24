# Indicator functions for mizer package

# Copyright 2012 Finlay Scott and Julia Blanchard.
# Copyright 2018 Gustav Delius and Richard Southwell.
# Distributed under the GPL 3 or later
# Maintainer: Gustav Delius, University of York, <gustav.delius@york.ac.uk>

#' Description of indicator functions
#'
#' Mizer provides a range of functions to calculate indicators
#' from a MizerSim or MizerParams object.
#'
#' When called with a `MizerSim` object, these functions return a time series
#' of values. When called with a `MizerParams` object, they return a single
#' value calculated from the initial abundances stored in the params object.
#'
#' A list of available indicator functions is given in the table below
#' \tabular{lll}{
#'   Function \tab Returns \tab Description \cr
#'   [getProportionOfLargeFish()] \tab A vector with values at each time step (or a single value for MizerParams). \tab Calculates the proportion of large fish through time. The threshold value can be specified. It is possible to calculation the proportion of large fish based on either length or weight. \cr
#'   [getMeanWeight()] \tab A vector with values at each saved time step (or a single value for MizerParams). \tab The mean weight of the community through time. This is calculated as the total biomass of the community divided by the total abundance. \cr
#'   [getMeanMaxWeight()] \tab Depends on the measure argument. If measure = “both” then you get a matrix with two columns, one with values by numbers, the other with values by biomass at each saved time step (or a named vector for MizerParams). If measure = “numbers” or “biomass” you get a vector of the respective values at each saved time step (or a single value for MizerParams). \tab The mean maximum weight of the community through time. This can be calculated by numbers or by biomass. See the help file for more details. \cr
#'   [getCommunitySlope()] \tab A data.frame with four columns: time step, slope, intercept and the coefficient of determination (or a single-row data.frame for MizerParams). \tab Calculates the slope of the community abundance spectrum through time by performing a linear regression on the logged total numerical abundance and logged body size. \cr
#' }
#'
#' @seealso [summary_functions], [plotting_functions]
#' @name indicator_functions
NULL


#' Calculate the proportion of large fish
#'
#' Calculates the proportion of large fish in a `MizerSim` or `MizerParams`
#' object within user defined size limits. The default option is to use the
#' whole size range. You can specify minimum and maximum size ranges for the
#' species and also the threshold size for large fish. Sizes can be expressed
#' as weight or length. Lengths take precedence over weights (i.e. if both
#' `min_l` and `min_w` are supplied, only `min_l` will be used, and if
#' `threshold_l` is supplied it takes precedence over `threshold_w`). You can
#' also specify the species to be used in the calculation. This function can be
#' used to calculate the Large Fish Index. The proportion is based on either
#' abundance or biomass.
#'
#' @inheritParams getMeanWeight
#' @param threshold_w The weight used as the cutoff between large and small
#'   fish.
#'   Default value is 100.
#' @param threshold_l The length used as the cutoff between large and small
#'   fish. If supplied, this takes precedence over `threshold_w`.
#' @param biomass_proportion A boolean value. If TRUE the proportion calculated
#'   is based on biomass, if FALSE it is based on numbers of individuals.
#'   Default is TRUE.
#' @inheritDotParams get_size_range_array -params
#'
#' @return A vector containing the proportion of large fish through time, or a
#'   single value if called with a `MizerParams` object.
#' @export
#' @family functions for calculating indicators
#' @concept summary_function
#' @examples
#' lfi <- getProportionOfLargeFish(NS_sim, min_w = 10, max_w = 5000,
#'                                 threshold_w = 500)
#' years <- c("1972", "2010")
#' lfi[years]
#' getProportionOfLargeFish(NS_sim)[years]
#' getProportionOfLargeFish(NS_sim, species=c("Herring","Sprat","N.pout"))[years]
#' getProportionOfLargeFish(NS_sim, min_w = 10, max_w = 5000)[years]
#' getProportionOfLargeFish(NS_sim, min_w = 10, max_w = 5000,
#'     threshold_w = 500, biomass_proportion = FALSE)[years]
#' getProportionOfLargeFish(NS_params)
getProportionOfLargeFish <- function(object,
                                     species = NULL,
                                     threshold_w = 100, threshold_l = NULL,
                                     biomass_proportion = TRUE, ...) {
    UseMethod("getProportionOfLargeFish")
}
#' @export
getProportionOfLargeFish.MizerSim <- function(object,
                                              species = NULL,
                                              threshold_w = 100,
                                              threshold_l = NULL,
                                              biomass_proportion = TRUE, ...) {
    sim <- object
    species <- valid_species_arg(sim, species)

    total_size_range <- get_size_range_array(sim@params, ...)
    # This args stuff is pretty ugly - couldn't work out another way of using ...
    args <- list(...)
    args[["params"]] <- sim@params
    args[["max_w"]] <- threshold_w
    args[["max_l"]] <- threshold_l
    large_size_range <- do.call("get_size_range_array", args = args)

    total_size_range <- total_size_range[species, , drop = FALSE]
    large_size_range <- large_size_range[species, , drop = FALSE]

    w <- sim@params@w
    if (!biomass_proportion) { # based on abundance numbers
        w[] <- 1
    }

    n <- sim@n[, species, , drop = FALSE]
    total_measure <-
        apply(sweep(sweep(n, c(2, 3), total_size_range, "*"),
                    3, w * sim@params@dw, "*"),
              1, sum)
    upto_threshold_measure <-
        apply(sweep(sweep(n, c(2, 3), large_size_range, "*"),
                    3, w * sim@params@dw, "*"),
              1, sum)

    1 - (upto_threshold_measure / total_measure)
}
#' @export
getProportionOfLargeFish.MizerParams <- function(object,
                                                 species = NULL,
                                                 threshold_w = 100,
                                                 threshold_l = NULL,
                                                 biomass_proportion = TRUE,
                                                 ...) {
    params <- object
    species <- valid_species_arg(params, species)

    total_size_range <- get_size_range_array(params, ...)
    args <- list(...)
    args[["params"]] <- params
    args[["max_w"]] <- threshold_w
    args[["max_l"]] <- threshold_l
    large_size_range <- do.call("get_size_range_array", args = args)

    total_size_range <- total_size_range[species, , drop = FALSE]
    large_size_range <- large_size_range[species, , drop = FALSE]

    w <- params@w
    if (!biomass_proportion) { # based on abundance numbers
        w[] <- 1
    }

    n <- params@initial_n[species, , drop = FALSE]
    total_measure <- sum(n * total_size_range * w * params@dw)
    upto_threshold_measure <- sum(n * large_size_range * w * params@dw)

    1 - (upto_threshold_measure / total_measure)
}


#' Calculate the mean weight of the community
#'
#' Calculates the mean weight of the community. This is simply the total
#' biomass of the community divided by the abundance in numbers. You can
#' specify minimum and maximum weight or length range for the species. Lengths
#' take precedence over weights (i.e. if both min_l and min_w are supplied, only
#' min_l will be used). You can also specify the species to be used in the
#' calculation.
#'
#' @param object A \linkS4class{MizerSim} or \linkS4class{MizerParams} object
#' @inheritParams valid_species_arg
#' @inheritDotParams get_size_range_array -params
#'
#' @return A vector containing the mean weight of the community through time,
#'   or a single value if called with a `MizerParams` object.
#' @export
#' @family functions for calculating indicators
#' @concept summary_function
#' @examples
#' mean_weight <- getMeanWeight(NS_sim)
#' years <- c("1967", "2010")
#' mean_weight[years]
#' getMeanWeight(NS_sim, species = c("Herring", "Sprat", "N.pout"))[years]
#' getMeanWeight(NS_sim, min_w = 10, max_w = 5000)[years]
#' getMeanWeight(NS_params)
getMeanWeight <- function(object, species = NULL, ...) {
    UseMethod("getMeanWeight")
}
#' @export
getMeanWeight.MizerSim <- function(object, species = NULL, ...) {
    sim <- object
    assert_that(is(sim, "MizerSim"))
    species <- valid_species_arg(sim, species)
    n_species <- getN(sim, ...)
    biomass_species <- getBiomass(sim, ...)
    n_total <- apply(n_species[, species, drop = FALSE], 1, sum)
    biomass_total <- apply(biomass_species[, species, drop = FALSE], 1, sum)
    biomass_total / n_total
}
#' @export
getMeanWeight.MizerParams <- function(object, species = NULL, ...) {
    params <- object
    species <- valid_species_arg(params, species)
    n_total <- sum(getN(params, ...)[species])
    biomass_total <- sum(getBiomass(params, ...)[species])
    biomass_total / n_total
}


#' Calculate the mean maximum weight of the community
#'
#' Calculates the mean maximum weight of the community. This can be calculated
#' by numbers or biomass. The calculation is the sum of the `w_inf` * abundance
#' of each species, divided by the total abundance community, where abundance is
#' either in biomass or numbers. You can specify minimum and maximum weight or
#' length range for the species. Lengths take precedence over weights (i.e. if
#' both min_l and min_w are supplied, only min_l will be used). You can also
#' specify the species to be used in the calculation.
#'
#' @inheritParams getMeanWeight
#' @param measure The measure to return. Can be 'numbers', 'biomass' or 'both'
#' @inheritDotParams get_size_range_array -params
#'
#' @return Depends on the `measure` argument. If \code{measure = “both”}
#'   then you get a matrix with two columns, one with values by numbers,
#'   the other with values by biomass at each saved time step (or a named
#'   vector with two entries for `MizerParams`). If \code{measure =
#'   “numbers”} or \code{“biomass”} you get a vector of the respective values
#'   at each saved time step (or a single value for `MizerParams`).
#' @export
#' @family functions for calculating indicators
#' @concept summary_function
#' @examples
#' mmw <- getMeanMaxWeight(NS_sim)
#' years <- c("1967", "2010")
#' mmw[years, ]
#' getMeanMaxWeight(NS_sim, species=c("Herring","Sprat","N.pout"))[years, ]
#' getMeanMaxWeight(NS_sim, min_w = 10, max_w = 5000)[years, ]
#' getMeanMaxWeight(NS_params)
getMeanMaxWeight <- function(object, species = NULL,
                             measure = "both", ...) {
    UseMethod("getMeanMaxWeight")
}
#' @export
getMeanMaxWeight.MizerSim <- function(object, species = NULL,
                                      measure = "both", ...) {
    sim <- object
    assert_that(is(sim, "MizerSim"))
    if (!(measure %in% c("both", "numbers", "biomass"))) {
        stop("measure must be one of 'both', 'numbers' or 'biomass'")
    }
    species <- valid_species_arg(sim, species)
    n_species <- getN(sim, ...)
    biomass_species <- getBiomass(sim, ...)
    n_winf <- apply(
        sweep(n_species, 2, sim@params@species_params$w_inf, "*")[,
            species, drop = FALSE
        ],
        1, sum
    )
    biomass_winf <- apply(
        sweep(biomass_species, 2, sim@params@species_params$w_inf, "*")[,
            species, drop = FALSE
        ],
        1, sum
    )
    mmw_numbers <- n_winf / apply(n_species[, species, drop = FALSE], 1, sum)
    mmw_biomass <- biomass_winf /
        apply(biomass_species[, species, drop = FALSE], 1, sum)
    if (measure == "numbers") {
        return(mmw_numbers)
    }
    if (measure == "biomass") {
        return(mmw_biomass)
    }
    cbind(mmw_numbers, mmw_biomass)
}
#' @export
getMeanMaxWeight.MizerParams <- function(object, species = NULL,
                                         measure = "both", ...) {
    params <- object
    if (!(measure %in% c("both", "numbers", "biomass"))) {
        stop("measure must be one of 'both', 'numbers' or 'biomass'")
    }
    species <- valid_species_arg(params, species)
    n_species <- getN(params, ...)[species]
    biomass_species <- getBiomass(params, ...)[species]
    w_inf <- params@species_params$w_inf[
        match(species, params@species_params$species)
    ]
    n_winf <- sum(n_species * w_inf)
    biomass_winf <- sum(biomass_species * w_inf)
    mmw_numbers <- n_winf / sum(n_species)
    mmw_biomass <- biomass_winf / sum(biomass_species)
    if (measure == "numbers") {
        return(mmw_numbers)
    }
    if (measure == "biomass") {
        return(mmw_biomass)
    }
    c(mmw_numbers = mmw_numbers, mmw_biomass = mmw_biomass)
}


#' Calculate the slope of the community abundance
#'
#' Calculates the slope of the community abundance by performing a linear
#' regression on the logged total numerical abundance at weight and logged
#' weights (natural logs, not log to base 10, are used). You can specify
#' minimum and maximum weight or length range for the species. Lengths take
#' precedence over weights (i.e. if both min_l and min_w are supplied, only
#' min_l will be used). You can also specify the species to be used in the
#' calculation.
#'
#' @inheritParams getMeanWeight
#' @param biomass Boolean. If TRUE (default), the abundance is based on biomass,
#'   if FALSE the abundance is based on numbers.
#' @inheritDotParams get_size_range_array -params
#'
#' @return A data.frame with columns slope, intercept and the coefficient of
#'   determination R^2 (and a time step column when called with a `MizerSim`
#'   object).
#' @export
#' @family functions for calculating indicators
#' @concept summary_function
#' @examples
#' # Slope based on biomass, using all species and sizes
#' slope_biomass <- getCommunitySlope(NS_sim)
#' slope_biomass[1, ] # in 1976
#' slope_biomass[idxFinalT(NS_sim), ] # in 2010
#'
#' # Slope based on numbers, using all species and sizes
#' slope_numbers <- getCommunitySlope(NS_sim, biomass = FALSE)
#' slope_numbers[1, ] # in 1976
#'
#' # Slope based on biomass, using all species and sizes between 10g and 1000g
#' slope_biomass <- getCommunitySlope(NS_sim, min_w = 10, max_w = 1000)
#' slope_biomass[1, ] # in 1976
#'
#' # Slope based on biomass, using only demersal species and
#' # sizes between 10g and 1000g
#' dem_species <- c("Dab","Whiting", "Sole", "Gurnard", "Plaice",
#'                  "Haddock", "Cod", "Saithe")
#' slope_biomass <- getCommunitySlope(NS_sim, species = dem_species,
#'                                    min_w = 10, max_w = 1000)
#' slope_biomass[1, ] # in 1976
#'
#' getCommunitySlope(NS_params)
getCommunitySlope <- function(object, species = NULL,
                              biomass = TRUE, ...) {
    UseMethod("getCommunitySlope")
}
#' @export
getCommunitySlope.MizerSim <- function(object, species = NULL,
                                       biomass = TRUE, ...) {
    sim <- object
    assert_that(is(sim, "MizerSim"))
    species <- valid_species_arg(sim, species)
    size_range <- get_size_range_array(sim@params, ...)
    # set entries for unwanted sizes to zero and sum over wanted species, giving
    # array (time x size)
    total_n <-
        apply(sweep(sim@n, c(2, 3), size_range, "*")[, species, , drop = FALSE],
              c(1, 3), sum)
    # numbers or biomass?
    if (biomass) {
        total_n <- sweep(total_n, 2, sim@params@w, "*")
    }
    # previously unwanted entries were set to zero, now set them to NA
    # so that they will be ignored when fitting the linear model
    total_n[total_n <= 0] <- NA
    # fit linear model at every time and put result in data frame
    slope <- plyr::adply(total_n, 1, function(x, w) {
        summary_fit <- summary(lm(log(x) ~ log(w)))
        data.frame(
            slope = summary_fit$coefficients[2, 1],
            intercept = summary_fit$coefficients[1, 1],
            r2 = summary_fit$r.squared
        )
    }, w = sim@params@w)
    dimnames(slope)[[1]] <- slope[, 1]
    slope[, -1]
}
#' @export
getCommunitySlope.MizerParams <- function(object, species = NULL,
                                          biomass = TRUE, ...) {
    params <- object
    species <- valid_species_arg(params, species)
    size_range <- get_size_range_array(params, ...)
    # set entries for unwanted sizes to zero and sum over wanted species,
    # giving a vector (size)
    total_n <- colSums(params@initial_n[species, , drop = FALSE] *
                           size_range[species, , drop = FALSE])
    # numbers or biomass?
    if (biomass) {
        total_n <- total_n * params@w
    }
    # set entries that are zero to NA so they are ignored in the linear model
    total_n[total_n <= 0] <- NA
    summary_fit <- summary(lm(log(total_n) ~ log(params@w)))
    data.frame(
        slope = summary_fit$coefficients[2, 1],
        intercept = summary_fit$coefficients[1, 1],
        r2 = summary_fit$r.squared
    )
}
