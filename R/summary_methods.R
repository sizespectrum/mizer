# Summary and indicator functions for mizer package

# Copyright 2012 Finlay Scott and Julia Blanchard.
# Copyright 2018 Gustav Delius and Richard Southwell.
# Development has received funding from the European Commission's Horizon 2020 
# Research and Innovation Programme under Grant Agreement No. 634495 
# for the project MINOUW (http://minouw-project.eu/).
# Distributed under the GPL 3 or later 
# Maintainer: Gustav Delius, University of York, <gustav.delius@york.ac.uk>

# Soundtrack: The Definitive Lead Belly

# Summary functions ####

#' Description of summary functions
#' 
#' Mizer provides a range of functions to summarise the results of a simulation.
#'
#' A list of available summary functions is given in the table below.
#' \tabular{lll}{
#'   Function \tab Returns \tab Description \cr
#'   [getDiet()] \tab Three dimensional array (predator x size x prey) \tab Diet of predator at size, resolved by prey species \cr
#'   [getSSB()] \tab Two dimensional array (time x species) \tab Total Spawning Stock Biomass (SSB) of each species through time where SSB is calculated as the sum of weight of all mature individuals. \cr
#'   [getBiomass()] \tab Two dimensional array (time x species) \tab Total biomass of each species through time. \cr
#'   [getN()] \tab Two dimensional array (time x species) \tab Total abundance of each species through time. \cr
#'   [getFeedingLevel()] \tab Three dimensional array (time x species x size) \tab Feeding level of each species by size through time. \cr
#'   \code{\link{getM2}} \tab Three dimensional array (time x species x size) \tab The predation mortality imposed on each species by size through time. \cr
#'   [getFMort()] \tab Three dimensional array (time x species x size) \tab Total fishing mortality on each species by size through time. \cr
#'   [getFMortGear()] \tab Four dimensional array (time x gear x species x size) \tab Fishing mortality on each species by each gear at size through time. \cr
#'   [getYieldGear()] \tab Three dimensional array (time x gear x species) \tab Total yield by gear and species through time. \cr
#'   [getYield()] \tab Two dimensional array (time x species) \tab Total yield of each species across all gears through time. \cr
#' }
#'
#' @seealso [indicator_functions], [plotting_functions]
#' @name summary_functions
NULL

#' Get diet of predator at size, resolved by prey species
#'
#' Calculates the rate at which a predator of a particular species and size
#' consumes biomass of each prey species, resource, and other components of the
#' ecosystem. Returns either the rates in grams/year or the proportion of the
#' total consumption rate.
#'
#' The rates \eqn{D_{ij}(w)} at which a predator of species \eqn{i}
#' and size \eqn{w} consumes biomass from prey species \eqn{j} are
#' calculated from the predation kernel \eqn{\phi_i(w, w_p)},
#' the search volume \eqn{\gamma_i(w)}, the feeding level \eqn{f_i(w)}, the
#' species interaction matrix \eqn{\theta_{ij}} and the prey abundance density
#' \eqn{N_j(w_p)}:
#' \deqn{
#' D_{ij}(w, w_p) = (1-f_i(w)) \gamma_i(w) \theta_{ij}
#' \int N_j(w_p) \phi_i(w, w_p) w_p dw_p.
#' }
#' The prey index \eqn{j} runs over all species and the resource. 
#' 
#' Extra columns are added for the external encounter rate and for any extra
#' ecosystem components in your model for which you have defined an encounter
#' rate function. These encounter rates are multiplied by \eqn{1-f_i(w)} to give
#' the rate of consumption of biomass from these extra components.
#' 
#' This function performs the same integration as [getEncounter()] but does not
#' aggregate over prey species, and multiplies by \eqn{1-f_i(w)} to get the
#' consumed biomass rather than the available biomass. Outside the range of
#' sizes for a predator species the returned rate is zero.
#'
#' @inheritParams getEncounter
#' @param proportion If TRUE (default) the function returns the diet as a
#'   proportion of the total consumption rate. If FALSE it returns the 
#'   consumption rate in grams per year.
#' 
#' @return An array (predator species  x predator size x 
#'   (prey species + resource + other components). Dimnames are "prey", "w",
#'   and "predator".
#'   
#' @export
#' @family summary functions
#' @concept summary_function
#' @seealso [plotDiet()]
#' @examples
#' diet <- getDiet(NS_params)
#' str(diet)
getDiet <- function(params,
                    n = initialN(params), 
                    n_pp = initialNResource(params),
                    n_other = initialNOther(params),
                    proportion = TRUE) {
    # The code is based on that for getEncounter()
    params <- validParams(params)
    species <- params@species_params$species
    no_sp <- length(species)
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)
    no_other <- length(params@other_encounter)
    other_names <- names(params@other_encounter)
    assert_that(identical(dim(n), c(no_sp, no_w)),
                length(n_pp) == no_w_full)
    diet <- array(0, dim = c(no_sp, no_w, no_sp + 2 + no_other),
                  dimnames = list("predator" = species,
                                  "w" = dimnames(params@initial_n)$w,
                                  "prey" = c(species, 
                                             "Resource", "External",
                                             other_names)))
    # idx_sp are the index values of object@w_full such that
    # object@w_full[idx_sp] = object@w
    idx_sp <- (no_w_full - no_w + 1):no_w_full
    
    # If the user has set a custom kernel we can not use fft.
    if (!is.null(comment(params@pred_kernel))) {
        # pred_kernel is predator species x predator size x prey size
        # We want to multiply this by the prey abundance, which is
        # prey species by prey size, sum over prey size. We use matrix
        # multiplication for this. Then we multiply 1st and 3rd 
        ae <- matrix(params@pred_kernel[, , idx_sp, drop = FALSE],
                     ncol = no_w) %*%
            t(sweep(n, 2, params@w * params@dw, "*"))
        diet[, , 1:no_sp] <- ae
        # Eating the resource
        diet[, , no_sp + 1] <- rowSums(sweep(
            params@pred_kernel, 3, params@dw_full * params@w_full * n_pp, "*"), 
            dims = 2)
    } else {
        prey <- matrix(0, nrow = no_sp + 1, ncol = no_w_full)
        prey[1:no_sp, idx_sp] <- sweep(n, 2, params@w * params@dw, "*")
        prey[no_sp + 1, ] <- n_pp * params@w_full * params@dw_full
        ft <- array(rep(params@ft_pred_kernel_e, times = no_sp + 1) *
                        rep(mvfft(t(prey)), each = no_sp),
                    dim = c(no_sp, no_w_full, no_sp + 1))
        # We now have an array predator x wave number x prey
        # To Fourier transform back we need a matrix of wave number x everything
        ft <- matrix(aperm(ft, c(2, 1, 3)), nrow = no_w_full)
        ae <- array(Re(mvfft(ft, inverse = TRUE) / no_w_full), 
                    dim = c(no_w_full, no_sp, no_sp + 1))
        ae <- ae[idx_sp, , , drop = FALSE]
        ae <- aperm(ae, c(2, 1, 3))
        # Due to numerical errors we might get negative or very small entries that
        # should be 0
        ae[ae < 1e-18] <- 0
        diet[, , 1:(no_sp + 1)] <- ae
    }
    # Multiply by interaction matrix, including resource, and then by 
    # search volume
    inter <- cbind(params@interaction, params@species_params$interaction_resource)
    diet[, , 1:(no_sp + 1)] <- sweep(sweep(diet[, , 1:(no_sp + 1), drop = FALSE],
                                           c(1, 3), inter, "*"), 
                                     c(1, 2), params@search_vol, "*")
    # Add diet from external sources
    diet[, , no_sp + 2] <- params@ext_encounter
    # Add diet from other components
    for (i in seq_along(params@other_encounter)) {
        diet[, , no_sp + 2 + i] <-
            do.call(params@other_encounter[[i]], 
                    list(params = params,
                         n = n, n_pp = n_pp, n_other = n_other,
                         component = names(params@other_encounter)[[i]]))
    }
    
    # Correct for satiation and keep only entries corresponding to fish sizes
    f <- getFeedingLevel(params, n, n_pp)
    fish_mask <- n > 0
    diet <- sweep(diet, c(1, 2), (1 - f) * fish_mask, "*")
    if (proportion) {
        total <- rowSums(diet, dims = 2)
        diet <- sweep(diet, c(1, 2), total, "/")
        diet[is.nan(diet)] <- 0
    }
    return(diet)
}


#' Calculate the SSB of species
#' 
#' Calculates the spawning stock biomass (SSB) through time of the species in
#' the `MizerSim` class. SSB is calculated as the total mass of all mature
#' individuals.
#' 
#' @param object An object of class `MizerParams` or MizerSim`.
#'
#' @return If called with a MizerParams object, a vector with the SSB in
#'   grams for each species in the model. If called with a MizerSim object, an
#'   array (time x species) containing the SSB in grams at each time step
#'   for all species.
#' @export
#' @family summary functions
#' @concept summary_function
#' @examples
#' ssb <- getSSB(NS_sim)
#' ssb[c("1972", "2010"), c("Herring", "Cod")]
getSSB <- function(object) {
    if (is(object, "MizerSim")) {
        sim <- object
        return(apply(sweep(sweep(sim@n, c(2, 3), sim@params@maturity, "*"), 3, 
                           sim@params@w * sim@params@dw, "*"), c(1, 2), sum) )
    }
    if (is(object, "MizerParams")) {
        params <- object
        return(((params@initial_n * params@maturity) %*% 
                    (params@w * params@dw))[, , drop = TRUE])
    }
    stop("'object' should be a MizerParams or a MizerSim object")
}


#' Calculate the total biomass of each species within a size range at each time 
#' step.
#' 
#' Calculates the total biomass through time within user defined size limits.
#' The default option is to use the size range starting at the size specified
#' by the `biomass_cutoff` species parameter, if it is set, or else the full
#' size range of each species. You can specify minimum
#' and maximum weight or length range for the species. Lengths take precedence
#' over weights (i.e. if both min_l and min_w are supplied, only min_l will be
#' used).
#' 
#' @details
#' When no size range arguments are provided, the function checks if the
#' `biomass_cutoff` column exists in the species parameters. If it does,
#' those values are used as the minimum weight for each species. For species
#' with NA values in `biomass_cutoff`, the default minimum weight (smallest
#' weight in the model) is used.
#' 
#' @param object An object of class `MizerParams` or `MizerSim`.
#' @param use_cutoff If TRUE, the `biomass_cutoff` column in the
#'   species parameters is used as the minimum weight for each species (ignoring any
#'   size range arguments in `...`). If FALSE (default), the specified size range
#'   arguments are used, if provided, or the full size range of the species is used.
#' @inheritDotParams get_size_range_array -params
#'
#' @return If called with a MizerParams object, a vector with the biomass in
#'   grams for each species in the model. If called with a MizerSim object, an
#'   array (time x species) containing the biomass in grams at each time step
#'   for all species.
#' @export
#' @family summary functions
#' @concept summary_function
#' @examples
#' biomass <- getBiomass(NS_sim)
#' biomass["1972", "Herring"]
#' biomass <- getBiomass(NS_sim, min_w = 10, max_w = 1000)
#' biomass["1972", "Herring"]
#' 
#' # If species_params contains a biomass_cutoff column, it can be used
#' # as the minimum weight when use_cutoff = TRUE
#' species_params(params)$biomass_cutoff <- species_params(params)$w_mat
#' biomass <- getBiomass(NS_sim, use_cutoff = TRUE)  # Uses biomass_cutoff as min_w
#' biomass["1972", "Herring"]
getBiomass <- function(object, use_cutoff = FALSE, ...) {
    if (is(object, "MizerSim")) {
        sim <- object
        
        if (use_cutoff && "biomass_cutoff" %in% names(sim@params@species_params)) {
            # Use biomass_cutoff as min_w for each species
            biomass_cutoff <- sim@params@species_params$biomass_cutoff
            # Replace NA values with the default minimum weight
            biomass_cutoff[is.na(biomass_cutoff)] <- min(sim@params@w)
            size_range <- get_size_range_array(sim@params, min_w = biomass_cutoff)
        } else {
            size_range <- get_size_range_array(sim@params, ...)
        }
        return(apply(sweep(sweep(sim@n, c(2, 3), size_range, "*"), 3,
                           sim@params@w * sim@params@dw, "*"), c(1, 2), sum))
    }
    if (is(object, "MizerParams")) {
        params <- object
        
        if (use_cutoff && "biomass_cutoff" %in% names(params@species_params)) {
            # Use biomass_cutoff as min_w for each species
            biomass_cutoff <- params@species_params$biomass_cutoff
            # Replace NA values with the default minimum weight
            biomass_cutoff[is.na(biomass_cutoff)] <- min(params@w)
            size_range <- get_size_range_array(params, min_w = biomass_cutoff)
        } else {
            size_range <- get_size_range_array(params, ...)
        }
        return(((params@initial_n * size_range) %*% 
                    (params@w * params@dw))[, , drop = TRUE])
    }
    stop("'object' should be a MizerParams or a MizerSim object")
}


#' Calculate the number of individuals within a size range
#'
#' Calculates the number of individuals within user-defined size limits. The
#' default option is to use the whole size range. You can specify minimum and
#' maximum weight or lengths for the species. Lengths take precedence over
#' weights (i.e. if both min_l and min_w are supplied, only min_l will be used)
#' 
#' @param object An object of class `MizerParams` or `MizerSim`.
#' @inheritDotParams get_size_range_array -params
#'
#' @return If called with a MizerParams object, a vector with the numbers for
#'   each species in the model. If called with a MizerSim object, an array (time
#'   x species) containing the numbers at each time step for all species.
#' @export
#' @family summary functions
#' @concept summary_function
#' @examples
#' numbers <- getN(NS_sim)
#' numbers["1972", "Herring"]
#' # The above gave a huge number, because that included all the larvae.
#' # The number of Herrings between 10g and 1kg is much smaller.
#' numbers <- getN(NS_sim, min_w = 10, max_w = 1000)
#' numbers["1972", "Herring"]
getN <- function(object, ...) {
    if (is(object, "MizerSim")) {
        sim <- object
        size_range <- get_size_range_array(sim@params, ...)
        return(apply(sweep(sweep(sim@n, c(2, 3), size_range, "*"), 3,
                           sim@params@dw, "*"), c(1, 2), sum))
    }
    if (is(object, "MizerParams")) {
        params <- object
        size_range <- get_size_range_array(params, ...)
        return(((params@initial_n * size_range) %*% params@dw)[, , drop = TRUE])
    }
    stop("'object' should be a MizerParams or a MizerSim object")
}


#' Calculate the rate at which biomass of each species is fished by each gear
#'
#' This yield rate is given in grams per year. It is calculated at each time
#' step saved in the MizerSim object. 
#' 
#' For details of how the yield rate is defined see the help page of
#' [getYield()].
#'
#' @param object An object of class `MizerParams` or `MizerSim`.
#'
#' @return If called with a MizerParams object, an array (gear x species) with
#'   the yield rate in grams per year from each gear for each species in the
#'   model. If called with a MizerSim object, an array (time x gear x species)
#'   containing the yield rate at each time step.
#' @export
#' @family summary functions
#' @concept summary_function
#' @seealso [getYield()]
#' @examples
#' yield <- getYieldGear(NS_sim)
#' yield["1972", "Herring", "Herring"]
#' # (In this example MizerSim object each species was set up with its own gear)
getYieldGear <- function(object) {
    if (is(object, "MizerSim")) {
        sim <- object
        biomass <- sweep(sim@n, 3, sim@params@w * sim@params@dw, "*")
        f_gear <- getFMortGear(sim)
        return(apply(sweep(f_gear, c(1, 3, 4), biomass, "*"), c(1, 2, 3), sum))
    }
    if (is(object, "MizerParams")) {
        params <- object
        biomass <- sweep(params@initial_n, 2, params@w * params@dw, "*")
        f_gear <- getFMortGear(params)
        return(apply(sweep(f_gear, c(2, 3), biomass, "*"), c(1, 2), sum))
    }
    stop("'object' should be a MizerParams or a MizerSim object")
}


#' Calculate the rate at which biomass of each species is fished
#'
#' This yield rate is given in grams per year. It is calculated at each time
#' step saved in the MizerSim object.
#' 
#' @details
#' The yield rate \eqn{y_i(t)} for species \eqn{i} at time \eqn{t} is defined as
#' \deqn{y_i(t)=\int\mu_{f.i}(w, t)N_i(w, t)w dw}
#' where \eqn{\mu_{f.i}(w, t)} is the fishing mortality of an individual of
#' species \eqn{i} and weight \eqn{w} at time \eqn{t} and \eqn{N_i(w, t)} is the
#' abundance density of such individuals.  The factor of \eqn{w} converts the
#' abundance density into a biomass density and the integral aggregates the
#' contribution from all sizes.
#' 
#' The total catch in a time period from \eqn{t_1} to  \eqn{t_2} is the integral
#' of the yield rate over that period:
#' \deqn{C = \int_{t_1}^{t2}y_i(t)dt}
#' In practice, as the yield rate is only available
#' at the saved times, one can only approximate this integral by averaging over
#' the available yield rates during the time period and multiplying by the time
#' period. The less the yield changes between the saved values, the more
#' accurate this approximation is. So the approximation can be improved by
#' saving simulation results at smaller intervals, using the `t_save` argument
#' to [project()]. But this is only a concern if abundances change quickly
#' during the time period of interest.
#'
#' @param object An object of class `MizerParams` or `MizerSim`.
#'
#' @return If called with a MizerParams object, a vector with the yield rate in
#'   grams per year for each species in the model. If called with a MizerSim
#'   object, an array (time x species) containing the yield rate at each time
#'   step for all species.
#' @export
#' @family summary functions
#' @concept summary_function
#' @seealso [getYieldGear()]
#' @examples
#' yield <- getYield(NS_sim)
#' yield[c("1972", "2010"), c("Herring", "Cod")]
#' 
#' # Running simulation for another year, saving intermediate time steps
#' params <- setInitialValues(getParams(NS_sim), NS_sim)
#' sim <- project(params, t_save = 0.1, t_max = 1, 
#'                t_start = 2010, progress_bar = FALSE)
#' # The yield rate for Herring decreases during the year
#' getYield(sim)[, "Herring"]
#' # We get the total catch in the year by averaging over the year
#' sum(getYield(sim)[1:10, "Herring"] / 10)
getYield <- function(object) {
    if (is(object, "MizerSim")) {
        sim <- object
        biomass <- sweep(sim@n, 3, sim@params@w * sim@params@dw, "*")
        f <- getFMort(sim, drop = FALSE)
        return(apply(f * biomass, c(1, 2), sum))
    }
    if (is(object, "MizerParams")) {
        params <- object
        biomass <- sweep(params@initial_n, 2, params@w * params@dw, "*")
        f <- getFMort(params, drop = FALSE)
        return(apply(f * biomass, 1, sum))
    }
    stop("'object' should be a MizerParams or a MizerSim object")
}


#' Get growth curves giving weight as a function of age
#' 
#' @param object MizerSim or MizerParams object. If given a
#'   \linkS4class{MizerSim} object, uses the growth rates at the final time of a
#'   simulation to calculate the size at age. If given a
#'   \linkS4class{MizerParams} object, uses the initial growth rates instead.
#' @inheritParams valid_species_arg
#' @param max_age The age up to which to run the growth curve. Default is 20.
#' @param percentage Boolean value. If TRUE, the size is given as a percentage
#'   of the maximal size.
#'
#' @return An array (species x age) containing the weight in grams.
#' @export
#' @family summary functions
#' @examples
#' growth_curves <- getGrowthCurves(NS_params, species = c("Cod", "Haddock"))
#' str(growth_curves)
#' 
#' library(ggplot2)
#' ggplot(melt(growth_curves)) +
#'   geom_line(aes(Age, value)) +
#'   facet_wrap(~ Species, scales = "free") +
#'   ylab("Size[g]") + xlab("Age[years]")
getGrowthCurves <- function(object, 
                            species = NULL,
                            max_age = 20,
                            percentage = FALSE) {
    if (is(object, "MizerSim")) {
        params <- object@params
        params <- setInitialValues(params, object)
    } else if (is(object, "MizerParams")) {
        params <- validParams(object)
    } else {
        stop("The first argument to `getGrowthCurves()` must be a ",
             "MizerParams or a MizerSim object.")
    }
    species <- valid_species_arg(params, species)
    # reorder list of species to coincide with order in params
    idx <- which(params@species_params$species %in% species)
    species <- params@species_params$species[idx]
    age <- seq(0, max_age, length.out = 50)
    ws <- array(dim = c(length(species), length(age)),
                dimnames = list(Species = species, Age = age))
    g <- getEGrowth(params)
    for (j in seq_along(species)) {
        i <- idx[j]
        g_fn <- stats::approxfun(c(params@w, params@species_params$w_max[[i]]),
                                 c(g[i, ], 0))
        myodefun <- function(t, state, parameters) {
            return(list(g_fn(state)))
        }
        ws[j, ] <- deSolve::ode(y = params@w[params@w_min_idx[i]], 
                                times = age, func = myodefun)[, 2]
        if (percentage) {
            ws[j, ] <- ws[j, ] / params@species_params$w_max[i] * 100
        }
    }
    return(ws)
}

#' Get size range array
#' 
#' Helper function that returns an array (species x size) of boolean values
#' indicating whether that size bin is within the size limits specified by the
#' arguments. Either the size limits can be the same for all species or they
#' can be specified as vectors with one value for each species in the model.
#' 
#' @param params MizerParams object
#' @param min_w Smallest weight in size range. Defaults to smallest weight in
#'   the model.
#' @param max_w Largest weight in size range. Defaults to largest weight in the
#'   model.
#' @param min_l Smallest length in size range. If supplied, this takes
#'   precedence over `min_w`.
#' @param max_l Largest length in size range. If supplied, this takes precedence
#'   over `max_w`.
#' @param ... Unused
#'   
#' @return Boolean array (species x size)
#' 
#' @section Length to weight conversion:
#' If `min_l` is specified there is no need to specify `min_w` and so on.
#' However, if a length is specified (minimum or maximum) then it is necessary
#' for the species parameter data.frame to include the parameters `a` and `b`
#' that determine the relation between length \eqn{l} and weight \eqn{w} by
#' \deqn{w = a l^b.}
#' 
#' It is possible to mix length and weight constraints, e.g. by supplying a
#' minimum weight and a maximum length, but this must be done the same for
#' all species. The default values are the minimum and
#' maximum weights of the spectrum, i.e., the full range of the size spectrum is
#' used.
#' @export
#' @concept helper
get_size_range_array <- function(params, min_w = min(params@w), 
                                 max_w = max(params@w), 
                                 min_l = NULL, max_l = NULL, ...) {
    no_sp <- nrow(params@species_params)
    if (!is.null(min_l) || !is.null(max_l)) {
        if (any(!c("a", "b") %in% names(params@species_params))) {
            stop("species_params slot must have columns 'a' and 'b' for ",
                 "length-weight conversion")
        }
        if (anyNA(params@species_params[["a"]]) ||
            anyNA(params@species_params[["a"]])) {
            stop("There must be no NAs in the species_params columns 'a' and 'b'.")
        }
    }
    if (!is.null(min_l)) {
        if (length(min_l) != 1 && length(min_l) != no_sp) {
            stop("min_l must be a single number or a vector with one value for each species.")
        }
        min_w <- params@species_params[["a"]] * 
            min_l ^ params@species_params[["b"]]
    }
    if (!is.null(max_l)) {
        if (length(max_l) != 1 && length(max_l) != no_sp) {
            stop("max_l must be a single number or a vector with one value for each species.")
        }
        max_w <- params@species_params[["a"]] *
            max_l ^ params@species_params[["b"]]
    }
    if (length(min_w) == 1) {
        min_w <- rep(min_w, no_sp)
    }
    if (length(max_w) == 1) {
        max_w <- rep(max_w, no_sp)
    }
    if (length(min_w) != no_sp || length(max_w) != no_sp) {
        stop("min_w and max_w must be a single number of a vector with one
             value for each species.")
    }
    if (!all(min_w < max_w)) stop("min_w must be less than max_w")
    min_n <- t(sapply(min_w, function(x) params@w >= x))
    max_n <- t(sapply(max_w, function(x) params@w <= x))
    size_n <- min_n & max_n
    dimnames(size_n) <- list(sp = params@species_params$species, 
                             w = signif(params@w, 3)) 
    return(size_n)
}


#### summary for MizerParams ####
#' Summarize MizerParams object 
#'
#' Outputs a general summary of the structure and content of the object
#' @param object A `MizerParams` object.
#' @param ... Other arguments (currently not used).
#' @return The MizerParams object, invisibly
#' @export
#' @concept summary_function
#' @examples
#' summary(NS_params)
setMethod("summary", signature(object = "MizerParams"), function(object, ...) {
    params <- validParams(object)
    cat("An object of class \"", as.character(class(params)), "\" \n", sep = "")
    cat("Consumer size spectrum:\n")
    cat("\tminimum size:\t", signif(min(params@w)), "\n", sep = "")
    cat("\tmaximum size:\t", signif(max(params@w)), "\n", sep = "")
    cat("\tno. size bins:\t", length(params@w), "\n", sep = "")
    # Length of background? 
    cat("Resource size spectrum:\n")
    cat("\tminimum size:\t", signif(min(params@w_full)), "\n", sep = "")
    cat("\tmaximum size:\t", signif(max(params@w_full[params@initial_n_pp > 0])), 
        "\n", sep = "")
    cat("\tno. size bins:\t", length(params@w_full[params@initial_n_pp > 0]), 
        "\t(", length(params@w_full), " size bins in total)\n", sep = "")
    cat("Species details:\n")
    sel_params <- intersect(c("species", "w_max", "w_mat", "w_min", "f0", "fc", 
                              "age_mat", "beta", "sigma"),
                            names(params@species_params))
    sp <- params@species_params[, sel_params]
    rownames(sp) <- NULL
    print(sp)
    cat("\nFishing gear details:\n")
    cat(sprintf("%-13s %s %s", "Gear", "Effort", " Target species"), "\n",
        "----------------------------------\n")
    gears <- dimnames(params@catchability)$gear
    for (i in seq_len(dim(params@catchability)[1])) {
        splist <- dimnames(params@catchability)$sp[params@catchability[i, ] > 0]
        cat(sprintf("%-14s %1.2f   %s",
                gears[i], params@initial_effort[[gears[i]]],
                toString(splist)), "\n")
    }
    invisible(params)
})


#### summary for MizerSim ####
#' Summarize MizerSim object 
#'
#' Outputs a general summary of the structure and content of the object
#' @param object A `MizerSim` object.
#' @param ... Other arguments (currently not used).
#' @return The MizerSim object, invisibly
#' @export
#' @concept summary_function
#' @examples
#' summary(NS_sim)
setMethod("summary", signature(object = "MizerSim"), function(object, ...) {
    cat("An object of class \"", as.character(class(object)), "\" \n", sep = "")
    cat("Parameters:\n")
    summary(object@params)
    cat("Simulation parameters:\n")
    # Need to store t_max and dt in a description slot? Or just in simulation 
    # time parameters? Like a list?
    cat("\tTime period: ", 
        min(as.numeric(dimnames(object@n)$time)), " to ",
        max(as.numeric(dimnames(object@n)$time)), 
        "\n", sep = "")
    cat("\tOutput stored every ", 
        as.numeric(dimnames(object@n)$time)[2] - 
            as.numeric(dimnames(object@n)$time)[1], " years\n", sep = "")
    invisible(object)
})

# Indicator functions ####
#' Description of indicator functions
#' 
#' Mizer provides a range of functions to calculate indicators 
#' from a MizerSim object.
#' 
#' A list of available indicator functions for MizerSim objects is given in the table below
#' \tabular{lll}{
#'   Function \tab Returns \tab Description \cr
#'   [getProportionOfLargeFish()] \tab A vector with values at each time step. \tab Calculates the proportion of large fish through time. The threshold value can be specified. It is possible to calculation the proportion of large fish based on either length or weight. \cr
#'   [getMeanWeight()] \tab A vector with values at each saved time step. \tab The mean weight of the community through time. This is calculated as the total biomass of the community divided by the total abundance. \cr
#'   [getMeanMaxWeight()] \tab Depends on the measure argument. If measure = “both” then you get a matrix with two columns, one with values by numbers, the other with values by biomass at each saved time step. If measure = “numbers” or “biomass” you get a vector of the respective values at each saved time step \tab The mean maximum weight of the community through time. This can be calculated by numbers or by biomass. See the help file for more details. \cr
#'   [getCommunitySlope()] \tab A data.frame with four columns: time step, slope, intercept and the coefficient of determination. \tab Calculates the slope of the community abundance spectrum through time by performing a linear regression on the logged total numerical abundance and logged body size. \cr
#' }
#'
#' @seealso [summary_functions], [plotting_functions]
#' @name indicator_functions
NULL


#' Calculate the proportion of large fish
#' 
#' Calculates the proportion of large fish through time in the `MizerSim`
#' class within user defined size limits. The default option is to use the whole
#' size range. You can specify minimum and maximum size ranges for the species
#' and also the threshold size for large fish. Sizes can be expressed as weight
#' or size. Lengths take precedence over weights (i.e. if both min_l and min_w
#' are supplied, only min_l will be used). You can also specify the species to
#' be used in the calculation. This function can be used to calculate the Large
#' Fish Index. The proportion is based on either abundance or biomass.
#' 
#' @inheritParams getMeanWeight
#' @param threshold_w the size used as the cutoff between large and small fish.
#'   Default value is 100.
#' @param threshold_l the size used as the cutoff between large and small fish.
#' @param biomass_proportion a boolean value. If TRUE the proportion calculated
#'   is based on biomass, if FALSE it is based on numbers of individuals.
#'   Default is TRUE.
#' @inheritDotParams get_size_range_array -params
#'   
#' @return A vector containing the proportion of large fish through time
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
getProportionOfLargeFish <- function(sim, 
                                     species = NULL, 
                                     threshold_w = 100, threshold_l = NULL, 
                                     biomass_proportion = TRUE, ...) {
    assert_that(is(sim, "MizerSim"))
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

    #lfi = data.frame(time = as.numeric(dimnames(sim@n)$time), 
    #                 proportion = 1 - (upto_threshold_measure / total_measure))
    #return(lfi)
    return(1 - (upto_threshold_measure / total_measure))
}


#' Calculate the mean weight of the community
#'
#' Calculates the mean weight of the community through time. This is simply the
#' total biomass of the community divided by the abundance in numbers. You can
#' specify minimum and maximum weight or length range for the species. Lengths
#' take precedence over weights (i.e. if both min_l and min_w are supplied, only
#' min_l will be used). You can also specify the species to be used in the
#' calculation.
#'
#' @param sim A \linkS4class{MizerSim} object
#' @inheritParams valid_species_arg
#' @inheritDotParams get_size_range_array -params
#'
#' @return A vector containing the mean weight of the community through time
#' @export
#' @family functions for calculating indicators
#' @concept summary_function
#' @examples
#' mean_weight <- getMeanWeight(NS_sim)
#' years <- c("1967", "2010")
#' mean_weight[years]
#' getMeanWeight(NS_sim, species = c("Herring", "Sprat", "N.pout"))[years]
#' getMeanWeight(NS_sim, min_w = 10, max_w = 5000)[years]
getMeanWeight <- function(sim, species = NULL, ...) {
    assert_that(is(sim, "MizerSim"))
    species <- valid_species_arg(sim, species)
    n_species <- getN(sim, ...)
    biomass_species <- getBiomass(sim, ...)
    n_total <- apply(n_species[, species, drop = FALSE], 1, sum)
    biomass_total <- apply(biomass_species[, species, drop = FALSE], 1, sum)
    return(biomass_total / n_total)
}


#' Calculate the mean maximum weight of the community
#'
#' Calculates the mean maximum weight of the community through time. This can be
#' calculated by numbers or biomass. The calculation is the sum of the w_max *
#' abundance of each species, divided by the total abundance community, where
#' abundance is either in biomass or numbers. You can specify minimum and
#' maximum weight or length range for the species. Lengths take precedence over
#' weights (i.e. if both min_l and min_w are supplied, only min_l will be used).
#' You can also specify the species to be used in the calculation.
#'
#' @inheritParams getMeanWeight
#' @param measure The measure to return. Can be 'numbers', 'biomass' or 'both'
#' @inheritDotParams get_size_range_array -params
#'
#' @return Depends on the `measure` argument. If \code{measure = “both”}
#'   then you get a matrix with two columns, one with values by numbers, the
#'   other with values by biomass at each saved time step. If \code{measure =
#'   “numbers”} or \code{“biomass”} you get a vector of the respective values at
#'   each saved time step.
#' @export
#' @family functions for calculating indicators
#' @concept summary_function
#' @examples
#' mmw <- getMeanMaxWeight(NS_sim)
#' years <- c("1967", "2010")
#' mmw[years, ]
#' getMeanMaxWeight(NS_sim, species=c("Herring","Sprat","N.pout"))[years, ]
#' getMeanMaxWeight(NS_sim, min_w = 10, max_w = 5000)[years, ]
getMeanMaxWeight <- function(sim, species = NULL, 
                             measure = "both", ...) {
    assert_that(is(sim, "MizerSim"))
    if (!(measure %in% c("both", "numbers", "biomass"))) {
        stop("measure must be one of 'both', 'numbers' or 'biomass'")
    }
    species <- valid_species_arg(sim, species)
    n_species <- getN(sim, ...)
    biomass_species <- getBiomass(sim, ...)
    n_winf <- apply(sweep(n_species, 2, sim@params@species_params$w_max, "*")[,species, drop = FALSE], 1, sum)
    biomass_winf <- apply(sweep(biomass_species, 2, sim@params@species_params$w_max,"*")[, species, drop = FALSE], 1, sum)
    mmw_numbers <- n_winf / apply(n_species, 1, sum)
    mmw_biomass <- biomass_winf / apply(biomass_species, 1, sum)
    if (measure == "numbers")
        return(mmw_numbers)
    if (measure == "biomass")
        return(mmw_biomass)
    if (measure == "both")
        return(cbind(mmw_numbers, mmw_biomass)) 
}


#' Calculate the slope of the community abundance
#'
#' Calculates the slope of the community abundance through time by performing a
#' linear regression on the logged total numerical abundance at weight and
#' logged weights (natural logs, not log to base 10, are used). You can specify
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
#' @return A data.frame with four columns: time step, slope, intercept and the
#'   coefficient of determination R^2.
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
getCommunitySlope <- function(sim, species = NULL,
                              biomass = TRUE, ...) {
    assert_that(is(sim, "MizerSim"))
    species <- valid_species_arg(sim, species)
    size_range <- get_size_range_array(sim@params, ...)
    # set entries for unwanted sizes to zero and sum over wanted species, giving
    # array (time x size)
    total_n <-
        apply(sweep(sim@n, c(2, 3), size_range, "*")[, species, , drop = FALSE],
              c(1, 3), sum)
    # numbers or biomass?
    if (biomass)
        total_n <- sweep(total_n, 2, sim@params@w, "*")
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
    slope <- slope[, -1]
    return(slope)
}
