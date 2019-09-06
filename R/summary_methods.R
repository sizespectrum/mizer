# Summary functions for mizer package

# Copyright 2012 Finlay Scott and Julia Blanchard.
# Copyright 2018 Gustav Delius and Richard Southwell.
# Development has received funding from the European Commission's Horizon 2020 
# Research and Innovation Programme under Grant Agreement No. 634495 
# for the project MINOUW (http://minouw-project.eu/).
# Distributed under the GPL 3 or later 
# Maintainer: Gustav Delius, University of York, <gustav.delius@york.ac.uk>

# Soundtrack: The Definitive Lead Belly

#' Description of summary functions
#' 
#' Mizer provides a range of functions to extract information about a simulation
#' from a MizerSim object.
#'
#' A list of available summary functions is given in the table below.
#' \tabular{lll}{
#'   Function \tab Returns \tab Description \cr
#'   \code{\link{getDiet}} \tab Three dimensional array (predator x size x prey) \tab Diet of predator at size, resolved by prey species \cr
#'   \code{\link{getSSB}} \tab Two dimensional array (time x species) \tab Total Spawning Stock Biomass (SSB) of each species through time where SSB is calculated as the sum of weight of all mature individuals. \cr
#'   \code{\link{getBiomass}} \tab Two dimensional array (time x species) \tab Total biomass of each species through time. \cr
#'   \code{\link{getN}} \tab Two dimensional array (time x species) \tab Total abundance of each species through time. \cr
#'   \code{\link{getFeedingLevel}} \tab Three dimensional array (time x species x size) \tab Feeding level of each species by size through time. \cr
#'   \code{\link{getM2}} \tab Three dimensional array (time x species x size) \tab The predation mortality imposed on each species by size through time. \cr
#'   \code{\link{getFMort}} \tab Three dimensional array (time x species x size) \tab Total fishing mortality on each species by size through time. \cr
#'   \code{\link{getFMortGear}} \tab Four dimensional array (time x gear x species x size) \tab Fishing mortality on each species by each gear at size through time. \cr
#'   \code{\link{getYieldGear}} \tab Three dimensional array (time x gear x species) \tab Total yield by gear and species through time. \cr
#'   \code{\link{getYield}} \tab Two dimensional array (time x species) \tab Total yield of each species across all gears through time. \cr
#' }
#' 
#' A list of available indicator functions for MizerSim objects is given in the table below
#' \tabular{lll}{
#'   Function \tab Returns \tab Description \cr
#'   \code{\link{getProportionOfLargeFish}} \tab A vector with values at each time step. \tab Calculates the proportion of large fish through time. The threshold value can be specified. It is possible to calculation the proportion of large fish based on either length or weight. \cr
#'   \code{\link{getMeanWeight}} \tab A vector with values at each saved time step. \tab The mean weight of the community through time. This is calculated as the total biomass of the community divided by the total abundance. \cr
#'   \code{\link{getMeanMaxWeight}} \tab Depends on the measure argument. If measure = “both” then you get a matrix with two columns, one with values by numbers, the other with values by biomass at each saved time step. If measure = “numbers” or “biomass” you get a vector of the respective values at each saved time step \tab The mean maximum weight of the community through time. This can be calculated by numbers or by biomass. See the help file for more details. \cr
#'   \code{\link{getCommunitySlope}} \tab A data.frame with four columns: time step, slope, intercept and the coefficient of determination. \tab Calculates the slope of the community abundance spectrum through time by performing a linear regression on the logged total numerical abundance and logged body size. \cr
#' }
#'
#' @seealso \code{\link{plotting_functions}}
#' @name summary_functions
NULL

#' Get diet of predator at size, resolved by prey species
#'
#' Calculates the rate at which a predator of a particular species and size
#' consumes biomass of each prey species, plankton and resources.
#' 
#' This function performs the same integration as
#' \code{\link{getEncounter}} but does not aggregate over prey species, and
#' multiplies by (1-feeding_level) to get the consumed biomass rather than the
#' available biomass. Outside the range of sizes for a predator species the
#' returned rate is zero.
#'
#' @inheritParams getEncounter
#' @param proportion If TRUE (default) the function returns the diet as a
#'   proportion of the total consumption rate. If FALSE it returns the 
#'   consumption rate in grams.
#' 
#' @return An array (predator species  x predator size x 
#'   (prey species + plankton + resources) )
#' @export
#' @family summary functions
#' @concept summary_function
getDiet <- function(params, 
                    n = params@initial_n, 
                    n_pp = params@initial_n_pp,
                    B = params@initial_B,
                    proportion = TRUE) {
    # The code is based on that for getEncounter()
    assert_that(is(params, "MizerParams"),
                is.array(n),
                is.vector(n_pp),
                is.vector(B))
    species <- params@species_params$species
    no_sp <- length(species)
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)
    no_res <- length(params@resource_dynamics)
    resource_names <- names(params@resource_dynamics)
    assert_that(identical(dim(n), c(no_sp, no_w)),
                length(n_pp) == no_w_full,
                length(B) == no_res)
    diet <- array(0, dim = c(no_sp, no_w, no_sp + 1 + no_res),
                  dimnames = list("predator" = species,
                                  "w" = dimnames(params@initial_n)$w,
                                  "prey" = c(as.character(species), 
                                             "Plankton", 
                                             resource_names)))
    # idx_sp are the index values of object@w_full such that
    # object@w_full[idx_sp] = object@w
    idx_sp <- (no_w_full - no_w + 1):no_w_full
    
    # If the feeding kernel does not have a fixed predator/prey mass ratio
    # then the integral is not a convolution integral and we can not use fft.
    if (length(params@ft_pred_kernel_e) == 1) {
        # pred_kernel is predator species x predator size x prey size
        # We want to multiply this by the prey abundance, which is
        # prey species by prey size, sum over prey size. We use matrix
        # multiplication for this. Then we multiply 1st and 3rd 
        ae <- matrix(params@pred_kernel[, , idx_sp, drop = FALSE],
                     ncol = no_w) %*%
            t(sweep(n, 2, params@w * params@dw, "*"))
        diet[, , 1:no_sp] <- ae
        # Eating the plankton
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
    # Multiply by interaction matrix, including plankton, and then by 
    # search volume
    inter <- cbind(params@interaction, params@species_params$interaction_p)
    diet[, , 1:(no_sp + 1)] <- sweep(sweep(diet[, , 1:(no_sp + 1), drop = FALSE],
                                           c(1, 3), inter, "*"), 
                                     c(1, 2), params@search_vol, "*")
    # Add diet from resources
    if (no_res > 0) {
        diet[, , (no_sp + 2):(no_sp + 1 + no_res)] <- 
            aperm(sweep(params@rho, 2, B, "*"), c(1, 3, 2))
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
#' the \code{MizerSim} class. SSB is calculated as the total mass of all mature
#' individuals.
#' 
#' @param sim An object of class \code{MizerSim}.
#'   
#' @return An array containing the SSB (time x species)
#' @export
#' @family summary functions
#' @concept summary_function
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- set_multispecies_model(NS_species_params_gears, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' getSSB(sim)
#' }
getSSB <- function(sim) {
    assert_that(is(sim, "MizerSim"))
    ssb <- apply(sweep(sweep(sim@n, c(2,3), sim@params@maturity,"*"), 3, 
                       sim@params@w * sim@params@dw, "*"), c(1, 2), sum) 
    return(ssb)
}


#' Calculate the total biomass of each species within a size range at each time 
#' step.
#' 
#' Calculates the total biomass through time of the species in the
#' \code{MizerSim} class within user defined size limits. The default option is
#' to use the whole size range. You can specify minimum and maximum weight or
#' length range for the species. Lengths take precedence over weights (i.e. if
#' both min_l and min_w are supplied, only min_l will be used).
#' 
#' @param sim An object of class \code{MizerSim}.
#' @inheritDotParams get_size_range_array -params
#'
#' @return An array containing the biomass (time x species)
#' @export
#' @family summary functions
#' @concept summary_function
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- set_multispecies_model(NS_species_params_gears, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' getBiomass(sim)
#' getBiomass(sim, min_w = 10, max_w = 1000)
#' }
getBiomass <- function(sim, ...) {
    assert_that(is(sim, "MizerSim"))
    size_range <- get_size_range_array(sim@params, ...)
    biomass <- apply(sweep(sweep(sim@n, c(2, 3), size_range, "*"), 3,
                           sim@params@w * sim@params@dw, "*"), c(1, 2), sum)
    return(biomass)
}


#' Calculate the number of individuals within a size range
#'
#' Calculates the number of individuals within user-defined size limits, for
#' each time and each species in the \code{MizerSim} object. The default option
#' is to use the whole size range. You can specify minimum and maximum weight or
#' lengths for the species. Lengths take precedence over weights (i.e. if both
#' min_l and min_w are supplied, only min_l will be used)
#' 
#' @param sim An object of class \code{MizerSim}.
#' @inheritDotParams get_size_range_array -params
#'
#' @return An array containing the total numbers (time x species)
#' @export
#' @family summary functions
#' @concept summary_function
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- set_multispecies_model(NS_species_params_gears, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' getN(sim)
#' getN(sim, min_w = 10, max_w = 1000)
#' }
getN <- function(sim, ...) {
    assert_that(is(sim, "MizerSim"))
    size_range <- get_size_range_array(sim@params, ...)
    n <- apply(sweep(sweep(sim@n, c(2, 3), size_range, "*"), 3,
                     sim@params@dw, "*"), c(1, 2), sum)
    return(n)
}


#' Calculate the total yield per gear and species
#'
#' Calculates the total yield per gear and species at each simulation
#' time step.
#'
#' @param sim An object of class \code{MizerSim}.
#'
#' @return An array containing the total yield (time x gear x species)
#' @export
#' @family summary functions
#' @concept summary_function
#' @seealso \code{\link{getYield}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- set_multispecies_model(NS_species_params_gears, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' getYieldGear(sim)
#' }
getYieldGear <- function(sim) {
    biomass <- sweep(sim@n, 3, sim@params@w * sim@params@dw, "*")
    f_gear <- getFMortGear(sim)
    yield_species_gear <- apply(sweep(f_gear, c(1, 3, 4), biomass, "*"),
                                c(1, 2, 3), sum)
    return(yield_species_gear)
}


#' Calculate the total yield of each species
#'
#' Calculates the total yield of each species across all gears at each
#' simulation time step.
#'
#' @param sim An object of class \code{MizerSim}.
#'
#' @return An array containing the total yield (time x species)
#' @export
#' @family summary functions
#' @concept summary_function
#' @seealso \code{\link{getYieldGear}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- set_multispecies_model(NS_species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=10)
#' y <- getYield(sim)
#' }
getYield <- function(sim) {
    # biomass less the first time step
    yield_gear_species <- getYieldGear(sim)
    return(apply(yield_gear_species, c(1, 3), sum))
}


#' Get growth curves giving weight as a function of age
#' 
#' If given a \linkS4class{MizerSim} object, uses the growth rates at the final
#' time of a simulation to calculate the size at age. If given a
#' \linkS4class{MizerParams} object, uses the initial growth rates instead.
#' 
#' @param object MizerSim or MizerParams object
#' @param species Name or vector of names of the species to be included. By
#'   default all species are included.
#' @param max_age The age up to which to run the growth curve. Default is 20.
#' @param percentage Boolean value. If TRUE, the size is given as a percentage
#'   of the maximal size.
#'
#' @return An array (species x age) containing the weight in grams.
#' @export
#' @family summary functions
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- suppressMessages(set_multispecies_model(NS_species_params_gears, inter))
#' getGrowthCurves(params)
#' sim <- project(params, effort=1, t_max = 20, t_save = 2, progress_bar = FALSE)
#' getGrowthCurves(sim, max_age = 24)
#' }
getGrowthCurves <- function(object, 
                            species,
                            max_age = 20,
                            percentage = FALSE) {
    if (is(object, "MizerSim")) {
        params <- object@params
        t <- dim(object@n)[1]
        n <- object@n[t, , ]
        n_pp <- object@n_pp[t, ]
    } else if (is(object, "MizerParams")) {
        params <- object
        n <- object@initial_n
        n_pp <- object@initial_n_pp
    }
    if (missing(species)) {
        species <- dimnames(n)$sp
    }
    # reorder list of species to coincide with order in params
    idx <- which(dimnames(n)$sp %in% species)
    species <- dimnames(n)$sp[idx]
    age <- seq(0, max_age, length.out = 50)
    ws <- array(dim = c(length(species), length(age)),
                dimnames = list(Species = species, Age = age))
    g <- getEGrowth(params, n, n_pp)
    for (j in 1:length(species)) {
        i <- idx[j]
        g_fn <- stats::approxfun(params@w, g[i, ])
        myodefun <- function(t, state, parameters){
            return(list(g_fn(state)))
        }
        ws[j, ] <- deSolve::ode(y = params@w[params@w_min_idx[i]], 
                                times = age, func = myodefun)[, 2]
        if (percentage) {
            ws[j, ] <- ws[j, ] / params@species_params$w_inf[i] * 100
        }
    }
    return(ws)
}

#' Get size range array
#' 
#' Helper function that returns an array (species x size) of boolean values
#' indicating whether that size bin is within the size limits specified by the
#' arguments.
#' 
#' @param params MizerParams object
#' @param min_w Smallest weight in size range. Defaults to smallest weight in
#'   the model.
#' @param max_w Largest weight in size range. Defaults to largest weight in the
#'   model.
#' @param min_l Smallest length in size range. If supplied, this takes
#'   precedence over \code{min_w}.
#' @param max_l Largest length in size range. If supplied, this takes precedence
#'   over \code{max_w}.
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
#' minimum weight and a maximum length. The default values are the minimum and
#' maximum weights of the spectrum, i.e. the full range of the size spectrum is
#' used.
#' @export
#' @keywords internal
#' @concept helper
get_size_range_array <- function(params, min_w = min(params@w), 
                                 max_w = max(params@w), 
                                 min_l = NULL, max_l = NULL, ...) {
    no_sp <- nrow(params@species_params)
    if (!is.null(min_l) | !is.null(max_l))
        if (any(!c("a","b") %in% names(params@species_params)))
            stop("species_params slot must have columns 'a' and 'b' for length-weight conversion")
    if (!is.null(min_l))
        min_w <- params@species_params$a * min_l ^ params@species_params$b
    else min_w <- rep(min_w,no_sp)
    if (!is.null(max_l))
        max_w <- params@species_params$a * max_l ^ params@species_params$b
    else max_w <- rep(max_w,no_sp)
    if (!all(min_w < max_w))
        stop("min_w must be less than max_w")
    min_n <- plyr::aaply(min_w, 1, function(x) params@w >= x, .drop = FALSE)
    max_n <- plyr::aaply(max_w, 1, function(x) params@w <= x, .drop = FALSE)
    size_n <- min_n & max_n
    # Add dimnames?
    dimnames(size_n) <- list(sp = params@species_params$species, w = signif(params@w,3)) 
    return(size_n)
}

# TODO: Check documentation for summary
#### summary for MizerParams ####
#' Summarize MizerParams object 
#'
#' Outputs a general summary of the structure and content of the object
#' @param object A \code{MizerParams} object.
#' @param ... Other arguments (currently not used).
#'
#' @export
#' @concept summary_function
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- set_multispecies_model(NS_species_params_gears,inter)
#' summary(params)
#' }
setMethod("summary", signature(object = "MizerParams"), function(object, ...) {
    cat("An object of class \"", as.character(class(object)), "\" \n", sep = "")
    cat("Consumer size spectrum:\n")
    cat("\tminimum size:\t", signif(min(object@w)), "\n", sep = "")
    cat("\tmaximum size:\t", signif(max(object@w)), "\n", sep = "")
    cat("\tno. size bins:\t", length(object@w), "\n", sep = "")
    # Length of background? 
    cat("Plankton size spectrum:\n")
    cat("\tminimum size:\t", signif(min(object@w_full)), "\n", sep = "")
    cat("\tmaximum size:\t", signif(max(object@w_full)), "\n", sep = "")
    cat("\tno. size bins:\t", length(object@w_full), "\n", sep = "")
    # w range - min, max, number of w
    # w background min max
    # no species and names and wInf,  - not all these wMat, beta, sigma
    # no gears, gear names catching what
    cat("Species details:\n")
    #cat("\tSpecies\t\tw_inf\n")
    #	for (i in 1:nrow(object@species_params))
    #	    cat("\t",as.character(object@species_params$species)[i], "\t\t ",signif(object@species_params$w_inf[i],3), "\n", sep = "")
    print(object@species_params[,c("species","w_inf","w_mat","beta","sigma")])
    cat("Fishing gear details:\n")
    cat("\tGear\t\t\tTarget species\n")
    for (i in 1:dim(object@catchability)[1]){
        cat("\t",dimnames(object@catchability)$gear[i], "\t\t",dimnames(object@catchability)$sp[object@catchability[i,]>0], "\n", sep=" ") 
    }
})


#### summary for MizerSim ####
#' Summarize MizerSim object 
#'
#' Outputs a general summary of the structure and content of the object
#' @param object A \code{MizerSim} object.
#' @param ... Other arguments (currently not used).
#' @export
#' @concept summary_function
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- set_multispecies_model(NS_species_params_gears,inter)
#' sim <- project(params, effort=1, t_max=5)
#' summary(sim)
#' }
setMethod("summary", signature(object = "MizerSim"), function(object, ...){
    cat("An object of class \"", as.character(class(object)), "\" \n", sep = "")
    cat("Parameters:\n")
    summary(object@params)
    cat("Simulation parameters:\n")
    # Need to store t_max and dt in a description slot? Or just in simulation time parameters? Like a list?
    cat("\tFinal time step: ", max(as.numeric(dimnames(object@n)$time)), "\n", sep = "")
    cat("\tOutput stored every ", as.numeric(dimnames(object@n)$time)[2] - as.numeric(dimnames(object@n)$time)[1], " time units\n", sep = "")
})


#' Calculate the proportion of large fish
#' 
#' Calculates the proportion of large fish through time in the \code{MizerSim}
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
#' @return An array containing the proportion of large fish through time
#' @export
#' @family functions for calculating indicators
#' @concept summary_function
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- set_multispecies_model(NS_species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=10)
#' getProportionOfLargeFish(sim)
#' getProportionOfLargeFish(sim, species=c("Herring","Sprat","N.pout"))
#' getProportionOfLargeFish(sim, min_w = 10, max_w = 5000)
#' getProportionOfLargeFish(sim, min_w = 10, max_w = 5000, threshold_w = 500)
#' getProportionOfLargeFish(sim, min_w = 10, max_w = 5000,
#'     threshold_w = 500, biomass_proportion=FALSE)
#' }
getProportionOfLargeFish <- function(sim, 
                                     species = 1:nrow(sim@params@species_params), 
                                     threshold_w = 100, threshold_l = NULL, 
                                     biomass_proportion=TRUE, ...) {
    check_species(sim,species)
    # This args stuff is pretty ugly - couldn't work out another way of using ...
    args <- list(...)
    args[["params"]] <- sim@params
    total_size_range <- do.call("get_size_range_array", args = args)
    args[["max_w"]] <- threshold_w
    args[["max_l"]] <- threshold_l
    large_size_range <- do.call("get_size_range_array", args = args)
    w <- sim@params@w
    if (!biomass_proportion) # based on abundance numbers
        w[] <- 1
    total_measure <- apply(sweep(sweep(sim@n[,species,,drop=FALSE],c(2,3),total_size_range[species,,drop=FALSE],"*"),3,w * sim@params@dw, "*"),1,sum)
    upto_threshold_measure <- apply(sweep(sweep(sim@n[,species,,drop=FALSE],c(2,3),large_size_range[species,,drop=FALSE],"*"),3,w * sim@params@dw, "*"),1,sum)
    #lfi = data.frame(time = as.numeric(dimnames(sim@n)$time), proportion = 1-(upto_threshold_measure / total_measure))
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
#' @param species Numeric or character vector of species to include in the
#'   calculation.
#' @inheritDotParams get_size_range_array -params
#'
#' @return A vector containing the mean weight of the community through time
#' @export
#' @family functions for calculating indicators
#' @concept summary_function
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- set_multispecies_model(NS_species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=10)
#' getMeanWeight(sim)
#' getMeanWeight(sim, species=c("Herring","Sprat","N.pout"))
#' getMeanWeight(sim, min_w = 10, max_w = 5000)
#' }
getMeanWeight <- function(sim, species = 1:nrow(sim@params@species_params), ...){
    check_species(sim, species)
    n_species <- getN(sim, ...)
    biomass_species <- getBiomass(sim, ...)
    n_total <- apply(n_species[, species, drop = FALSE], 1, sum)
    biomass_total <- apply(biomass_species[, species, drop = FALSE], 1, sum)
    return(biomass_total / n_total)
}


#' Calculate the mean maximum weight of the community
#'
#' Calculates the mean maximum weight of the community through time. This can be
#' calculated by numbers or biomass. The calculation is the sum of the w_inf *
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
#' @return A matrix or vector containing the mean maximum weight of the
#'   community through time
#' @export
#' @family functions for calculating indicators
#' @concept summary_function
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- set_multispecies_model(NS_species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=10)
#' getMeanMaxWeight(sim)
#' getMeanMaxWeight(sim, species=c("Herring","Sprat","N.pout"))
#' getMeanMaxWeight(sim, min_w = 10, max_w = 5000)
#' }
getMeanMaxWeight <- function(sim, species = 1:nrow(sim@params@species_params), 
                             measure = "both", ...) {
    if (!(measure %in% c("both","numbers","biomass"))) {
        stop("measure must be one of 'both', 'numbers' or 'biomass'")
    }
    check_species(sim, species)
    n_species <- getN(sim, ...)
    biomass_species <- getBiomass(sim, ...)
    n_winf <- apply(sweep(n_species, 2, sim@params@species_params$w_inf,"*")[,species,drop=FALSE], 1, sum)
    biomass_winf <- apply(sweep(biomass_species, 2, sim@params@species_params$w_inf,"*")[,species,drop=FALSE], 1, sum)
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
#' Calculates the slope of the community abundance through time by performing a linear regression on the logged total numerical abundance at weight and logged weights (natural logs, not log to base 10, are used).
#' You can specify minimum and maximum weight or length range for the species. Lengths take precedence over weights (i.e. if both min_l and min_w are supplied, only min_l will be used).
#' You can also specify the species to be used in the calculation.
#'
#' @inheritParams getMeanWeight
#' @param biomass Boolean. If TRUE (default), the abundance is based on biomass,
#'   if FALSE the abundance is based on numbers.
#' @inheritDotParams get_size_range_array -params
#'
#' @return A data frame with slope, intercept and R2 values.
#' @export
#' @family functions for calculating indicators
#' @concept summary_function
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- set_multispecies_model(NS_species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=40, dt = 1, t_save = 1)
#' # Slope based on biomass, using all species and sizes
#' slope_biomass <- getCommunitySlope(sim)
#' # Slope based on numbers, using all species and sizes
#' slope_numbers <- getCommunitySlope(sim, biomass=FALSE)
#' # Slope based on biomass, using all species and sizes between 10g and 1000g
#' slope_biomass <- getCommunitySlope(sim, min_w = 10, max_w = 1000)
#' # Slope based on biomass, using only demersal species and sizes between 10g and 1000g
#' dem_species <- c("Dab","Whiting","Sole","Gurnard","Plaice","Haddock", "Cod","Saithe")
#' slope_biomass <- getCommunitySlope(sim, species = dem_species, min_w = 10, max_w = 1000)
#' }
getCommunitySlope <- function(sim, species = 1:nrow(sim@params@species_params),
                              biomass = TRUE, ...) {
    check_species(sim, species)
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
        out_df <- data.frame(
            slope = summary_fit$coefficients[2, 1],
            intercept = summary_fit$coefficients[1, 1],
            r2 = summary_fit$r.squared
        )
    }, w = sim@params@w)
    dimnames(slope)[[1]] <- slope[, 1]
    slope <- slope[, -1]
    return(slope)
}


# internal
check_species <- function(object, species){
    if (!(is(species,"character") | is(species,"numeric")))
        stop("species argument must be either a numeric or character vector")
    if (is(species,"character"))
        check <- all(species %in% dimnames(object@n)$sp)  
    if (is(species,"numeric"))
        check <- all(species %in% 1:dim(object@n)[2])
    if (!check)
        stop("species argument not in the model species. species must be a character vector of names in the model, or a numeric vector referencing the species")
    return(check)
}

