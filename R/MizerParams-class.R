# Class specification and constructors for mizer base parameters class
# Class has members to store parameters of size based model

# Copyright 2012 Finlay Scott and Julia Blanchard.
# Copyright 2018 Gustav Delius and Richard Southwell.
# Development has received funding from the European Commission's Horizon 2020 
# Research and Innovation Programme under Grant Agreement No. 634495 
# for the project MINOUW (http://minouw-project.eu/).
# Distributed under the GPL 3 or later 
# Maintainer: Gustav Delius, University of York, <gustav.delius@york.ac.uk>

#Naming conventions:
#S4 classes and constructors: AClass
#functions aFunction


# Validity function ---------------------------------------------------------
# Not documented as removed later on
validMizerParams <- function(object) {
    
    errors <- character()
    # grab some dims
    no_w <- length(object@w)
    no_w_full <- length(object@w_full)
    w_idx <- (no_w_full - no_w + 1):no_w_full
    
    # Check dw and dw_full are correct length
    if (length(object@dw) != no_w) {
        msg <- paste("dw is length ", length(object@dw),
                     " and w is length ", no_w,
                     ". These should be the same length", sep = "")
        errors <- c(errors, msg)
    }
    
    if (length(object@dw_full) != no_w_full) {
        msg <- paste("dw_full is length ", length(object@dw_full),
                     " and w_full is length ", no_w_full,
                     ". These should be the same length", sep = "")
        errors <- c(errors, msg)
    }
    
    # Check that the last entries of w_full and dw_full agree with w and dw
    if (!identical(object@w, object@w_full[w_idx])) {
        msg <- "The later entries of w_full should be equal to those of w."
        errors <- c(errors, msg)
    }
    if (!identical(object@dw, object@dw_full[w_idx])) {
        msg <- "The later entries of dw_full should be equal to those of dw."
        errors <- c(errors, msg)
    }

    # Check the array dimensions are good
    # 2D arrays
    if (!all(c(length(dim(object@psi)),
               length(dim(object@intake_max)),
               length(dim(object@search_vol)),
               length(dim(object@metab)),
               length(dim(object@mu_b)),
               length(dim(object@interaction)),
               length(dim(object@catchability))) == 2)) {
        msg <- "psi, intake_max, search_vol, metab, mu_b, interaction and catchability must all be two dimensions"
        errors <- c(errors, msg)
    }
    # 3D arrays
    if (length(dim(object@selectivity)) != 3) {
        msg <- "selectivity must be three dimensions"
        errors <- c(errors, msg)
    }
    # Check number of species is equal across relevant slots
    if (!all(c(
        dim(object@psi)[1],
        dim(object@intake_max)[1],
        dim(object@search_vol)[1],
        dim(object@metab)[1],
        dim(object@mu_b)[1],
        dim(object@selectivity)[2],
        dim(object@catchability)[2],
        dim(object@interaction)[1],
        dim(object@interaction)[2]) == 
        dim(object@species_params)[1])) {
        msg <- "The number of species in the model must be consistent across the species_params, psi, intake_max, search_vol, mu_b, interaction (dim 1), selectivity, catchability and interaction (dim 2) slots"
        errors <- c(errors, msg)
    }
    # Check number of size groups
    if (!all(c(
        dim(object@psi)[2],
        dim(object@intake_max)[2],
        dim(object@search_vol)[2],
        dim(object@metab)[2],
        dim(object@selectivity)[3]) ==
        no_w)) {
        msg <- "The number of size bins in the model must be consistent across the w, psi, intake_max, search_vol, and selectivity (dim 3) slots"
        errors <- c(errors, msg)
    }
    # Check numbe of gears
    if (!isTRUE(all.equal(dim(object@selectivity)[1], dim(object@catchability)[1]))) {
        msg <- "The number of fishing gears must be consistent across the catchability and selectivity (dim 1) slots"
        errors <- c(errors, msg)
    }
    # Check names of dimnames of arrays
    # sp dimension
    if (!all(c(
        names(dimnames(object@psi))[1],
        names(dimnames(object@intake_max))[1],
        names(dimnames(object@search_vol))[1],
        names(dimnames(object@metab))[1],
        names(dimnames(object@mu_b))[1],
        names(dimnames(object@selectivity))[2],
        names(dimnames(object@catchability))[2]) == "sp")) {
        msg <- "Name of first dimension of psi, intake_max, search_vol, metab, mu_b, and the second dimension of selectivity and catchability must be 'sp'"
        errors <- c(errors, msg)
    }
    #interaction dimension names
    if (names(dimnames(object@interaction))[1] != "predator") {
        msg <- "The first dimension of interaction must be called 'predator'"
        errors <- c(errors, msg)
    }
    if (names(dimnames(object@interaction))[2] != "prey") {
        msg <- "The first dimension of interaction must be called 'prey'"
        errors <- c(errors, msg)
    }
    # w dimension
    if (!all(c(
        names(dimnames(object@psi))[2],
        names(dimnames(object@intake_max))[2],
        names(dimnames(object@search_vol))[2],
        names(dimnames(object@metab))[2],
        names(dimnames(object@selectivity))[3]) == "w")) {
        msg <- "Name of second dimension of psi, intake_max, search_vol, metab and third dimension of selectivity must be 'w'"
        errors <- c(errors, msg)
    }
    if (!all(c(
        names(dimnames(object@selectivity))[1],
        names(dimnames(object@catchability))[1]) == "gear")) {
        msg <- "Name of first dimension of selectivity and catchability must be 'gear'"
        errors <- c(errors, msg)
    }
    
    # Check dimnames of species are identical
    # Bit tricky this one as I don't know of a way to compare lots of vectors 
    # at the same time. Just use == and the recycling rule
    if (!all(c(
        dimnames(object@psi)[[1]],
        dimnames(object@intake_max)[[1]],
        dimnames(object@search_vol)[[1]],
        dimnames(object@metab)[[1]],
        dimnames(object@mu_b)[[1]],
        dimnames(object@selectivity)[[2]],
        dimnames(object@catchability)[[2]],
        dimnames(object@interaction)[[1]],
        dimnames(object@interaction)[[2]]) ==
        object@species_params$species)) {
        msg <- "The species names of species_params, psi, intake_max, search_vol, metab, mu_b, selectivity, catchability and interaction must all be the same"
        errors <- c(errors, msg)
    }
    # Check dimnames of w
    if (!all(c(
        dimnames(object@psi)[[2]],
        dimnames(object@intake_max)[[2]],
        dimnames(object@search_vol)[[2]],
        dimnames(object@metab)[[2]]) == 
        dimnames(object@selectivity)[[3]])) {
        msg <- "The size names of psi, intake_max, search_vol, metab and selectivity must all be the same"
        errors <- c(errors, msg)
    }
    # Check dimnames of gear
    if (!isTRUE(all.equal(
        dimnames(object@catchability)[[1]],
        dimnames(object@selectivity)[[1]]))) {
        msg <- "The gear names of selectivity and catchability must all be the same"
        errors <- c(errors, msg)
    }
    # Check the vector slots
    if (length(object@rr_pp) != length(object@w_full)) {
        msg <- "rr_pp must be the same length as w_full"
        errors <- c(errors, msg)
    }
    if (length(object@cc_pp) != length(object@w_full)) {
        msg <- "cc_pp must be the same length as w_full"
        errors <- c(errors, msg)
    }
    
    # SRR
    # Must have two arguments: rdi amd species_params
    if (!isTRUE(all.equal(names(formals(object@srr)), c("rdi", "species_params")))) {
        msg <- "Arguments of srr function must be 'rdi' and 'species_params'"
        errors <- c(errors, msg)
    }
    
    # # species_params data.frame must have columns: 
    # # species, z0, alpha, eRepro
    # species_params_cols <- c("species","z0","alpha","erepro")
    # if (!all(species_params_cols %in% names(object@species_params))) {
    #     msg <- "species_params data.frame must have 'species', 'z0', 'alpha' and 'erepro' columms"
    #     errors <- c(errors,msg)
    # }
    # must also have SRR params but sorted out yet
    
    # species_params
    # Column check done in constructor
    # If everything is OK
    if (length(errors) == 0) TRUE else errors
}


#### Class definition ####
#' A class to hold the parameters for a size based model. 
#' 
#' MizerParams objects can be created using a range of constructor functions.
#' 
#' Dynamic simulations are performed using the \code{\link{project}} function on
#' objects of this class.
#' 
#' @slot w The size grid for the fish part of the spectrum. An increasing
#'   vector of weights (in grams) running from the smallest egg size to the
#'   largest asymptotic size.
#' @slot dw The spacing in the size grid. So dw[i] = w[i+1] - w[i]. A vector 
#'   the same length as the w_full slot. The last entry is not determined by
#'   the w slot but represents the size of the last size bin.
#' @slot w_full The size grid for the full size range including the plankton
#'   spectrum. An increasing vector of weights (in grams) running from the
#'   smallest plankton size to the largest asymptotic size of fish. The
#'   last entries of the vector have to be equal to the content of the w slot.
#' @slot dw_full The spacing in the full size grid. 
#'   So dw_full[i] = w_full[i+1] - w_full[i]. The last entries have to be
#'   equal to the content of the dw slot.
#' @slot w_min_idx A vector holding the index of the weight of the egg size
#'   of each species
#' @slot maturity An array (species x size) that holds the proportion of
#'   individuals of each species at size that are mature. This enters in the
#'   calculation of the spawning stock biomass with \code{\link{getSSB}}. Set 
#'   with \code{\link{setReproduction}}.
#' @slot psi An array (species x size) that holds the allocation to reproduction
#'   for each species at size, \eqn{\psi_i(w)}. Changed with 
#'   \code{\link{setReproduction}}.
#' @slot intake_max An array (species x size) that holds the maximum intake for
#'   each species at size. Changed with \code{\link{setIntakeMax}}.
#' @slot search_vol An array (species x size) that holds the search volume for
#'   each species at size. Changed with \code{\link{setSearchVolume}}.
#' @slot rho A 3-dim array (species x resource x size) holding the encounter
#'   rates for unstructured resources. Changed with 
#'   \code{\link{setResourceEncounter}}.
#' @slot metab An array (species x size) that holds the metabolism
#'   for each species at size. Changed with \code{\link{setMetab}}.
#' @slot mu_b An array (species x size) that holds the background death 
#'   \eqn{\mu_{b.i}(w)}. Changed with \code{\link{setBMort}}.
#' @slot pred_kernel An array (species x predator size x prey size) that holds
#'   the predation coefficient of each predator at size on each prey size. If
#'   this is NA then the following two slots will be used. Changed with 
#'   \code{\link{setPredKernel}}.
#' @slot ft_pred_kernel_e An array (species x log of predator/prey size ratio)
#'   that holds the Fourier transform of the feeding kernel in a form
#'   appropriate for evaluating the encounter rate integral. If this is NA
#'   then the \code{pred_kernel} will be used to calculate the available 
#'   energy integral. Changed with \code{\link{setPredKernel}}.
#' @slot ft_pred_kernel_p An array (species x log of predator/prey size ratio)
#'   that holds the Fourier transform of the feeding kernel in a form
#'   appropriate for evaluating the predation mortality integral. If this is NA
#'   then the \code{pred_kernel} will be used to calculate the integral.
#'   Changed with \code{\link{setPredKernel}}.
#' @slot rr_pp A vector the same length as the w_full slot. The size specific
#'   growth rate of the plankton spectrum. Changed with \code{\link{setPlankton}}.
#' @slot cc_pp A vector the same length as the w_full slot. The size specific
#'   carrying capacity of the plankton spectrum. Changed with 
#'   \code{\link{setPlankton}}.
#' @slot plankton_dynamics A function for projecting the plankton abundance
#'   density by one timestep. The default is 
#'   \code{\link{plankton_semichemostat}}. 
#'   Changed with \code{\link{setPlankton}}.
#' @slot resource_dynamics A named list of functions for projecting the
#'   biomasses in the unstructured resource components by one timestep. The
#'   names of the list entries are the resource names. Changed with 
#'   \code{\link{setResourceDynamics}}.
#' @slot resource_params A list containing the parameters needed by the
#'   \code{resource_dynamics} functions.  Changed with 
#'   \code{\link{setResourceDynamics}}.
#' @slot sc The community abundance of the scaling community
#' @slot species_params A data.frame to hold the species specific parameters.
#'   See \code{\link{set_multispecies_model}} for details.
#' @slot interaction The species specific interaction matrix, \eqn{\theta_{ij}}.
#'   Changed with \code{\link{setInteraction}}.
#' @slot srr Function to calculate the realised (density dependent) recruitment.
#'   Has two arguments which are rdi and species_params.
#' @slot selectivity An array (gear x species x w) that holds the selectivity of
#'   each gear for species and size, \eqn{S_{g,i,w}}. Changed with 
#'   \code{\link{setFishing}}.
#' @slot catchability An array (gear x species) that holds the catchability of
#'   each species by each gear, \eqn{Q_{g,i}}. Changed with 
#'   \code{\link{setFishing}}.
#' @slot initial_n An array (species x size) that holds abundance of each species
#'   at each weight at our candidate steady state solution.
#' @slot initial_n_pp A vector the same length as the w_full slot that describes
#'   the plankton abundance at each weight.
#' @slot initial_B A vector containing the biomasses of the unstructured
#'   resource components.
#'   Has length zero if there are no unstructured resources.
#' @slot n Exponent of maximum intake rate. Changed with \code{\link{setIntakeMax}}.
#' @slot p Exponent of metabolic cost. Changed with \code{\link{setMetab}}.
#' @slot lambda Exponent of plankton spectrum. Changed with \code{\link{setPlankton}}.
#' @slot kappa Magnitude of plankton spectrum. Changed with \code{\link{setPlankton}}.
#' @slot q Exponent for volumetric search rate. 
#'   Changed with \code{\link{setSearchVolume}}.
#' @slot f0 Initial feeding level.
#' @slot A Abundance multipliers.
#' @slot linecolour A named vector of colour values, named by species. Used 
#'   to give consistent colours to species in plots.
#' @slot linetype A named vector of linetypes, named by species. Used 
#'   to give consistent line types to species in plots.

#' @note The \linkS4class{MizerParams} class is fairly complex with a large number of
#'   slots, many of which are multidimensional arrays. The dimensions of these
#'   arrays is strictly enforced so that \code{MizerParams} objects are
#'   consistent in terms of number of species and number of size classes.
#'   
#'   Although it is possible to build a \code{MizerParams} object by hand it is
#'   not recommended and several constructors are available.
#'   
#'   The \code{MizerParams} class does not hold any dynamic information, e.g.
#'   abundances or harvest effort through time. These are held in
#'   \linkS4class{MizerSim} objects.
#' @seealso \code{\link{project}} \code{\link{MizerSim}}
#'   \code{\link{emptyParams}} \code{\link{set_multispecies_model}}
#'   \code{\link{set_community_model}}
#'   \code{\link{set_trait_model}} \code{\link{set_scaling_model}}
#' @export
setClass(
    "MizerParams",
    slots = c(
        w = "numeric",
        dw = "numeric",
        w_full = "numeric",
        dw_full = "numeric",
        w_min_idx = "numeric",
        maturity = "array",
        psi = "array",
        initial_n = "array",
        intake_max = "array",
        search_vol = "array",
        rho = "array",
        metab = "array",
        pred_kernel = "array",
        ft_pred_kernel_e = "array",
        ft_pred_kernel_p = "array",
        mu_b = "array",
        rr_pp = "numeric",
        cc_pp = "numeric",
        resource_dynamics = "list",
        plankton_dynamics = "function",
        resource_params = "list",
        sc = "numeric",
        initial_n_pp = "numeric",
        initial_B = "numeric",
        species_params = "data.frame",
        interaction = "array",
        srr  = "function",
        selectivity = "array",
        catchability = "array",
        n = "numeric",
        p = "numeric",
        lambda = "numeric",
        q = "numeric",
        f0 = "numeric",
        kappa = "numeric",
        A = "numeric",
        linecolour = "character",
        linetype = "character"
    ),
)

setValidity("MizerParams", validMizerParams)
remove(validMizerParams)


#' Create empty MizerParams object of the right size
#' 
#' Sets up a valid \linkS4class{MizerParams} object with all the slots
#' initialised and given dimension names, but with some slots left empty. This
#' function is to be used by other functions to set up full parameter objects.
#' 
#' See \code{\link{set_multispecies_model}} for a function that fills the
#' slots left empty by this function.
#' 
# Some code is commented out that would allow the user to 
# specify a grid with a non-constant log spacing. But we comment this out
# for now because of the fft.
# #' When the `w_full` argument is not given, then 
#' A size grid is created so that
#' the log-sizes are equally spaced. The spacing is chosen so that there will be
#' `no_w` fish size bins, with the smallest starting at `min_w` and the largest
#' starting at `max_w`. For `w_full` additional size bins are added below
#' `min_w`, with the same log size. The number of extra bins is such that
#' `min_w_pp` comes to lie within the smallest bin. 
#' 
#' The \code{species_params} slot of the returned MizerParams object may differ
#' slightly from the data frame supplied as argument to this function in the
#' following ways:
#' \itemize{
#'   \item Default values are set for \code{w_min, w_inf, alpha, gear, interaction_p}.
#'   \item The egg sizes in \code{w_min} are rounded down to lie on a grid point.
#' }
#' Note that the other characteristic sizes of the species, like \code{w_mat} and
#' \code{w_inf}, are not modified to lie on grid points.
#' 
#' @param species_params A data frame of species-specific parameter values.
#' @param no_w The number of size bins in the consumer spectrum.
#' @param min_w Sets the size of the eggs of all species for which this is not
#'   given in the \code{w_min} column of the \code{species_params} dataframe.
# #' @param w_full Increasing vector of weights giving the boundaries of size
# #'   classes. Must include the value min_w. Has one more entry than the number
# #'   of size bins. The last entry is the upper end of the largest size class. It
# #'   be used to calculate the sizes of the size bins but will not be stored in
# #'   the w_full slot of the returned MizerParams object. If this argument is not
# #'   provided then size classes are set by the other arguments as described in
# #'   the Details.
#' @param max_w The largest size of the consumer spectrum. By default this is
#'   set to the largest w_inf specified in the species_params data frame.
#' @param min_w_pp The smallest size of the plankton spectrum.
# #'   Ignored if w_full is specified.
#' @param no_w_pp  No longer used
#' 
#' @return An empty but valid MizerParams object
#' 
#' @export
emptyParams <- function(species_params,
                        no_w = 100,
                        min_w = 0.001,
                        # w_full = NA,
                        max_w = NA,
                        min_w_pp = NA,
                        no_w_pp = NA) {
    if (!is.na(no_w_pp)) {
        warning("New mizer code does not support the parameter no_w_pp")
    }
    assert_that(is.data.frame(species_params))
    assert_that(no_w > 10)
    
    if (!("species" %in% colnames(species_params))) {
        stop("The species params dataframe needs a column 'species' with the species names")
    }
    species_names <- species_params$species
    row.names(species_params) <- species_names
    no_sp <- nrow(species_params)
    
    ## Set defaults ----
    species_params <- set_species_param_default(species_params, "w_min", min_w)
    min_w <- min(species_params$w_min)
    
    if (!("w_inf" %in% colnames(species_params))) {
        species_params$w_inf <- rep(NA, no_sp)
    }
    missing <- is.na(species_params$w_inf)
    if (any(missing)) {
        stop("You need to specify maximum sizes for all species.")
    }
    if (is.na(max_w)) {
        max_w <- max(species_params$w_inf)
    } else {
        if (max(species_params$w_inf) > max_w * (1 + 1e-9)) { # The fudge factor
            # is there to avoid false alerts due to rounding errors.
            too_large <- species_params$species[max_w < species_params$w_inf]
            stop(paste0("Some of your species have an maximum size larger than max_w: ",
                        toString(too_large)))
        }
    }
    
    # If no gear_name column in species_params, then named after species
    if (!("gear" %in% colnames(species_params))) {
        species_params$gear <- species_params$species
    }
    gear_names <- unique(species_params$gear)
    
    # If no alpha (conversion efficiency), then set to 0.6
    species_params <- set_species_param_default(species_params,
                                                "alpha", 0.6)
    
    # Set up grids ----
    # The following code anticipates that in future we might allow the user to 
    # specify a grid with a non-constant log spacing. But we comment this out
    # for now because of the fft.
    # if (missing(w_full)) {
        # set up logarithmic grids
        dx <- log10(max_w / min_w) / (no_w - 1)
        # Community grid
        w <- 10^(seq(from = log10(min_w), by = dx, length.out = no_w))
        # dw[i] = w[i+1] - w[i]. Following formula works also for last entry dw[no_w]
        dw <- (10^dx - 1) * w
        # To avoid issues due to numerical imprecission
        min_w <- w[1]
        
        # If not provided, set min_w_pp so that all fish have their full feeding
        # kernel inside plankton spectrum
        ppmr <- 10^(seq(from = 0, by = dx, length.out = 3 * no_w))
        phis <- get_phi(species_params, ppmr)
        max_ppmr <- apply(phis, 1, function(x) ppmr[max(which(x != 0)) + 1])
        min_w_feeding <- species_params$w_min / max_ppmr
        species_params <- set_species_param_default(species_params,
                                                    "interaction_p", 1)
        if (any(species_params$interaction_p > 0)) {
            min_w_feeding <- min_w_feeding[species_params$interaction_p > 0]
        } else {
            min_w_feeding <- min_w
        }
        if (is.na(min_w_pp)) {
            min_w_pp <- min(min_w_feeding)
        } else {
            assert_that(min_w_pp <= min_w)
            hungry_sp <- species_params$species[min_w_feeding < min_w_pp]
            if (length(hungry_sp) > 0) {
                message(paste(
                    "Note: The following species have feeding kernels that extend",
                    "below the smallest plankton size specified by min_w_pp:",
                    toString(hungry_sp)))
            }
        }
        # For fft methods we need a constant log bin size throughout. 
        # Therefore we use as many steps as are necessary so that the first size
        # class includes min_w_pp.
        x_pp <- rev(seq(from = log10(min_w),
                        to = log10(min_w_pp),
                        by = -dx)) - dx
        w_full <- c(10^x_pp, w)
        # If min_w_pp happened to lie exactly on a grid point, we now added
        # one grid point too much which we need to remove again
        if (w_full[2] == min_w_pp) {
            w_full <- w_full[2:length(w_full)]
        }
        no_w_full <- length(w_full)
        dw_full <- (10^dx - 1) * w_full	
    # } else {
    #     # use supplied w_full
    #     no_w_full <- length(w_full) - 1
    #     dw_full <- diff(w_full)
    #     w_full <- w_full[seq_along(dw_full)]
    #     # Check that sizes are increasing
    #     if (any(dw_full <= 0)) {
    #         stop("w_full must be increasing.")
    #     }
    #     w_min_idx <- match(min_w, w_full)
    #     if (is.na(w_min_idx)) {
    #         stop("w_min must be contained in w_full.")
    #     }
    #     w <- w_full[w_min_idx:no_w_full]
    #     dw <- dw_full[w_min_idx:no_w_full]
    #     no_w <- length(w)
    #     min_w_pp <- w_full[1]
    # }
    
    # Basic arrays for templates ----
    mat1 <- array(NA, dim = c(no_sp, no_w), 
                  dimnames = list(sp = species_names, w = signif(w,3)))
    ft_pred_kernel <- array(NA, dim = c(no_sp, no_w_full),
                            dimnames = list(sp = species_names, k = 1:no_w_full))
    
    selectivity <- array(0, dim = c(length(gear_names), no_sp, no_w), 
                         dimnames = list(gear = gear_names, sp = species_names, 
                                         w = signif(w, 3)))
    catchability <- array(0, dim = c(length(gear_names), no_sp), 
                          dimnames = list(gear = gear_names, sp = species_names))
    interaction <- array(1, dim = c(no_sp, no_sp),
                         dimnames = list(predator = species_names,
                                         prey = species_names))
    
    vec1 <- as.numeric(rep(NA, no_w_full))
    names(vec1) <- signif(w_full, 3)
    
    # Round down w_min to lie on grid points and store the indices of these
    # grid points in w_min_idx
    w_min_idx <- as.vector(suppressWarnings(
        tapply(species_params$w_min, 1:no_sp,
               function(w_min, wx) max(which(wx <= w_min)), wx = w)))
    # Due to rounding errors this might happen:
    w_min_idx[w_min_idx == -Inf] <- 1
    names(w_min_idx) = species_names
    species_params$w_min <- w[w_min_idx]
    
    # Colour and linetype scales ----
    # for use in plots
    # Colour-blind-friendly palettes
    # From http://dr-k-lo.blogspot.co.uk/2013/07/a-color-blind-friendly-palette-for-r.html
    # cbbPalette <- c("#000000", "#009E73", "#e79f00", "#9ad0f3", "#0072B2", "#D55E00", 
    #                 "#CC79A7", "#F0E442")
    # From http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
    cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", 
                    "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    
    if ("linecolour" %in% names(species_params)) {
        linecolour <- species_params$linecolour
    } else {
        linecolour <- rep(cbbPalette, length.out = no_sp)
    }
    names(linecolour) <- as.character(species_names)
    linecolour <- c(linecolour, "Total" = "black", "Plankton" = "green",
                    "Background" = "grey")
    
    if ("linetype" %in% names(species_params)) {
        linetype <- species_params$linetype
    } else {
        linetype <- rep(c("solid", "dashed", "dotted", "dotdash", "longdash", 
                          "twodash"), length.out = no_sp)
    }
    names(linetype) <- as.character(species_names)
    linetype <- c(linetype, "Total" = "solid", "Plankton" = "solid",
                  "Background" = "solid")
    
    # Make object ----
    # Should Z0, rrPP and ccPP have names (species names etc)?
    params <- new(
        "MizerParams",
        w = w,
        dw = dw,
        w_full = w_full,
        dw_full = dw_full,
        w_min_idx = w_min_idx,
        maturity = mat1,
        psi = mat1,
        initial_n = mat1,
        intake_max = mat1,
        search_vol = mat1,
        rho = array(),
        metab = mat1,
        mu_b = mat1,
        ft_pred_kernel_e = ft_pred_kernel,
        ft_pred_kernel_p = ft_pred_kernel,
        pred_kernel = array(),
        selectivity = selectivity,
        catchability = catchability,
        rr_pp = vec1,
        cc_pp = vec1,
        sc = w,
        initial_n_pp = vec1,
        species_params = species_params,
        interaction = interaction,
        srr = srrBevertonHolt,
        resource_dynamics = list(),
        plankton_dynamics = plankton_semichemostat,
        resource_params = list(),
        initial_B = vector(mode = "numeric"),
        A = as.numeric(rep(NA, no_sp)),
        linecolour = linecolour,
        linetype = linetype
    )
    
    return(params)
}

#' Construct MizerParams object for general multispecies model
#'
#' Sets up a multi-species size spectrum model by filling all slots in the
#' \linkS4class{MizerParams} object based on user-provided or default
#' parameters. It does this by creating an empty MizerParams object with
#' \code{\link{emptyParams}} and then filling the slots by passing its arguments
#' to \code{\link{setParams}}. There is a long list of arguments, but almost
#' all of them have sensible default values. All arguments are described in more
#' details in the sections below the list.
#' 
#' @inheritParams emptyParams
#' @param min_w_pp The smallest size of the plankton spectrum. By default this
#'   is set to the smallest value at which any of the consumers can feed.
#' @param no_w_pp Obsolete argument that is no longer used because the number
#'    of plankton size bins is determined because all size bins have to
#'    be logarithmically equally spaced.
#' @inheritParams setParams
#'
#' @return An object of type \linkS4class{MizerParams}
#' 
#' @section Species parameters:
#' The only essential argument is a data frame that contains the species
#' parameters. The data frame is arranged species by parameter, so each column
#' of the parameter data frame is a parameter and each row has the values of the
#' parameters for one of the species in the model.
#'
#' There are two essential columns that must be included in the species
#' parameter data.frame and that do not have default values: the 
#' \code{species} column that should hold strings with the names of the
#' species and the \code{w_inf} column with the asymptotic sizes of the species. 
#' 
#' The species_params dataframe also needs to contain the parameters needed
#' by any predation kernel function or size selectivity function. This will
#' be mentioned in the appropriate sections below.
#' 
#' For all other species parameters, mizer will calculate default values if they
#' are not included in the species parameter data frame. They will be
#' automatically added when the \code{MizerParams} object is created. For these
#' parameters you can also specify values for only some species and leave the
#' other entries as NA and the missing values will be set to the defaults.
#' 
#' All the parameters will be mentioned in the following sections.
#' 
#' @inheritSection setParams Units in mizer
#' @inheritSection setInteraction Setting interactions
#' @inheritSection setPredKernel Setting predation kernel
#' @inheritSection setSearchVolume Setting search volume
#' @inheritSection setIntakeMax Setting maximum intake rate
#' @inheritSection setMetab Setting metabolic rate
#' @inheritSection setBMort Setting background mortality rate
#' @inheritSection setReproduction Setting reproduction
#' @inheritSection setFishing Setting fishing
#' @inheritSection setPlankton Setting plankton dynamics
#' @inheritSection setResourceDynamics Setting resource dynamics
#' @inheritSection setResourceEncounter Setting resource encounter rate
#'   
#' @seealso \code{\link{project}}, \linkS4class{MizerSim},
#'   \code{\link{set_community_model}}, \code{\link{set_trait_model}}
#' @export
#' @family functions for setting up models
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- set_multispecies_model(NS_species_params_gears, inter)
#' }
set_multispecies_model <- function(species_params,
                                   interaction = matrix(1,
                                                        nrow = nrow(species_params),
                                                        ncol = nrow(species_params)),
                                   no_w = 100,
                                   min_w = 0.001,
                                   max_w = NA,
                                   min_w_pp = NA,
                                   no_w_pp = NA,
                                   n = 2 / 3,
                                   p = 0.7,
                                   q = 0.8,
                                   kappa = 1e11,
                                   lambda = (2 + q - n),
                                   f0 = 0.6,
                                   # setPredKernel()
                                   pred_kernel = NULL,
                                   # setSearchVolume()
                                   search_vol = NULL,
                                   # setIntakeMax()
                                   intake_max = NULL,
                                   # setMetab()
                                   metab = NULL,
                                   # setBMort
                                   mu_b = NULL,
                                   z0pre = 0.6,
                                   z0exp = n - 1,
                                   # setReproduction
                                   maturity = NULL,
                                   repro_prop = NULL,
                                   # setPlankton
                                   r_pp = 10,
                                   w_pp_cutoff = 10,
                                   plankton_dynamics = plankton_semichemostat,
                                   # setResourceDynamics
                                   resource_dynamics = list(),
                                   resource_params = list(),
                                   # setResourceEncounter
                                   rho = NULL,
                                   srr = srrBevertonHolt) {
    
    ## Create MizerParams object ----
    params <- emptyParams(species_params,
                          no_w = no_w, 
                          min_w = min_w,  
                          max_w = max_w, 
                          min_w_pp = min_w_pp, 
                          no_w_pp = NA)
    
    ## Fill the slots ----
    params@n <- n
    params@lambda <- lambda
    params@f0 <- f0
    params@kappa <- kappa
    
    params <- setParams(params,
                        # setInteraction
                        interaction = interaction,
                        # setPredKernel()
                        pred_kernel = pred_kernel,
                        # setSearchVolume()
                        search_vol = search_vol,
                        q = q,
                        f0 = f0,
                        # setIntakeMax()
                        intake_max = intake_max,
                        n = n,
                        # setMetab()
                        metab = metab,
                        p = p,
                        # setBMort
                        mu_b = mu_b,
                        z0pre = z0pre,
                        z0exp = z0exp,
                        # setReproduction
                        maturity = maturity,
                        repro_prop = repro_prop,
                        srr = srr,
                        # setPlankton
                        kappa = kappa,
                        lambda = lambda,
                        r_pp = r_pp,
                        w_pp_cutoff = w_pp_cutoff,
                        plankton_dynamics = plankton_dynamics,
                        # setResourceDynamics
                        resource_dynamics = resource_dynamics,
                        resource_params = resource_params,
                        # setResourceEncounter
                        rho = rho)
    
    params@initial_n <- get_initial_n(params)
    params@initial_n_pp <- params@cc_pp
    params@A <- rep(1, nrow(species_params))
    
    return(params)
}

#' Alias for set_multispecies_model
#' 
#' An alias provided for backward compatibility with mizer version <= 1.0
#' @inherit set_multispecies_model
#' @export
MizerParams <- set_multispecies_model

#' Set or change any model parameters
#' 
#' This is a convenient wrapper function calling each of the following
#' functions
#' \itemize{
#' \item \code{\link{setPredKernel}}
#' \item \code{\link{setSearchVolume}}
#' \item \code{\link{setInteraction}}
#' \item \code{\link{setIntakeMax}}
#' \item \code{\link{setMetab}}
#' \item \code{\link{setBMort}}
#' \item \code{\link{setReproduction}}
#' \item \code{\link{setFishing}}
#' \item \code{\link{setPlankton}}
#' \item \code{\link{setResourceDynamics}}
#' \item \code{\link{setResourceEncounter}}
#' }
#' See the Details section below for a discussion of how to use this function.
#' 
#' @param params A \linkS4class{MizerParams} object
#' @param f0 Average feeding level. Used to calculated \code{h} and \code{gamma}
#'   if those are not columns in the species data frame. Also requires
#'   \code{k_vb} (the von Bertalanffy K parameter) to be a column in the species
#'   data frame. 
#' @inheritParams setInteraction
#' @inheritParams setPredKernel
#' @inheritParams setSearchVolume
#' @inheritParams setIntakeMax
#' @inheritParams setMetab
#' @inheritParams setBMort
#' @inheritParams setReproduction
#' @inheritParams setFishing
#' @inheritParams setPlankton
#' @inheritParams setResourceDynamics
#' @inheritParams setResourceEncounter
#' 
#' @return A \linkS4class{MizerParams} object
#' 
#' @details 
#' Usually, if you are happy with the way mizer calculates its model functions
#' from the species parameters and only want to change the values of some
#' species parameters, you would make those changes in the `species_params` data
#' frame contained in the `params` object and then call the `setParams()`
#' function to effect the change. Note that just changing the species parameters
#' by themselves is not changing the model until you call `setParams()` or the
#' appropriate one of its sub-functions. Here is an example which assumes that
#' you have have a MizerParams object `params` in which you just want to change
#' one parameter of the third species:
#' ```
#' params@species_params$gamma[3] <- 1000
#' params <- setParams(params)
#' ```
#' Because of the way the R language works, `setParams` does not make the
#' changes to the `params` object that you pass to it but instead returns a new
#' params object. So to affect the change you call the function in the form
#' `params <- setParams(params, ...)`.
#' 
#' If you are not happy with the assumptions that mizer makes by default about
#' the shape of the model functions, for example if you want to change one of
#' the allometric scaling assumptions, you can do this by providing your
#' choice as an array in the appropriate argument to `setParams()`. The
#' sections below discuss all the model functions that you can change this way.
#' 
#' This function will use the species parameters in the `params` object to reset
#' the values of all the model functions that you do not specify explictly when
#' calling this function. If you have changed any of the model functions in the
#' `params` object previously and now want to make changes to a different slot,
#' you will want to call the appropriate change function individually. So in the
#' above example you would have used `params <- setSearchVolume(params)`
#' instead of `params <- setParams(params)`.
#' 
#' @section Units in mizer:
#' Mizer uses grams to measure weight, centimetres to measure lengths, and
#' years to measure time.
#' 
#' Mizer is agnostic about whether abundances are given as 
#' 1. numbers per area, 
#' 2. numbers per volume or
#' 3. total numbers for the entire study area. 
#' 
#' You should make the choice most convenient for your application and then
#' stick with it. If you make choice 1 or 2 you will also have to choose a unit
#' for area or volume. Your choice will then determine the units for some of
#' the parameters. This will be mentioned when the parameters are discussed in
#' the sections below.
#' 
#' You choice will also affect the units of the quantities you may want to
#' calculate with the model. For example, the yield will be in grams/year/m^2 in
#' case 1 if you choose m^2 as your measure of area, in grams/year/m^3 in case 2
#' if you choose m^3 as your unit of volume, or simply grams/year in case 3. The
#' same comment applies for other measures, like total biomass, which will be
#' grams/area in case 1, grams/volume in case 2 or simply grams in case 3. When
#' mizer puts units on axes, for example in \code{plotBiomass}, it will simply
#' put grams, as appropriate for case 3.
#' 
#' You can convert between these choices. For example, if you use case 1, you
#' need to multiply with the area of the ecosystem to get the total quantity. 
#' If you work with case 2, you need to multiply by both area and the thickness 
#' of the productive layer. In that respect, case 2 is a bit cumbersome.
#' 
#' @inheritSection setInteraction Setting interactions
#' @inheritSection setPredKernel Setting predation kernel
#' @inheritSection setSearchVolume Setting search volume
#' @inheritSection setIntakeMax Setting maximum intake rate
#' @inheritSection setMetab Setting metabolic rate
#' @inheritSection setBMort Setting background mortality rate
#' @inheritSection setReproduction Setting reproduction
#' @inheritSection setFishing Setting fishing
#' @inheritSection setPlankton Setting plankton dynamics
#' @inheritSection setResourceDynamics Setting resource dynamics
#' @inheritSection setResourceEncounter Setting resource encounter rate
#' @md
#' @export
#' @family functions for setting parameters
#' @examples
#' \dontrun{
#' params <- set_trait_model()
#' params@species_params$gamma[3] <- 1000
#' params <- setParams(params)
#' }
setParams <- function(params,
                      # setInteraction()
                      interaction = NULL,
                      # setPredKernel()
                      pred_kernel = NULL,
                      # setSearchVolume()
                      search_vol = NULL,
                      q = params@q,
                      f0 = params@f0,
                      # setIntakeMax()
                      intake_max = NULL,
                      n = params@n,
                      # setMetab()
                      metab = NULL,
                      p = params@p,
                      # setBMort
                      mu_b = NULL,
                      z0pre = 0.6,
                      z0exp = n - 1,
                      # setReproduction
                      maturity = NULL,
                      repro_prop = NULL,
                      srr = params@srr,
                      # setPlankton
                      kappa = params@kappa,
                      lambda = params@lambda,
                      r_pp = 10,
                      w_pp_cutoff = 10,
                      plankton_dynamics = NULL,
                      # setResourceDynamics
                      resource_dynamics = params@resource_dynamics,
                      resource_params = params@resource_params,
                      # setResourceEncounter
                      rho = NULL) {
    validObject(params)
    params <- setInteraction(params,
                             interaction = interaction)
    params <- setFishing(params)
    params <- setPredKernel(params,
                            pred_kernel = pred_kernel)
    params <- setIntakeMax(params,
                           intake_max = intake_max,
                           n = n)
    params <- setMetab(params,
                       metab = metab,
                       p = p)
    params <- setBMort(params,
                       mu_b = mu_b,
                       z0pre = z0pre,
                       z0exp = z0exp)
    params <- setSearchVolume(params,
                              search_vol = search_vol,
                              q = q)
    params <- setReproduction(params,
                              maturity = maturity,
                              repro_prop = repro_prop,
                              srr = srr)
    params <- setPlankton(params,
                          kappa = kappa,
                          lambda = lambda,
                          r_pp = r_pp,
                          w_pp_cutoff = w_pp_cutoff,
                          plankton_dynamics = plankton_dynamics)
    params <- setResourceDynamics(params,
                                  resource_dynamics = resource_dynamics,
                                  resource_params = resource_params)
    params <- setResourceEncounter(params,
                                   rho = rho,
                                   n = params@n)
    return(params)
}


#' Set species interation matrix
#' 
#' @section Setting interactions:
#' 
#' The species interaction matrix \eqn{\theta_{ij}}, is used when calculating the
#' food encounter rate in \code{\link{getEncounter}} and the predation mortality rate in
#' \code{\link{getPredMort}}. Its entries are dimensionless numbers between
#' 0 and 1 that characterise the strength at which predator species \eqn{i}
#' predates on prey species \eqn{j}. 
#' 
#' This function checks that the supplied interaction
#' matrix is valid and then stores it in the \code{interaction} slot of the
#' params object before returning that object.
#' 
#' The order of the columns and rows of the \code{interaction} argument should be the 
#' same as the order in the species params dataframe in the \code{params} object.
#' If you supply a named array then the function will check the order and warn 
#' if it is different.
#' 
#' The interaction of the species with the plankton are set via a column
#' \code{interaction_p} in the \code{species_params} data frame. Again the entries
#' have to be numbers between 0 and 1. By default this column is set to all
#' 1s.
#' 
#' @param params MizerParams object
#' @param interaction Interaction matrix of the species (predator by prey).
#'   Entries should be numbers between 0 and 1. See "Setting interactions"
#'   section below.
#' 
#' @return MizerParams object
#' @export
#' @family functions for setting parameters
#' @examples
#' \dontrun{
#' params <- set_trait_model()
#' interaction <- params@interaction
#' interaction[1, 3] <- 0
#' params <- setInteraction(params, interaction)
#' }
setInteraction <- function(params,
                           interaction = NULL) {
    assert_that(is(params, "MizerParams"))
    if (is.null(interaction)) {
        interaction <- params@interaction
    }
    # Check dims of interaction argument
    if (!identical(dim(params@interaction), dim(interaction))) {
        stop("interaction matrix is not of the right dimensions. Must be number of species x number of species.")
    }
    # Check that all values of interaction matrix are 0 - 1.
    if (!all((interaction >= 0) & (interaction <= 1))) {
        stop("Values in the interaction matrix must be between 0 and 1")
    }
    # In case user has supplied names to interaction matrix, check them.
    if (!is.null(dimnames(interaction))) {
        if (!is.null(names(dimnames(interaction)))) {
            if (!identical(names(dimnames(interaction)),
                           names(dimnames(params@interaction)))) {
                message(paste0("Note: Your interaction matrix has dimensions called: `",
                               toString(names(dimnames(interaction))),
                               "`. I expected 'predator, prey'. I will now ignore your names."))
            }
        }
        names(dimnames(interaction)) <- names(dimnames(params@interaction))
        if (!identical(dimnames(params@interaction),
                       dimnames(interaction))) {
            message("Note: Dimnames of interaction matrix do not match the order of species names in the species data.frame. I am now ignoring your dimnames so your interaction matrix may be in the wrong order.")
        }
    }
    params@interaction[] <- interaction
    
    # Check the interaction_p column in species_params
    message <- "Note: No interaction_p column in species data frame so assuming all species feed on plankton."
    species_params <- set_species_param_default(params@species_params,
                                                "interaction_p", 1,
                                                message = message)
    # Check that all values of interaction vector are 0 - 1.
    if (!all((species_params$interaction_p >= 0) & 
             (species_params$interaction_p <= 1))) {
        stop("Values in the plantkon interaction vector should be between 0 and 1")
    }
    params@species_params$interaction_p <- species_params$interaction_p
    
    return(params)
}


#' Set predation kernel
#' 
#' The predation kernel determines the distribution of prey sizes that a
#' predator feeds on. It is used in \code{\link{getEncounter}} when calculating
#' the rate at which food is encountered and in \code{\link{getPredRate}} when
#' calculating the rate at which a prey is predated upon. The predation kernel
#' can be a function of the predator/prey size ratio or it can be a function of
#' the predator size and the prey size separately. Both types can be set up with
#' this function.
#' 
#' @section Setting predation kernel:
#' \strong{Kernel dependent on predator to prey size ratio}
#' 
#' If the \code{pred_kernel} argument is not supplied, then this function sets a
#' predation kernel that depends only on the ratio of predator mass to prey
#' mass, not on the two masses independently. The shape of that kernel is then
#' determined by the \code{pred_kernel_type} column in species_params.
#'
#' The default pred_kernel_type is "lognormal". This will call the function
#' \code{\link{lognormal_pred_kernel}} to calculate the predation kernel.
#' An alternative pred_kernel type is "box", implemented by the functions
#' \code{\link{box_pred_kernel}}. These functions require certain species
#' parameters in the species_params data frame. For the lognormal kernel these
#' are \code{beta} and \code{sigma}, for the box kernely they are
#' \code{ppmr_min} and \code{ppmr_max}. They are explained in the help pages
#' for the kernel functions. No defaults are set for these parameters. If they
#' are missing from the species_params data frame then mizer will issue an
#' error message.
#'
#' You can use any other string as the type. If for example you choose "my" then
#' you need to define a function \code{my_pred_kernel} that you can model on the
#' existing functions like \code{\link{lognormal_pred_kernel}}.
#' 
#' When using a kernel that depends on the predator/prey size ratio only, mizer
#' does not need to store the entire three dimensional array in the MizerParams
#' object. Such an array can be very big when there is a large number of size
#' bins. Instead, mizer only needs to store two two-dimensional arrays that hold
#' Fourier transforms of the feeding kernel function that allow the encounter
#' rate and the predation rate to be calculated very efficiently. However, if
#' you need the full three-dimensional array you can calculate it with the
#' \code{\link{getPredKernel}} function.
#' 
#' \strong{Kernel dependent on both predator and prey size}
#' 
#' If you want to work with a feeding kernel that depends on predator mass and
#' prey mass independently, you can specify the full feeding kernel as a
#' three-dimensional array (predator species x predator size x prey size).
#' The dimensions are thus (no_sp, no_w, no_w_full). 
#'
#' You should use this option only if a kernel dependent only on the
#' predator/prey mass ratio is not appropriate. Using a kernel dependent on
#' predator/prey mass ratio only allows mizer to use fast Fourier transform
#' methods to significantly reduce the running time of simulations.
#'
#' The order of the predator species in \code{pred_kernel} should be the same
#' as the order in the species params dataframe in the `params` object. If you
#' supply a named array then the function will check the order and warn if it is
#' different.
#' 
#' @param params A MizerParams object
#' @param pred_kernel Optional. An array (species x predator size x prey size)
#'   that holds the predation coefficient of each predator at size on each prey
#'   size. If not supplied, a default is set as described in section "Setting
#'   predation kernel".
#' 
#' @return A MizerParams object
#' @export
#' @family functions for setting parameters
#' @examples
#' \dontrun{
#' ## Set up a MizerParams object
#' data(NS_species_params_gears)
#' data(inter)
#' params <- set_multispecies_model(NS_species_params_gears, inter)
#' 
#' ## If you change predation kernel parameters after setting up a model, you
#' # need to call setPredKernel
#' params@species_params["Cod", "beta"] <- 200
#' params <- setPredKernel(params)
#' 
#' ## You can change to a different predation kernel type
#' params@species_params$pred_kernel_type <- "box"
#' params@species_params$ppmr_min <- 2
#' params@species_params$ppmr_max <- 4
#' params <- setPredKernel(params)
#' 
#' ## If you need a kernel that depends also on prey size you need to define
#' # it yourself.
#' pred_kernel <- getPredKernel(params)
#' pred_kernel["Herring", , ] <- sweep(pred_kernel["Herring", , ], 2, 
#'                                     params@w_full, "*")
#' params<- setPredKernel(params, pred_kernel = pred_kernel)
#' }
setPredKernel <- function(params,
                          pred_kernel = NULL) {
    assert_that(is(params, "MizerParams"))
    if (!is.null(pred_kernel)) {
        # A pred kernel was supplied, so check it and store it
        assert_that(is.array(pred_kernel))
        # psi is used in the next line just because it has the right dimension
        assert_that(identical(dim(pred_kernel), 
                              c(dim(params@psi), length(params@w_full))))
        if (!is.null(dimnames(pred_kernel)) && 
            !all(dimnames(pred_kernel)[[1]] == params@species_params$species)) {
            stop(paste0("You need to use the same ordering of species as in the ",
                        "params object: ", toString(params@species_params$species)))
        }
        assert_that(all(pred_kernel >= 0))
        dimnames(pred_kernel) <- 
            list(sp = params@species_params$species,
                 w_pred = signif(params@w, 3),
                 w_prey = signif(params@w_full, 3))
        params@pred_kernel <- pred_kernel
        # Empty the Fourier transforms of kernel, to ensure that the FFT is not
        # used by model
        params@ft_pred_kernel_e <- array()
        params@ft_pred_kernel_p <- array()
        return(params)
    }
    
    ## Set a pred kernel dependent on predator/prey size ratio only
    # If pred_kernel_type is not supplied use "lognormal"
    params@species_params <- set_species_param_default(params@species_params,
                                                       "pred_kernel_type",
                                                       "lognormal")
    species_params <- params@species_params
    pred_kernel_type <- species_params$pred_kernel_type
    no_sp <- nrow(species_params)
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)
    # Vector of predator/prey mass ratios
    # The smallest predator/prey mass ratio is 1
    ppmr <- params@w_full / params@w_full[1]
    phis <- get_phi(species_params, ppmr)
    for (i in 1:no_sp) {
        phi <- phis[i, ]
        # Fourier transform of feeding kernel for evaluating available energy
        params@ft_pred_kernel_e[i, ] <- fft(phi)
        # Fourier transform of feeding kernel for evaluating predation rate
        ri <- max(which(phi > 0))  # index of largest ppmr
        phi_p <- rep(0, no_w_full)
        phi_p[(no_w_full - ri + 1):no_w_full] <- phi[(ri + 1):2]
        params@ft_pred_kernel_p[i, ] <- fft(phi_p)
    }
    
    return(params)
}


#' Get predation kernel
#' 
#' If no explicit predation kernel \eqn{\phi_i(w, w_p)} is stored in the params
#' object, then this function calculates it from the information in the species
#' parameter data frame in the params object.
#' 
#' For more detail about the predation kernel see \code{\link{setPredKernel}}.
#' 
#' @param params A MizerParams object
#' @return An array (predator x predator_size x prey_size)
#' @export
getPredKernel <- function(params) {
    assert_that(is(params, "MizerParams"))
    if (length(dim(params@pred_kernel)) > 1) {
        return(params@pred_kernel)
    }
    species_params <- set_species_param_default(params@species_params,
                                                "pred_kernel_type",
                                                "lognormal")
    pred_kernel_type <- species_params$pred_kernel_type
    no_sp <- nrow(species_params)
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)
    # Vector of predator/prey mass ratios
    # The smallest predator/prey mass ratio is 1
    ppmr <- params@w_full / params@w_full[1]
    phis <- get_phi(species_params, ppmr)
    pred_kernel <- 
        array(0,
              dim = c(no_sp, no_w, no_w_full),
              dimnames = list(sp = species_params$species,
                              w_pred = signif(params@w, 3),
                              w_prey = signif(params@w_full, 3)))
    for (i in 1:no_sp) {
        min_w_idx <- no_w_full - no_w + 1
        for (k in seq_len(no_w)) {
            pred_kernel[i, k, (min_w_idx - 1 + k):1] <-
                phis[i, 1:(min_w_idx - 1 + k)]
        }
    }
    return(pred_kernel)
}


#' Set search volume
#' 
#' @section Setting search volume:
#' The search volume \eqn{\gamma_i(w)} of an individual of species \eqn{i}
#' and weight \eqn{w} multiplies the predation kernel when
#' calculating the encounter rate in \code{\link{getEncounter}} and the 
#' predation rate in \code{\link{getPredRate}}.
#' 
#' The name "search volume" is a bit misleading, because \eqn{\gamma_i(w)} does
#' not have units of volume. It is simply a parameter that determines the rate
#' of predation. Its units depend on your choice, see section "Units in mizer".
#' If you have chose to work with total abundances, then it is a rate with units
#' 1/year. If you have chosen to work with abundances per m^2 then it has units
#' of m^2/year. If you have chosen to work with abundances per m^3 then it has
#' units of m^3/year.
#' 
#' If the \code{search_vol} argument is not supplied, then the search volume is set to
#' \deqn{\gamma_i(w) = \gamma_i w^q.} The values of \eqn{\gamma_i} are taken from
#' the \code{gamma} column in the species parameter dataframe. If the \code{gamma}
#' column is not supplied in the species parameter dataframe, a default is
#' calculated by the \code{\link{get_gamma_default}} function. Note that only
#' for predators of size \eqn{w = 1} gram is the value of the species parameter
#' \eqn{\gamma_i} the same as the value of the search volume \eqn{\gamma_i(w)}.
#' 
#' @param params MizerParams
#' @param search_vol Optional. An array (species x size) holding the search volume
#'   for each species at size. If not supplied, a default is set as described in
#'   the section "Setting search volume". 
#' @param q Exponent of the allometric search volume. Not needed if a
#'   \code{search_vol} array is specified.
#' 
#' @return MizerParams
#' @export
#' @family functions for setting parameters
#' @examples
#' \dontrun{
#' params <- set_trait_model()
#' params@species_params$gamma[3] <- 1000
#' params <- setSearchVolume(params)
#' }
setSearchVolume <- function(params, 
                            search_vol = NULL,
                            q = params@q) {
    assert_that(is(params, "MizerParams"),
                is.number(q))
    species_params <- params@species_params
    params@q <- q
    # If search_vol array is supplied, check it, store it and return
    if (!is.null(search_vol)) {
        assert_that(is.array(search_vol))
        assert_that(identical(dim(search_vol), dim(params@search_vol)))
        if (!is.null(dimnames(search_vol)) && 
            !all(dimnames(search_vol)[[1]] == species_params$species)) {
            stop(paste0("You need to use the same ordering of species in the search_vol array as in the ",
                        "params object: ", toString(species_params$species)))
        }
        assert_that(all(search_vol > 0))
        params@search_vol[] <- search_vol
        return(params)
    }
    # Calculate default for any missing gammas
    params@species_params$gamma <- get_gamma_default(params)
    
    params@search_vol[] <- outer(params@species_params$gamma, params@w^q)
    return(params)
}


#' Set maximum intake rate
#'
#' @section Setting maximum intake rate:
#' The maximum intake rate \eqn{h_i(w)} of an individual of species \eqn{i} and
#' weight \eqn{w} determines the feeding level, calculated with
#' \code{\link{getFeedingLevel}}. It is measured in grams/year.
#'
#' If the \code{intake_max} argument is not supplied, then the maximum intake
#' rate is set to \deqn{h_i(w) = h_i w^n.} 
#' The values of \eqn{h_i} (the maximum intake rate of an individual of size
#' 1 gram) are taken from the \code{h} column in the
#' species parameter dataframe. If the \code{h} column is not supplied in the
#' species parameter dataframe, it is calculated by the 
#' \code{\link{get_h_default}} function, using \code{f0} and the \code{k_vb}
#' column, if they are supplied.
#' 
#' If \eqn{h_i} is set to \code{Inf}, fish will consume all encountered food.
#'
#' @param params MizerParams
#' @param intake_max Optional. An array (species x size) holding the maximum
#'   intake rate for each species at size. If not supplied, a default is set as
#'   described in the section "Setting maximum intake rate".
#' @param n Scaling exponent of the intake rate.
#' @return MizerParams
#' @export
#' @family functions for setting parameters
setIntakeMax <- function(params, 
                         intake_max = NULL, 
                         n = params@n) {
    assert_that(is(params, "MizerParams"),
                is.number(n))
    species_params <- params@species_params
    params@n <- n
    
    # If intake_max array is supplied, check it, store it and return
    if (!is.null(intake_max)) {
        assert_that(is.array(intake_max),
                    identical(dim(intake_max), dim(params@intake_max)))
        if (!is.null(dimnames(intake_max)) && 
            !all(dimnames(intake_max)[[1]] == species_params$species)) {
            stop(paste0("You need to use the same ordering of species in the intake_max array as in the ",
                        "params object: ", toString(species_params$species)))
        }
        assert_that(all(intake_max > 0))
        params@intake_max[] <- intake_max
        return(params)
    }
    
    params@species_params$h <- get_h_default(params)
    
    params@intake_max[] <- outer(params@species_params$h, params@w^params@n)
    return(params)
}


#' Set metabolic rate
#' 
#' Sets the rate at which energy is used for metabolism and activity
#' 
#' @section Setting metabolic rate:
#' The metabolic rate is subtracted from the energy income rate to calculate
#' the rate at which energy is available for growth and reproduction, see
#' \code{\link{getEReproAndGrowth}}. It is measured in grams/year.
#' 
#' If the \code{metab} argument is not supplied, then the metabolic
#' rate \eqn{k_i(w)} for an individual of species \eqn{i} and size \eqn{w}
#' is set to \deqn{k_i(w) = ks_i\, w^p + k_i\, w,}{k_i(w) = ks_i w^p + k_i w}
#' where \eqn{ks_i w^p} represents the rate of standard metabolism and 
#' \eqn{k_i w} is the rate at which energy is expended on activity and movement.
#' The values of \eqn{ks_i} and \eqn{k_i} are taken from the \code{ks} and
#' \code{k} columns in the
#' species parameter dataframe. If these parameters are not supplied, the
#' defaults are \eqn{ks_i = 0.2 h_i} and \eqn{k_i = 0}, where \eqn{h_i} is the
#' coefficient of the maximum intake rate and is taken from the species
#' parameter data frame in \code{params}.
#' 
#' @param params MizerParams
#' @param metab Optional. An array (species x size) holding the metabolic rate
#'   for each species at size. If not supplied, a default is set as described in
#'   the section "Setting metabolic rate".
#' @param p Scaling exponent of the standard metabolic rate.
#' 
#' @return MizerParams
#' @export
#' @family functions for setting parameters
setMetab <- function(params, 
                     metab = NULL, 
                     p = params@p) {
    assert_that(is(params, "MizerParams"),
                is.number(p))
    species_params <- params@species_params
    params@p <- p
    if (!is.null(metab)) {
        assert_that(is.array(metab),
                    identical(dim(metab), dim(params@metab)))
        if (!is.null(dimnames(metab)) && 
            !all(dimnames(metab)[[1]] == species_params$species)) {
            stop(paste0("You need to use the same ordering of species as in the ",
                        "params object: ", toString(species_params$species)))
        }
        assert_that(all(metab > 0))
        params@metab[] <- metab
        return(params)
    }
    
    params <- set_species_param_default(params, "k", 0)
    params@species_params$ks <- get_ks_default(params)
    params@metab[] <- 
        outer(params@species_params$ks, params@w^p) +
        outer(params@species_params$k, params@w)
    return(params)
}


#' Set background mortality rate
#' 
#' @section Setting background mortality rate:
#' The background mortality is all the mortality that is not due to either
#' predation or fishing. It is a rate with units 1/year.
#' 
#' The \code{mu_b} argument allows you to specify a background mortality rate
#' that depends on species and body size. You can see an example of this in
#' the Examples section of the help page for \code{\link{setBMort}}.
#' 
#' If the \code{mu_b} argument is not supplied, then the background mortality
#' is assumed to depend only on the asymptotic size of the species, not on the
#' size of the individual: \eqn{\mu_{b.i}(w) = z_{0.i}}. The value of the
#' constant \eqn{z_0} for each species is taken from the \code{z0} column of the
#' species_params data frame, if that column exists. Otherwise it is calculated
#' as \deqn{z_{0.i} = {\tt z0pre}_i\, w_{inf}^{\tt z0exp}.}{z_{0.i} = z0pre_i w_{inf}^{z0exp}.}
#' 
#' @param params MizerParams
#' @param mu_b Optional. An array (species x size) holding the background
#'   mortality rate.
#' @param z0pre If \code{z0}, the mortality from other sources, is not a column
#'   in the species data frame, it is calculated as z0pre * w_inf ^ z0exp.
#'   Default value is 0.6.
#' @param z0exp If \code{z0}, the mortality from other sources, is not a column
#'   in the species data frame, it is calculated as \code{z0pre * w_inf ^ z0exp}.
#'   Default value is \code{n-1}.
#' 
#' @return MizerParams
#' @export
#' @family functions for setting parameters
#' @examples
#' \dontrun{
#' data("NS_species_params")
#' params <- set_multispecies_model(NS_species_params)
#'
#' #### Setting allometric death rate #######################
#' 
#' # Set coefficient for each species. Here we choose 0.1 for each species
#' mu_b0 <- rep(0.1, nrow(params@species_params))
#' 
#' # Multiply by power of size with exponent, here chosen to be -1/4
#' # The outer() function makes it an array species x size
#' mu_b <- outer(mu_b0, params@w^(-1/4))
#' 
#' # Change the background death in the params object
#' params <- setBMort(params, mu_b = mu_b)
#' }
setBMort <- function(params, mu_b = NULL, z0pre = 0.6, z0exp = params@n - 1) {
    assert_that(is(params, "MizerParams"))
    if (!is.null(mu_b)) {
        assert_that(is.array(mu_b),
                    identical(dim(mu_b), dim(params@mu_b)))
        params@mu_b[] <- mu_b
        return(params)
    }
    assert_that(is.number(z0pre), z0pre >= 0,
                is.number(z0exp))
    species_params <- params@species_params
    assert_that(noNA(species_params$w_inf))
    # Sort out z0 (background mortality)
    message <- ("Note: Using z0 = z0pre * w_inf ^ z0exp for missing z0 values.")
    params <- set_species_param_default(params, "z0",
                                        z0pre * species_params$w_inf^z0exp,
                                        message)
    params@mu_b[] <- params@species_params$z0
    return(params)
}


#' Set proportion of energy that is invested into reproduction
#' 
#' Sets the proportion of the total energy available for reproduction and growth
#' that is invested into reproduction as a function of the size of the
#' individual and sets the reproductive efficiency.
#' 
#' @section Setting reproduction:
#' For each species and at each size, the proportion of the available energy 
#' that is invested into reproduction is the product of two factors: the
#' proportion \code{maturity} of individuals that are mature and the proportion
#' \code{repro_prop} of the energy available to a mature individual that is 
#' invested into reproduction.
#' 
#' If the \code{maturity} argument is not supplied, then it is set to a sigmoidal 
#' maturity ogive that changes from 0 to 1 at around the maturity size:
#' \deqn{{\tt maturity}(w) = \left[1+\left(\frac{w}{w_{mat}}\right)^{-U}\right]^{-1}.}{
#'   maturity(w) = [1+(w/w_mat)^(-U)]^(-1)}
#' (To avoid clutter, we are not showing the species index in the equations.)
#' The maturity weights are taken from the \code{w_mat} column of the 
#' species_params data frame. Any missing maturity weights are set to 1/4 of the
#' asymptotic weight in the \code{w_inf} column.
#' The exponent \eqn{U} determines the steepness of the maturity ogive. By
#' default it is chosen as \eqn{U = 10}, however this can be overridden by
#' including a column \code{w_mat25} in the species parameter dataframe that
#' specifies the weight at which 25\% of individuals are mature, which sets
#' \eqn{U = \log(3) / \log(w_{mat} / w_{25}).}{U = log(3) / log(w_mat / w_25).}
#' 
#' The sigmoidal function given above would strictly reach 1 only asymptotically.
#' Mizer instead sets the function equal to 1 already at the species' 
#' maximum size, taken from the compulsory \code{w_inf} column in the
#' \code{species_params} data frame.
#' 
#' If the \code{repro_prop} argument is not supplied, it is set to the
#' allometric form
#' \deqn{{\tt repro\_prop}(w) = \left(\frac{w}{w_{inf}}\right)^{m-n}.}{
#'   repro_prop = (w/w_inf)^(m - n).}
#' Here \eqn{n} is the scaling exponent of the energy income rate. Hence
#' the exponent \eqn{m} determines the scaling of the investment into
#' reproduction for mature individuals. By default it is chosen to be 
#' \eqn{m = 1} so that the rate at which energy is invested into reproduction 
#' scales linearly with the size. This default can be overridden by including a 
#' column \code{m} in the species parameter dataframe. The asymptotic sizes
#' are taken from the compulsory \code{w_inf} column in the species_params
#' data frame.
#' 
#' The reproductive efficiency, i.e., the proportion of energy allocated to
#' reproduction that results in offspring biomass, is set from the 
#' \code{erepro} column in the species_params data frame. If that is not
#' provided the default is set to 1 (which you will want to override).
#' The offspring biomass divided by the egg biomass gives the
#' density-independent rate of egg production, returned by \code{\link{getRDI}}.
#' 
#' Mizer allows some density dependence in the production of eggs by putting
#' the density-independent rate of egg production through a stock-recruitment
#' function. The result is returned by \code{\link{getRDD}}. The
#' stock-recruitment function is specified by the \code{srr} argument. The default
#' is the BevertonHolt function \code{\link{srrBevertonHolt}}, which requires
#' an \code{r_max} column in the species_params data frame giving the maximum
#' egg production rate. If this column does not exist, it is initialised to 
#' \code{Inf}, leading to no density-dependence.
#' 
#' @param params A MizerParams object
#' @param maturity Optional. An array (species x size) that holds the proportion
#'   of individuals of each species at size that are mature. If not supplied, a
#'   default is set as described in the section "Setting reproduction".
#' @param repro_prop Optional. An array (species x size) that holds the
#'   proportion of consumed energy that a mature individual allocates to
#'   reproduction for each species at size. If not supplied, a default is set as
#'   described in the section "Setting reproduction".
#' @param srr Optional. The stock recruitment function. 
#' 
#' @return The MizerParams object.
#' @export
#' @family functions for setting parameters
setReproduction <- function(params, maturity = NULL, repro_prop = NULL,
                            srr = params@srr) {
    assert_that(is(params, "MizerParams"),
                is.function(srr))
    species_params <- params@species_params

    # Check maximum sizes
    if (!("w_inf" %in% colnames(species_params))) {
        stop("The maximum sizes of the species must be specified in the w_inf column of the species parameter data frame.")
    }
    missing <- is.na(species_params$w_inf)
    if (any(missing)) {
        stop(paste("The following species are missing data for their maximum size w_inf:"),
             toString(species_params$species[missing]))
    }
    if (any(species_params$w_inf <= species_params$w_min)) {
        stop("Some of the asymptotic sizes are smaller than the egg sizes.")
    }
    # # Round maximum sizes to nearest grid point
    # for (i in seq_along(species_params$w_inf)) {
    #     idx <- which.min(abs(species_params$w_inf[i] - params@w))
    #     params@species_params$w_inf[i] < params@w[idx]
    # }
    
    if (!is.null(maturity)) {
        assert_that(is.array(maturity),
                    identical(dim(maturity), dim(params@psi)))
        if (!is.null(dimnames(maturity)) && 
            !all(dimnames(maturity)[[1]] == species_params$species)) {
            stop(paste0("You need to use the same ordering of species as in the ",
                        "params object: ", toString(species_params$species)))
        }
    } else {
        # Check maturity sizes
        if (!("w_mat" %in% colnames(species_params))) {
            species_params$w_mat <- rep(NA, nrow(species_params))
        }
        missing <- is.na(species_params$w_mat)
        if (any(missing)) {
            message("Note: The following species were missing data for ",
                    "their maturity size w_mat: ",
                    toString(species_params$species[missing]),
                    ". These have been set to 1/4 w_inf.")
            species_params$w_mat[missing] <- species_params$w_inf[missing] / 4
        }
        assert_that(all(species_params$w_mat > species_params$w_min))
        
        # Set defaults for w_mat25
        if (!("w_mat25" %in% colnames(species_params))) {
            species_params$w_mat25 <- species_params$w_mat/(3^(1/10))
        }
        missing <- is.na(species_params$w_mat25)
        if (any(missing)) {
            species_params$w_mat25[missing] <- species_params$w_mat[missing]/(3^(1/10))
        }
        # Check w_mat25
        assert_that(all(species_params$w_mat25 > species_params$w_min))
        assert_that(all(species_params$w_mat25 < species_params$w_mat))
        params@species_params$w_mat25 <- species_params$w_mat25
        
        maturity <- 
            unlist(
                tapply(params@w, 1:length(params@w),
                       function(wx, w_inf, w_mat, w_mat25) {
                           U <- log(3) / log(w_mat / w_mat25)
                           return((1 + (wx / w_mat)^-U)^-1)
                       },
                       w_inf = species_params$w_inf,
                       w_mat = species_params$w_mat,
                       w_mat25 = species_params$w_mat25
                )
            )
    }
    assert_that(all(maturity >= 0 & maturity <= 1))
    params@maturity[] <- maturity
    
    if (!is.null(repro_prop)) {
        assert_that(is.array(repro_prop),
                    identical(dim(repro_prop), dim(params@psi)))
        if (!is.null(dimnames(repro_prop)) && 
            !all(dimnames(repro_prop)[[1]] == species_params$species)) {
            stop(paste0("You need to use the same ordering of species as in the ",
                        "params object: ", toString(species_params$species)))
        }
    } else {
        # Set defaults for m
        params <- set_species_param_default(params, "m", 1)
        
        repro_prop <- array(
            unlist(
                tapply(params@w, 1:length(params@w),
                       function(wx, w_inf, mn) (wx / w_inf)^(mn),
                       w_inf = params@species_params$w_inf,
                       mn = params@species_params$m - params@n
                )
            ), dim = c(nrow(species_params), length(params@w)))
    }
    
    params@psi[] <- params@maturity * repro_prop
    
    # psi should never be larger than 1
    params@psi[params@psi > 1] <- 1
    # For reasons of efficiency we next set all very small values to 0 
    # Set w < 10% of w_mat to 0
    params@psi[outer(species_params$w_mat * 0.1, params@w, ">")] <- 0
    # Set all w > w_inf to 1
    params@psi[outer(species_params$w_inf, params@w, "<")] <- 1
    assert_that(all(params@psi >= 0 & params@psi <= 1))
    
    # If no erepro (reproductive efficiency), then set to 1
    params <- set_species_param_default(params, "erepro", 1)
    assert_that(all(params@species_params$erepro > 0))
    
    # srr must have two arguments: rdi amd species_params
    if (!isTRUE(all.equal(names(formals(srr)), c("rdi", "species_params")))) {
        stop("Arguments of srr function must be 'rdi' and 'species_params'")
    }
    params@srr <- srr
    if (identical(params@srr, srrBevertonHolt)) {
        set_species_param_default(params, "r_max", Inf)
    }
    
    return(params)
}

#' Set up plankton
#' 
#' @section Setting plankton dynamics:
#' Still need to document
#' 
#' @param params A MizerParams object
#' @param kappa Carrying capacity of the plankton spectrum.
#' @param lambda Exponent of the plankton spectrum.
#' @param r_pp Growth rate of the primary productivity. Default is 10 g/year.
#' @param w_pp_cutoff The upper cut off size of the plankton spectrum. 
#'   Default is 10 g.
#' @param plankton_dynamics Function that determines plankton dynamics by
#'   calculating the plankton spectrum at the next time step from the current
#'   state.
#' 
#' @return A MizerParams object
#' @export
#' @family functions for setting parameters
setPlankton <- function(params,
                        kappa = params@kappa,
                        lambda = params@lambda,
                        r_pp = 10, 
                        w_pp_cutoff = 10,
                        plankton_dynamics = NULL) {
    assert_that(is(params, "MizerParams"),
                is.number(kappa), kappa > 0,
                is.number(lambda),
                is.number(r_pp), r_pp > 0,
                is.number(w_pp_cutoff),
                is.number(params@n))
    params@kappa <- kappa
    params@lambda <- lambda
    # weight specific plankton growth rate
    params@rr_pp[] <- r_pp * params@w_full^(params@n - 1)
    # the plankton carrying capacity
    params@cc_pp[] <- kappa*params@w_full^(-lambda)
    params@cc_pp[params@w_full > w_pp_cutoff] <- 0
    if (!is.null(plankton_dynamics)) {
        assert_that(is.function(plankton_dynamics))
        params@plankton_dynamics <- plankton_dynamics
    }
    
    return(params)
}

#' Set resource dynamics
#' 
#' @section Setting resource dynamics:
#' Besides the size-structured planktonic resource, mizer can also model any number
#' of unstructured resource components. Such unstructured components are
#' appropriate whenever the predation on these componente is not size based.
#' Examples include detritus as a resource for detritovores, carrion as a
#' resource for scavengers, or macroflora on which fish can graze. 
#' 
#' During a simulation using [project()], the biomasses of the resources are
#' updated at each time step by calling the functions specified in the
#' `resource_dynamics` list, which has one named entry for each unstructured
#' resource component. By default a model is set up without unstructured
#' resource components, so you do not need to provide this list. But if you
#' do want to use it, you can see an example of how to set up a
#' `resource_dynamics` list in the Examples section of
#' `setResourceDynamics()`.
#' 
#' Mizer provides two functions that you can use to model resource dynamics:
#' [detritus_dynamics()] and [carrion_dynamics()], but you can easily implement
#' others by following those templates.
#' As you can see in the documentation of these functions, their arguments are:
#' the `MizerParams` object `params`, the current fish size spectra `n`, the
#' current plankton spectrum `n_pp`, the current resource biomasses `B` and the
#' current rates calculated by the [getRates()] function.
#' 
#' The other arguments to the resource dynamics functions are model parameters,
#' like for example growth rates. These need to be provided in the
#' `resource_params` argument which is a named list. One model parameter that
#' should always be present in this list is the rate of change due to external
#' causes. This should be given a name of the form `resource_external` where
#' `resource` should be replaced by the name of the resource, see for example
#' `detritus_external` in [detritus_dynamics()].
#' 
#' When writing your own resource dynamics functions, you can choose any names
#' for your other model parameters, but you must make sure not to use the same
#' name in the function for another resource component. One way to ensure this
#' is to prefix all parameter namse with your resource name.
#' 
#' The dynamics for a resource should always have a loss term accounting for
#' the consumption of the resource. This should always have the form used in the
#' example function [detritus_dynamics()], in order to be in agreement with the
#' feeding by consumers that is set with \code{\link{setResourceEncounter}}.
#' 
#' @param params A MizerParams object
#' @param resource_dynamics A named list of functions that determine the
#'   dynamics of the unstructured resources by calculating their biomasses at
#'   the next time step from the current state. An empty list if the model
#'   does not have unstructured resources.
#' @param resource_params A named list of parameters needed by the
#'   \code{resource_dynamics} functions. An empty list if no parameters are
#'   needed.
#' 
#' @return A MizerParams object
#' @export
#' @md
#' @family functions for setting parameters
setResourceDynamics <- function(params,
                         resource_dynamics = NULL,
                         resource_params = NULL) {
    assert_that(is(params, "MizerParams"))
    if (!is.null(resource_dynamics)) {
        assert_that(is.list(resource_dynamics))
        no_res <- length(resource_dynamics)
        resource_names = names(resource_dynamics)
            if (no_res > 0  && is.null(resource_names)) {
                stop("The resource_dynamics list must be a named list.")
            }
            params@resource_dynamics <- resource_dynamics
            params@rho <- 
                array(NA,
                      dim = c(nrow(params@species_params),
                              no_res,
                              length(params@w)),
                      dimnames = list(sp = params@species_params$species,
                                      res = resource_names,
                                      w = signif(params@w, 3)))
            if (length(params@initial_B) != no_res) {
                params@initial_B <- rep(0, no_res)
                names(params@initial_B) <- resource_names
            }
    }
    if (!is.null(resource_params)) {
        assert_that(is.list(resource_params))
        params@resource_params <- resource_params
    }
    return(params)
}

#' Set resource encounter rate
#' 
#' @section Setting resource encounter rate:
#' The resource encounter rate \eqn{\rho_{id}(w)} (units 1/year) determines the
#' rate at which an individual of species \eqn{i} encounters biomass of resource
#' \eqn{d}, so that the contribution from all unstructured resources to the
#' total encounter rate is
#' \deqn{E_{u.i}(w) = \sum_d\rho_{id}(w) B_d,} 
#' where \eqn{B_d} is the biomass of the d-th unstructured resource component.
#' 
#' Resource consumption is subject to satiation in the same way as other food,
#' so that a consumer only consumes a fraction \eqn{1-f_i(w)} of the encountered
#' resource biomass, where \eqn{f_i(w)} is the feeding level.
#' 
#' If the \code{rho} array is not supplied, then the resource encounter rate is
#' set to a power law
#' \deqn{\rho_{id}(w) = \rho_{id} w^n.}
#' The coefficients \eqn{\rho_{id}} are parameters in the species_params dataframe.
#' For example if there is a resource called "detritus" then the species_params
#' data frame needs to have a column called \code{rho_detritus} and similarly
#' for each other resource.
#' 
#' If the \code{rho} array is supplied, the ordering of the entries in the array
#' is important. The order of the species in the first array dimension needs to
#' be the same as that in the species parameter dataframe. The order of the
#' resources in the second array dimension must be the same as in the list of
#' resource dynamics. The third dimension is the size dimension.
#' 
#' @param params A MizerParams object
#' @param rho Optional. An array (species x resource x size)
#'   holding the rate at which a consumer of a particular size feeds on each
#'   resource. Described in the section "Setting resource encounter rate".
#' @param n Scaling exponent of the intake rate.
#' 
#' @return A MizerParams object
#' @export
#' @family functions for setting parameters
setResourceEncounter <- function(params, rho = NULL, n = params@n) {
    assert_that(is(params, "MizerParams"),
                is.number(n))
    params@n <- n
    if (is.null(rho)) {
        # Use columns in species_params
        for (res in names(params@resource_dynamics)) {
            rho_var <- paste0("rho_", res)
            if (!rho_var %in% names(params@species_params)) {
                stop(paste("The species_params data frame needs a column ",
                           rho_var))
            }
            params@rho[, res, ] <- 
                outer(params@species_params[[rho_var]], params@w^n)
        }
        return(params)
    }
    # Check validity of arguments
    assert_that(is.array(rho),
                length(dim(rho)) == 3)
    if (nrow(params@species_params) != dim(rho)[1]) {
        stop("The first dimension of the rho argument should equal the number of species.")
    }
    no_res <- dim(rho)[2]
    if (length(params@resource_dynamics) != no_res) {
        stop("The second dimension of the rho argument should equal the number of resources.")
    }
    if (is.character(dimnames(rho)["res"])) {
        assert_that(are_equal(dimnames(rho)["res"],
                              names(params@resource_dynamics)))
    }
    if (length(params@w) != dim(rho)[3]) {
            stop("The third dimension of the rho array should have one entry for every consumer size.")
    }
    assert_that(all(rho > 0))
    params@rho[] <- rho
    
    return(params)
}


#' Set fishing parameters
#' 
#' @section Setting fishing:
#' In `mizer`, fishing mortality is imposed on species by fishing gears. The
#' total fishing mortality is obtained by summing over the mortality from all
#' gears,
#' \deqn{\mu_{f.i}(w) = \sum_g F_{g,i}(w),}
#' where the fishing mortality \eqn{F_{g,i}(w)} imposed by gear \eqn{g} on
#' species \eqn{i} at size \eqn{w} is calculated as:
#' \deqn{F_{g,i}(w) = S_{g,i}(w) Q_{g,i} E_{g},}
#' where \eqn{S} is the selectivity by species, gear and size, \eqn{Q} is the 
#' catchability by species and gear and \eqn{E} is the fishing effort by gear.
#' At the moment a species can only be selected by one fishing gear, although 
#' each gear can select more than one species (this is a limitation with the 
#' current package that will be developed in future releases).
#' 
#' \strong{Selectivity}
#' 
#' The selectivity at size of each gear has a range between 0 (not selected at
#' that size) to 1 (fully selected at that size). It is given by a selectivity
#' function. The name of the selectivity function is given by the `sel_func`
#' column in the species parameters data frame. Some selectivity functions are
#' included in the package: `knife_edge()` and `sigmoid_length()`. New functions
#' can be defined by the user. Each gear has the same selectivity function for
#' all the species it selects, but the parameter values for each species may be
#' different, e.g. the lengths of species that a gear selects may be different.
#' 
#' Each selectivity function has a range of arguments. Values for these
#' arguments must be included as columns in the species parameters data.frame.
#' The names of the columns must exactly match the names of the arguments. For
#' example, the default selectivity function is `knife_edge()` which has sudden
#' change of selectivity from 0 to 1 at a certain size.
#' In its help page you can see that the `knife_edge()` function has arguments `w` and
#' `knife_edge_size` The first argument, `w`, is size (the function calculates
#' selectivity at size). All selectivity functions must have `w` as the first
#' argument. The values for the other arguments must be found in the species
#' parameters data.frame. So for the `knife_edge()` function there should be a
#' `knife_edge_size` column. Because `knife_edge()` is the default
#' selectivity function, the `knife_edge_size` argument has a default
#' value = `w_mat`.
#' 
#' \strong{Catchability}
#' 
#' Catchability is used as an additional scalar to make the link between gear 
#' selectivity, fishing effort and fishing mortality. For example, it can be set 
#' so that an effort of 1 gives a desired fishing mortality.
#' In this way effort can then be specified relative to a 'base effort', e.g. 
#' the effort in a particular year.
#' 
#' Because currently mizer only allows one gear to select each species, the
#' catchability matrix \eqn{Q_{g,i}} reduces to a catchability vector 
#' \eqn{Q_{i}}, and this is given as a column `catchability` in the species
#' parameter data frame. If it is not specified, it defaults to 1.
#' 
#' Fishing effort is not stored in the `MizerParams` object.
#' Instead, effort is set when the simulation is run and can vary through time 
#' with `project()`.
#'         
#' 
#' @param params A MizerParams object
#' 
#' @return MizerParams object
#' @export
#' @md
#' @family functions for setting parameters
setFishing <- function(params) {
    assert_that(is(params, "MizerParams"))
    species_params <- params@species_params
    no_sp <- nrow(species_params)
    
    # If no gear_name column in species_params, then named after species
    if (!("gear" %in% colnames(species_params))) {
        species_params$gear <- species_params$species
    }
    
    # If no sel_func column in species_params, set to 'knife_edge'
    if (!("sel_func" %in% colnames(species_params))) {
        message("Note: No sel_func column in species data frame. Setting selectivity to be 'knife_edge' for all species.")
        species_params$sel_func <- 'knife_edge'
        # Set default selectivity size
        if (!("knife_edge_size" %in% colnames(species_params))) {
            message("Note: No knife_edge_size column in species data frame. Setting knife edge selectivity equal to w_mat.")
            species_params$knife_edge_size <- species_params$w_mat
        }
    }
    
    # If no catchability column in species_params, set to 1
    species_params <- set_species_param_default(species_params,
                                                "catchability", 1)
    
    # At the moment, each species is only caught by 1 gear so in species_params
    # there are the columns: gear_name and sel_func.
    # BEWARE! This routine assumes that each species has only one gear operating on it
    # So we can just go row by row through the species parameters
    # However, I really hope we can do something better soon
    for (g in 1:nrow(species_params)) {
        # Do selectivity first
        # get args
        # These as.characters are annoying - but factors everywhere
        arg <- names(formals(as.character(species_params[g, 'sel_func'])))
        # lop off w as that is always the first argument of the selectivity functions
        arg <- arg[!(arg %in% "w")]
        if (!all(arg %in% colnames(species_params))) {
            stop("Some arguments needed for the selectivity function are ",
                 "missing in the parameter dataframe.")
        }
        # Check that there are no missing values for selectivity parameters
        if (any(is.na(as.list(species_params[g, arg])))) {
            stop("Some selectivity parameters are NA.")
        }
        # Call selectivity function with selectivity parameters
        par <- c(w = list(params@w), as.list(species_params[g, arg]))
        sel <- do.call(as.character(species_params[g, 'sel_func']), args = par)
        # Dump Sel in the right place
        params@selectivity[as.character(species_params[g, 'gear']), g, ] <- sel
        # Now do catchability
        params@catchability[as.character(species_params[g,'gear']), g] <- 
            species_params[g, "catchability"]
    }
    
    return(params)
}

#' Set initial abundances
#'
#' Sets the slots in the \code{MizerParams} object holding the initial
#' abundances, \code{initial_n}, \code{initial_n_pp} and \code{initial_B}.
#'
#' @param params A \code{MizerParams} object
#' @param initial_n The initial abundances of species. A matrix with dimensions
#'   species x size. The order of species must be the same as in the MizerParams
#'   argument. Optional. Ignored if `sim` is supplied.
#' @param initial_n_pp The initial abundances of plankton. A numeric vector.
#'   Optional. Ignored if `sim` is supplied.
#' @param initial_B The initial biomasses of the unstructured resources. A named
#'   vector with one entry for each resource. Optional. Ignored if `sim` is
#'   supplied.
#' @param sim A \code{MizerSim} object. If supplied, the `initial_n`, 
#'   `initial_n_pp` and `initial_B` arguments are ignored and the information
#'   is taken from the last timestep of the simulation in `sim`.
#' @export
#' @family functions for setting parameters
setInitial <- function(params,
                          initial_n = params@initial_n,
                          initial_n_pp = params@initial_n_pp,
                          initial_B = params@initial_B,
                          sim) {
    if (!missing(sim)) {
        assert_that(is(sim, "MizerSim"))
        no_t <- dim(sim@n)[1]
        initial_n < sim@params@initial_n # Needed to get the right dimensions
        initial_n[] <- sim@n[no_t, , ]
        initial_n_pp < sim@params@initial_n_pp # Needed to get the right dimensions
        initial_n_pp[] <- sim@n_pp[no_t, ]
        initial_B < sim@params@initial_B # Needed to get the right dimensions
        initial_B[] <- sim@B[no_t, ]
    }
    assert_that(identical(dim(initial_n), dim(params@initial_n)),
                all(initial_n >= 0),
                identical(dim(initial_n_pp), dim(params@initial_n_pp)),
                all(initial_n_pp >= 0),
                identical(dim(initial_B), dim(params@initial_B)),
                all(initial_B >= 0))
    if (!is.null(dimnames(initial_n)) &&
        !identical(dimnames(initial_n), dimnames(params@initial_n))) {
        warning("The dimnames of initial_n are not as expected. I will ignore them.")
    }
    if (!is.null(dimnames(initial_n_pp)) &&
        !identical(dimnames(initial_n_pp), dimnames(params@initial_n_pp))) {
        warning("The dimnames of initial_n_pp are not as expected. I will ignore them.")
    }
    if (!is.null(dimnames(initial_B)) &&
        !identical(dimnames(initial_B), dimnames(params@initial_B))) {
        warning("The dimnames of initial_B are not as expected. I will ignore them.")
    }
    params@initial_n[] <- initial_n
    params@initial_n_pp[] <- initial_n_pp
    params@initial_B[] <- initial_B
    return(params)
}


#' Beverton Holt stock-recruitment function
#' 
#' @param rdi x
#' @param species_params x
#' 
#' @return rdd
#' @export
srrBevertonHolt <- function(rdi, species_params) {
    return(rdi / (1 + rdi/species_params$r_max))
}

#' Set a species parameter to a default value
#'
#' If the species parameter does not yet exist in the species parameter data
#' frame, then create it and fill it with the default. Otherwise use the default
#' only to fill in any NAs. Optionally gives a message if the parameter
#' did not already exist.
#' @param object Either a MizerParams object or a species parameter data frame
#' @param parname A string with the name of the species parameter to set
#' @param default A single defaul value or a vector with one default value for
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
        species_params <- cbind(species_params, default)
        colnames(species_params)[[ncol(species_params)]] <- parname
    } else {
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


#' Get values from feeding kernel function
#' 
#' This involves finding the feeding kernel function for each species, using the
#' pred_kernel_type parameter in the species_params data frame, checking that it
#' is valid and all its arguments are contained in the species_params data
#' frame, and then calling this function with the ppmr vector.
#' 
#' @param species_params A species parameter data frame
#' @param ppmr Values of the predator/prey mass ratio at which to evaluate the
#'   predation kernel function
#' @return An array (species x ppmr) with the values of the predation kernel
#'   function
#' @export
#' @keywords internal
#' @concept helper
get_phi <- function(species_params, ppmr) {
    assert_that(is.data.frame(species_params))
    no_sp <- nrow(species_params)
    species_params <- set_species_param_default(species_params, 
                                                "pred_kernel_type",
                                                "lognormal")
    phis <- array(dim = c(no_sp, length(ppmr)))
    for (i in 1:no_sp) {
        pred_kernel_func_name <- paste0(species_params$pred_kernel_type[i],
                                        "_pred_kernel")
        pred_kernel_func <- get0(pred_kernel_func_name)
        assert_that(is.function(pred_kernel_func))
        args <- names(formals(pred_kernel_func))
        if (!("ppmr" %in% args)) {
            stop(paste("The predation kernel function",
                       pred_kernel_func_name,
                       "needs the argument 'ppmr'."))
        }
        # lop off the compulsory arg
        args <- args[!(args %in% "ppmr")]
        missing <- !(args %in% colnames(species_params))
        if (any(missing)) {
            stop(paste("The following arguments for the predation kernel function",
                       pred_kernel_func_name,
                       "are missing from the parameter dataframe:",
                       toString(args[missing])))
        }
        pars <- c(ppmr = list(ppmr), as.list(species_params[i, args]))
        phi <- do.call(pred_kernel_func_name, args = pars)
        if (any(is.na(phi)) && 
            (species_params$interaction_p[i] > 0 ||
             any(interaction[i, ] > 0))) {
            stop(paste("The function", pred_kernel_func_name,
                       "returned NA. Did you correctly specify all required",
                       "parameters in the species parameter dataframe?"))
        }
        phis[i, ] <- phi
    }
    return(phis)
}

#' Get default value for h
#' 
#' Fills in any missing values for h according to the formula
#' \deqn{h = 3 k_{vb} w_{inf}^{1/3}/ (\alpha f_0)}.
#' @param params A MizerParams object
#' @return A vector with the values of h for all species
#' @export
#' @keywords internal
#' @concept helper
get_h_default <- function(params) {
    species_params <- params@species_params
    if (!("h" %in% colnames(species_params))) {
        species_params$h <- rep(NA, nrow(species_params))
    }
    missing <- is.na(species_params$h)
    if (any(missing)) {
        assert_that(is.number(params@f0),
                    noNA(species_params$alpha))
        message("Note: No h provided for some species, so using f0 and k_vb to calculate it.")
        if (!("k_vb" %in% colnames(species_params))) {
            stop("\tExcept I can't because there is no k_vb column in the species data frame")
        }
        h <- ((3 * species_params$k_vb) / (species_params$alpha * params@f0)) * 
            (species_params$w_inf ^ (1/3))
        if (any(is.na(h[missing]))) {
            stop("Could not calculate h, perhaps k_vb is missing?")
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
#' Currently this is implemented only for the lognormal predation kernel and uses
#' an analytic expression for the encounter rate, but in future we can do this
#' numerically for any predation kernel. Also currently feeding on unstructured
#' resources is not taken into account.
#' @param params A MizerParams object
#' @return A vector with the values of gamma for all species
#' @export
#' @keywords internal
#' @concept helper
get_gamma_default <- function(params) {
    species_params <- params@species_params
    if (!("gamma" %in% colnames(species_params))) {
        species_params$gamma <- rep(NA, nrow(species_params))
    }
    missing <- is.na(species_params$gamma)
    if (any(missing)) {
        if (!all(params@species_params$pred_kernel_type == "lognormal")) {
            stop("Calculation of a default for gamma has not yet been implemented for kernels that are not lognormal, so you will have to supply a gamma column in the species_params data frame.")
        }
        assert_that(is.number(params@lambda),
                    is.number(params@kappa),
                    is.number(params@f0))
        message("Note: Using f0, h, beta, sigma, lambda and kappa to calculate gamma.")
        lm2 <- params@lambda - 2
        ae <- sqrt(2 * pi) * species_params$sigma * species_params$beta^lm2 *
            exp(lm2^2 * species_params$sigma^2 / 2) *
            # The factor on the following lines takes into account the cutoff
            # of the integral at 0 and at beta + 3 sigma
            (pnorm(3 - lm2 * species_params$sigma) + 
                 pnorm(log(species_params$beta)/species_params$sigma + 
                           lm2 * species_params$sigma) - 1)
        if (!"h" %in% names(params@species_params) || 
            any(is.na(species_params$h))) {
            species_params$h <- get_h_default(params)
        }
        gamma_default <- (species_params$h / (params@kappa * ae)) * 
            (params@f0 / (1 - params@f0))
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
#' Fills in any missing values for ks according to the formula
#' \eqn{ks_i = 0.2 h_i}.
#' 
#' @param params A MizerParams object
#' @return A vector with the values of ks for all species
#' @export
#' @keywords internal
#' @concept helper
get_ks_default <- function(params) {
    species_params <- params@species_params
    if (!"h" %in% names(params@species_params) ||
        any(is.na(params@species_params$h))) {
        params@species_params$h <- get_h_default(params)
    }
    message <- ("Note: No ks column in species data frame so using ks = h * 0.2.")
    params <- set_species_param_default(params, "ks",
                                        params@species_params$h * 0.2,
                                        message)
    if (any(is.na(params@species_params$ks) | 
            is.infinite(params@species_params$ks))) {
        stop(paste("Could not calculate default values for the missing species",
             "parameter ks. Got:", params@species_params$ks))
    }
    return(params@species_params$ks)
}