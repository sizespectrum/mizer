# Class specification and constructors for mizer base parameters class
# Class has members to store parameters of size based model
# Copyright 2012 Finlay Scott and Julia Blanchard.
# Copyright 2018 Gustav Delius and Richard Southwell.
# Development has received funding from the European Commission's Horizon 2020 
# Research and Innovation Programme under Grant Agreement No. 634495 
# for the project MINOUW (http://minouw-project.eu/).
# Distributed under the GPL 3 or later 
# Maintainer: Gustav Delius, University of York, <gustav.delius@york.ac.uk>

# Naming conventions:
# S4 classes and constructors: AClass
# Functions: aFunction
# Variables: a_variable

# Validity function ---------------------------------------------------------
# Not documented as removed later on
validMizerParams <- function(object) {
    
    errors <- character()
    # grab some dims
    no_w <- length(object@w)
    no_w_full <- length(object@w_full)
    w_idx <- (no_w_full - no_w + 1):no_w_full
    
    # Check weight grids ----
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
    if (any(object@w[] != object@w_full[w_idx])) {
        msg <- "The later entries of w_full should be equal to those of w."
        errors <- c(errors, msg)
    }
    if (any(object@dw[] != object@dw_full[w_idx])) {
        msg <- "The later entries of dw_full should be equal to those of dw."
        errors <- c(errors, msg)
    }

    # Check the array dimensions are good ----
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
    # Check names of dimnames of arrays ----
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
    # Check the vector slots ----
    if (length(object@rr_pp) != length(object@w_full)) {
        msg <- "rr_pp must be the same length as w_full"
        errors <- c(errors, msg)
    }
    if (length(object@cc_pp) != length(object@w_full)) {
        msg <- "cc_pp must be the same length as w_full"
        errors <- c(errors, msg)
    }
    
    # TODO: Rewrite the following into a test of the @rates_funcs slot ----
    # SRR 
    # if (!is.string(object@srr)) {
    #     msg <- "srr needs to be specified as a string giving the name of the function"
    #     errors <- c(errors, msg)
    # } else {
    #     if (!exists(object@srr)) {
    #         msg <- paste0("The stock-recruitment function ",
    #                       object@srr,
    #                       "does not exist.")
    #         errors <- c(errors, msg)
    #     } else {
    #         srr <- get(object@srr)
    #         if (!is.function(srr)) {
    #             msg <- "The specified srr is not a function."
    #             errors <- c(errors, msg)
    #         } else {
    #             # Must have two arguments: rdi amd species_params
    #             if (!isTRUE(all.equal(names(formals(srr)), c("rdi", "species_params")))) {
    #                 msg <- "Arguments of srr function must be 'rdi' and 'species_params'"
    #                 errors <- c(errors, msg)
    #             }
    #         }
    #     }
    # }
    
    # Should not have legacy r_max column (has been renamed to R_max)
    if ("r_max" %in% names(object@species_params)) {
        msg <- "The 'r_max' column in species_params should be called 'R_max'. You can use 'upgradeParams()' to upgrade your params object."
        errors <- c(errors, msg)
    }
    # # species_params data.frame must have columns: 
    # # species, z0, alpha, eRepro
    # species_params_cols <- c("species","z0","alpha","erepro")
    # if (!all(species_params_cols %in% names(object@species_params))) {
    #     msg <- "species_params data.frame must have 'species', 'z0', 'alpha' and 'erepro' columms"
    #     errors <- c(errors,msg)
    # }
    # must also have SRR params but not sorted out yet
    
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
#'   each species at size. Changed with \code{\link{setMaxIntakeRate}}.
#' @slot search_vol An array (species x size) that holds the search volume for
#'   each species at size. Changed with \code{\link{setSearchVolume}}.
#' @slot metab An array (species x size) that holds the metabolism
#'   for each species at size. Changed with \code{\link{setMetabolicRate}}.
#' @slot mu_b An array (species x size) that holds the external mortality rate
#'   \eqn{\mu_{b.i}(w)}. Changed with \code{\link{setExtMortality}}.
#' @slot pred_kernel An array (species x predator size x prey size) that holds
#'   the predation coefficient of each predator at size on each prey size. If
#'   this is NA then the following two slots will be used. Changed with 
#'   \code{\link{setPredationKernel}}.
#' @slot ft_pred_kernel_e An array (species x log of predator/prey size ratio)
#'   that holds the Fourier transform of the feeding kernel in a form
#'   appropriate for evaluating the encounter rate integral. If this is NA
#'   then the \code{pred_kernel} will be used to calculate the available 
#'   energy integral. Changed with \code{\link{setPredationKernel}}.
#' @slot ft_pred_kernel_p An array (species x log of predator/prey size ratio)
#'   that holds the Fourier transform of the feeding kernel in a form
#'   appropriate for evaluating the predation mortality integral. If this is NA
#'   then the \code{pred_kernel} will be used to calculate the integral.
#'   Changed with \code{\link{setPredationKernel}}.
#' @slot rr_pp A vector the same length as the w_full slot. The size specific
#'   growth rate of the plankton spectrum. Changed with \code{\link{setPlankton}}.
#' @slot cc_pp A vector the same length as the w_full slot. The size specific
#'   carrying capacity of the plankton spectrum. Changed with 
#'   \code{\link{setPlankton}}.
#' @slot plankton_dynamics Name of the function for projecting the plankton abundance
#'   density by one timestep. The default is 
#'   \code{\link{plankton_semichemostat}}. 
#'   Changed with \code{\link{setPlankton}}.
#' @slot other_dynamics A named list of functions for projecting the
#'   values of other dynamical components of the ecosystem that may be modelled
#'   by a mizer extensions you have installed. The names of the list entries
#'   are the names of those components.
#' @slot other_encounter A named list of functions for calculating the 
#'   contribution to the encounter rate from each other dynamical component.
#' @slot other_mort A named list of functions for calculating the 
#'   contribution to the mortality rate from each other dynamical components.
#' @slot other_params A list containing the parameters needed by any mizer
#'   extensions you may have installed to model other dynamical components of
#'   the ecosystem.
#' @slot rates_funcs A named list with the names of the functions that should be
#'   used to calculate the rates needed by `project()`. By default this will be
#'   set to the names of the built-in rate functions.
#' @slot sc The community abundance of the scaling community
#' @slot species_params A data.frame to hold the species specific parameters.
#'   See \code{\link{newMultispeciesParams}} for details.
#' @slot interaction The species specific interaction matrix, \eqn{\theta_{ij}}.
#'   Changed with \code{\link{setInteraction}}.
#' @slot selectivity An array (gear x species x w) that holds the selectivity of
#'   each gear for species and size, \eqn{S_{g,i,w}}. Changed with 
#'   \code{\link{setFishing}}.
#' @slot catchability An array (gear x species) that holds the catchability of
#'   each species by each gear, \eqn{Q_{g,i}}. Changed with 
#'   \code{\link{setFishing}}.
#' @slot initial_effort A vector containing the initial fishing effort for each
#'   gear. Changed with \code{\link{setFishing}}.
#' @slot initial_n An array (species x size) that holds the initial abundance of
#'   each species at each weight.
#' @slot initial_n_pp A vector the same length as the w_full slot that describes
#'   the initial plankton abundance at each weight.
#' @slot initial_n_other A list with the initial abundances of all other
#'   ecosystem components. Has length zero if there are no other components.
#' @slot plankton_params Parameters for plankton. See \code{\link{setPlankton}}.
#' @slot A Abundance multipliers.
#' @slot linecolour A named vector of colour values, named by species.
#'   Used to give consistent colours in plots.
#' @slot linetype A named vector of linetypes, named by species. 
#'   Used to give consistent line types in plots.

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
#'   \code{\link{emptyParams}} \code{\link{newMultispeciesParams}}
#'   \code{\link{newCommunityParams}}
#'   \code{\link{newTraitParams}}
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
        metab = "array",
        pred_kernel = "array",
        ft_pred_kernel_e = "array",
        ft_pred_kernel_p = "array",
        mu_b = "array",
        rr_pp = "numeric",
        cc_pp = "numeric",
        plankton_dynamics = "character",
        plankton_params = "list",
        other_dynamics = "list",
        other_params = "list",
        other_encounter = "list",
        other_mort = "list",
        rates_funcs = "list",
        sc = "numeric",
        initial_n_pp = "numeric",
        initial_n_other = "list",
        species_params = "data.frame",
        interaction = "array",
        selectivity = "array",
        catchability = "array",
        initial_effort = "numeric",
        A = "numeric",
        linecolour = "character",
        linetype = "character"
    ),
)

setValidity("MizerParams", validMizerParams)
remove(validMizerParams)


#' Create empty MizerParams object of the right size
#' 
#' An internal function.
#' Sets up a valid \linkS4class{MizerParams} object with all the slots
#' initialised and given dimension names, but with some slots left empty. This
#' function is to be used by other functions to set up full parameter objects.
#' 
#' @section Size grid:
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
#' @section Changes to species params:
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
#'   set to the largest \code{w_inf} specified in the \code{species_params} data
#'   frame.
#' @param min_w_pp The smallest size of the plankton spectrum.
# #'   Ignored if w_full is specified.
#' 
#' @return An empty but valid MizerParams object
#' @seealso See \code{\link{newMultispeciesParams}} for a function that fills
#'   the slots left empty by this function.
#' @md
#' @export
emptyParams <- function(species_params,
                        no_w = 100,
                        min_w = 0.001,
                        # w_full = NA,
                        max_w = NA,
                        min_w_pp = NA) {
    assert_that(no_w > 10)
    
    ## Set defaults ----
    species_params <- set_species_param_default(species_params, "w_min", min_w)
    min_w <- min(species_params$w_min)
    
    species_params <- validSpeciesParams(species_params)
    
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
    no_sp <- nrow(species_params)
    species_names <- as.character(species_params$species)
    gear_names <- unique(species_params$gear)
    mat1 <- array(NA, dim = c(no_sp, no_w), 
                  dimnames = list(sp = species_names, w = signif(w,3)))
    ft_pred_kernel <- array(NA, dim = c(no_sp, no_w_full),
                            dimnames = list(sp = species_names, k = 1:no_w_full))
    
    selectivity <- array(0, dim = c(length(gear_names), no_sp, no_w), 
                         dimnames = list(gear = gear_names, sp = species_names, 
                                         w = signif(w, 3)))
    catchability <- array(0, dim = c(length(gear_names), no_sp), 
                          dimnames = list(gear = gear_names, sp = species_names))
    initial_effort <- rep(0, length(gear_names))
    names(initial_effort) <- gear_names
    
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
    # cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", 
    #                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    # Random palette gemerated pm https://medialab.github.io/iwanthue/
    colour_palette <- c("#815f00",
                        "#6237e2",
                        "#8da600",
                        "#de53ff",
                        "#0e4300",
                        "#430079",
                        "#6caa72",
                        "#ee0053",
                        "#007957",
                        "#b42979",
                        "#142300",
                        "#a08dfb",
                        "#644500",
                        "#04004c",
                        "#b79955",
                        "#0060a8",
                        "#dc8852",
                        "#007ca9",
                        "#ab003c",
                        "#9796d9",
                        "#472c00",
                        "#b492b0",
                        "#140000",
                        "#dc8488",
                        "#005c67",
                        "#5c585a")
    # type_palette <- c("solid", "dashed", "dotdash", "longdash", 
    #                   "twodash")
    type_palette <- c("solid")
    
    if ("linecolour" %in% names(species_params)) {
        linecolour <- species_params$linecolour
        # If any NA's first fill them with unused colours
        linecolour[is.na(linecolour)] <- 
            setdiff(colour_palette, linecolour)[1:sum(is.na(linecolour))]
        # if there are still NAs, start from beginning of palette again
        linecolour[is.na(linecolour)] <- 
            colour_palette[1:sum(is.na(linecolour))]
    } else {
        linecolour <- rep(colour_palette, length.out = no_sp)
    }
    names(linecolour) <- as.character(species_names)
    linecolour <- c(linecolour, "Total" = "black", "Plankton" = "green",
                    "Background" = "grey", "Fishing" = "red")
    
    if ("linetype" %in% names(species_params)) {
        linetype <- species_params$linetype
        linetype[is.na(linetype)] <- "solid"
    } else {
        linetype <- rep(type_palette, length.out = no_sp)
    }
    names(linetype) <- as.character(species_names)
    linetype <- c(linetype, "Total" = "solid", "Plankton" = "solid",
                  "Background" = "solid", "Fishing" = "solid")
    
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
        metab = mat1,
        mu_b = mat1,
        ft_pred_kernel_e = ft_pred_kernel,
        ft_pred_kernel_p = ft_pred_kernel,
        pred_kernel = array(),
        selectivity = selectivity,
        catchability = catchability,
        initial_effort = initial_effort,
        rr_pp = vec1,
        cc_pp = vec1,
        sc = w,
        initial_n_pp = vec1,
        species_params = species_params,
        interaction = interaction,
        other_dynamics = list(),
        other_encounter = list(),
        other_mort = list(),
        rates_funcs = list(
            Rates = "mizerRates",
            Encounter = "mizerEncounter",
            FeedingLevel = "mizerFeedingLevel",
            EReproAndGrowth = "mizerEReproAndGrowth",
            PredRate = "mizerPredRate",
            PredMort = "mizerPredMort",
            FMort = "mizerFMort",
            Mort = "mizerMort",
            ERepro = "mizerERepro",
            EGrowth = "mizerEGrowth",
            PlanktonMort = "mizerPlanktonMort",
            RDI = "mizerRDI",
            RDD = "BevertonHoltRDD"),
        plankton_dynamics = "plankton_semichemostat",
        other_params = list(),
        initial_n_other = list(),
        A = as.numeric(rep(NA, no_sp)),
        linecolour = linecolour,
        linetype = linetype
    )
    
    return(params)
}

#' Set up parameters for a general multispecies model
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
#' @inheritParams setParams
#' @param n The allometric growth exponent. This can be overruled for individual
#'   species by including a \code{n} column in the \code{species_params}. 
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
#' @inheritSection emptyParams Changes to species params
#' @inheritSection emptyParams Size grid
#' @inheritSection setParams Units in mizer
#' @inheritSection setInteraction Setting interactions
#' @inheritSection setPredationKernel Setting predation kernel
#' @inheritSection setSearchVolume Setting search volume
#' @inheritSection setMaxIntakeRate Setting maximum intake rate
#' @inheritSection setMetabolicRate Setting metabolic rate
#' @inheritSection setExtMortality Setting external mortality rate
#' @inheritSection setReproduction Setting reproduction
#' @inheritSection setFishing Setting fishing
#' @inheritSection setPlankton Setting plankton dynamics
#'   
#' @seealso \code{\link{project}}, \linkS4class{MizerSim},
#'   \code{\link{newCommunityParams}}, \code{\link{newTraitParams}}
#' @export
#' @family functions for setting up models
#' @examples
#' \dontrun{
#' params <- newMultispeciesParams(NS_species_params_gears, inter)
#' }
newMultispeciesParams <- function(
    species_params,
    interaction = NULL,
    no_w = 100,
    min_w = 0.001,
    max_w = NA,
    min_w_pp = NA,
    # setPredationKernel()
    pred_kernel = NULL,
    # setSearchVolume()
    search_vol = NULL,
    # setMaxIntakeRate()
    intake_max = NULL,
    # setMetabolicRate()
    metab = NULL,
    p = 0.7,
    # setExtMortality
    z0 = NULL,
    z0pre = 0.6,
    z0exp = n - 1,
    # setReproduction
    maturity = NULL,
    repro_prop = NULL,
    RDD = "BevertonHoltRDD",
    # setPlankton
    rate = NULL,
    capacity = NULL,
    n = 2 / 3,
    r_pp = 10,
    kappa = 1e11,
    lambda = 2.05,
    w_pp_cutoff = 10,
    plankton_dynamics = "plankton_semichemostat",
    # setFishing
    initial_effort = NULL) {
    no_sp <- nrow(species_params)
    
    ## For backwards compatibility, allow r_max instead of R_max
    if (!("R_max" %in% names(species_params)) &&
        "r_max" %in% names(species_params)) {
        names(species_params)[names(species_params) == "r_max"] <- "R_max"
    }
    
    ## Create MizerParams object ----
    params <- emptyParams(species_params,
                          no_w = no_w, 
                          min_w = min_w,  
                          max_w = max_w, 
                          min_w_pp = min_w_pp)
    
    ## Fill the slots ----
    params <- params %>% 
        set_species_param_default("n", n) %>% 
        set_species_param_default("p", p)
    params <- set_species_param_default(params, "q", 
                                        lambda - 2 + params@species_params$n)
    if (is.null(interaction)) {
        interaction <- matrix(1, nrow = no_sp, ncol = no_sp)
    }
    params <-
        setParams(params,
                  # setInteraction
                  interaction = interaction,
                  # setPredationKernel()
                  pred_kernel = pred_kernel,
                  # setSearchVolume()
                  search_vol = search_vol,
                  # setMaxIntakeRate()
                  intake_max = intake_max,
                  # setMetabolicRate()
                  metab = metab,
                  # setExtMortality
                  z0 = z0,
                  z0pre = z0pre,
                  z0exp = z0exp,
                  # setReproduction
                  maturity = maturity,
                  repro_prop = repro_prop,
                  RDD = RDD,
                  # setPlankton
                  rate = rate,
                  capacity = capacity,
                  r_pp = r_pp,
                  kappa = kappa,
                  lambda = lambda,
                  n = n,
                  w_pp_cutoff = w_pp_cutoff,
                  plankton_dynamics = plankton_dynamics,
                  # setFishing
                  initial_effort = initial_effort)
    
    params@initial_n <- get_initial_n(params)
    params@initial_n_pp <- params@cc_pp
    params@A <- rep(1, nrow(species_params))
    
    return(params)
}

#' Set or change any model parameters
#' 
#' This is a convenient wrapper function calling each of the following
#' functions
#' \itemize{
#' \item \code{\link{setPredationKernel}}
#' \item \code{\link{setSearchVolume}}
#' \item \code{\link{setInteraction}}
#' \item \code{\link{setMaxIntakeRate}}
#' \item \code{\link{setMetabolicRate}}
#' \item \code{\link{setExtMortality}}
#' \item \code{\link{setReproduction}}
#' \item \code{\link{setFishing}}
#' \item \code{\link{setPlankton}}
#' }
#' See the Details section below for a discussion of how to use this function.
#' 
#' @param params A \linkS4class{MizerParams} object
#' @inheritParams setInteraction
#' @inheritParams setPredationKernel
#' @inheritParams setSearchVolume
#' @inheritParams setMaxIntakeRate
#' @inheritParams setMetabolicRate
#' @inheritParams setExtMortality
#' @inheritParams setReproduction
#' @inheritParams setFishing
#' @inheritParams setPlankton
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
#' the values of all the model functions that you do not specify explicitly when
#' calling this function, unless you have protected the corresponding slots with
#' a comment. If you have changed any of the model functions in the
#' `params` object previously and now want to make changes to a different slot,
#' you will want to call the appropriate change function individually. So in the
#' above example you would have used `params <- setSearchVolume(params)`
#' instead of `params <- setParams(params)`. 
#' 
#' If you have added a comment to a slot of the params object, then setParams()
#' and its subfunctions will not recalculate the value for that slot from the
#' species parameters. For example after 
#' ```
#' comment(params@search_vol) <- "This should not change"
#' params@species_params$gamma <- 10
#' params <- setParams(params)
#' ```
#' will just issue a warning "The search volume has been commented and therefore
#' will not be recalculated from the species parameters". You can remove the
#' comment, and therefore allow recalculation of the slot, with
#' `comment(params@search_vol) <- NULL`.
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
#' @inheritSection setPredationKernel Setting predation kernel
#' @inheritSection setSearchVolume Setting search volume
#' @inheritSection setMaxIntakeRate Setting maximum intake rate
#' @inheritSection setMetabolicRate Setting metabolic rate
#' @inheritSection setExtMortality Setting external mortality rate
#' @inheritSection setReproduction Setting reproduction
#' @inheritSection setFishing Setting fishing
#' @inheritSection setPlankton Setting plankton dynamics
#' @md
#' @export
#' @family functions for setting parameters
#' @examples
#' \dontrun{
#' params <- newTraitParams()
#' params@species_params$gamma[3] <- 1000
#' params <- setParams(params)
#' }
setParams <- function(params,
                      # setPlankton
                      rate = NULL,
                      capacity = NULL,
                      r_pp = params@plankton_params[["r_pp"]],
                      kappa = params@plankton_params[["kappa"]],
                      lambda = params@plankton_params[["lambda"]],
                      n = params@plankton_params[["n"]],
                      w_pp_cutoff = params@plankton_params[["w_pp_cutoff"]],
                      plankton_dynamics = NULL,
                      # setInteraction()
                      interaction = NULL,
                      # setPredationKernel()
                      pred_kernel = NULL,
                      # setSearchVolume()
                      search_vol = NULL,
                      # setMaxIntakeRate()
                      intake_max = NULL,
                      # setMetabolicRate()
                      metab = NULL,
                      p = NULL,
                      # setExtMortality
                      z0 = NULL,
                      z0pre = 0.6,
                      z0exp = n - 1,
                      # setReproduction
                      maturity = NULL,
                      repro_prop = NULL,
                      RDD = params@rates_funcs$RDD,
                      # setFishing
                      initial_effort = NULL) {
    validObject(params)
    params <- setPlankton(params,
                          rate = rate,
                          capacity = capacity,
                          r_pp = r_pp,
                          kappa = kappa,
                          lambda = lambda,
                          n = n,
                          w_pp_cutoff = w_pp_cutoff,
                          plankton_dynamics = plankton_dynamics)
    params <- setInteraction(params,
                             interaction = interaction)
    params <- setPredationKernel(params,
                            pred_kernel = pred_kernel)
    params <- setMaxIntakeRate(params,
                           intake_max = intake_max)
    params <- setMetabolicRate(params,
                       metab = metab, p = p)
    params <- setExtMortality(params,
                       z0 = z0,
                       z0pre = z0pre,
                       z0exp = z0exp)
    # setSearchVolume() should be called only after 
    # setMaxIntakeRate() and setPredationKernel()
    params <- setSearchVolume(params,
                              search_vol = search_vol)
    params <- setReproduction(params,
                              maturity = maturity,
                              repro_prop = repro_prop,
                              RDD = RDD)
    params <- setFishing(params, initial_effort = initial_effort)
    return(params)
}


#' Set species interaction matrix
#'
#' @section Setting interactions:
#'
#'   The species interaction matrix \eqn{\theta_{ij}}, is used when calculating
#'   the food encounter rate in \code{\link{getEncounter}} and the predation
#'   mortality rate in \code{\link{getPredMort}}. Its entries are dimensionless
#'   numbers between 0 and 1 that characterise the strength at which predator
#'   species \eqn{i} predates on prey species \eqn{j}.
#'
#'   This function checks that the supplied interaction matrix is valid and then
#'   stores it in the \code{interaction} slot of the params object before
#'   returning that object.
#'
#'   The order of the columns and rows of the \code{interaction} argument should
#'   be the same as the order in the species params dataframe in the
#'   \code{params} object. If you supply a named array then the function will
#'   check the order and warn if it is different.
#'
#'   The interaction of the species with the plankton are set via a column
#'   \code{interaction_p} in the \code{species_params} data frame. Again the
#'   entries have to be numbers between 0 and 1. By default this column is set
#'   to all 1s.
#'
#' @param params MizerParams object
#' @param interaction Optional interaction matrix of the species (predator
#'   species x prey species). Entries should be numbers between 0 and 1. By
#'   default all entries are 1. See "Setting interactions" section below.
#'
#' @return MizerParams object with updated interaction matrix. Because of the
#'   way the R language works, `setInteraction()` does not make the changes to
#'   the params object that you pass to it but instead returns a new params
#'   object. So to affect the change you call the function in the form
#'   `params <- setInteraction(params, ...)`.
#' @export
#' @md
#' @family functions for setting parameters
#' @examples
#' \dontrun{
#' params <- newTraitParams()
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
        stop("Values in the plankton interaction vector should be between 0 and 1")
    }
    params@species_params$interaction_p <- species_params$interaction_p
    
    return(params)
}

#' @rdname setInteraction
#' @export
getInteraction <- function(params) {
    params@interaction
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
#' An alternative pred_kernel type is "box", implemented by the function
#' \code{\link{box_pred_kernel}}, and "power_law", implemented by the function
#' \code{\link{power_law_pred_kernel}}. These functions require certain species
#' parameters in the species_params data frame. For the lognormal kernel these
#' are \code{beta} and \code{sigma}, for the box kernel they are \code{ppmr_min}
#' and \code{ppmr_max}. They are explained in the help pages for the kernel
#' functions. Except for \code{beta} and \code{sigma}, no defaults are set for
#' these parameters. If they are missing from the species_params data frame then
#' mizer will issue an error message.
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
#' \code{\link{getPredationKernel}} function.
#' 
#' \strong{Kernel dependent on both predator and prey size}
#' 
#' If you want to work with a feeding kernel that depends on predator mass and
#' prey mass independently, you can specify the full feeding kernel as a
#' three-dimensional array (predator species x predator size x prey size).
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
#' @return A MizerParams object with updated predation kernel. Because of the
#'   way the R language works, `setPredationKernel()` does not make the changes
#'   to the params object that you pass to it but instead returns a new params
#'   object. So to affect the change you call the function in the form
#'   `params <- setPredationKernel(params, ...)`.
#' @export
#' @family functions for setting parameters
#' @examples
#' \dontrun{
#' ## Set up a MizerParams object
#' params <- newMultispeciesParams(NS_species_params_gears, inter)
#' 
#' ## If you change predation kernel parameters after setting up a model, you
#' # need to call setPredationKernel
#' params@species_params["Cod", "beta"] <- 200
#' params <- setPredationKernel(params)
#' 
#' ## You can change to a different predation kernel type
#' params@species_params$pred_kernel_type <- "box"
#' params@species_params$ppmr_min <- 2
#' params@species_params$ppmr_max <- 4
#' params <- setPredationKernel(params)
#' 
#' ## If you need a kernel that depends also on prey size you need to define
#' # it yourself.
#' pred_kernel <- getPredationKernel(params)
#' pred_kernel["Herring", , ] <- sweep(pred_kernel["Herring", , ], 2, 
#'                                     params@w_full, "*")
#' params<- setPredationKernel(params, pred_kernel = pred_kernel)
#' }
setPredationKernel <- function(params,
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
    
    if (!is.null(comment(params@pred_kernel))) {
        message("The predation kernel has been commented and therefore will ",
                "not be recalculated from the species parameters.")
        return(params)
    }
    
    ## Set a pred kernel dependent on predator/prey size ratio only
    
    # If pred_kernel_type is not supplied use "lognormal"
    params <- default_pred_kernel_params(params)
    
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
        ri <- min(max(which(phi > 0)), no_w_full - 1)  # index of largest ppmr
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
#' For more detail about the predation kernel see \code{\link{setPredationKernel}}.
#' 
#' @param params A MizerParams object
#' @return An array (predator species x predator_size x prey_size)
#' @export
getPredationKernel <- function(params) {
    assert_that(is(params, "MizerParams"))
    if (length(dim(params@pred_kernel)) > 1) {
        return(params@pred_kernel)
    }
    species_params <- default_pred_kernel_params(params@species_params)
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
#' If the \code{search_vol} argument is not supplied, then the search volume is 
#' set to
#' \deqn{\gamma_i(w) = \gamma_i w^q_i.} 
#' The values of \eqn{\gamma_i} (the search volume at 1g) and \eqn{q_i} (the
#' allometric exponent of the search volume) are taken from the \code{gamma} and
#' \code{q} columns in the species parameter dataframe. If the \code{gamma}
#' column is not supplied in the species parameter dataframe, a default is
#' calculated by the \code{\link{get_gamma_default}} function. Note that only
#' for predators of size \eqn{w = 1} gram is the value of the species parameter
#' \eqn{\gamma_i} the same as the value of the search volume \eqn{\gamma_i(w)}.
#' 
#' @param params MizerParams
#' @param search_vol Optional. An array (species x size) holding the search volume
#'   for each species at size. If not supplied, a default is set as described in
#'   the section "Setting search volume". 
#' 
#' @return MizerParams with updated search volume. Because of the way the R
#'   language works, `setSearchVolume()` does not make the changes to the params
#'   object that you pass to it but instead returns a new params object. So to
#'   affect the change you call the function in the form
#'   `params <- setSearchVolume(params, ...)`.
#' @export
#' @family functions for setting parameters
#' @examples
#' \dontrun{
#' params <- newTraitParams()
#' params@species_params$gamma[3] <- 1000
#' params <- setSearchVolume(params)
#' }
setSearchVolume <- function(params, 
                            search_vol = NULL) {
    assert_that(is(params, "MizerParams"))
    species_params <- params@species_params
    # If search_vol array is supplied, check it, store it and return
    if (!is.null(search_vol)) {
        assert_that(is.array(search_vol))
        assert_that(identical(dim(search_vol), dim(params@search_vol)))
        if (!is.null(dimnames(search_vol)) && 
            !all(dimnames(search_vol)[[1]] == species_params$species)) {
            stop(paste0("You need to use the same ordering of species in the search_vol array as in the ",
                        "params object: ", toString(species_params$species)))
        }
        assert_that(all(search_vol >= 0))
        params@search_vol[] <- search_vol
        comment(params@search_vol) <- comment(search_vol)
        return(params)
    }
    
    if (!is.null(comment(params@search_vol))) {
        message("The search volume has been commented and therefore will ",
                "not be recalculated from the species parameters.")
        return(params)
    }
    
    # Calculate default for any missing gammas
    params@species_params$gamma <- get_gamma_default(params)
    
    params@search_vol[] <- 
        sweep(outer(params@species_params[["q"]], params@w,
                    function(x, y) y ^ x),
              1, params@species_params$gamma, "*")
    
    return(params)
}

#' @rdname setSearchVolume
#' @export
getSearchVolume <- function(params) {
    params@search_vol
}


#' Set maximum intake rate
#'
#' @section Setting maximum intake rate:
#' The maximum intake rate \eqn{h_i(w)} of an individual of species \eqn{i} and
#' weight \eqn{w} determines the feeding level, calculated with
#' \code{\link{getFeedingLevel}}. It is measured in grams/year.
#'
#' If the \code{intake_max} argument is not supplied, then the maximum intake
#' rate is set to \deqn{h_i(w) = h_i w^n_i.} 
#' The values of \eqn{h_i} (the maximum intake rate of an individual of size 1
#' gram) and \eqn{n_i} (the allometric exponent for the intake rate) are taken
#' from the \code{h} and \code{n} columns in the species parameter dataframe. If
#' the \code{h} column is not supplied in the species parameter dataframe, it is
#' calculated by the \code{\link{get_h_default}} function, using \code{f0} and
#' the \code{k_vb} column, if they are supplied.
#' 
#' If \eqn{h_i} is set to \code{Inf}, fish will consume all encountered food.
#'
#' @param params MizerParams
#' @param intake_max Optional. An array (species x size) holding the maximum
#'   intake rate for each species at size. If not supplied, a default is set as
#'   described in the section "Setting maximum intake rate".
#' @return A \code{MizerParams} object with updated maximum intake rate. Because
#'   of the way the R language works, `setMaxIntakeRate()` does not make the
#'   changes to the params object that you pass to it but instead returns a new
#'   params object. So to affect the change you call the function in the form
#'   `params <- setMaxIntakeRate(params, ...)`.
#' @export
#' @family functions for setting parameters
#' @examples
#' \dontrun{
#' params <- NS_params
#' params@species_params$h[3] <- 35
#' params <- setMaxIntakeRate(params)
#' }
setMaxIntakeRate <- function(params, 
                         intake_max = NULL) {
    assert_that(is(params, "MizerParams"))
    species_params <- params@species_params
    
    # If intake_max array is supplied, check it, store it and return
    if (!is.null(intake_max)) {
        assert_that(is.array(intake_max),
                    identical(dim(intake_max), dim(params@intake_max)))
        if (!is.null(dimnames(intake_max)) && 
            !all(dimnames(intake_max)[[1]] == species_params$species)) {
            stop(paste0("You need to use the same ordering of species in the intake_max array as in the ",
                        "params object: ", toString(species_params$species)))
        }
        assert_that(all(intake_max >= 0))
        params@intake_max[] <- intake_max
        comment(params@intake_max) <- comment(intake_max)
        return(params)
    }
    
    if (!is.null(comment(params@intake_max))) {
        message("The max intake rate has been commented and therefore will ",
                "not be recalculated from the species parameters.")
        return(params)
    }
    
    params@species_params$h <- get_h_default(params)
    
    params@intake_max[] <- sweep(outer(params@species_params[["n"]], 
                                       params@w, function(x, y) y^x),
                                 1, params@species_params[["h"]], "*") 
    
    return(params)
}

#' @rdname setMaxIntakeRate
#' @export
getMaxIntakeRate <- function(params) {
    params@intake_max
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
#' If the \code{metab} argument is not supplied, then for each species the
#' metabolic rate \eqn{k(w)} for an individual of size \eqn{w} is set to
#' \deqn{k(w) = ks w^p + k w,}
#' where \eqn{ks w^p} represents the rate of standard metabolism and \eqn{k w}
#' is the rate at which energy is expended on activity and movement. The values
#' of \eqn{ks}, \eqn{p} and \eqn{k} are taken from the \code{ks}, \code{p} and
#' \code{k} columns in the species parameter dataframe. If any of these
#' parameters are not supplied, the defaults are \eqn{k = 0}, \eqn{p = n} and
#' \deqn{ks = f_c h \alpha w_{mat}^{n-p},}{ks = f_c * h * alpha * w_mat^(n - p),}
#' where \eqn{f_c} is the critical feeding level taken from the \code{fc} column
#' in the species parameter data frame. If the critical feeding level is not
#' specified, a default of \eqn{f_c = 0.2} is used.
#' 
#' @param params MizerParams
#' @param metab Optional. An array (species x size) holding the metabolic rate
#'   for each species at size. If not supplied, a default is set as described in
#'   the section "Setting metabolic rate".
#' @param p The allometric metabolic exponent. This is only used if \code{metab}
#'   is not given explicitly and if the exponent is not specified in a \code{p}
#'   column in the \code{species_params}.
#' 
#' @return MizerParams object with updated metabolic rate. Because of the way
#'   the R language works, `setMetabolicRate()` does not make the changes to the
#'   params object that you pass to it but instead returns a new params object.
#'   So to affect the change you call the function in the form
#'   `params <- setMetabolicRate(params, ...)`.
#' @export
#' @family functions for setting parameters
#' @examples
#' \dontrun{
#' params <- NS_params
#' # Change activity coefficient for species 3
#' params@species_params$k[3] <- 8
#' params <- setMetabolicRate(params)
#' }
setMetabolicRate <- function(params, 
                     metab = NULL, p = NULL) {
    assert_that(is(params, "MizerParams"))
    if (!is.null(p)) {
        assert_that(is.numeric(p))
        params <- set_species_param_default(params, "p", p)
    }
    species_params <- params@species_params
    if (!is.null(metab)) {
        assert_that(is.array(metab),
                    identical(dim(metab), dim(params@metab)))
        if (!is.null(dimnames(metab)) && 
            !all(dimnames(metab)[[1]] == species_params$species)) {
            stop(paste0("You need to use the same ordering of species as in the ",
                        "params object: ", toString(species_params$species)))
        }
        assert_that(all(metab >= 0))
        params@metab[] <- metab
        comment(params@metab) <- comment(metab)
        return(params)
    }
    
    if (!is.null(comment(params@metab))) {
        message("The metabolic rate has been commented and therefore will ",
                "not be recalculated from the species parameters.")
        return(params)
    }
    
    params <- set_species_param_default(params, "k", 0)
    params@species_params$ks <- get_ks_default(params)
    params@metab[] <- 
        sweep(outer(params@species_params$p, params@w,
                    function(x, y) y ^ x),
              1, params@species_params$ks, "*") +
        outer(params@species_params$k, params@w)
    return(params)
}

#' @rdname setMetabolicRate
#' @export
getMetabolicRate <- function(params) {
    params@metab
}

#' Set external mortality rate
#' 
#' @section Setting external mortality rate:
#' The external mortality is all the mortality that is not due to fishing or
#' predation by predators included in the model. The external mortality could be
#' due to predation by predators that are not explicitly included in the model
#' (e.g. mammals or seabirds) or due to other causes like illness. It is a rate
#' with units 1/year.
#' 
#' The \code{z0} argument allows you to specify an external mortality rate
#' that depends on species and body size. You can see an example of this in
#' the Examples section of the help page for \code{\link{setExtMortality}}.
#' 
#' If the \code{z0} argument is not supplied, then the external mortality
#' is assumed to depend only on the species, not on the
#' size of the individual: \eqn{\mu_{b.i}(w) = z_{0.i}}. The value of the
#' constant \eqn{z_0} for each species is taken from the \code{z0} column of the
#' species_params data frame, if that column exists. Otherwise it is calculated
#' as 
#' \deqn{z_{0.i} = {\tt z0pre}_i\, w_{inf}^{\tt z0exp}.}{z_{0.i} = z0pre_i w_{inf}^{z0exp}.}
#' 
#' @param params MizerParams
#' @param z0 Optional. An array (species x size) holding the external
#'   mortality rate.
#' @param z0pre If \code{z0}, the mortality from other sources, is not a column
#'   in the species data frame, it is calculated as z0pre * w_inf ^ z0exp.
#'   Default value is 0.6.
#' @param z0exp If \code{z0}, the mortality from other sources, is not a column
#'   in the species data frame, it is calculated as \code{z0pre * w_inf ^ z0exp}.
#'   Default value is \code{n-1}.
#' 
#' @return MizerParams object with updated external mortality rate. Because of
#'   the way the R language works, `setExtMortality()` does not make the changes
#'   to the params object that you pass to it but instead returns a new params
#'   object. So to affect the change you call the function in the form
#'   `params <- setExtMortality(params, ...)`.
#' @export
#' @family functions for setting parameters
#' @examples
#' \dontrun{
#' params <- newMultispeciesParams(NS_species_params)
#'
#' #### Setting allometric death rate #######################
#' 
#' # Set coefficient for each species. Here we choose 0.1 for each species
#' z0pre <- rep(0.1, nrow(params@species_params))
#' 
#' # Multiply by power of size with exponent, here chosen to be -1/4
#' # The outer() function makes it an array species x size
#' z0 <- outer(z0pre, params@w^(-1/4))
#' 
#' # Change the external mortality rate in the params object
#' params <- setExtMortality(params, z0 = z0)
#' }
setExtMortality <- function(params, z0 = NULL, z0pre = 0.6, z0exp = -1/4) {
    assert_that(is(params, "MizerParams"))
    if (!is.null(z0)) {
        assert_that(is.array(z0),
                    identical(dim(z0), dim(params@mu_b)))
        params@mu_b[] <- z0
        comment(params@mu_b) <- comment(z0)
        return(params)
    }
    
    if (!is.null(comment(params@mu_b))) {
        message("The external mortality rate has been commented and therefore will ",
                "not be recalculated from the species parameters.")
        return(params)
    }
    
    assert_that(is.number(z0pre), z0pre >= 0,
                is.number(z0exp))
    species_params <- params@species_params
    assert_that(noNA(species_params$w_inf))
    # Sort out z0 (external mortality)
    message <- ("Note: Using z0 = z0pre * w_inf ^ z0exp for missing z0 values.")
    params <- set_species_param_default(params, "z0",
                                        z0pre * species_params$w_inf^z0exp,
                                        message)
    params@mu_b[] <- params@species_params$z0
    return(params)
}

#' @rdname setExtMortality
#' @export
getExtMortality <- function(params) {
    params@mu_b
}


#' Set reproduction parameters
#' 
#' Sets the proportion of the total energy available for reproduction and growth
#' that is invested into reproduction as a function of the size of the
#' individual and sets the reproductive efficiency.
#' 
#' @section Setting reproduction:
#' 
#' \subsection{Investment}{
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
#' }
#' 
#' \subsection{Efficiency}{
#' The reproductive efficiency, i.e., the proportion of energy allocated to
#' reproduction that results in egg biomass, is set from the \code{erepro}
#' column in the species_params data frame. If that is not provided, the default
#' is set to 1 (which you will want to override). The offspring biomass divided
#' by the egg biomass gives the rate of egg production, returned by
#' \code{\link{getRDI}}.
#' }
#' 
#' \subsection{Density dependence}{
#' The stock-recruitment relationship is an emergent phenomenon in mizer, with
#' several sources of density dependence. Firstly, the amount of energy invested
#' into reproduction depends on the energy income of the spawners, which is
#' density-dependent due to competition for prey. Secondly, the proportion of
#' larvae that grow up to recruitment size depends on the larval mortality,
#' which depends on the density of predators, and on larval growth rate, which
#' depends on density of prey.
#' 
#' Finally, the proportion of eggs that are viable and hatch to larvae can be
#' density dependent. Somewhat misleadingly, mizer refers to this relationship
#' between the number of eggs and the number of hatched larvae as the
#' stock-recruitment relationship, even though it is only one part of the full
#' stock-recruitment relationship. However it is the only part that can be set
#' independently, while the other parts are already determined by the predation
#' parameters and other model parameters. Thus in practice this part of the
#' density dependence is used to encode all the density dependence that is not
#' already included in the other two sources of density dependence.
#' 
#' To calculate the density-dependent rate of larvae production, mizer puts the
#' the density-independent rate of egg production through a "stock-recruitment"
#' function. The result is returned by \code{\link{getRDD}}. The name of the
#' stock-recruitment function is specified by the \code{RDD} argument. The
#' default is the Beverton-Holt function \code{\link{BevertonHoltRDD}}, which
#' requires an \code{R_max} column in the species_params data frame giving the
#' maximum egg production rate. If this column does not exist, it is initialised
#' to \code{Inf}, leading to no density-dependence. Other functions provided by
#' mizer are \code{\link{RickerRDD}} and \code{\link{SheperdRDD}} and you can
#' easily use these as models for writing your own functions.
#' }
#' @param params A MizerParams object
#' @param maturity Optional. An array (species x size) that holds the proportion
#'   of individuals of each species at size that are mature. If not supplied, a
#'   default is set as described in the section "Setting reproduction".
#' @param repro_prop Optional. An array (species x size) that holds the
#'   proportion of consumed energy that a mature individual allocates to
#'   reproduction for each species at size. If not supplied, a default is set as
#'   described in the section "Setting reproduction".
#' @param RDD The name of the stock recruitment function. Defaults to 
#'   "\code{\link{BevertonHoltRDD}}".
#' 
#' @return The updated MizerParams object. Because of the way the R language
#'   works, `setReproduction()` does not make the changes to the params object
#'   that you pass to it but instead returns a new params object. So to affect
#'   the change you call the function in the form
#'   `params <- setReproduction(params, ...)`.
#' @export
#' @family functions for setting parameters
#' @examples
#' \dontrun{
#' params <- NS_params
#' # Change maturity size for species 3
#' params@species_params$w_mat[3] <- 24
#' params <- setReproduction(params)
#' }
setReproduction <- function(params, maturity = NULL, repro_prop = NULL,
                            RDD = params@rates_funcs$RDD) {
    assert_that(is(params, "MizerParams"),
                is.string(RDD),
                exists(RDD),
                is.function(get(RDD)))
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
        
    } else if (!is.null(comment(params@maturity))) {
        message("The maturity ogive has been commented and therefore will ",
                "not be recalculated from the species parameters.")
        maturity <- params@maturity
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
    comment(params@maturity) <- comment(maturity)
    
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
                       mn = params@species_params$m - params@species_params$n
                )
            ), dim = c(nrow(species_params), length(params@w)))
    }
    
    if (!is.null(comment(params@psi))) {
        message("The reproductive proportion has been commented and therefore will ",
                "not be recalculated from the species parameters.")
    } else {
        params@psi[] <- params@maturity * repro_prop
        comment(params@psi) <- comment(repro_prop)
    }
    
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
    
    # RDD function is currently called only with three arguments
    RDD_fn <- get(RDD)
    if (!all(names(formals(RDD)) %in%  c("rdi", "species_params", "t", "..."))) {
        stop("Arguments of RDD function can only contain 'rdi', 'species_params' and `t`.")
    }
    if (!all(c("rdi", "...") %in% names(formals(RDD)))) {
        stop("The RDD function needs to have at least arguments `rdi` and `...`.")
    }
    params@rates_funcs$RDD <- RDD
    if (identical(params@rates_funcs$RDD, "BevertonHoltRDD")) {
        
        # for legacy reasons (R_max used to be called r_max):
        if ("r_max" %in% names(params@species_params)) {
            params@species_params$R_max <- params@species_params$r_max
            params@species_params$r_max <- NULL
            message("The 'r_max' column has been renamed to 'R_max'.")
        }
        
        params <- set_species_param_default(params, "R_max", Inf)
    }
    
    return(params)
}

#' @rdname setReproduction
#' @export
getMaturityProportion <- function(params) {
    params@maturity
}

#' @rdname setReproduction
#' @export
getReproductionProportion <- function(params) {
    repro_prop <- params@psi / params@maturity
    repro_prop[is.nan(repro_prop)] <- 0
    comment(repro_prop) <- comment(params@psi)
    repro_prop
}

#' Set up plankton
#' 
#' Sets the intrinsic plankton growth rate and the intrinsic plankton carrying
#' capacity as well as the name of the function used to simulate the plankton
#' dynamics
#' 
#' @section Setting plankton dynamics:
#' By default, mizer uses a semichemostat model to describe the plankton
#' dynamics in each size class independently. This semichemostat dynamics is implemented
#' by the function \code{\link{plankton_semichemostat}}. You can change the
#' plankton dynamics by writing your own function, modelled on
#' \code{\link{plankton_semichemostat}}, and then passing the name of your
#' function in the \code{plankton_dynamics} argument.
#' 
#' The \code{rate} argument is a vector specifying the intrinsic plankton
#' growth rate for each size class. If it is not supplied, then the intrinsic growth
#' rate \eqn{r_p(w)} at size \eqn{w}
#' is set to \deqn{r_p(w) = r_p\, w^{n-1}.}{r_p(w) = r_p w^{n-1}}
#' The values of \eqn{r_p} and \eqn{n} are taken from the \code{r_pp}
#' and \code{n} arguments.
#' 
#' The \code{capacity} argument is a vector specifying the intrinsic plankton
#' carrying capacity for each size class. If it is not supplied, then the intrinsic carrying
#' capacity \eqn{c_p(w)} at size \eqn{w}
#' is set to \deqn{c_p(w) = \kappa\, w^{-\lambda}}{c_p(w) = \kappa w^{-\lambda}}
#' for all \eqn{w} less than \code{w_pp_cutoff} and zero for larger sizes.
#' The values of \eqn{\kappa} and \eqn{\lambda} are taken from the \code{kappa}
#' and \code{lambda} arguments.
#' 
#' @param params A MizerParams object
#' @param rate Optional. Vector of plankton intrinsic birth rates
#' @param capacity Optional. Vector of plankton intrinsic carrying capacity
#' @param r_pp Coefficient of the intrinsic plankton birth rate
#' @param n Allometric growth exponent for plankton
#' @param kappa Coefficient of the intrinsic plankton carrying capacity
#' @param lambda Scaling exponent of the intrinsic plankton carrying capacity
#' @param w_pp_cutoff The upper cut off size of the plankton spectrum. 
#'   Default is 10 g.
#' @param plankton_dynamics Function that determines plankton dynamics by
#'   calculating the plankton spectrum at the next time step from the current
#'   state.
#' 
#' @return A MizerParams object with updated plankton parameters. Because of the
#'   way the R language works, `setPlankton()` does not make the changes to the
#'   params object that you pass to it but instead returns a new params object.
#'   So to affect the change you call the function in the form
#'   `params <- setPlankton(params, ...)`.
#' @export
#' @family functions for setting parameters
setPlankton <- function(params,
                        rate = NULL,
                        capacity = NULL,
                        r_pp = params@plankton_params[["r_pp"]],
                        kappa = params@plankton_params[["kappa"]],
                        lambda = params@plankton_params[["lambda"]],
                        n = params@plankton_params[["n"]],
                        w_pp_cutoff = params@plankton_params[["w_pp_cutoff"]],
                        plankton_dynamics = NULL) {
    assert_that(is(params, "MizerParams"),
                is.number(kappa), kappa > 0,
                is.number(lambda),
                is.number(r_pp), r_pp > 0,
                is.number(w_pp_cutoff),
                is.number(n))
    params@plankton_params[["kappa"]] <- kappa
    params@plankton_params[["lambda"]] <- lambda
    params@plankton_params[["r_pp"]] <- r_pp
    params@plankton_params[["n"]] <- n
    params@plankton_params[["w_pp_cutoff"]] <- w_pp_cutoff
    # weight specific plankton growth rate
    if (!is.null(rate)) {
        assert_that(is.numeric(rate),
                    identical(length(rate), length(params@rr_pp)))
        params@rr_pp[] <- rate
        comment(params@rr_pp) <- comment(rate)
    } else {
        if (!is.null(comment(params@rr_pp))) {
            message("The plankton intrinsic growth rate has been commented and therefore will ",
                    "not be recalculated from the species parameters.")
        } else {
            params@rr_pp[] <- r_pp * params@w_full^(n - 1)
        }
    }
    # the plankton carrying capacity
    if (!is.null(capacity)) {
        assert_that(is.numeric(capacity),
                    identical(length(capacity), length(params@cc_pp)))
        params@cc_pp[] <- capacity
        comment(params@cc_pp) <- comment(capacity)
    } else {
        if (!is.null(comment(params@cc_pp))) {
            message("The plankton intrinsic growth rate has been commented and therefore will ",
                    "not be recalculated from the species parameters.")
        } else {
            params@cc_pp[] <- kappa*params@w_full^(-lambda)
            params@cc_pp[params@w_full > w_pp_cutoff] <- 0
        }
    }
    if (!is.null(plankton_dynamics)) {
        assert_that(is.character(plankton_dynamics))
        if (!is.function(get0(plankton_dynamics))) {
            stop("The function ", plankton_dynamics, "is not defined.")
        }
        params@plankton_dynamics <- plankton_dynamics
    }
    
    return(params)
}

#' @rdname setPlankton
#' @export
getPlanktonRate <- function(params) {
    params@rr_pp
}

#' @rdname setPlankton
#' @export
getPlanktonCapacity <- function(params) {
    params@cc_pp
}

#' @rdname setPlankton
#' @export
getPlanktonParams <- function(params) {
    params@plankton_params
}

#' @rdname setPlankton
#' @export
getPlanktonDynamics <- function(params) {
    params@plankton_dynamics
}


#' Set fishing parameters
#' 
#' @section Setting fishing:
#' 
#' \strong{Gears}
#' 
#' In `mizer`, fishing mortality is imposed on species by fishing gears. The
#' total fishing mortality is obtained by summing over the mortality from all
#' gears,
#' \deqn{\mu_{f.i}(w) = \sum_g F_{g,i}(w),}
#' where the fishing mortality \eqn{F_{g,i}(w)} imposed by gear \eqn{g} on
#' species \eqn{i} at size \eqn{w} is calculated as:
#' \deqn{F_{g,i}(w) = S_{g,i}(w) Q_{g,i} E_{g},}
#' where \eqn{S} is the selectivity by species, gear and size, \eqn{Q} is the 
#' catchability by species and gear and \eqn{E} is the fishing effort by gear.
#' 
#' At the moment a species can only be selected by one fishing gear, although 
#' each gear can select more than one species (this is a limitation with the 
#' current package that will be developed in future releases). The gear
#' selecting each species can be specified in the `gear` column in the
#' species_params data frame. If no gear is specified, the default gear is
#' "knife_edge_gear".
#' 
#' \strong{Selectivity}
#' 
#' The selectivity at size of each gear has a range between 0 (not selected at
#' that size) to 1 (fully selected at that size). It is given by a selectivity
#' function. The name of the selectivity function is given by the `sel_func`
#' column in the species parameters data frame. Some selectivity functions are
#' included in the package: `knife_edge()`, `sigmoid_length()`,
#' `double_sigmoid_length()`, and `sigmoid_weight()`. New functions can be defined 
#' by the user. Each
#' gear has the same selectivity function for all the species it selects, but
#' the parameter values for each species may be different, e.g. the lengths of
#' species that a gear selects may be different.
#' 
#' Each selectivity function has a range of parameters. Values for these
#' parameters must be included as columns in the species parameters data.frame.
#' The names of the columns must exactly match the names of the corresponding
#' arguments of the selectivity function. For example, the default selectivity
#' function is `knife_edge()` which has sudden change of selectivity from 0 to 1
#' at a certain size. In its help page you can see that the `knife_edge()`
#' function has arguments `w` and `knife_edge_size` The first argument, `w`, is
#' size (the function calculates selectivity at size). All selectivity functions
#' must have `w` as the first argument. The values for the other arguments must
#' be found in the species parameters data.frame. So for the `knife_edge()`
#' function there should be a `knife_edge_size` column. Because `knife_edge()`
#' is the default selectivity function, the `knife_edge_size` argument has a
#' default value = `w_mat`.
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
#' \strong{Effort}
#' 
#' The initial fishing effort is stored in the `MizerParams` object. If it is
#' not supplied, it is set to zero. The initial effort can be overruled when
#' the simulation is run with `project()`, where it is also possible to specify
#' an effort that varies through time.
#' 
#' @param params A MizerParams object
#' @param initial_effort Optional. A number or a named numeric vector specifying
#'   the fishing effort. If a number, the same effort is used for all gears. If
#'   a vector, must be named by gear.
#'   
#' @return MizerParams object with updated catchability and selectivity. Because
#'   of the way the R language works, `setFishing()` does not make the changes
#'   to the params object that you pass to it but instead returns a new params
#'   object. So to affect the change you call the function in the form
#'   `params <- setFishing(params, ...)`.
#' @export
#' @md
#' @family functions for setting parameters
#' @examples
#' \dontrun{
#' params <- NS_params
#' # Change knife edge size for species 1
#' params@species_params$knife_edge_size[1] <- 15
#' params <- setFishing(params)
#' }
setFishing <- function(params, initial_effort = NULL) {
    assert_that(is(params, "MizerParams"))
    species_params <- params@species_params
    no_sp <- nrow(species_params)
    
    # If no gear specified in species_params, then use `knife_edge_gear`
    species_params <- set_species_param_default(
        species_params, "gear", default = "knife_edge_gear")
    
    # If no sel_func column in species_params, set to 'knife_edge'
    species_params <- set_species_param_default(
        species_params, "sel_func", default = "knife_edge",
        message = "Note: Setting missing selectivity function to be 'knife_edge'.")

    # Provide default for knife_edge_size if needed
    if ("knife_edge" %in% species_params$sel_func) {
        species_params <- set_species_param_default(
            species_params, "knife_edge_size", 
            default = species_params$w_mat,
            message = "Note: Setting missing knife edge selectivity equal to w_mat.")
    }
    
    # If no catchability column in species_params, set to 1
    species_params <- set_species_param_default(species_params,
                                                "catchability", 1)
    
    if (!is.null(initial_effort)) {
        validate_effort_vector(params, initial_effort)
        params@initial_effort[] <- initial_effort
        comment(params@initial_effort) <- comment(initial_effort)
    }
    
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
        if (!is.null(comment(params@selectivity))) {
            message("The selectivity has been commented and therefore will ",
                    "not be recalculated from the species parameters.")
        } else {
            params@selectivity[as.character(species_params[g, 'gear']), g, ] <- sel
        }
        # Now do catchability
        if (!is.null(comment(params@catchability))) {
            message("The catchability has been commented and therefore will ",
                    "not be recalculated from the species parameters.")
        } else {
            params@catchability[as.character(species_params[g,'gear']), g] <- 
                species_params[g, "catchability"]
        }
    }
    params@species_params <- species_params
    return(params)
}

#' @rdname setFishing
#' @export
getCatchability <- function(params) {
    params@catchability
}

#' @rdname setFishing
#' @export
getSelectivity <- function(params) {
    params@selectivity
}


#' Set initial values
#'
#' @param params A \code{\link{MizerParams}} object
#' @param sim A \code{MizerSim} object. If supplied, the `initial_n`, 
#'   `initial_n_pp` and `initial_n_pp` arguments are ignored and the information
#'   is taken from the last timestep of the simulation in `sim`.
#' @param initial_n The initial abundances of species. A matrix with dimensions
#'   species x size. The order of species must be the same as in the MizerParams
#'   argument. Optional. Ignored if `sim` is supplied.
#' @param initial_n_pp The initial abundances of plankton. A numeric vector.
#'   Optional. Ignored if `sim` is supplied.
#' @param initial_n_other The initial abundances of the other dynamic ecosystem
#'   components. A named list with one entry for each component. Optional.
#'   Ignored if `sim` is supplied.
#'   
#' @return A MizerParams object with updated initial values. Because of the
#'   way the R language works, `setInitialValues()` does not make the changes to the
#'   params object that you pass to it but instead returns a new params object.
#'   So to affect the change you call the function in the form
#'   `params <- setInitialValues(params, ...)`.
#' @export
#' @family functions for setting parameters
#' @examples
#' \dontrun{
#' params <- NS_params
#' sim <- project(params, t_max = 20, effort = 0.5)
#' params <- setInitialValues(params, sim)
#' }
setInitialValues <- function(params, sim,
                             initial_n = getInitial_n(params),
                             initial_n_pp = getInitial_n_pp(params),
                             initial_n_other = getInitial_n_other(params)) {
    if (!missing(sim)) {
        assert_that(is(sim, "MizerSim"))
        no_t <- dim(sim@n)[1]
        initial_n < sim@params@initial_n # Needed to get the right dimensions
        initial_n[] <- sim@n[no_t, , ]
        initial_n_pp < sim@params@initial_n_pp # Needed to get the right dimensions
        initial_n_pp[] <- sim@n_pp[no_t, ]
        initial_n_other < sim@n_other[[no_t]]
    }
    assert_that(identical(dim(initial_n), dim(params@initial_n)),
                all(initial_n >= 0),
                identical(dim(initial_n_pp), dim(params@initial_n_pp)),
                all(initial_n_pp >= 0),
                identical(length(initial_n_other), length(params@initial_n_other)))
    if (!is.null(dimnames(initial_n)) &&
        !identical(dimnames(initial_n), dimnames(params@initial_n))) {
        warning("The dimnames of initial_n are not as expected. I will ignore them.")
    }
    if (!is.null(dimnames(initial_n_pp)) &&
        !identical(dimnames(initial_n_pp), dimnames(params@initial_n_pp))) {
        warning("The dimnames of initial_n_pp are not as expected. I will ignore them.")
    }
    if (!identical(names(initial_n_other), names(params@initial_n_other))) {
        stop("The names of initial_n_other do not match those in params.")
    }
    params@initial_n[] <- initial_n
    params@initial_n_pp[] <- initial_n_pp
    params@initial_n_other <- initial_n_other
    return(params)
}

#' @rdname setInitialValues
#' @export
getInitial_n <- function(params) {
    params@initial_n
}

#' @rdname setInitialValues
#' @export
getInitial_n_pp <- function(params) {
    params@initial_n_pp
}
#' @rdname setInitialValues
#' @export
getInitial_n_other <- function(params) {
    params@initial_n_other
}

#' Set line colours to be used in mizer plots
#' 
#' @param params A MizerParams object
#' @param colours A named list or named vector of line colours.
#' 
#' @return The MizerParams object with updated line colours
#' @export
#' @examples
#' params <- NS_params
#' params <- setColours(params, list("Cod" = "red", "Haddock" = "#00ff00"))
#' plotSpectra(params)
#' getColours(params)
setColours <- function(params, colours) {
    assert_that(is(params, "MizerParams"),
                all(validColour(colours)))
    params@linecolour <- unlist(
        modifyList(as.list(params@linecolour), as.list(colours)))
    params
}

#' @rdname setColours
#' @export
getColours <- function(params) {
    as.list(params@linecolour)
}

validColour <- function(colour) {
    sapply(colour, function(X) {
        tryCatch(is.matrix(col2rgb(X)), 
                 error = function(e) FALSE)
    })
}

#' Set linetypes to be used in mizer plots
#' 
#' @param params A MizerParams object
#' @param linetypes A named list or named vector of linetypes.
#' 
#' @return The MizerParams object with updated linetypes
#' @export
#' @examples
#' params <- NS_params
#' params <- setLinetypes(params, list("Cod" = "solid"))
#' plotSpectra(params)
#' getLinetypes(params)
setLinetypes <- function(params, linetypes) {
    assert_that(is(params, "MizerParams"))
    params@linetype <- unlist(
        modifyList(as.list(params@linetype), as.list(linetypes)))
    params
}

#' @rdname setLinetypes
#' @export
getLinetypes <- function(params) {
    as.list(params@linetype)
}

#' Set own rate function to replace mizer rate function
#' 
#' At each time step during a simulation with the [project()] function, mizer
#' needs to calculate the instantaneous values of the various rates. By
#' default it calls the [mizerRates()] function which creates a list with the
#' following components:
#' * `encounter` from [mizerEncounter()]
#' * `feeding_level` from [mizerFeedingLevel()]
#' * `pred_rate` from [mizerPredRate()]
#' * `pred_mort` from [mizerPredMort()]
#' * `fishing_mort` from [mizerFMort()]
#' * `mort` from [mizerMort()]
#' * `plankton_mort` from [mizerPlanktonMort()]
#' * `e` from [mizerEReproAndGrowth()]
#' * `e_repro` from [mizerERepro()]
#' * `e_growth` from [mizerEGrowth()]
#' * `rdi` from [mizerRDI()]
#' * `rdd` from [BervertonHoltRDD()]
#' 
#' You can modify these in two ways.
#' 
#' @param params A `MizerParams` object
#' @param rate Name of the rate for which a new function is to be set.
#' @md
#' @export
setRateFunction <- function(params, rate = "Rates", fun) {
    assert_that(is(params, "MizerParams"),
                is.string(rate),
                is.string(fun),
                is.function(get(fun)))
    if (!(rate %in% names(params@rates_funcs))) {
        stop("The `rate` argument must be one of ", 
             toString(names(params@rates_funcs)), ".")
    }
    f <- get0(fun, mode = "function")
    if (is.null(f)) {
        stop(fun, " should be a function")
    }
    # TODO: put some code to test that the function has the right kind of
    # arguments
    params@rates_funcs[[rate]] <- fun
    
    validObject(params)
    params
}

#' @rdname setRateFunction
#' @export
getRateFunction <- function(params, rate = "Rates") {
    assert_that(is(params, "MizerParams"),
                is.string(rate))
    validObject(params)
    if (rate == "All") {
        return(params@rates_funcs)
    }
    if (!(rate %in% names(params@rates_funcs))) {
        stop("The `rate` argument must be one of ", 
             toString(names(params@rates_funcs)), ".")
    }
    params@rates_funcs[[rate]]
}


#' Update the initial values
#' 
#' Recalculates the steady-state abundances in a fixed background
#' given by the current abundances, keeping the abundances fixed in the
#' smallest size class for each species. Then readjusts the \code{erepro}
#' values.
#' 
#' @param params A MizerParams object
#'   
#' @return The MizerParams object with updated \code{initial_n} and 
#'   \code{initial_n_pp} slots.
#' @export
updateInitialValues <- function(params) {
    assert_that(is(params, "MizerParams"))
    # Calculate the rates in the current background
    plankton_mort <- getPlanktonMort(params)
    mumu <- getMort(params)
    gg <- getEGrowth(params)
    # Recompute plankton
    params@initial_n_pp <- params@rr_pp * params@cc_pp / 
        (params@rr_pp + plankton_mort)
    # Recompute all species
    for (sp in 1:length(params@species_params$species)) {
        w_inf_idx <- min(sum(params@w < params@species_params[sp, "w_inf"]) + 1,
                         length(params@w))
        idx <- params@w_min_idx[sp]:(w_inf_idx - 1)
        if (any(gg[sp, idx] == 0)) {
            stop("Can not compute steady state due to zero growth rates")
        }
        n0 <- params@initial_n[sp, params@w_min_idx[sp]]
        params@initial_n[sp, ] <- 0
        params@initial_n[sp, params@w_min_idx[sp]:w_inf_idx] <- 
            c(1, cumprod(gg[sp, idx] / ((gg[sp, ] + mumu[sp, ] * params@dw)[idx + 1]))) *
            n0
    }
    
    # Retune the values of erepro so that we get the correct level of
    # recruitment
    mumu <- getMort(params)
    gg <- getEGrowth(params)
    rdd <- getRDD(params)
    # TODO: vectorise this
    for (i in (1:length(params@species_params$species))) {
        gg0 <- gg[i, params@w_min_idx[i]]
        mumu0 <- mumu[i, params@w_min_idx[i]]
        DW <- params@dw[params@w_min_idx[i]]
        params@species_params$erepro[i] <- params@species_params$erepro[i] *
            params@initial_n[i, params@w_min_idx[i]] *
            (gg0 + DW * mumu0) / rdd[i]
    }
    return(params)
}

#' Upgrade MizerParams object from earlier mizer versions
#' 
#' Occasionally during the development of new features for mizer, the
#' MizerParams object gains extra slots. MizerParams objects created in older
#' versions of mizer are then no longer valid in the new version because of
#' the missing slots. This function adds the missing slots and fills them
#' with default values.
#' 
#' Uses newMultispeciesParams() to create a new MizerParams object using the
#' parameters extracted from the old MizerParams object.
#' 
#' If you only have a serialised version of the old object, for example
#' created via `saveRDS()`, and you get an error when trying to read it in
#' with `readRDS()` then unfortunately you will need to install the old version
#' of mizer first to read the params object into your workspace, then switch
#' to the current version and then call `upgradeParams()`. You can then save
#' the new version again with `saveRDS()`.
#' 
#' @param params An old MizerParams object to be upgraded
#' 
#' @return The upgraded MizerParams object
#' @export
upgradeParams <- function(params) {
    
    if (.hasSlot(params, "srr")) {
        if (is.function(params@srr)) {
            RDD <- "BevertonHoltRDD"
            message('The density-dependent reproduction rate function has been set to "BevertonHoltRDD".')
        } else {
            RDD <- params@srr
        }
    } else {
        RDD <- "BevertonHoltRDD"
    }
    
    if (is.function(params@plankton_dynamics)) {
        params@plankton_dynamics <- "plankton_semichemostat"
        message('The plankton dynamics function has been set to "plankton_semichemostat".')
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
    } else {
        maturity <- NULL
        repro_prop <- NULL
    }
    
    if ("r_max" %in% names(params@species_params)) {
        params@species_params$R_max <- params@species_params$r_max
        params@species_params$r_max <- NULL
        message("The 'r_max' column has been renamed to 'R_max'.")
    }
    
    if (.hasSlot(params, "p")) {
        params@species_params[["p"]] <- params@p
    }
    if (.hasSlot(params, "q")) {
        params@species_params[["q"]] <- params@q
    }
    if (.hasSlot(params, "n")) {
        params@species_params[["n"]] <- params@n
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
        rate = params@rr_pp,
        capacity = params@cc_pp,
        pred_kernel = pred_kernel,
        search_vol = params@search_vol,
        intake_max = params@intake_max,
        metab = metab,
        z0 = params@mu_b,
        maturity = maturity,
        repro_prop = repro_prop,
        RDD = RDD,
        initial_effort = initial_effort)
    
    pnew@linecolour <- params@linecolour
    pnew@linetype <- params@linetype
    pnew@initial_n <- params@initial_n
    pnew@initial_n_pp <- params@initial_n_pp
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
        pnew@plankton_dynamics <- params@plankton_dynamics
    }
    if (.hasSlot(params, "plankton_params")) {
        pnew@plankton_params <- params@plankton_params
    }
    if (.hasSlot(params, "lambda")) {
        pnew@plankton_params[["lambda"]] <- params@lambda
        pnew@plankton_params[["kappa"]] <- params@kappa
        pnew@plankton_params[["n"]] <- params@n
        pnew@plankton_params[["w_pp_cutoff"]] <- max(pnew@w_full[pnew@cc_pp > 0])
    }
    
    if (.hasSlot(params, "A")) {
        pnew@A <- params@A
    }
    
    if (.hasSlot(params, "rates_func")) {
        pnew@rates_funcs$getRates <- params@rates_func
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

#' Set defaults for predation kernel parameters
#'
#' If the predation kernel type has not been specified for a species, then it
#' is set to "lognormal" and the default values are set for the parameters
#' `beta` and `sigma`.
#' @param object Either a MizerParams object or a species parameter data frame
#' @return The `object` with updated columns in the species params data frame.
#' @export
#' @keywords internal
#' @concept helper
default_pred_kernel_params <- function(object) {
    if (is(object, "MizerParams")) {
        # Nothing to do if full pred kernel has been specified
        if (length(dim(object@pred_kernel)) > 1) {
            return(object)
        }
        species_params <- object@species_params
    } else {
        species_params <- object
    }
    
    species_params <- set_species_param_default(species_params,
                                                "pred_kernel_type",
                                                "lognormal")
    # For species where the pred_kernel_type is lognormal, set defaults for
    # sigma and beta if none are supplied
    if (any(species_params$pred_kernel_type == "lognormal")) {
        species_params <- set_species_param_default(species_params,
                                                    "beta", 30)
        species_params <- set_species_param_default(species_params,
                                                    "sigma", 2)
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
    species_params <- default_pred_kernel_params(species_params)
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
#' Sets \code{h} so that the species reaches maturity 
#' size at the age predicted by the von Bertalanffy growth curve parameters
#' \code{k_vb} and (optionally \code{t0}) taken from the species parameter
#' data frame. Also needs the exponent \code{b} from the length-weight
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
        if (anyNA(species_params$k_vb[missing])) {
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
        assert_that(is.number(params@plankton_params$lambda),
                    is.number(params@plankton_params$kappa),
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
        params@initial_n_pp[] <- params@plankton_params$kappa * 
            params@w_full^(-params@plankton_params$lambda)
        avail_energy <- getEncounter(params)[, length(params@w)] /
            params@w[length(params@w)] ^ 
            (2 + params@species_params$q - params@plankton_params$lambda)
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
        stop(paste("Could not calculate default values for the missing species",
             "parameter ks. Got:", params@species_params$ks))
    }
    return(params@species_params$ks)
}

#' Check that an effort vector is specified correctly
#' 
#' Throws an error with an explanatory message when the supplied \code{effort}
#' vector is not valid for the model described by \code{params}.
#' 
#' @param params A MizerParams object
#' @param effort An effort vector
#' 
#' @return TRUE if \code{effort} is valid. Throws an error otherwise.
#' @export
#' @keywords internal
#' @concept helper
validate_effort_vector <- function(params, effort) {
    assert_that(is(params, "MizerParams"),
                is.numeric(effort))
    no_gears <- dim(params@catchability)[1]
    if ((length(effort) > 1) & (length(effort) != no_gears)) {
        stop("Effort vector must be the same length as the number of fishing gears\n")
    }
    # If more than 1 gear need to check that gear names match
    gear_names <- dimnames(params@catchability)[[1]]
    effort_gear_names <- names(effort)
    if (length(effort) == 1 & is.null(effort_gear_names)) {
        effort_gear_names <- gear_names
    }
    if (!all(gear_names %in% effort_gear_names)) {
        stop(paste0("Gear names in the MizerParams object (", 
                    paste(gear_names, collapse = ", "), 
                    ") do not match those in the effort vector."))
    }
    return(TRUE)
}

#' Check validity of species parameters and set defaults for missing but
#' required parameters
#' 
#' @param species_params The user-supplied species parameter data frame
#' @return A valid species parameter data frame
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
        set_species_param_default("gear", "knife_edge_gear") %>% 
        set_species_param_default("alpha", 0.6) %>% 
        set_species_param_default("interaction_p", 1)
}

#' Get critical feeding level
#' 
#' The critical feeding level is the feeding level at which the food intake is
#' just high enough to cover the metabolic costs, with nothing left over for
#' growth or reproduction. 
#' 
#' @param params A MizerParams object
#' @return A matrix (species x size) with the critical feeding level
#' @export
getCriticalFeedingLevel <- function(params) {
    validObject(params)
    params@metab/params@intake_max/params@species_params$alpha
}
