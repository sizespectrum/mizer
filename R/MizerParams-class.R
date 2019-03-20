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
#' Dynamic simulations are performed using the \code{\link{project}} method on
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
#' @slot psi An array (species x size) that holds the allocation to reproduction
#'   for each species at size, \eqn{\psi_i(w)}
#' @slot intake_max An array (species x size) that holds the maximum intake for
#'   each species at size.
#' @slot search_vol An array (species x size) that holds the search volume for
#'   each species at size.
#' @slot rho A 3-dim array (species x resource x size) holding the encounter
#'   rates for unstructured resources. See [resource_dynamics] for details.
#' @slot metab An array (species x size) that holds the metabolism
#'   for each species at size.
#' @slot mu_b An array (species x size) that holds the background death 
#'   \eqn{\mu_{b.i}(w)}
#' @slot pred_kernel An array (species x predator size x prey size) that holds
#'   the predation coefficient of each predator at size on each prey size. If
#'   this is NA then the following two slots will be used.
#' @slot ft_pred_kernel_e An array (species x log of predator/prey size ratio)
#'   that holds the Fourier transform of the feeding kernel in a form
#'   appropriate for evaluating the encounter rate integral. If this is NA
#'   then the \code{pred_kernel} will be used to calculate the available 
#'   energy integral.
#' @slot ft_pred_kernel_p An array (species x log of predator/prey size ratio)
#'   that holds the Fourier transform of the feeding kernel in a form
#'   appropriate for evaluating the predation mortality integral. If this is NA
#'   then the \code{pred_kernel} will be used to calculate the integral.
#' @slot rr_pp A vector the same length as the w_full slot. The size specific
#'   growth rate of the plankton spectrum.
#' @slot cc_pp A vector the same length as the w_full slot. The size specific
#'   carrying capacity of the plankton spectrum.
#' @slot resource_dynamics A named list of functions for projecting the
#'   biomasses in the unstructured resource components by one timestep. The
#'   names of the list entries are the resource names. See
#'   \code{\link{resource_dynamics}} for details.
#' @slot plankton_dynamics A function for projecting the plankton abundance
#'   density by one timestep. See \code{\link{plankton_semichemostat}} for 
#'   an example.
#' @slot resource_params A list containing the parameters needed by the
#'   `resource_dynamics` functions, see \code{\link{resource_dynamics}} for
#'   details.
#' @slot sc The community abundance of the scaling community
#' @slot species_params A data.frame to hold the species specific parameters
#'   (see the mizer vignette, Table 2, for details)
#' @slot interaction The species specific interaction matrix, \eqn{\theta_{ij}}
#' @slot interaction_p The species specific interaction with plankton,
#'   \eqn{\theta_{ip}}
#' @slot srr Function to calculate the realised (density dependent) recruitment.
#'   Has two arguments which are rdi and species_params
#' @slot selectivity An array (gear x species x w) that holds the selectivity of
#'   each gear for species and size, \eqn{S_{g,i,w}}
#' @slot catchability An array (gear x species) that holds the catchability of
#'   each species by each gear, \eqn{Q_{g,i}}
#' @slot initial_n An array (species x size) that holds abundance of each species
#'   at each weight at our candidate steady state solution.
#' @slot initial_n_pp A vector the same length as the w_full slot that describes
#'   the plankton abundance at each weight.
#' @slot initial_B A vector containing the biomasses of the unstructured
#'   resource components, see \code{\link{resource_dynamics}} for details.
#' @slot n Exponent of maximum intake rate.
#' @slot p Exponent of metabolic cost.
#' @slot lambda Exponent of plankton spectrum.
#' @slot q Exponent for volumetric search rate.
#' @slot f0 Initial feeding level.
#' @slot kappa Magnitude of plankton spectrum.
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
    representation(
        w = "numeric",
        dw = "numeric",
        w_full = "numeric",
        dw_full = "numeric",
        w_min_idx = "numeric",
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
        interaction_p = "numeric",
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
    prototype = prototype(
        w = NA_real_,
        dw = NA_real_,
        w_full = NA_real_,
        dw_full = NA_real_,
        w_min_idx = NA_real_,
        n = NA_real_,
        p = NA_real_,
        lambda = NA_real_,
        q = NA_real_,
        f0 = NA_real_,
        kappa = NA_real_,
        psi = array(NA,dim = c(1,1), dimnames = list(sp = NULL,w = NULL)),
        initial_n = array(NA,dim = c(1,1), dimnames = list(sp = NULL,w = NULL)),
        intake_max = array(NA,dim = c(1,1), dimnames = list(sp = NULL,w = NULL)),
        search_vol = array(NA,dim = c(1,1), dimnames = list(sp = NULL,w = NULL)),
        rho = array(NA,dim = c(1,1), dimnames = list(sp = NULL,w = NULL)),
        metab = array(NA,dim = c(1,1), dimnames = list(sp = NULL,w = NULL)),
        pred_kernel = array(
            NA, dim = c(1,1,1), dimnames = list(
                sp = NULL, w_pred = NULL, w_prey = NULL
            )
        ),
        ft_pred_kernel_e = array(NA,dim = c(1,1), dimnames = list(sp = NULL,k = NULL)),
        ft_pred_kernel_p = array(NA,dim = c(1,1), dimnames = list(sp = NULL,k = NULL)),
        mu_b = array(NA,dim = c(1,1), dimnames = list(sp = NULL,w = NULL)),
        rr_pp = NA_real_,
        cc_pp = NA_real_,
        sc = NA_real_,
        initial_n_pp = NA_real_,
        initial_B = NA_real_,
        A = NA_real_,
        linecolour = NA_character_,
        linetype = NA_character_,
        #speciesParams = data.frame(),
        interaction = array(
            NA,dim = c(1,1), dimnames = list(predator = NULL, prey = NULL)
        ),
        interaction_p = NA_real_,
        selectivity = array(
            NA, dim = c(1,1,1), dimnames = list(gear = NULL, sp = NULL, w = NULL)
        ),
        catchability = array(
            NA, dim = c(1,1), dimnames = list(gear = NULL, sp = NULL)
        )
    ),
    validity = validMizerParams
)


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
#' @param min_w_pp The smallest size of the plankton spectrum. By default this
#'   is set to the smallest value at which any of the consumers can feed.
# #'   Ignored if w_full is specified.
#' @param no_w_pp  No longer used
#' @param resource_dynamics A named list of functions that determine the
#'   dynamics of the unstructured resources by calculating their biomasses at
#'   the next time step from the current state. See
#'   \code{\link{resource_dynamics}} for details. An empty list if the model
#'   does not have unstructured resources.
#' @param resource_params A list of parameters needed by the
#'   \code{resource_dynamics} functions. An empty list if no parameters are
#'   needed.
#' @inheritParams setInteraction
#' @param srr The stock recruitment function. Default is
#'   \code{\link{srrBevertonHolt}}.
#' 
#' @return An empty but valid MizerParams object
#' 
#' @export
emptyParams <- function(species_params,
                        interaction = matrix(1,
                                             nrow = nrow(species_params),
                                             ncol = nrow(species_params)),
                        no_w = 100,
                        min_w = 0.001,
                        # w_full = NA,
                        max_w = NA,
                        min_w_pp = 1e-10,
                        no_w_pp = NA,
                        resource_dynamics = list(),
                        resource_params = list(),
                        srr = srrBevertonHolt) {
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
    if (!("w_min" %in% colnames(species_params))) {
        species_params$w_min <- rep(NA, no_sp)
    }
    missing <- is.na(species_params$w_min)
    if (any(missing)) {
        species_params$w_min[missing] <- min_w
    }
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
        if (max(species_params$w_inf) < max_w) {
            too_large <- species_params$species[max_w > species_params$w_inf]
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
    if (!("alpha" %in% colnames(species_params))) {
        species_params$alpha <- rep(NA, no_sp)
    }
    missing <- is.na(species_params$alpha)
    if (any(missing)) {
        species_params$alpha[missing] <- 0.6
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
        
        # For fft methods we need a constant log bin size throughout. 
        # Therefore we use as many steps as are necessary so that the first size
        # class includes min_w_pp.
        x_pp <- rev(seq(from = log10(min_w),
                        to = log10(min_w_pp),
                        by = -dx)) - dx
        # If min_w_pp happened to lie exactly on a grid point, we now added
        # one grid point to much which we need to remove again
        if (x_pp[2] == log10(min_w_pp)) {
            x_pp <- x_pp[2:length(x_pp)]
        }
        w_full <- c(10^x_pp, w)
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
    # mat2 <- array(NA, dim = c(no_sp, no_w, no_w_full), 
    #               dimnames = list(sp = species_names, w_pred = signif(w,3), 
    #                               w_prey = signif(w_full,3)))
        
    mat3 <- array(1, dim = c(no_sp, no_sp),
                  dimnames = list(predator = species_names, 
                                  prey = species_names))
    ft_pred_kernel <- array(NA, dim = c(no_sp, no_w_full),
                            dimnames = list(sp = species_names, k = 1:no_w_full))
    
    selectivity <- array(0, dim = c(length(gear_names), no_sp, no_w), 
                         dimnames = list(gear = gear_names, sp = species_names, 
                                         w = signif(w, 3)))
    catchability <- array(0, dim = c(length(gear_names), no_sp), 
                          dimnames = list(gear = gear_names, sp = species_names))
    interaction_p <- rep(1, no_sp)
    names(interaction_p) <- species_names
    
    vec1 <- as.numeric(rep(NA, no_w_full))
    names(vec1) <- signif(w_full,3)
    
    w_min_idx <- as.vector(
        tapply(species_params$w_min, 1:no_sp,
               function(w_min, wx) max(which(wx <= w_min)), wx = w))
    names(w_min_idx) = species_names
    
    ## Resources  set up----
    assert_that(is.list(resource_dynamics))
    assert_that(is.list(resource_params))
    no_res <- length(resource_dynamics)
    resource_names = names(resource_dynamics)
    if (no_res == 0) {
        rho <- array(0, dim = 0)
        resource_dynamics <- list()
        initial_B <- 0
    } else {
        rho <- array(NA, dim = c(no_sp, no_res, no_w), 
                          dimnames = list(sp = species_names, 
                                          res = resource_names,
                                          w = signif(w,3)))
        initial_B <- rep(0, no_res)
        names(initial_B) <- resource_names
    }
    
    # Colour and linetype scales ----
    # for use in plots
    # Colour-blind-friendly palettes
    # From http://dr-k-lo.blogspot.co.uk/2013/07/a-color-blind-friendly-palette-for-r.html
    # cbbPalette <- c("#000000", "#009E73", "#e79f00", "#9ad0f3", "#0072B2", "#D55E00", 
    #                 "#CC79A7", "#F0E442")
    # From http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
    cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", 
                    "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    linecolour <- rep(cbbPalette, length.out = no_sp)
    names(linecolour) <- as.character(species_names)
    linecolour <- c(linecolour, "Total" = "black", "Plankton" = "green",
                    "Background" = "grey")
    linetype <- rep(c("solid", "dashed", "dotted", "dotdash", "longdash", 
                      "twodash"), length.out = no_sp)
    names(linetype) <- as.character(species_names)
    linetype <- c(linetype, "Total" = "solid", "Plankton" = "solid",
                  "Background" = "solid")
    # Override default if colours or linetypes are contained in species parameters
    if ("linetype" %in% names(species_params)) {
        linetype[!is.na(species_params$linetype)] <- 
            species_params$linetype[!is.na(species_params$linetype)]
    }
    if ("linecolour" %in% names(species_params)) {
        linecolour[!is.na(species_params$linecolour)] <- 
            species_params$linecolour[!is.na(species_params$linecolour)]
    }
    
    # Make object ----
    # Should Z0, rrPP and ccPP have names (species names etc)?
    params <- new(
        "MizerParams",
        w = w,
        dw = dw,
        w_full = w_full,
        dw_full = dw_full,
        w_min_idx = w_min_idx,
        psi = mat1,
        initial_n = mat1,
        intake_max = mat1,
        search_vol = mat1,
        rho = rho,
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
        interaction = mat3,
        interaction_p = interaction_p,
        srr = srr,
        resource_dynamics = resource_dynamics,
        plankton_dynamics = plankton_semichemostat,
        resource_params = resource_params,
        initial_B = initial_B,
        A = as.numeric(rep(NA, dim(interaction)[1])),
        linecolour = linecolour,
        linetype = linetype
    )
    
    params <- setInteraction(params, interaction = interaction)
    params <- setFishing(params)
    
    # If no erepro (reproductive efficiency), then set to 1
    if (!("erepro" %in% colnames(species_params))) {
        species_params$erepro <- rep(NA, no_sp)
    }
    missing <- is.na(species_params$erepro)
    if (any(missing)) {
        species_params$erepro[missing] <- 1
        params@species_params$erepro <- species_params$erepro
    }
    return(params)
}


#' Set species interation matrix
#' 
#' Checks that the supplied interaction matrix is valid and then stores it in
#' the \code{interaction} slot of the params object before returning that 
#' object.
#' 
#' Any dimnames of the interaction matrix argument are ignored by the
#' constructor. The dimnames of the interaction matrix in the returned
#' \code{MizerParams} object are taken from the species names in the
#' \code{species_params} slot. This means that the order of the columns and rows
#' of the interaction matrix argument should be the same as the species name in
#' the \code{species_params} slot.
#' 
#' @param params MizerParams object
#' @param interaction Interaction matrix of the species (predator by prey). By
#'   default all interactions between species are set to 1.
#' 
#' @return MizerParams object
#' @export
setInteraction <- function(params, interaction) {
    
    # Check dims of interaction argument - make sure it's right
    if (!isTRUE(all.equal(dim(params@interaction), dim(interaction))))
        stop("interaction matrix is not of the right dimensions. Must be number of species x number of species")
    # Check that all values of interaction matrix are 0 - 1. Issue warning if not
    if (!all((interaction >= 0) & (interaction <= 1))) {
        warning("Values in the interaction matrix should be between 0 and 1")
    }
    # In case user has supplied names to interaction matrix which are wrong order
    for (dim_check in 1:length(dimnames(params@interaction))) {
        if (!is.null(dimnames(interaction)[[dim_check]]) & 
            (!(isTRUE(all.equal(dimnames(params@interaction)[[dim_check]],dimnames(interaction)[[dim_check]]))))) {
            warning("Dimnames of interaction matrix do not match the order of species names in the species data.frame. I am now ignoring your dimnames so your interaction matrix may be in the wrong order.")
        }
    }
    params@interaction[] <- interaction
    return(params)
}

#' Set fishing parameters
#' 
#' Needs to be documented
#' 
#' @param params A MizerParams object
#' 
#' @return MizerParams object
#' @export
setFishing <- function(params) {
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
    if (!("catchability" %in% colnames(species_params))) {
        species_params$catchability <- rep(NA, no_sp)
    }
    missing <- is.na(species_params$catchability)
    if (any(missing)) {
        species_params$catchability[missing] <- 1
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
            stop("All of the arguments needed for the selectivity function are not in the parameter dataframe")
        }
        # Check that there is only one column in species_params with the same name
        # Check that column of arguments exists
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

srrBevertonHolt <- function(rdi, species_params) {
    return(rdi / (1 + rdi/species_params$r_max))
}