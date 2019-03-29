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
#' @param min_w_pp The smallest size of the plankton spectrum.
# #'   Ignored if w_full is specified.
#' @param no_w_pp  No longer used
#' @param srr The stock recruitment function. Default is
#'   \code{\link{srrBevertonHolt}}.
#' 
#' @return An empty but valid MizerParams object
#' 
#' @export
emptyParams <- function(species_params,
                        no_w = 100,
                        min_w = 0.001,
                        # w_full = NA,
                        max_w = NA,
                        min_w_pp = 1e-10,
                        no_w_pp = NA,
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
    interaction <- array(1, dim = c(no_sp, no_sp),
                         dimnames = list(predator = species_names,
                                         prey = species_names))
    interaction_p <- rep(1, no_sp)
    names(interaction_p) <- species_names
    
    vec1 <- as.numeric(rep(NA, no_w_full))
    names(vec1) <- signif(w_full,3)
    
    w_min_idx <- as.vector(
        tapply(species_params$w_min, 1:no_sp,
               function(w_min, wx) max(which(wx <= w_min)), wx = w))
    names(w_min_idx) = species_names
    
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
        interaction_p = interaction_p,
        srr = srr,
        resource_dynamics = list(),
        plankton_dynamics = plankton_semichemostat,
        resource_params = list(),
        initial_B = 0,
        A = as.numeric(rep(NA, no_sp)),
        linecolour = linecolour,
        linetype = linetype
    )
    
    return(params)
}


#' Set species interation matrix and plankton interaction vector
#' 
#' @section Setting interactions:
#' 
#' The species interaction matrix \eqn{\theta_{ij}} and the plankton interaction
#' vector \eqn{\theta_{ip}} are used when calculating the food encounter rate in
#' \code{\link{getEncounter}} and the predation rate in
#' \code{\link{getPredRate}}. This function checks that the supplied arguments
#' are valid and then stores them in the \code{interaction} and
#' \code{interaction_p} slot of the params object before returning that object.
#' Arguments that are not supplied or NULL will not change the corresponding
#' slot.
#' 
#' Any dimnames of the interaction matrix argument are ignored by the
#' constructor. The dimnames of the interaction matrix in the returned
#' \code{MizerParams} object are taken from the species names in the
#' \code{species_params} slot. This means that the order of the columns and rows
#' of the interaction matrix argument should be the same as the species name in
#' the \code{species_params} slot. The same applies to the names of the
#' \code{interaction_p} argument.
#' 
#' @param params MizerParams object
#' @param interaction Interaction matrix of the species (predator by prey).
#'   Entries should be numbers between 0 and 1. See "Setting interactions"
#'   section below.
#' @param interaction_p Vector specifying for each species its interaction with
#'   plankton, similar to what the interaction matrix does for the interaction
#'   with other species. See "Setting interactions" section below.
#' 
#' @return MizerParams object
#' @export
setInteraction <- function(params,
                           interaction = NULL,
                           interaction_p = NULL) {
    if (!is.null(interaction)) {
        # Check dims of interaction argument
        if (!identical(dim(params@interaction), dim(interaction))) {
            stop( "interaction matrix is not of the right dimensions. Must be number of species x number of species.")
        }
        # Check that all values of interaction matrix are 0 - 1.
        if (!all((interaction >= 0) & (interaction <= 1))) {
            warning("Values in the interaction matrix should be between 0 and 1")
        }
        # In case user has supplied names to interaction matrix, check them.
        if (!is.null(dimnames(interaction))) {
            if (!is.null(names(dimnames(interaction)))) {
                if (!identical(names(dimnames(interaction)),
                               names(dimnames(params@interaction)))) {
                    warning(paste0("Your interaction matrix has dimensions called: ",
                                   toString(names(dimnames(interaction))),
                                   ". I expected 'predator, prey'. I will now ignore your names."))
                }
            }
            names(dimnames(interaction)) <- names(dimnames(params@interaction))
            if (!identical(dimnames(params@interaction),
                           dimnames(interaction))) {
                warning("Dimnames of interaction matrix do not match the order of species names in the species data.frame. I am now ignoring your dimnames so your interaction matrix may be in the wrong order.")
            }
        }
        params@interaction[] <- interaction
    }
    if (!is.null(interaction_p)) {
        # Check dims of interaction_p argument
        if (!identical(length(params@interaction_p), length(interaction_p))) {
            stop("Plankton interaction vector has wrong length. Must be equal to number of species.")
        }
        # Check that all values of interaction vector are 0 - 1.
        if (!all((interaction_p >= 0) & (interaction_p <= 1))) {
            warning("Values in the plantkon interaction vector should be between 0 and 1")
        }
        # In case user has supplied names to interaction vector, check them.
        if (!is.null(names(interaction_p)) &
            (!identical(names(params@interaction_p), names(interaction_p)))) {
            warning("Names in the plankton interaction vector do not match the order of species names in the species data.frame. I am now ignoring your names so your interaction vector may be in the wrong order.")
        }
        params@interaction_p[] <- interaction_p
    }
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
#' \subsection{Kernel dependent on predator to prey size ratio}{
#' If the \code{pred_kernel} argument is not supplied, then this function sets a
#' predation kernel that depends only on the ratio of predator mass to prey
#' mass, not on the two masses independently. The shape of that kernel is then
#' determined by the \code{pred_kernel_type} argument.
#'
#' The default pred_kernel_type is "lognormal". This will call the function
#' \code{\link{lognormal_pred_kernel}} to calculate the predation kernel and the
#' function \code{\link{lognormal_max_ppmr}} to return the maximal predator/prey
#' mass ratio for each species.
#'
#' An alternative pred_kernel type is "box", implemented by the functions
#' \code{\link{box_pred_kernel}} and \code{\link{box_max_ppmr}}.
#'
#' You can use any other string as the type. If for example you choose "my" then
#' you need to define a function \code{my_pred_kernel} that you can model on the
#' existing functions like \code{\link{lognormal_pred_kernel}}. You can also
#' define a function \code{my_max_ppmr} if you want the smallest plankton size
#' to be calculated automatically by \code{\link{set_multispecies_model}},
#' otherwise you need to specify the \code{min_w_pp} argument explicitly.
#' 
#' When using a kernel dependent on the predator/prey size ratio only, mizer
#' does not need to store the entire three dimensional array in the MizerParams
#' object. Such an array can be very big when there is a large number of size
#' bins. Instead mizer only needs to store two two-dimensional arrays that hold
#' Fourier transforms of the feeding kernel function that allow the encounter
#' rate and the predation rate to be calculated very efficiently. However mizer
#' gives you the option of storing the full three-dimensional array anyway by
#' setting \code{store_kernel = TRUE}. This might be useful if you have code for
#' analysing the results of a mizer simulation that relies on the full array.
#' }
#' 
#' \subsection{Kernel dependent on both predator and prey size}{
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
#' as the order in the species params dataframe. If you supply a named array
#' then the function will check the order and warn if it is different.
#' }
#' @param params A MizerParams object
#' @param pred_kernel Optional. An array (species x predator size x prey size)
#'   that holds the predation coefficient of each predator at size on each prey
#'   size. The dimensions are thus no_sp, no_w, no_w_full.
#' @param pred_kernel_type Only used if \code{pred_kernel} is not supplied. A
#'   string or vector of strings that determines which predation kernel function
#'   is used for each species. If not supplied, is taken from species parameter
#'   dataframe if possible, or defaults to "lognormal".
#' @param store_kernel Only used if \code{pred_kernel} is not supplied. A
#'   boolean flag that determines whether the full three dimensional predation
#'   kernel array is calculated and stored, even though it is not needed by
#'   mizer. Only useful if you have your own code that relies on this array.
#'   Default is FALSE.
#' 
#' @return A MizerParams object
#' @export
#' @examples
#' \dontrun{
#' ## Set up a MizerParams object
#' data(NS_species_params_gears)
#' data(inter)
#' params <- set_multispecies_model(NS_species_params_gears, inter)
#' 
#' ## If you change predation kernel parameters after setting up a model, you
#' # need to call setPredKernel
#' params@species_params$beta["Cod"] <- 200
#' params <- setPredKernel(params)
#' 
#' ## You can change to a different predation kernel type
#' params <- setPredKernel(params, pred_kernel_type = "box")
#
# ## If only some species need a kernel that depends on predator and prey
# # size separately, you can modify the default kernel.
# pred_kernel <- getPredKernel(params)
# pred_kernel["Herring"] <-
# params<- setPredKernel(params, pred_kernel = pred_kernel)
#' }
setPredKernel <- function(params,
                          pred_kernel = NULL,
                          pred_kernel_type = "lognormal",
                          store_kernel = FALSE) {
    if (!is.null(pred_kernel)) {
        # A pred kernel was supplied, so check it and store it
        if (!identical(dim(pred_kernel), c(dim(params@psi), length(params@w_full)))) {
            stop("The pred_kerel has the wrong dimensions")
        }
        if (!is.null(dimnames(pred_kernel)) && 
            !all(dimnames(pred_kernel)[[1]] == params@species_params$species)) {
            stop(paste0("You need to use the same ordering of species as in the ",
                        "params object: ", toString(params@species_params$species)))
        }
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
    # If pred_kernel_type is not supplied use the one from species_params
    # if exists or "lognormal"
    if (missing(pred_kernel_type) &
        !"pred_kernel_type" %in% names(params@species_params)) {
        params@species_params$pred_kernel_type <- "lognormal"
    } else {
        params@species_params$pred_kernel_type <- pred_kernel_type
    }
    species_params <- params@species_params
    pred_kernel_type <- species_params$pred_kernel_type
    no_sp <- nrow(species_params)
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)
    if (store_kernel) {
        params@pred_kernel <- 
            array(0,
                  dim = c(no_sp, no_w, no_w_full),
                  dimnames = list(sp = species_params$species,
                                  w_pred = signif(params@w, 3),
                                  w_prey = signif(params@w_full, 3))
            )
    }
    # Vector of predator/prey mass ratios
    # The smallest predator/prey mass ratio is 1
    ppmr <- params@w_full / params@w_full[1]
    
    for (i in 1:no_sp) {
        pred_kernel_func <- get0(paste0(pred_kernel_type[i], "_pred_kernel"))
        assert_that(is.function(pred_kernel_func))
        # TODO: check that the arguments of pred_kernel_func are contained in
        # species params
        phi <- pred_kernel_func(ppmr, i, params)
        if (anyNA(phi)) {
            stop(paste0("The function ", pred_kernel_func,
                        "returned NA. Did you correctly specify all required",
                        "parameters in the species parameter dataframe?"))
        }
        # Fourier transform of feeding kernel for evaluating available energy
        params@ft_pred_kernel_e[i, ] <- fft(phi)
        # Fourier transform of feeding kernel for evaluating predation rate
        ri <- max(which(phi > 0))  # index of largest ppmr
        phi_p <- rep(0, no_w_full)
        phi_p[(no_w_full - ri + 1):no_w_full] <- phi[(ri + 1):2]
        params@ft_pred_kernel_p[i, ] <- fft(phi_p)
        # Full feeding kernel array
        if (store_kernel) {
            min_w_idx <- no_w_full - no_w + 1
            for (k in seq_len(no_w)) {
                params@pred_kernel[i, k, (min_w_idx - 1 + k):1] <-
                    phi[1:(min_w_idx - 1 + k)]
            }
        }
    }
    
    return(params)
}


#' Set search volume
#' 
#' @section Setting search volume:
#' The search volume \eqn{\gamma_i(w)} of an individual of species \eqn{i}
#' and weight \eqn{w} multiplies the predation kernel when
#' calculating the encounter rate and the predation rate.
#' 
#' If the \code{gamma} argument is not supplied, then it is set to
#' \deqn{\gamma_i(w) = \gamma_i w^q.} The values of \eqn{gamma_i} are taken from
#' the \code{gamma} column in the species parameter dataframe. If the \code{gamma}
#' column is not supplied in the species parameter dataframe, it is calculated
#' from \code{f0, h, beta, sigma, lambda} and \code{kappa}.
#' 
#' @param params MizerParams
#' @param gamma Optional. An array (species x size) holding the search volume
#'   for each species at size or a vector or number giving the coefficient in
#'   an allometric search volume. If not supplied, a default is set as described in
#'   the section "Setting search volume". 
#' @param q Exponent of the allometric search volume. Only needed if 
#'   \code{gamma} is not an array.
#' 
#' @return MizerParams
#' @export
setSearchVolume <- function(params, 
                            gamma = NULL,
                            q) {
    species_params <- params@species_params
    params@q <- q
    if (!is.null(dim(gamma))) {
        if (!identical(dim(gamma), dim(params@search_vol))) {
            stop("The gamma array has the wrong dimensions.")
        }
        if (!is.null(dimnames(gamma)) && 
            !all(dimnames(gamma)[[1]] == species_params$species)) {
            stop(paste0("You need to use the same ordering of species as in the ",
                        "params object: ", toString(species_params$species)))
        }
        params@search_vol[] <- gamma
        return(params)
    }
    if (!is.null(gamma)) {
        if (length(gamma) == nrow(params@species_params) || length(gamma) == 1) {
            species_params$gamma <- gamma
        } else {
            stop("The gamma argument has the wrong length")
        }
    }
    # Sorting out gamma column in species_params
    if (!("gamma" %in% colnames(species_params))) {
        species_params$gamma <- rep(NA, nrow(species_params))
    }
    missing <- is.na(species_params$gamma)
    if (any(missing)) {
        message("Note: No gamma provided for some species, so using f0, h, beta, sigma, lambda and kappa to calculate it.")
        lm2 <- params@lambda - 2
        ae <- sqrt(2 * pi) * species_params$sigma * species_params$beta^lm2 *
            exp(lm2^2 * species_params$sigma^2 / 2) *
            # The factor on the following lines takes into account the cutoff
            # of the integral at 0 and at beta + 3 sigma
            (pnorm(3 - lm2 * species_params$sigma) + 
                 pnorm(log(species_params$beta)/species_params$sigma + 
                           lm2 * species_params$sigma) - 1)
        gamma_default <- (species_params$h / (params@kappa * ae)) * 
            (params@f0 / (1 - params@f0))
        # Only overwrite missing gammas with calculated values
        if (any(is.na(gamma_default[missing]))) {
            stop("Could not calculate gamma.")
        }
        species_params$gamma[missing] <- gamma_default[missing]
    }
    params@species_params$gamma <- species_params$gamma
    params@search_vol[] <- outer(species_params$gamma, params@w^q)
    return(params)
}


#' Set maximum intake rate
#'
#' @section Setting maximum intake rate:
#' The maximum intake rate \eqn{h_i(w)} of an individual of species \eqn{i} and
#' weight \eqn{w} determines the feeding level, calculated with
#' \code{\link{getFeedingLevel}}.
#'
#' If the \code{h} argument is not supplied, then it is set to \deqn{h_i(w) =
#' h_i w^n.} The values of \eqn{h_i} are taken from the \code{h} column in the
#' species parameter dataframe. If the \code{h} column is not supplied in the
#' species parameter dataframe, it is calculated from \code{f0} and the
#' \code{k_vb} column, if they are supplied.
#' 
#' If \eqn{h_i} is set to \code{Inf}, fish will consume all encountered food.
#'
#' @param params MizerParams
#' @param h Optional. An array (species x size) holding the maximum intake rate
#'   for each species at size. If not supplied, a default is set as described in
#'   the section "Setting maximum intake rate".
#' @param n Scaling exponent of the intake rate.
#' @return MizerParams
#' @export
setIntakeMax <- function(params, h = NULL, n) {
    species_params <- params@species_params
    params@n <- n
    if (!is.null(h)) {
        if (!identical(dim(h), dim(params@intake_max))) {
            stop("The h array has the wrong dimensions.")
        }
        if (!is.null(dimnames(h)) && 
            !all(dimnames(h)[[1]] == species_params$species)) {
            stop(paste0("You need to use the same ordering of species as in the ",
                        "params object: ", toString(species_params$species)))
        }
        params@intake_max[] <- h
        return(params)
    }
    # 
    if (!("h" %in% colnames(species_params))) {
        species_params$h <- rep(NA, nrow(species_params))
    }
    if (any(is.na(species_params$h))) {
        message("Note: No h provided for some species, so using f0 and k_vb to calculate it.")
        if (!("k_vb" %in% colnames(species_params))) {
            stop("\tExcept I can't because there is no k_vb column in the species data frame")
        }
        h <- ((3 * species_params$k_vb) / (species_params$alpha * params@f0)) * 
            (species_params$w_inf ^ (1/3))
        # Only overwrite missing h with calculated values
        missing <- is.na(species_params$h)
        if (any(is.na(h[missing]))) {
            stop("Could not calculate h, perhaps k_vb is missing?")
        }
        species_params$h[missing] <- h[missing]
        params@species_params$h <- species_params$h
    }
    params@intake_max[] <- outer(species_params$h, params@w^params@n)
    return(params)
}


#' Set metabolic rate
#' 
#' Sets the rate at which energy is used for metabolism and activity
#' @section Setting metabolic rate:
#' To be documented
#' 
#' @param params MizerParams
#' @param metab Optional. An array (species x size) holding the metabolic rate
#'   for each species at size. If not supplied, a default is set as described in
#'   the section "Setting metabolic rate".
#' @param p Scaling exponent of the standard metabolic rate.
#' 
#' @return MizerParams
#' @export
setMetab <- function(params, metab = NULL, p) {
    species_params <- params@species_params
    params@p <- p
    if (!is.null(metab)) {
        if (!identical(dim(metab), dim(params@metab))) {
            stop("The metab array has the wrong dimensions.")
        }
        if (!is.null(dimnames(metab)) && 
            !all(dimnames(metab)[[1]] == species_params$species)) {
            stop(paste0("You need to use the same ordering of species as in the ",
                        "params object: ", toString(species_params$species)))
        }
        params@metab[] <- metab
        return(params)
    }
    
    # If no k (activity coefficient), then set to 0
    if (!("k" %in% colnames(species_params))) {
        species_params$k <- rep(NA, nrow(species_params))
    }
    missing <- is.na(species_params$k)
    if (any(missing)) {
        species_params$k[missing] <- 0
    }
    
    # Sort out ks column
    if (!("ks" %in% colnames(species_params))) {
        message("Note: No ks column in species data frame so using ks = h * 0.2.")
        species_params$ks <- rep(NA, nrow(species_params))
    }
    missing <- is.na(species_params$ks)
    if (any(missing)) {
        species_params$ks[missing] <- species_params$h[missing] * 0.2
        params@species_params$ks <- species_params$ks
    }
    
    params@metab[] <- 
        outer(species_params$ks, params@w^p) +
        outer(species_params$k, params@w)
    return(params)
}


#' Set background mortality rate
#' 
#' @section Setting background mortality rate:
#' Still need to document
#' 
#' @param params MizerParams
#' @param z0pre If \code{z0}, the mortality from other sources, is not a column
#'   in the species data frame, it is calculated as z0pre * w_inf ^ z0exp.
#'   Default value is 0.6.
#' @param z0exp If \code{z0}, the mortality from other sources, is not a column
#'   in the species data frame, it is calculated as \code{z0pre * w_inf ^ z0exp}.
#'   Default value is \code{n-1}.
#' 
#' @return MizerParams
#' @export
setBMort <- function(params, z0pre = 0.6, z0exp = params@n - 1) {
    species_params <- params@species_params
    
    # Sort out z0 (background mortality)
    if (!("z0" %in% colnames(species_params))) {
        species_params$z0 <- rep(NA, nrow(species_params))
    }
    missing <- is.na(species_params$z0)
    if (any(missing)) {
        message("Note: Using z0 = z0pre * w_inf ^ z0exp for missing z0 values.")
        species_params$z0[missing] <- z0pre * species_params$w_inf[missing]^z0exp
        params@species_params$z0 <- species_params$z0
    }
    
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
#' If the \code{psi} argument is not supplied,
#' the proportion is set to the product of a sigmoidal maturity ogive that 
#' gives the proportion of individuals of a given species and size that are
#' mature, and a factor that describes how investment into reproduction by mature
#' individuals scales with size. In formulas:
#' \deqn{\psi(w) = \left[1+\left(\frac{w}{w_{mat}}\right)^{-U}\right]^{-1}
#'   \left(\frac{w}{w_{inf}}\right)^{m-n}.}{
#'   [1+(w/w_mat)^(-U)]^(-1) * (w/w_inf)^(m - n)}
#' Here \eqn{n} is the scaling exponent of the energy income rate.
#' The exponent \eqn{m} determines the scaling of the investment into
#' reproduction for mature individuals. By default it is chosen to be \eqn{m =
#' 1} so that the rate at which energy is invested into reproduction scales
#' linearly with the size. This default can be overridden by including a column
#' \code{m} in the species parameter dataframe.
#' 
#' The exponent \eqn{U} determines the steepness of the maturity ogive. By default it is
#' chosen as \eqn{U = 10}, however this can be overridden by including a 
#' column \code{w_mat25} in the species parameter dataframe that specifies the weight
#' at which 25\% of individuals are mature, which sets 
#' \deqn{U = \frac{\log(3)}{\log(w_{mat} / w_{25})}}{U = log(3)/ log(w_mat / w_25)}
#' 
#' @param params A MizerParams object
#' @param psi Optional. An array (species x size) that holds the allocation to
#'   reproduction for each species at size. If not supplied, a default is set
#'   as described in the section "Setting reproduction".
#' 
#' @return The MizerParams object with the updated \code{psi} slot.
#' @export
setReproduction <- function(params, psi = NULL) {
    species_params <- params@species_params
    
    # If no erepro (reproductive efficiency), then set to 1
    if (!("erepro" %in% colnames(species_params))) {
        species_params$erepro <- rep(NA, nrow(species_params))
    }
    missing <- is.na(species_params$erepro)
    if (any(missing)) {
        species_params$erepro[missing] <- 1
        params@species_params$erepro <- species_params$erepro
    }
    
    if (!is.null(psi)) {
        if (!identical(dim(psi), dim(params@psi))) {
            stop("The psi array has the wrong dimensions.")
        }
        if (!is.null(dimnames(psi)) && 
            !all(dimnames(psi)[[1]] == species_params$species)) {
            stop(paste0("You need to use the same ordering of species as in the ",
                        "params object: ", toString(species_params$species)))
        }
        params@psi[] <- psi
        return(params)
    }
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
    # Check maturity sizes
    if (!("w_mat" %in% colnames(species_params))) {
        stop("The maturity sizes of the species must be specified in the w_mat column of the species parameter data frame.")
    }
    missing <- is.na(species_params$w_mat)
    if (any(missing)) {
        stop(paste("The following species are missing data for their maturity size w_mat:"),
             toString(species_params$species[missing]))
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
    
    # Set defaults for m
    if (!("m" %in% colnames(species_params))) {
        species_params$m <- 1 - params@n
    }
    missing <- is.na(species_params$m)
    if (any(missing)) {
        species_params$m[missing] <- 1 - params@n
    }
    # Check m
    assert_that(is.numeric(species_params$m), 
                all(species_params$m > 0 & species_params$m < 2))
    
    params@psi[] <- 
        unlist(
            tapply(params@w, 1:length(params@w),
                   function(wx, w_inf, w_mat, m, w_mat25) {
                       U <- log(3) / log(w_mat / w_mat25)
                       return((1 + (wx / w_mat)^-U)^-1 * (wx / w_inf)^m)
                   },
                   w_inf = species_params$w_inf, 
                   w_mat = species_params$w_mat, 
                   m = species_params$m,
                   w_mat25 = species_params$w_mat25
            )
        )
    # For reasons of efficiency we next set all very small values to 0 
    # Set w < 10% of w_mat to 0
    params@psi[outer(species_params$w_mat * 0.1, params@w, ">")] <- 0
    # Set all w > w_inf to 1
    params@psi[outer(species_params$w_inf, params@w, "<")] <- 1
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
#'   state. The default is \code{"\link{plankton_semichemostat}"}.
#' 
#' @return A MizerParams object
#' @export
setPlankton <- function(params,
                        kappa,
                        lambda,
                        r_pp = 10, 
                        w_pp_cutoff = 10,
                        plankton_dynamics = plankton_semichemostat) {
    # TODO: check arguments
    params@kappa <- kappa
    params@lambda <- lambda
    # weight specific plankton growth rate
    params@rr_pp[] <- r_pp * params@w_full^(params@n - 1)
    # the plankton carrying capacity
    params@cc_pp[] <- kappa*params@w_full^(-lambda)
    params@cc_pp[params@w_full > w_pp_cutoff] <- 0
    params@plankton_dynamics <- plankton_dynamics
    
    return(params)
}

#' Set resource dynamics
#' 
#' @section Setting resource dynamics:
#' Still need to document
#' 
#' @param params A MizerParams object
#' @param resource_dynamics A named list of functions that determine the
#'   dynamics of the unstructured resources by calculating their biomasses at
#'   the next time step from the current state. See
#'   \code{\link{resource_dynamics}} for details. An empty list if the model
#'   does not have unstructured resources.
#' @param resource_params A list of parameters needed by the
#'   \code{resource_dynamics} functions. An empty list if no parameters are
#'   needed.
setResources <- function(params,
                         resource_dynamics,
                         resource_params) {
    if (!missing(resource_dynamics)) {
        assert_that(is.list(resource_dynamics))
        no_res <- length(resource_dynamics)
        resource_names = names(resource_dynamics)
        if (no_res == 0) {
            params@rho <- array(0, dim = 0)
            params@resource_dynamics <- list()
            params@initial_B <- 0
        } else {
            if (is.null(resource_names)) {
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
            params@initial_B <- rep(0, no_res)
            names(params@initial_B) <- resource_names
        }
    }
    if (!missing(resource_params)) {
        assert_that(is.list(resource_params))
        params@resource_params <- resource_params
    }
    return(params)
}

#' Set resource encounter rate
#' 
#' @section Setting resource encounter rate:
#' The resource encounter rate \eqn{\rho_{id}(w)} determines the rate at which
#' an individual of species \eqn{i} encounters biomass of resource \eqn{d},
#' \deqn{\sum_d\rho_{id}(w) B_d,} where \eqn{B_d} is the biomass of the d-th
#' unstructured resource component. 
#' See \code{\link{resource_dynamics}} for more details.
#' 
#' If only a two-dimensional array (species x resource) is given for the
#' \code{rho} argument, then the size dependence is set to a power-law:
#' \deqn{\rho_{id}(w) = \rho_{id} w^n.}
#' 
#' The ordering of the entries in the array \code{rho} is important. The order
#' of the species in the first array dimension needs to be the same as that in
#' the species parameter dataframe. The order of the resources in the second
#' array dimension must be the same as in the list of resource dynamics. The
#' third dimension, if given, is the size dimension.
#' 
#' @param params A MizerParams object
#' @param rho Either a 3-dim array (predator species x resource x predator size)
#'   or a 2-dim array (predator species x resource). In the latter case the
#'   size dependence is given by allometry.
#' @param n Scaling exponent of the intake rate.
#' 
#' @return A MizerParams object
#' @export
setResourceEncounter <- function(params, rho = NULL, n) {
    params@n <- n
    if (is.null(rho)) {
        return(params)
    }
    # Check validity of arguments
    if (!length(dim(rho)) %in% c(2, 3)) {
        stop("The rho argument must be a two- or three-dimensional array.")
    }
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
    if (length(dim(rho)) == 2) {
        rho <- outer(rho, params@w^n)
    } else {
        if (length(params@w) != dim(rho)[3]) {
            stop("The third dimension of the rho array should have one entry for every consumer size.")
        }
    }
    params@rho[] <- rho
    params@initial_B[] <- rep(1, no_res)  # TODO: find better initial value
    
    return(params)
}


#' Set fishing parameters
#' 
#' @section Setting fishing:
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


#' Beverton Holt stock-recruitment function
#' 
#' @param rdi x
#' @param species_params x
#' 
#' @return rdd
srrBevertonHolt <- function(rdi, species_params) {
    return(rdi / (1 + rdi/species_params$r_max))
}