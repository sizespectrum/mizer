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
    length_w <- length(object@w)
    length_w_full <- length(object@w_full)
    
    # Check dw and dw_full are correct length
    if (length(object@dw) != length_w) {
        msg <- paste("dw is length ", length(object@dw),
                     " and w is length ", length_w,
                     ". These should be the same length", sep = "")
        errors <- c(errors, msg)
    }
    
    if (length(object@dw_full) != length_w_full) {
        msg <- paste("dw_full is length ", length(object@dw_full),
                     " and w_full is length ", length_w_full,
                     ". These should be the same length", sep = "")
        errors <- c(errors, msg)
    }
    # Check the array dimensions are good
    # 2D arrays
    if (!all(c(length(dim(object@psi)),
        length(dim(object@intake_max)),
        length(dim(object@search_vol)),
        length(dim(object@metab)),
        length(dim(object@mu_b)),
        length(dim(object@ft_pred_kernel_e)),
        length(dim(object@ft_pred_kernel_p)),
        length(dim(object@interaction)),
        length(dim(object@catchability))) == 2)) {
        msg <- "psi, intake_max, search_vol, metab, mu_b, ft_pred_kernel_e, ft_pred_kernel_p, interaction and catchability must all be two dimensions"
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
        dim(object@ft_pred_kernel_e)[1],
        dim(object@ft_pred_kernel_p)[1],
        dim(object@mu_b)[1],
        dim(object@selectivity)[2],
        dim(object@catchability)[2],
        dim(object@interaction)[1],
        dim(object@interaction)[2]) == 
        dim(object@species_params)[1])) {
        msg <- "The number of species in the model must be consistent across the species_params, psi, intake_max, search_vol, mu_b, interaction (dim 1), selectivity, ft_pred_kernel_e, ft_pred_kernel_p, catchability and interaction (dim 2) slots"
        errors <- c(errors, msg)
    }
    # Check number of size groups
    if (!all(c(
        dim(object@psi)[2],
        dim(object@intake_max)[2],
        dim(object@search_vol)[2],
        dim(object@metab)[2],
        dim(object@selectivity)[3]) ==
        length_w)) {
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
        names(dimnames(object@ft_pred_kernel_e))[1],
        names(dimnames(object@ft_pred_kernel_p))[1],
        names(dimnames(object@selectivity))[2],
        names(dimnames(object@catchability))[2]) == "sp")) {
        msg <- "Name of first dimension of psi, intake_max, search_vol, metab, mu_b, ft_pred_kernel_e, ft_pred_kernel_p and the second dimension of selectivity and catchability must be 'sp'"
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
        dimnames(object@ft_pred_kernel_e)[[1]],
        dimnames(object@ft_pred_kernel_p)[[1]],
        dimnames(object@mu_b)[[1]],
        dimnames(object@selectivity)[[2]],
        dimnames(object@catchability)[[2]],
        dimnames(object@interaction)[[1]],
        dimnames(object@interaction)[[2]]) ==
        object@species_params$species)) {
        msg <- "The species names of species_params, psi, intake_max, search_vol, metab, mu_b, ft_pred_kernel_e, ft_pred_kernel_p, selectivity, catchability and interaction must all be the same"
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
    
    # species_params data.frame must have columns: 
    # species, z0, alpha, eRepro
    species_params_cols <- c("species","z0","alpha","erepro")
    if (!all(species_params_cols %in% names(object@species_params))) {
        msg <- "species_params data.frame must have 'species', 'z0', 'alpha' and 'erepro' columms"
        errors <- c(errors,msg)
    }
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
#' @slot w A numeric vector of size bins used for the community (i.e. fish) part
#'   of the model. These are usually spaced on a log10 scale
#' @slot dw The absolute difference between the size bins specified in the w
#'   slot. A vector the same length as the w slot. The final value is the same
#'   as the second to last value
#' @slot w_full A numeric vector of size bins used for the whole model (i.e. the
#'   community and plankton spectra) . These are usually spaced on a log10
#'   scale
#' @slot dw_full The absolute difference between the size bins specified in the
#'   w_full slot. A vector the same length as the w_full slot. The final value
#'   is the same as the second to last value
#' @slot w_min_idx A vector holding the index of the weight of the egg size
#'   of each species
#' @slot psi An array (species x size) that holds the allocation to reproduction
#'   for each species at size, \eqn{\psi_i(w)}
#' @slot intake_max An array (species x size) that holds the maximum intake for
#'   each species at size. Default \eqn{h_i w^n}
#' @slot search_vol An array (species x size) that holds the search volume for
#'   each species at size. Default \eqn{\gamma_i w^q}
#' @slot metab An array (species x size) that holds the metabolism
#'   for each species at size. Default \eqn{k_{s.i} w^p + k_i w}
#' @slot mu_b An array (species x size) that holds the background death 
#'   \eqn{\mu_{b.i}(w)}
#' @slot ft_pred_kernel_e An array (species x log of predator/prey size ratio) that holds 
#'   the Fourier transform of the feeding kernel in a form appropriate for
#'   evaluating the available energy integral
#' @slot ft_pred_kernel_p An array (species x log of predator/prey size ratio) that holds 
#'   the Fourier transform of the feeding kernel in a form appropriate for
#'   evaluating the predation mortality integral
#' @slot rr_pp A vector the same length as the w_full slot. The size specific
#'   growth rate of the plankton spectrum. Default \eqn{r_0 w^{p-1}}
#' @slot cc_pp A vector the same length as the w_full slot. The size specific
#'   carrying capacity of the plankton spectrum. Default \eqn{\kappa w^{-\lambda}}
#' @slot sc The community abundance of the scaling community
#' @slot species_params A data.frame to hold the species specific parameters
#'   (see the mizer vignette, Table 2, for details)
#' @slot interaction The species specific interaction matrix, \eqn{\theta_{ij}}
#' @slot srr Function to calculate the realised (density dependent) recruitment.
#'   Has two arguments which are rdi and species_params
#' @slot selectivity An array (gear x species x w) that holds the selectivity of
#'   each gear for species and size, \eqn{S_{g,i,w}}
#' @slot catchability An array (gear x species) that holds the catchability of
#'   each species by each gear, \eqn{Q_{g,i}}
#' @slot initial_n An array (species x size) that holds abundance of each species
#'  at each weight at our candidate steady state solution.
#' @slot initial_n_pp A vector the same length as the w_full slot that describes
#'  the plankton abundance at each weight.
#' @slot n Exponent of maximum intake rate.
#' @slot p Exponent of metabolic cost.
#' @slot lambda Exponent of resource spectrum.
#' @slot q Exponent for volumetric search rate.
#' @slot f0 Initial feeding level.
#' @slot kappa Magnitude of resource spectrum.
#' @slot A Abundance multipliers.
#' @slot linecolour A named vector of colour values, named by species. Used 
#'   to give consistent colours to species in plots.
#' @slot linetype A named vector of linetypes, named by species. Used 
#'   to give consistent colours to species in plots.

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
        metab = "array",
        ft_pred_kernel_e = "array",
        ft_pred_kernel_p = "array",
        mu_b = "array",
        rr_pp = "numeric",
        cc_pp = "numeric", # was NinPP, carrying capacity of plankton
        sc = "numeric",
        initial_n_pp = "numeric",
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
        metab = array(NA,dim = c(1,1), dimnames = list(sp = NULL,w = NULL)),
        ft_pred_kernel_e = array(NA,dim = c(1,1), dimnames = list(sp = NULL,k = NULL)),
        ft_pred_kernel_p = array(NA,dim = c(1,1), dimnames = list(sp = NULL,k = NULL)),
        mu_b = array(NA,dim = c(1,1), dimnames = list(sp = NULL,w = NULL)),
        rr_pp = NA_real_,
        cc_pp = NA_real_,
        sc = NA_real_,
        initial_n_pp = NA_real_,
        A = NA_real_,
        linecolour = NA_character_,
        linetype = NA_character_,
        #speciesParams = data.frame(),
        interaction = array(
            NA,dim = c(1,1), dimnames = list(predator = NULL, prey = NULL)
        ), # which dimension is prey and which is prey?
        selectivity = array(
            NA, dim = c(1,1,1), dimnames = list(gear = NULL, sp = NULL, w = NULL)
        ),
        catchability = array(
            NA, dim = c(1,1), dimnames = list(gear = NULL, sp = NULL)
        )
    ),
    validity = validMizerParams
)

#' Basic constructor that creates empty MizerParams object of the right size
#' 
#' @param object Number of species
#' @param min_w  Smallest weight
#' @param max_w  Largest weight
#' @param no_w   Number of weight brackets
#' @param min_w_pp Smallest plankton weight
#' @param no_w_pp  No longer used
#' @param species_names Names of species
#' @param gear_names Names of gears
#' 
#' @return An empty but valid MizerParams object
#' 
emptyParams <- function(object, min_w = 0.001, max_w = 1000, no_w = 100,  
                        min_w_pp = 1e-10, no_w_pp = NA, 
                        species_names=1:object, gear_names=species_names) {
    if (!is.na(no_w_pp))
        warning("New mizer code does not support the parameter no_w_pp")
    # Some checks
    if (length(species_names) != object)
        stop("species_names must be the same length as the value of object argument")
    no_sp <- length(species_names)
    
    # Set up grids:
    # Community grid
    w <- 10^(seq(from = log10(min_w),to = log10(max_w), length.out = no_w))
    dw <- diff(w)
    # Correctly defined dw by using the proper ratio 
    # (successive dw's have a fixed ratio). 
    dw[no_w] <- dw[no_w - 1] * (dw[no_w - 1] / dw[no_w - 2])	
    
    # Set up full grid - plankton + community
    # ERROR if dw > w, nw must be at least... depends on minw, maxw and nw
    if (w[1] <= dw[1]) {
        stop("Your size bins are too close together. You should consider increasing the number of bins, or changing the size range")
    }
    
    # For fft methods we need a constant log step size throughout. 
    # Therefore we use as many steps as are necessary to almost reach min_w_pp. 
    x_pp <- rev(seq(from = log10(min_w), log10(min_w_pp), 
                    by = log10(min_w / max_w) / (no_w - 1))[-1])
    w_full <- c(10^x_pp, w)
    no_w_full <- length(w_full)
    dw_full <- diff(w_full)
    dw_full[no_w_full] <- dw_full[no_w_full - 1] *
        (dw_full[no_w_full - 1] / dw_full[no_w_full - 2])	
    
    # Basic arrays for templates
    mat1 <- array(NA, dim = c(no_sp, no_w), 
                  dimnames = list(sp = species_names, w = signif(w,3)))
    mat2 <- array(NA, dim = c(no_sp, no_w, no_w_full), 
                  dimnames = list(sp = species_names, w_pred = signif(w,3), 
                                  w_prey = signif(w_full,3)))
    
    ft_pred_kernel_e <- array(NA, dim = c(no_sp, no_w_full), 
                              dimnames = list(sp = species_names, k = 1:no_w_full))
    
    # We do not know the second dimension of ft_pred_kernel_p until the species
    # parameters determining the predation kernel are know. 
    # So for now we set it to 2
    ft_pred_kernel_p <- array(NA, dim = c(no_sp, 2), 
                              dimnames = list(sp = species_names, k = 1:2))
    
    selectivity <- array(0, dim = c(length(gear_names), no_sp, no_w), 
                         dimnames = list(gear = gear_names, sp = species_names, 
                                       w = signif(w, 3)))
    catchability <- array(0, dim = c(length(gear_names), no_sp), 
                          dimnames = list(gear = gear_names, sp = species_names))
    interaction <- array(1, dim = c(no_sp, no_sp), 
                         dimnames = list(predator = species_names, 
                                         prey = species_names))
    vec1 <- as.numeric(rep(NA, no_w_full))
    names(vec1) <- signif(w_full,3)
    w_min_idx <- rep(1, no_sp)
    names(w_min_idx) = species_names
    
    # Make an empty data.frame for species_params
    # This is just to pass validity check. 
    # The project method uses the columns species z0 alpha erepro
    # so these must be in there
    # There is also a seperate function to check the dataframe that is
    # passed in by users (not used in validity check)
    species_params <- data.frame(species = species_names,
                                 z0 = NA, alpha = NA, erepro = NA)
    
    # Make an empty srr function, just to pass validity check
    srr <- function(rdi, species_params) return(0)
    
    # Make colour and linetype scales for use in plots
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
    
    # Make the new object
    # Should Z0, rrPP and ccPP have names (species names etc)?
    res <- new("MizerParams",
               w = w, dw = dw, w_full = w_full, dw_full = dw_full, w_min_idx = w_min_idx,
               psi = mat1, initial_n = mat1, intake_max = mat1, search_vol = mat1,
               metab = mat1, mu_b = mat1, ft_pred_kernel_e = ft_pred_kernel_e, 
               ft_pred_kernel_p = ft_pred_kernel_p,
               selectivity = selectivity, catchability = catchability,
               rr_pp = vec1, cc_pp = vec1, sc = w, initial_n_pp = vec1, 
               species_params = species_params,
               interaction = interaction, srr = srr, 
               A = as.numeric(rep(NA, dim(interaction)[1])),
               linecolour = linecolour, linetype = linetype) 
    return(res)
}



#' Construct \code{MizerParams} cobject for multispecies model
#'
#' Provides default functional forms for all slots in the MizerParams object
#' based on user-provided species parameters.
#' 
#' @param object A data frame of species specific parameter values (see notes
#'   below).
#' @param interaction Optional argument to specify the interaction matrix of the
#'   species (predator by prey). If missing a default interaction is used where
#'   all interactions between species are set to 1. Note that any dimnames of
#'   the interaction matrix argument are ignored by the constructor. The
#'   dimnames of the interaction matrix in the returned \code{MizerParams}
#'   object are taken from the species names in the \code{species_params} slot.
#'   This means that the order of the columns and rows of the interaction matrix
#'   argument should be the same as the species name in the
#'   \code{species_params} slot.
#' @param min_w The smallest size of the community spectrum.
#' @param max_w The largest size of the community spectrum.
#'    Default value is the largest w_inf in the community x 1.1.
#' @param no_w The number of size bins in the community spectrum.
#' @param min_w_pp The smallest size of the plankton spectrum.
#' @param no_w_pp Obsolete argument that is no longer used because the number
#'    of plankton size bins is determined because all size bins have to
#'    be logarithmically equally spaced.
#' @param n Scaling of the intake. Default value is 2/3.
#' @param p Scaling of the standard metabolism. Default value is 0.7. 
#' @param q Exponent of the search volume. Default value is 0.8. 
#' @param r_pp Growth rate of the primary productivity. Default value is 10. 
#' @param kappa Carrying capacity of the resource spectrum. Default
#'       value is 1e11.
#' @param lambda Exponent of the resource spectrum. Default value is
#'       (2+q-n).
#' @param w_pp_cutoff The cut off size of the plankton spectrum.
#'       Default value is 10.
#' @param f0 Average feeding level. Used to calculated \code{h} and
#'       \code{gamma} if those are not columns in the species data frame. Also
#'       requires \code{k_vb} (the von Bertalanffy K parameter) to be a column
#'       in the species data frame. If \code{h} and \code{gamma} are supplied
#'       then this argument is ignored. Default is 0.6..
#' @param z0pre If \code{z0}, the mortality from other sources, is not
#'       a column in the species data frame, it is calculated as 
#'       z0pre * w_inf ^ z0exp. Default value is 0.6.
#' @param z0exp If \code{z0}, the mortality from other sources, is not
#'       a column in the species data frame, it is calculated as 
#'       z0pre * w_inf ^ z0exp. Default value is n-1.
#' @param species_names Names of the species. Generally not needed as normally
#'   taken from the \code{object} data.frame.
#' @param gear_names Names of the gears that catch each species. Generally not
#'   needed as normally taken from the \code{object} data.frame. Default is
#'   \code{species_names}.
#' @param ... Additional arguments.
#'
#' @return An object of type \code{MizerParams}
#' 
#' @note The only essential argument to the \code{MizerParams} constructor is a
#'   data frame which contains the species data. The data frame is arranged
#'   species by parameter, so each column of the parameter data frame is a
#'   parameter and each row has the parameters for one of the species in the
#'   model.
#'   
#'   There are some essential columns that must be included in the parameter
#'   data.frame and that do not have default values. Other columns do have
#'   default values, so that if they are not included in the species parameter
#'   data frame, they will be automatically added when the \code{MizerParams}
#'   object is created. See the accompanying vignette for details of these
#'   columns.
#'   
#' @seealso \code{\link{project}} \linkS4class{MizerSim}
#' @export
#' @examples
#' data(NS_species_params_gears)
#' data(inter)
#' params <- multispeciesParams(NS_species_params_gears, inter)
multispeciesParams <- function(object, interaction,
                    min_w = 0.001, max_w = max(object$w_inf) * 1.1, no_w = 100,
                    min_w_pp = 1e-10, no_w_pp = NA,
                    n = 2/3, p = 0.7, q = 0.8, r_pp = 10,
                    kappa = 1e11, lambda = (2 + q - n), w_pp_cutoff = 10,
                    f0 = 0.6, z0pre = 0.6, z0exp = n - 1) {
    
    row.names(object) <- object$species
    no_sp <- nrow(object)
    
    if (missing(interaction)) {
        interaction <- matrix(1, nrow = no_sp, ncol = no_sp)
    }
    
    ## Set default values for missing values in species params  --------------
    # If no gear_name column in object, then named after species
    if (!("gear" %in% colnames(object))) {
        object$gear <- object$species
    }
    
    # If no k (activity coefficient), then set to 0
    if (!("k" %in% colnames(object))) {
        object$k <- rep(NA, no_sp)
    }
    missing <- is.na(object$k)
    if (any(missing)) {
        object$k[missing] <- 0
    }
    
    # If no alpha (conversion efficiency), then set to 0.6
    if (!("alpha" %in% colnames(object))) {
        object$alpha <- rep(NA, no_sp)
    }
    missing <- is.na(object$alpha)
    if (any(missing)) {
        object$alpha[missing] <- 0.6
    }
    
    # If no erepro (reproductive efficiency), then set to 1
    if (!("erepro" %in% colnames(object))) {
        object$erepro <- rep(NA, no_sp)
    }
    missing <- is.na(object$erepro)
    if (any(missing)) {
        object$erepro[missing] <- 1
    }
    
    # If no sel_func column in species_params, set to 'knife_edge'
    if (!("sel_func" %in% colnames(object))) {
        message("\tNote: No sel_func column in species data frame. Setting selectivity to be 'knife_edge' for all species.")
        object$sel_func <- 'knife_edge'
        # Set default selectivity size
        if (!("knife_edge_size" %in% colnames(object))) {
            message("Note: \tNo knife_edge_size column in species data frame. Setting knife edge selectivity equal to w_mat.")
            object$knife_edge_size <- object$w_mat
        }
    }
    
    # If no catchability column in species_params, set to 1
    if (!("catchability" %in% colnames(object))) {
        object$catchability <- rep(NA, no_sp)
    }
    missing <- is.na(object$catchability)
    if (any(missing)) {
        object$catchability[missing] <- 1
    }
    
    # Sort out h column If not passed in directly, is calculated from f0 and
    # k_vb if they are also passed in
    if (!("h" %in% colnames(object))) {
        object$h <- rep(NA, no_sp)
    }
    if (any(is.na(object$h))) {
        message("Note: \tNo h provided for some species, so using f0 and k_vb to calculate it.")
        if (!("k_vb" %in% colnames(object))) {
            stop("\t\tExcept I can't because there is no k_vb column in the species data frame")
        }
        h <- ((3 * object$k_vb) / (object$alpha * f0)) * (object$w_inf ^ (1/3))
        # Only overwrite missing h with calculated values
        missing <- is.na(object$h)
        if (any(is.na(h[missing]))) {
            stop("Could not calculate h, perhaps k_vb is missing?")
        }
        object$h[missing] <- h[missing]
    }
    
    # Sorting out gamma column
    if (!("gamma" %in% colnames(object))) {
        object$gamma <- rep(NA, no_sp)
    }
    if (any(is.na(object$gamma))) {
        message("Note: \tNo gamma provided for some species, so using f0, h, beta, sigma, lambda and kappa to calculate it.")
        lm2 <- lambda - 2
        ae <- sqrt(2 * pi) * object$sigma * object$beta^lm2 *
            exp(lm2^2 * object$sigma^2 / 2) *
            # The factor on the following lines takes into account the cutoff
            # of the integral at 0 and at beta + 3 sigma
            (pnorm(3 - lm2 * object$sigma) + 
                 pnorm(log(object$beta)/object$sigma + lm2 * object$sigma) - 1)
        gamma <- (object$h / (kappa * ae)) * (f0 / (1 - f0))
        # Only overwrite missing gammas with calculated values
        missing <- is.na(object$gamma)
        if (any(is.na(gamma[missing]))) {
            stop("Could not calculate gamma.")
        }
        object$gamma[missing] <- gamma[missing]
    }
    
    # Sort out z0 (background mortality)
    if (!("z0" %in% colnames(object))) {
        object$z0 <- rep(NA, no_sp)
    }
    missing <- is.na(object$z0)
    if (any(missing)) {
        message("Note: \tUsing z0 = z0pre * w_inf ^ z0exp for missing z0 values.")
        object$z0[missing] <- z0pre * object$w_inf[missing]^z0exp
    }
    
    # Sort out ks column
    if (!("ks" %in% colnames(object))) {
        message("Note: \tNo ks column in species data frame so using ks = h * 0.2.")
        object$ks <- object$h * 0.2
    }
    missing <- is.na(object$ks)
    if (any(missing)) {
        object$ks[missing] <- object$h[missing] * 0.2
    }
    
    # Check essential columns: species (name), wInf, wMat, h, gamma,  ks, beta, sigma 
    check_species_params_dataframe(object)
    
    ## Make an empty object of the right dimensions -----------------------------
    res <- emptyParams(no_sp, min_w = min_w, max_w = max_w, no_w = no_w,  
                       min_w_pp = min_w_pp, no_w_pp = NA, 
                       species_names = object$species, 
                       gear_names = unique(object$gear))
    res@n <- n
    res@p <- p
    res@lambda <- lambda
    res@q <- q
    res@f0 <- f0
    res@kappa <- kappa
    
    # If not w_min column in species_params, set to w_min of community
    if (!("w_min" %in% colnames(object)))
        object$w_min <- min(res@w)
    # Check min_w argument is not > w_min in species_params
    if (any(object$w_min < min(res@w))) {
        stop("One or more of your w_min values is less than the smallest size of the community spectrum")
    }
    
    ## Start filling the slots ---------------------------------------------
    res@species_params <- object
    res@w_min_idx <- as.vector(
        tapply(object$w_min, 1:length(object$w_min),
               function(w_min, wx) max(which(wx <= w_min)), wx = res@w))
    # Check dims of interaction argument - make sure it's right
    if (!isTRUE(all.equal(dim(res@interaction), dim(interaction))))
        stop("interaction matrix is not of the right dimensions. Must be number of species x number of species")
    # Check that all values of interaction matrix are 0 - 1. Issue warning if not
    if (!all((interaction >= 0) & (interaction <= 1))) {
        warning("Values in the interaction matrix should be between 0 and 1")
    }
    # In case user has supplied names to interaction matrix which are wrong order
    for (dim_check in 1:length(dimnames(res@interaction))) {
        if (!is.null(dimnames(interaction)[[dim_check]]) & (!(isTRUE(all.equal(dimnames(res@interaction)[[dim_check]],dimnames(interaction)[[dim_check]]))))) {
            warning("Dimnames of interaction matrix do not match the order of species names in the species data.frame. I am now ignoring your dimnames so your interaction matrix may be in the wrong order.")
        }
    }
    res@interaction[] <- interaction
    
    # Now fill up the slots using default formulations:
    # psi - allocation to reproduction - from original Setup() function
    res@psi[] <- 
        unlist(
            tapply(res@w, 1:length(res@w),
                   function(wx, w_inf, w_mat,n) {
                       ((1 + (wx / (w_mat))^-10)^-1) * (wx / w_inf)^(1 - n)
                   },
                   w_inf = object$w_inf, w_mat = object$w_mat, n = n
            )
        )
    # Set w < 10% of w_mat to 0
    res@psi[unlist(
        tapply(res@w, 1:length(res@w),
               function(wx, w_mat) wx < (w_mat * 0.1),
               w_mat = object$w_mat))] <- 0
    # Set all w > w_inf to 1 # Check this is right...
    res@psi[unlist(
        tapply(res@w, 1:length(res@w),
               function(wx, w_inf) (wx/w_inf) > 1,
               w_inf = object$w_inf))] <- 1
    # note sure what a and n0_mult are in get_initial_n
    
    res@intake_max[] <- unlist(tapply(res@w,1:length(res@w),function(wx,h,n)h * wx^n,h=object$h,n=n))
    res@search_vol[] <- unlist(tapply(res@w,1:length(res@w),function(wx,gamma,q)gamma * wx^q, gamma=object$gamma, q=q))
    res@metab[] <-  unlist(tapply(res@w,1:length(res@w),function(wx,ks,k,p)ks * wx^p + k * wx, ks=object$ks,k=object$k,p=p))
    res@mu_b[] <- res@species_params$z0
    
    # Set up predation kernels ------------------------------------------------
    Beta <- log(res@species_params$beta)
    sigma <- res@species_params$sigma
    # w_full has the weights from the smallest relevant plankton, to the largest fish
    x_full <- log(res@w_full)
    # We choose the origin of the x axis to be at the smallest plankton size
    x_full <- x_full - x_full[1]
    dx <- x_full[2] - x_full[1]
    # The first choice makes the calculation agree with the old mizer
    # Dx <- res@w[2]/res@w[1] - 1  # dw = w Dx, 
    # The following gives a better agreement with analytic results
    Dx <- dx
    
    # rr is the log of the maximal predator/prey mass ratio
    # Here we use default rr = Beta + 3*sigma
    rr <- Beta + 3*sigma
    # Perturb rr so it falls on grid points
    rr <- dx*ceiling(rr/dx)
    
    # ft_pred_kernel_e is an array (no_sp x no_w_full) 
    # that holds the Fourier transform of the feeding kernel in a form 
    # appropriate for evaluating the available energy integral
    res@ft_pred_kernel_e <- matrix(0, nrow = no_sp, ncol = length(x_full))
    for (i in 1:no_sp) {
        # We compute the feeding kernel terms and their fft.
        phi <- exp(-(x_full - Beta[i])^2 / (2 * sigma[i]^2))
        phi[x_full > rr[i]] <- 0
        res@ft_pred_kernel_e[i, ] <- Dx * fft(phi)
    }
    
    # ft_pred_kernel_p is an array (no_sp x P (to be determined below)) 
    # that holds the Fourier transform of the feeding kernel in a form 
    # appropriate for evaluating the predation mortality rate integral
    # Determine period used
    P <- max(x_full[length(x_full)] + rr)
    # Determine number of x points used in period
    no_P <- 1 + ceiling(P/dx)  # P/dx should already be integer
    # vector of values for log predator/prey mass ratio
    x_P <- (-1:(no_P - 2)) * dx
    
    # The dimension of ft_pred_kernel_p was not know at the time the res object
    # was initialised. Hence we need to create it with the right dimension here.
    res@ft_pred_kernel_p <- matrix(0, nrow = no_sp, ncol = no_P)
    dimnames(res@ft_pred_kernel_p) <- list(sp = rownames(res@metab),
                                           k = (1:no_P))
    
    for (j in 1:no_sp) {
        phi <- rep(0, no_P)
        # Our phi is a periodic extension of the normal feeding kernel.
        # For 0<=x<=P we use phi[x-P] as our
        # value of the period P extension of phi, since support(phi)=[-rr,0]
        phi[x_P-P >= -rr[j]] <- exp(-(Beta[j]-P+x_P[x_P-P >= -rr[j]])^2/(2*sigma[j]^2))
        # We also save the fft of this vector, so we don't have to use too many fft s in the time evolution
        res@ft_pred_kernel_p[j, ] <- Dx*fft(phi)
    }
    
    # plankton spectrum -------------------------------------------------
    res@rr_pp[] <- r_pp * res@w_full^(n - 1) # weight specific plankton growth rate
    res@cc_pp[] <- kappa*res@w_full^(-lambda) # the resource carrying capacity - one for each mp and m (130 of them)
    res@cc_pp[res@w_full > w_pp_cutoff] <- 0  # set density of sizes < plankton cutoff size
    
    # Beverton Holt esque stock-recruitment relationship ----------------------
    # Can add more functional forms or user specifies own
    res@initial_n_pp <- res@cc_pp
    res@srr <- function(rdi, species_params){
        return(rdi / (1 + rdi/species_params$r_max))
    }
    
    # Set fishing parameters: selectivity and catchability -------------
    # At the moment, each species is only caught by 1 gear so in species_params
    # there are the columns: gear_name and sel_func.
    # BEWARE! This routine assumes that each species has only one gear operating on it
    # So we can just go row by row through the species parameters
    # However, I really hope we can do something better soon
    for (g in 1:nrow(object)) {
        # Do selectivity first
        # get args
        # These as.characters are annoying - but factors everywhere
        arg <- names(formals(as.character(object[g,'sel_func'])))
        # lop off w as that is always the first argument of the selectivity functions
        arg <- arg[!(arg %in% "w")]
        if (!all(arg %in% colnames(object))) {
            stop("All of the arguments needed for the selectivity function are not in the parameter dataframe")
        }
        # Check that there is only one column in object with the same name
        # Check that column of arguments exists
        par <- c(w = list(res@w), as.list(object[g,arg]))
        sel <- do.call(as.character(object[g, 'sel_func']), args = par)
        # Dump Sel in the right place
        res@selectivity[as.character(object[g,'gear']), g, ] <- sel
        # Now do catchability
        res@catchability[as.character(object[g,'gear']), g] <- object[g,"catchability"]
    }
    
    # Store colours and linetypes in slots if contained in species parameters
    if ("linetype" %in% names(object)) {
        linetype <- object$linetype[!is.na(object$linetype)]
        res@linetype[object$species[!is.na(object$linetype)]] <- linetype
    }
    if ("linecolour" %in% names(object)) {
        linecolour <- object$linecolour[!is.na(object$linecolour)]
        res@linecolour[object$species[!is.na(object$linecolour)]] <- linecolour
    }
    
    # Remove catchabiliy from species data.frame, now stored in slot
    #params@species_params[,names(params@species_params) != "catchability"]
    res@species_params <- res@species_params[, -which(names(res@species_params) == "catchability")]
    res@initial_n <- res@psi
    res@initial_n <- get_initial_n(res)
    res@A <- rep(1, no_sp)
    return(res)
}

#' Alias for multispeciesParams
#' 
#' An alias provided for backward compatibility with mizer version <= 1.0
#' @inherit multispeciesParams
#' @export
MizerParams <- multispeciesParams


# Check that the species_params dataset is OK
# internal only
check_species_params_dataframe <- function(species_params) {
    # Check species_params dataframe (with a function) for essential cols
    # Essential columns: species (name) # wInf # wMat # h # gamma - search Volume #  ks # beta # z0
    essential_cols <- c("species","w_inf","w_mat","h","gamma","ks","beta","sigma", "z0")
    missing_cols <- !(essential_cols %in% colnames(species_params))
    if (any(missing_cols)) {
        errors <- character()
        for (i in essential_cols[missing_cols]) {
            errors <- paste(errors, i, sep = " ")
        }
        stop("You are missing these columns from the input dataframe:\n", errors)
    }
    return(TRUE)
}

