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
#' Although it is possible to build a `MizerParams` object by hand it is
#' not recommended and several constructors are available. Dynamic simulations
#' are performed using [project()] function on objects of this class. As a 
#' user you should never need to access the slots inside a `MizerParams` object
#' directly. 
#' 
#' @slot w The size grid for the fish part of the spectrum. An increasing
#'   vector of weights (in grams) running from the smallest egg size to the
#'   largest asymptotic size.
#' @slot dw The widths (in grams) of the size bins
#' @slot w_full The size grid for the full size range including the resource
#'   spectrum. An increasing vector of weights (in grams) running from the
#'   smallest resource size to the largest asymptotic size of fish. The
#'   last entries of the vector have to be equal to the content of the w slot.
#' @slot dw_full The width of the size bins for the full spectrum. The last
#'   entries have to be equal to the content of the dw slot.
#' @slot w_min_idx A vector holding the index of the weight of the egg size
#'   of each species
#' @slot maturity An array (species x size) that holds the proportion of
#'   individuals of each species at size that are mature. This enters in the
#'   calculation of the spawning stock biomass with [getSSB()]. Set 
#'   with [setReproduction()].
#' @slot psi An array (species x size) that holds the allocation to reproduction
#'   for each species at size, \eqn{\psi_i(w)}. Changed with 
#'   [setReproduction()].
#' @slot intake_max An array (species x size) that holds the maximum intake for
#'   each species at size. Changed with [setMaxIntakeRate()].
#' @slot search_vol An array (species x size) that holds the search volume for
#'   each species at size. Changed with [setSearchVolume()].
#' @slot metab An array (species x size) that holds the metabolism
#'   for each species at size. Changed with [setMetabolicRate()].
#' @slot mu_b An array (species x size) that holds the external mortality rate
#'   \eqn{\mu_{b.i}(w)}. Changed with [setExtMort()].
#' @slot pred_kernel An array (species x predator size x prey size) that holds
#'   the predation coefficient of each predator at size on each prey size. If
#'   this is NA then the following two slots will be used. Changed with 
#'   [setPredKernel()].
#' @slot ft_pred_kernel_e An array (species x log of predator/prey size ratio)
#'   that holds the Fourier transform of the feeding kernel in a form
#'   appropriate for evaluating the encounter rate integral. If this is NA
#'   then the `pred_kernel` will be used to calculate the available 
#'   energy integral. Changed with [setPredKernel()].
#' @slot ft_pred_kernel_p An array (species x log of predator/prey size ratio)
#'   that holds the Fourier transform of the feeding kernel in a form
#'   appropriate for evaluating the predation mortality integral. If this is NA
#'   then the `pred_kernel` will be used to calculate the integral.
#'   Changed with [setPredKernel()].
#' @slot rr_pp A vector the same length as the w_full slot. The size specific
#'   growth rate of the resource spectrum. Changed with [setResource()].
#' @slot cc_pp A vector the same length as the w_full slot. The size specific
#'   carrying capacity of the resource spectrum. Changed with 
#'   [setResource()].
#' @slot resource_dynamics Name of the function for projecting the resource abundance
#'   density by one timestep. The default is 
#'   [resource_semichemostat()]. 
#'   Changed with [setResource()].
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
#'   See [newMultispeciesParams()] for details.
#' @slot gear_params Data frame with parameters for gear selectivity. See 
#'   [setFishing()] for details.
#' @slot interaction The species specific interaction matrix, \eqn{\theta_{ij}}.
#'   Changed with [setInteraction()].
#' @slot selectivity An array (gear x species x w) that holds the selectivity of
#'   each gear for species and size, \eqn{S_{g,i,w}}. Changed with 
#'   [setFishing()].
#' @slot catchability An array (gear x species) that holds the catchability of
#'   each species by each gear, \eqn{Q_{g,i}}. Changed with 
#'   [setFishing()].
#' @slot initial_effort A vector containing the initial fishing effort for each
#'   gear. Changed with [setFishing()].
#' @slot initial_n An array (species x size) that holds the initial abundance of
#'   each species at each weight.
#' @slot initial_n_pp A vector the same length as the w_full slot that describes
#'   the initial resource abundance at each weight.
#' @slot initial_n_other A list with the initial abundances of all other
#'   ecosystem components. Has length zero if there are no other components.
#' @slot resource_params List with parameters for resource. See [setResource()].
#' @slot A Abundance multipliers.
#' @slot linecolour A named vector of colour values, named by species.
#'   Used to give consistent colours in plots.
#' @slot linetype A named vector of linetypes, named by species. 
#'   Used to give consistent line types in plots.
#' @slot ft_mask An array (species x w_full) with zeros for weights larger than
#'   the asymptotic weight of each species. Used to efficiently minimize
#'   wrap-around errors in Fourier transform calculations.
#' 
#' The \linkS4class{MizerParams} class is fairly complex with a large number of
#' slots, many of which are multidimensional arrays. The dimensions of these
#' arrays is strictly enforced so that `MizerParams` objects are consistent
#' in terms of number of species and number of size classes.
#'   
#' The `MizerParams` class does not hold any dynamic information, e.g.
#' abundances or harvest effort through time. These are held in
#' \linkS4class{MizerSim} objects.
#' 
#' @seealso [project()] [MizerSim()]
#'   [emptyParams()] [newMultispeciesParams()]
#'   [newCommunityParams()]
#'   [newTraitParams()]
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
        resource_dynamics = "character",
        resource_params = "list",
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
        gear_params = "data.frame",
        selectivity = "array",
        catchability = "array",
        initial_effort = "numeric",
        A = "numeric",
        linecolour = "character",
        linetype = "character",
        ft_mask = "array"
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
#' The `species_params` slot of the returned MizerParams object may differ
#' slightly from the data frame supplied as argument to this function in the
#' following ways:
#' \itemize{
#'   \item Default values are set for \code{w_min, w_inf, alpha, gear, interaction_resource}.
#'   \item The egg sizes in `w_min` are rounded down to lie on a grid point.
#' }
#' Note that the other characteristic sizes of the species, like `w_mat` and
#' `w_inf`, are not modified to lie on grid points.
#' 
#' @param species_params A data frame of species-specific parameter values.
#' @param gear_params A data frame with gear-specific parameter values.
#' @param no_w The number of size bins in the consumer spectrum.
#' @param min_w Sets the size of the eggs of all species for which this is not
#'   given in the `w_min` column of the `species_params` dataframe.
# #' @param w_full Increasing vector of weights giving the boundaries of size
# #'   classes. Must include the value min_w. Has one more entry than the number
# #'   of size bins. The last entry is the upper end of the largest size class. It
# #'   be used to calculate the sizes of the size bins but will not be stored in
# #'   the w_full slot of the returned MizerParams object. If this argument is not
# #'   provided then size classes are set by the other arguments as described in
# #'   the Details.
#' @param max_w The largest size of the consumer spectrum. By default this is
#'   set to the largest `w_inf` specified in the `species_params` data
#'   frame.
#' @param min_w_pp The smallest size of the resource spectrum.
# #'   Ignored if w_full is specified.
#' 
#' @return An empty but valid MizerParams object
#' @seealso See [newMultispeciesParams()] for a function that fills
#'   the slots left empty by this function.
#' @export
emptyParams <- function(species_params,
                        gear_params = data.frame(),
                        no_w = 100,
                        min_w = 0.001,
                        # w_full = NA,
                        max_w = NA,
                        min_w_pp = 1e-12) {
    assert_that(is.data.frame(species_params),
                is.data.frame(gear_params),
                no_w > 10)
    
    ## Set defaults ----
    if (is.na(min_w_pp)) min_w_pp <- 1e-12
    species_params <- set_species_param_default(species_params, "w_min", min_w)
    min_w <- min(species_params$w_min)
    
    species_params <- validSpeciesParams(species_params)
    gear_params <- validGearParams(gear_params, species_params)
    
    if (is.na(max_w)) {
        max_w <- max(species_params$w_inf)
    } else {
        if (max(species_params$w_inf) > max_w * (1 + 1e-9)) { # The fudge factor
            # is there to avoid false alerts due to rounding errors.
            too_large <- species_params$species[max_w < species_params$w_inf]
            stop("Some of your species have an maximum size larger than max_w: ",
                 toString(too_large))
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
        # To avoid issues due to numerical imprecision
        min_w <- w[1]
        
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
    gear_names <- unique(gear_params$gear)
    mat1 <- array(NA, dim = c(no_sp, no_w), 
                  dimnames = list(sp = species_names, w = signif(w,3)))
    ft_pred_kernel <- array(NA, dim = c(no_sp, no_w_full),
                            dimnames = list(sp = species_names, k = 1:no_w_full))
    ft_mask <- plyr::aaply(species_params$w_inf, 1,
                           function(x) w_full < x, .drop = FALSE)
    
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
    names(w_min_idx) <- species_names
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
    linecolour <- c(linecolour, "Total" = "black", "Resource" = "green",
                    "Background" = "grey", "Fishing" = "red")
    
    if ("linetype" %in% names(species_params)) {
        linetype <- species_params$linetype
        linetype[is.na(linetype)] <- "solid"
    } else {
        linetype <- rep(type_palette, length.out = no_sp)
    }
    names(linetype) <- as.character(species_names)
    linetype <- c(linetype, "Total" = "solid", "Resource" = "solid",
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
        gear_params = gear_params,
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
            ResourceMort = "mizerResourceMort",
            RDI = "mizerRDI",
            RDD = "BevertonHoltRDD"),
        resource_dynamics = "resource_semichemostat",
        other_params = list(),
        initial_n_other = list(),
        A = as.numeric(rep(NA, no_sp)),
        linecolour = linecolour,
        linetype = linetype,
        ft_mask = ft_mask
    )
    
    return(params)
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
    params@linecolour
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


#' Size bins
#' 
#' This is a good place to explain how mizer discretises the size
#' 
#' @param params A MizerParams object
#' 
#' @export
w <- function(params) {
    params@w
}

#' @rdname w
#' @export
w_full <- function(params) {
    params@w_full
}

#' @rdname w
#' @export
dw <- function(params) {
    params@dw
}

#' @rdname w
#' @export
dw_full <- function(params) {
    params@dw_full
}
