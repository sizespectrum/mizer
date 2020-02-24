#' Functions used for setting up models
#'
#' The functions defined in the file wrapper functions set up MizerParams
#' objects for various kinds of size-spectrum models.
#'
#' @section List of Functions:
#' In this list we relate the functions in this file to the sections in
#' the mizer vignette where the corresponding model is described.
#' \tabular{llll}{
#'   Function name \tab Description \tab Section in vignette\cr
#'   \code{\link{newCommunityParams}} \tab Community model \tab 5 \cr
#'   \code{\link{newTraitParams}} \tab Trait-based model \tab 6 \cr
#' }
#'
#' The file also contains functions for adding, removing or rescaling species.
#'
#' @name wrapper_functions
NULL

# Copyright 2012 Finlay Scott and Julia Blanchard.
# Copyright 2018 Gustav Delius and Richard Southwell.
# Development has received funding from the European Commission's Horizon 2020 
# Research and Innovation Programme under Grant Agreement No. 634495 
# for the project MINOUW (http://minouw-project.eu/).
# Distributed under the GPL 3 or later 
# Maintainer: Gustav Delius, University of York, <gustav.delius@york.ac.uk>

#' Set up parameters for a community-type model
#' 
#' This functions creates a \code{\linkS4class{MizerParams}} object so that
#' community-type models can be easily set up and run. 
#' 
#' A community model has
#' several features that distinguish it from the food-web type models. Only one
#' 'species' is resolved, i.e. one 'species' is used to represent the whole
#' community. The plankton spectrum only extends to the start of the community
#' spectrum. Recruitment to the smallest size in the community spectrum is
#' constant and set by the user. As recruitment is constant, the proportion of
#' energy invested in reproduction (the slot \code{psi} of the returned 
#' \code{MizerParams} object) is set to 0. Standard metabolism has been turned 
#' off (the parameter \code{ks} is set to 0). Consequently, the growth rate is 
#' now determined solely by the assimilated food (see the package vignette for 
#' more details).
#' 
#' The function has many arguments, all of which have default values.
#' 
#' Fishing selectivity is modelled as a knife-edge function with one parameter, 
#' \code{knife_edge_size}, which determines the size at which species are 
#' selected.
#' 
#' The resulting \code{MizerParams} object can be projected forward using 
#' \code{project()} like any other \code{MizerParams} object. When projecting 
#' the community model it may be necessary to keep a small time step size
#' \code{dt} of around 0.1 to avoid any instabilities with the solver. You can
#' check for these numerical instabilities by plotting the biomass or abundance
#' through time after the projection.
#' 
#' @param z0 The background mortality of the community. Default value is 0.1.
#' @param alpha The assimilation efficiency of the community. Default value 0.2
#' @param f0 The average feeding level of individuals who feed on a power-law 
#'   spectrum. This value is used to calculate the search rate parameter 
#'   \code{gamma} (see the package vignette). Default value is 0.7.
#' @param h The maximum food intake rate. Default value is 10.
#' @param beta The preferred predator prey mass ratio. Default value is 100.
#' @param sigma The width of the prey preference. Default value is 2.0.
#' @param n The scaling exponent of the maximum intake rate. Default value is 2/3.
#' @param kappa The carrying capacity of the plankton spectrum. Default value
#'   is 1000.
#' @param lambda The exponent of the plankton spectrum. Default value is 2.05.
#' @param r_pp Growth rate parameter for the plankton spectrum. Default value is 10.
#' @param gamma Volumetric search rate. Estimated using \code{h}, \code{f0} and 
#'   \code{kappa} if not supplied.
#' @param recruitment The constant recruitment in the smallest size class of the
#'   community spectrum. By default this is set so that the community spectrum 
#'   continues the plankton spectrum.
#' @param knife_edge_size The size at the edge of the knife-selectivity 
#'   function. Default value is 1000.
#' @param max_w The maximum size of the community. The \code{w_inf} of the 
#'   species used to represent the community is set to this value. The 
#'   default value is 1e6.
#' @param min_w The minimum size of the community. Default value is 1e-3.
#' @inheritParams newTraitParams
#' @param ... Other arguments to pass to the \code{MizerParams} constructor.
#' @export
#' @return An object of type \code{\linkS4class{MizerParams}}
#' @references K. H. Andersen,J. E. Beyer and P. Lundberg, 2009, Trophic and 
#'   individual efficiencies of size-structured communities, Proceedings of the 
#'   Royal Society, 276, 109-114
#' @family functions for setting up models
#' @examples
#' \dontrun{
#' params <- newCommunityParams(f0=0.7, z0=0.2, recruitment=3e7)
#' sim <- project(params, effort = 0, t_max = 100, dt=0.1)
#' plotBiomass(sim)
#' plotSpectra(sim)
#' }
newCommunityParams <- function(max_w = 1e6,
                                min_w = 1e-3,
                                z0 = 0.1,
                                alpha = 0.2,
                                h = 10,
                                beta = 100,
                                sigma = 2.0,
                                n = 2/3,
                                kappa = 1000,
                                lambda = 2.05,
                                f0 = 0.7,
                                r_pp = 10,
                                gamma = NA,
                                knife_edge_size = 1000,
                                recruitment,
                                ...) {
    w_inf <- max_w
    w_pp_cutoff <- min_w
    ks <- 0 # Turn off standard metabolism
    p <- n # But not used as ks = 0
    
    # Make the species data.frame
    species_params <- data.frame(
        species = "Community",
        w_inf = w_inf,
        w_mat = 1e12, # Has no affect as psi set to 0 but we set it to something 
                      # to help the constructor
        h = h, # max food intake
        gamma = gamma, # vol. search rate,
        ks = ks, # standard metabolism coefficient,
        beta = beta,
        sigma = sigma,
        z0 = z0, # background mortality
        alpha = alpha,
        erepro = 1, # not used
        interaction_p = 1,
        sel_func = "knife_edge",
        knife_edge_size = knife_edge_size,
        stringsAsFactors = FALSE
    )
    params <- 
        newMultispeciesParams(species_params, f0 = f0,
                              p = p, n = n, lambda = lambda, 
                              kappa = kappa, min_w = min_w,
                              w_pp_cutoff = w_pp_cutoff, r_pp = r_pp, ...)
    
    initial_n <- array(kappa * params@w ^ (-lambda), 
                       dim = c(1, length(params@w)))
    params <- setInitial(params, initial_n = initial_n)

    params@srr <- "srrConstant"
    if (missing(recruitment)) {
        recruitment <- get_required_recruitment(params)
    }
    params@species_params$constant_recruitment <- recruitment
    params@psi[] <- 0 # Need to force to be 0. Can try setting w_mat but 
                          # due to slope still not 0
    # Set w_mat to NA for clarity - it is not actually being used
    params@species_params$w_mat[] <- NA
    return(params)
}

#' Set up parameters for a trait-based model
#' 
#' This functions creates a \code{MizerParams} object so that trait-based
#' models can be easily set up and run. A trait-based size spectrum model is
#' a simplification of the general size-based model used in
#' \code{mizer}. The species-specific parameters are the same for all species, 
#' except for
#' the asymptotic size, which is considered the most important trait
#' characterizing a species. Other parameters are related to the asymptotic
#' size. For example, the size at maturity is given by \code{w_inf * eta}, 
#' where \code{eta} is
#' the same for all species. For the trait-based model the number of species is
#' not important. For applications of the trait-based model see Andersen &
#' Pedersen (2010). See the \code{mizer} vignette for more details and examples
#' of the trait-based model.
#'
#' The function has many arguments, all of which have default values. Of
#' particular interest to the user are the number of species in the model and
#' the minimum and maximum asymptotic sizes.
#'
#' The characteristic weights of the smallest species are defined by
#' \code{min_w} (egg size), \code{min_w_mat} (maturity size) and 
#' \code{min_w_inf} (asymptotic size). The asymptotic sizes of 
#' the \code{no_sp} species
#' are logarithmically evenly spaced, ranging from \code{min_w_inf} to
#' \code{max_w_inf}. 
#' Similarly the maturity sizes of the species are logarithmically evenly
#' spaced, so that the ratio \code{eta} between maturity size and asymptotic
#' size is the same for all species. If \code{egg_size_scaling = TRUE} then also
#' the ratio between asymptotic size and egg size is the same for all species.
#' Otherwise all species have the same egg size.
#'
#' In addition to setting up the parameters, this function also sets up an
#' initial condition that is close to steady state.
#'
#' Although the trait based model's steady state is often
#' stable without imposing a stock recruitment relationship, the function can
#' set a Beverton-Holt type stock recruitment relationship that imposes a
#' maximal reproduction rate that is a multiple of the recruitment rate at
#' steady state. That multiple is set by the argument \code{rfac}.
#'
#' The search rate coefficient \code{gamma} is calculated using the expected
#' feeding level, \code{f0}.
#'
#' The option of including fishing is given, but the steady state may lose its
#' natural stability if too much fishing is included. In such a case the user
#' may wish to include stabilising effects (like \code{rfac}) to ensure the
#' steady state is stable. Fishing selectivity is modelled as a knife-edge
#' function with one parameter, \code{knife_edge_size}, which is the size at
#' which species are selected. Each species can either be fished by the same
#' gear (\code{knife_edge_size} has a length of 1) or by a different gear (the
#' length of \code{knife_edge_size} has the same length as the number of species
#' and the order of selectivity size is that of the asymptotic size).
#'
#' The resulting \code{MizerParams} object can be projected forward using
#' \code{project()} like any other \code{MizerParams} object. When projecting
#' the model it may be necessary to reduce \code{dt} below 0.1 to avoid any
#' instabilities with the solver. You can check this by plotting the biomass or
#' abundance through time after the projection.
#'
#' @param no_sp The number of species in the model. The default value is 11.
#' @param min_w_inf The asymptotic size of the smallest species in the
#'   community. Default value is 10. This will be rounded to lie on a grid
#'   point.
#' @param max_w_inf The asymptotic size of the largest species in the community.
#'   Default value is 1000. This will be rounded to lie on a grid point.
#' @param min_w The size of the the egg of the smallest species. Default value
#'   is 10^(-4). This also defines the start of the community size spectrum.
#' @param max_w The largest size in the model. By default this is set to
#'   the largest asymptotic size \code{max_w_inf}. Setting it to something
#'   larger only makes sense if you plan to add larger species to the model
#'   later.
#' @param eta Ratio between maturity size and asymptotic size of a species.
#'   Ignored if \code{min_w_mat} is supplied. Default is 10^(-0.6), 
#'   approximately 1/4.
#' @param min_w_mat The maturity size of the smallest species. Default value is
#'   \code{eta * min_w_inf}. This will be rounded to lie on a grid point.
#' @param no_w The number of size bins in the community spectrum. These bins 
#'   will be equally spaced on a logarithmic scale. Default value
#'   is such that there are 50 bins for each factor of 10 in weight.
#' @param min_w_pp The smallest size of the plankton spectrum. By default this
#'   is set to the smallest value at which any of the consumers can feed.
#' @param w_pp_cutoff The largest size of the plankton spectrum. Default
#'   value is max_w_inf unless \code{perfect_scaling = TRUE} when it is Inf.
#' @param n Scaling exponent of the maximum intake rate. Default value is 2/3.
#' @param p Scaling exponent of the standard metabolic rate. By default this is
#'   equal to the exponent \code{n}.
#' @param lambda Exponent of the abundance power law.
#' @param r_pp Growth rate parameter for the plankton spectrum. Default value is 0.1.
#' @param kappa Coefficient in abundance power law. Default value is
#'   0.005.
#' @param alpha The assimilation efficiency of the community. The default value
#'   is 0.4.
#' @param ks Standard metabolism coefficient. Default value is 4.
#' @param h Maximum food intake rate. Default value is 30.
#' @param beta Preferred predator prey mass ratio. Default value is 100.
#' @param sigma Width of prey size preference. Default value is 1.3.
#' @param f0 Expected average feeding level. Used to set \code{gamma}, the
#'   coefficient in the search rate. The default value is 0.6. Ignored if 
#'   \code{gamma} is given explicitly.
#' @param gamma Volumetric search rate. If not provided, default is determined
#'   by \code{\link{get_gamma_default}} using the value of \code{f0}.
#' @param bmort_prop The proportion of the total mortality that comes from
#'   background mortality, i.e., from sources other than predation or fishing. A
#'   number in the interval [0, 1). Default 0.
#' @param rfac The factor such that \code{R_max = rfac * R}, where \code{R_max}
#'   is the maximum recruitment allowed and \code{R} is the steady-state
#'   recruitment. Thus the larger \code{rfac} the less the impact of the
#'   non-linear stock-recruitment curve. The default is 4.
#' @param gear_names The names of the fishing gears. A character vector, the
#'   same length as the number of gears. Default is "knife_edge_gear".
#' @param knife_edge_size The minimum size at which the gear or gears select
#'   fish. A vector with the length equal to the number of gears.
#' @param egg_size_scaling Boolean. Default FALSE. If TRUE, the egg size is a
#'   constant fraction of the maximum size of each species. This fraction is
#'   \code{min_w / min_w_inf}. If FALSE, all species have the egg size
#'   \code{w_min}.
#' @param perfect_scaling Boolean. Default FALSE. If TRUE then parameters are set so
#'   that the community abundance, growth before reproduction and death are
#'   perfect power laws.
#' @param ... Other arguments to pass to the \code{MizerParams} constructor.
#' @export
#' @return An object of type \code{MizerParams}
#' @family functions for setting up models
#' @examples
#' \dontrun{
#' params <- newTraitParams()
#' sim <- project(params, t_max = 5, effort = 0)
#' plotSpectra(sim)
#' }
newTraitParams <- function(no_sp = 11,
                           min_w_inf = 10,
                           max_w_inf = 10 ^ 3,
                           min_w = 10 ^ (-4),
                           max_w = max_w_inf,
                           eta = 10^(-0.6),
                           min_w_mat = min_w_inf * eta,
                           no_w = log10(max_w_inf / min_w) * 50 + 1,
                           min_w_pp = NA,
                           w_pp_cutoff = min_w_inf,
                           n = 2 / 3,
                           p = n,
                           lambda = 2.05,
                           r_pp = 0.1,
                           kappa = 0.005,
                           alpha = 0.4,
                           ks = 4,
                           h = 30,
                           beta = 100,
                           sigma = 1.3,
                           f0 = 0.6,
                           gamma = NA,
                           bmort_prop = 0,
                           rfac = 4,
                           knife_edge_size = 1000,
                           gear_names = "knife_edge_gear",
                           egg_size_scaling = FALSE,
                           perfect_scaling = FALSE,
                           ...) {
    
    ## Check validity of parameters ----
    assert_that(is.logical(egg_size_scaling),
                is.logical(perfect_scaling))
    if (bmort_prop >= 1 || bmort_prop < 0) {
        stop("bmort_prop can not take the value ", bmort_prop,
             " because it should be the proportion of the total mortality",
             " coming from sources other than predation.")
    }
    if (rfac <= 1) {
        message("rfac needs to be larger than 1. Setting rfac=1.01")
        rfac <- 1.01
    }
    no_w <- round(no_w)
    if (no_w < 1) {
        stop("The number of size bins no_w must be a positive integer")
    }
    if (no_w < log10(max_w_inf/min_w)*5) {
        no_w <- round(log10(max_w_inf / min_w) * 5 + 1)
        message(paste("Increased no_w to", no_w, "so that there are 5 bins ",
                      "for an interval from w and 10w."))
    }
    if (no_w > 10000) {
        message("Running a simulation with ", no_w, 
                " size bins is going to be very slow.")
    }
    if (min_w <= 0) {
        stop("The smallest egg size min_w must be greater than zero.")
    }
    if (min_w_inf >= max_w_inf) {
        stop("The asymptotic size of the smallest species min_w_inf must be ",
             "smaller than the asymptotic size of the largest species max_w_inf")
    }
    if (min_w >= min_w_mat) {
        stop("The egg size of the smallest species min_w must be smaller than ",
             "its maturity size min_w_mat")
    }
    if (min_w_mat >= min_w_inf) {
        stop("The maturity size of the smallest species min_w_mat must be ",
             "smaller than its maximum size min_w_inf")
    }
    no_sp <- as.integer(no_sp)
    if (no_sp < 2) {
        stop("The number of species must be at least 2.")
    }
    if (!all(c(n, r_pp, lambda, kappa, alpha, h, beta, sigma, ks, f0) > 0)) {
        stop("The parameters n, lambda, r_pp, kappa, alpha, h, beta, sigma, ks ",
             "and f0, if supplied, need to be positive.")
    }    
    # Check gears
    if (length(knife_edge_size) > no_sp) {
        stop("There cannot be more gears than species in the model")
    }
    if ((length(knife_edge_size) > 1) & (length(knife_edge_size) != no_sp)) {
        warning("Number of gears is less than number of species so gear information is being recycled. Is this what you want?")
    }
    if ((length(gear_names) != 1) & (length(gear_names) != no_sp)) {
        stop("Length of gear_names argument must equal the number of species.")
    }
    
    if (perfect_scaling) {
        egg_size_scaling <- TRUE
        w_pp_cutoff <- Inf
        p <- n
    }
    
    ## Set grid points and characteristic sizes ----
    # in such a way that the sizes all line up with the grid and the species are
    # all equally spaced.
    
    # Divide the range from min_w to max_w into (no_w - 1) logarithmic bins of
    # log size dx so that the last bin starts at max_w
    min_x <- log10(min_w)
    max_x <- log10(max_w)
    dx <- (max_x - min_x) / (no_w - 1) 
    x <- seq(min_x, by = dx, length.out = no_w)
    w <- 10 ^ x
    
    # Find index of nearest grid point to min_w_inf that is an integer multiple
    # of the species spacing away from max_w
    min_x_inf <- log10(min_w_inf)
    max_x_inf <- log10(max_w_inf)
    # bins_per_sp is the number of bins separating species
    bins_per_sp <- round((max_x_inf - min_x_inf) / (dx * (no_sp - 1)))
    min_i_inf <- no_w - (no_sp - 1) * bins_per_sp
    # Maximum sizes for all species
    w_inf <- w[seq(min_i_inf, by = bins_per_sp, length.out = no_sp)]
    
    # Find index of nearest grid point to min_w_mat
    min_x_mat <- log10(min_w_mat)
    min_i_mat <- round((min_x_mat - min_x) / dx) + 1
    # Maturity sizes for all species
    w_mat <- w[seq(min_i_mat, by = bins_per_sp, length.out = no_sp)]
    
    if (egg_size_scaling) {
        # Determine egg weights w_min for all species
        w_min_idx <- seq(1, by = bins_per_sp, length.out = no_sp)
        w_min <- w[w_min_idx]
    } else {
        w_min <- rep(min_w, no_sp)
        w_min_idx <- rep(1, no_sp)
    }

    ## Build Params Object ----
    erepro <- 0.1  # Will be changed later to achieve coexistence
    species_params <- data.frame(
        species = as.factor(1:no_sp),
        w_min = w_min,
        w_inf = w_inf,
        w_mat = w_mat,
        w_min_idx = w_min_idx,
        h = h,
        ks = ks,
        beta = beta,
        sigma = sigma,
        z0 = 0,
        alpha = alpha,
        erepro = erepro,
        sel_func = "knife_edge",
        knife_edge_size = knife_edge_size,
        gear = gear_names,
        stringsAsFactors = FALSE
    )
    params <-
        newMultispeciesParams(
            species_params,
            min_w = min_w,
            no_w = no_w,
            max_w = max_w,
            w_pp_cutoff = max_w,
            lambda = lambda,
            kappa = kappa,
            n = n,
            p = p,
            min_w_pp = min_w_pp,
            r_pp = r_pp
        )
    
    w <- params@w
    dw <- params@dw
    
    ## Construct steady state solution ----
    
    # Get constants for steady-state solution
    # Predation mortality rate coefficient
    mu0 <- get_power_law_mort(params)
    # Add backgound mortality rate
    mu0 <- mu0 / (1 - bmort_prop)
    hbar <- alpha * h * f0 - ks
    if (hbar < 0) {
        stop("The feeding level is not sufficient to maintain the fish.")
    }
    pow <- mu0 / hbar / (1 - n)
    if (pow < 1) {
        message("The ratio of death rate to growth rate is too small, leading to
                an accumulation of fish at their largest size.")
    }
    
    initial_n <- params@psi  # get array with correct dimensions and names
    initial_n[, ] <- 0
    mumu <- mu0 * w^(n - 1)  # Death rate
    i_inf <- min_i_inf  # index of asymptotic size
    i_min <- 1  # index of natural egg size
    for (i in 1:no_sp) {
        gg <- hbar * w^n * (1 - params@psi[i, ])  # Growth rate
        idx <- w_min_idx[i]:(i_inf - 2)
        # Steady state solution of the upwind-difference scheme used in project
        n_exact <- c(1, cumprod(gg[idx] / ((gg + mumu * dw)[idx + 1])))
        # Use the first species for normalisation
        if (i == 1) {
            dist_sp <- bins_per_sp * dx
            mult <- kappa / 
                sum(n_exact * (w^(lambda - 1) * dw)[1:(min_i_inf - 1)]) *
                (10^(dist_sp*(1-lambda)/2) - 10^(-dist_sp*(1-lambda)/2)) / 
                (1-lambda)
        }
        if (!egg_size_scaling) {
            n_exact <- n_exact / n_exact[i_min]
        }
        idxs <- w_min_idx[i]:(i_inf - 1) 
        initial_n[i, idxs] <- n_exact * mult * 
            (w_inf[1] / w_inf[i]) ^ lambda
        i_inf <- i_inf + bins_per_sp
        i_min <- i_min + bins_per_sp
    }
    
    # Calculate the community spectrum
    sc <- colSums(initial_n)
    params@sc <- sc
    
    ##  Setup plankton ----
    plankton_vec <- (kappa * w ^ (-lambda)) - sc
    # Cut off plankton at w_pp_cutoff
    plankton_vec[w >= w_pp_cutoff] <- 0
    if (any(plankton_vec < 0)) {
        message("Note: Negative plankton abundances")
        if (!perfect_scaling) {
            # Do not allow negative plankton abundance
            message("Note: Negative plankton abundance values overwritten with zeros")
            plankton_vec[plankton_vec < 0] <- 0
        }
    }
    params@cc_pp[sum(params@w_full <= w[1]):length(params@cc_pp)] <-
        plankton_vec
    initial_n_pp <- params@cc_pp
    # The cc_pp factor needs to be higher than the desired steady state in
    # order to compensate for predation mortality
    m2_background <- getPlanktonMort(params, initial_n, initial_n_pp)
    params@cc_pp <- (params@rr_pp + m2_background ) * initial_n_pp/params@rr_pp
    
    ## Setup background death ----
    m2 <- getPredMort(params, initial_n, initial_n_pp)
    flag <- FALSE
    for (i in 1:no_sp) {
        # The steplike psi was only needed when we wanted to use the analytic
        # expression for the steady-state solution
        # params@psi[i,] <- (w / w_inf[i]) ^ (1 - n)
        # params@psi[i, w < (w_mat[i] - 1e-10)] <- 0
        # params@psi[i, w > (w_inf[i] - 1e-10)] <- 1
        params@mu_b[i,] <- mu0 * w ^ (n - 1) - m2[i, ]
        if (!perfect_scaling && any(params@mu_b[i,] < 0)) {
            params@mu_b[i, params@mu_b[i,] < 0] <- 0
            flag <- TRUE
        }
    }
    if (flag) {
        message("Note: Negative background mortality rates overwritten with zeros")
    }
    
    ## Set erepro to meet boundary condition ----
    rdi <- getRDI(params, initial_n, initial_n_pp)
    gg <- getEGrowth(params, initial_n, initial_n_pp)
    mumu <- getMort(params, initial_n, initial_n_pp)
    erepro_final <- 1:no_sp  # set up vector of right dimension
    for (i in (1:no_sp)) {
        gg0 <- gg[i, params@w_min_idx[i]]
        mumu0 <- mumu[i, params@w_min_idx[i]]
        DW <- params@dw[params@w_min_idx[i]]
        erepro_final[i] <- erepro * 
            (initial_n[i, params@w_min_idx[i]] *
                 (gg0 + DW * mumu0)) / rdi[i]
    }
    if (is.finite(rfac)) {
        # erepro has been multiplied by a factor of (rfac/(rfac-1)) to
        # compensate for using a stock recruitment relationship.
        erepro_final <- (rfac / (rfac - 1)) * erepro_final
    }
    params@species_params$erepro <- erepro_final
    # Record abundance of fish and plankton at steady state, as slots.
    params@initial_n <- initial_n
    params@initial_n_pp <- initial_n_pp
    # set rmax=fac*RDD
    # note that erepro has been multiplied by a factor of (rfac/(rfac-1)) to
    # compensate for using a stock recruitment relationship.
    params@species_params$R_max <-
        (rfac - 1) * getRDI(params, initial_n, initial_n_pp)

    return(params)
}

#' Set up parameters for a single-species in a Sheldon power-law background
#' 
#' This functions creates a \code{MizerParams} object with a single 
#' species. This species is embedded in a fixed power-law community spectrum 
#' \deqn{N_c(w) = \kappa w^{-\lambda}} 
#' This community provides the food income for the species. Cannibalism is
#' switched off. The predation mortality arises only from the predators in the
#' power-law community and it is assumed that the predators in the community
#' have the same feeding parameters as the foreground species. The function has
#' many arguments, all of which have default values.
#'
#' In addition to setting up the parameters, this function also sets up an
#' initial condition that is close to steady state, under the assumption of
#' no fishing.
#'
#' Although the steady state is often stable without imposing a stock
#' recruitment relationship, the function can set a Beverton-Holt type stock
#' recruitment relationship that imposes a maximal reproduction rate that is a
#' multiple of the recruitment rate at steady state. That multiple is set by the
#' argument \code{rfac}.
#'
#' @inheritParams newTraitParams
#' @param w_inf Asymptotic size of species
#' @param w_min Egg size of species
#' @param eta Ratio between maturity size \code{w_mat} and asymptotic size
#'   \code{w_inf}. Default is 10^(-0.6), approximately 1/4.. Ignored if
#'   \code{w_mat} is supplied explicitly.
#' @param w_mat Maturity size of species. Default value is 
#'   \code{eta * w_inf}.
#' @param ... Other arguments to pass to the \code{MizerParams} constructor.
#' @export
#' @return An object of type \code{MizerParams}
#' @family functions for setting up models
#' @examples
#' \dontrun{
#' params <- newTraitParams()
#' sim <- project(params, t_max = 5, effort = 0)
#' plotSpectra(sim)
#' }
newSheldonParams <- function(w_inf = 100,
                             w_min = 0.001,
                             eta = 10^(-0.6),
                             w_mat = w_inf * eta,
                             no_w = log10(w_inf / w_min) * 50 + 1,
                             n = 3/4,
                             p = n,
                             lambda = 2.05,
                             kappa = 0.005,
                             alpha = 0.4,
                             ks = 4,
                             k_vb = 1,
                             beta = 100,
                             sigma = 1.3,
                             f0 = 0.6,
                             gamma = NA,
                             bmort_prop = 0,
                             rfac = 4,
                             ...) {
    no_sp <- 1
    ## Much of the following code is copied from newTraitParams
    
    ## Check validity of parameters ----
    if (bmort_prop >= 1 || bmort_prop < 0) {
        stop("bmort_prop can not take the value ", bmort_prop,
             " because it should be the proportion of the total mortality",
             " coming from sources other than predation.")
    }
    if (rfac <= 1) {
        message("rfac needs to be larger than 1. Setting rfac=1.01")
        rfac <- 1.01
    }
    no_w <- round(no_w)
    if (no_w < 1) {
        stop("The number of size bins no_w must be a positive integer")
    }
    if (no_w < log10(w_inf/w_min)*5) {
        no_w <- round(log10(w_inf / w_min) * 5 + 1)
        message(paste("Increased no_w to", no_w, "so that there are 5 bins ",
                      "for an interval from w and 10w."))
    }
    if (no_w > 10000) {
        message("Running a simulation with ", no_w, 
                " size bins is going to be very slow.")
    }
    if (w_min <= 0) {
        stop("The smallest egg size w_min must be greater than zero.")
    }
    if (w_min >= w_mat) {
        stop("The egg size of the smallest species w_min must be smaller than ",
             "its maturity size w_mat")
    }
    if (w_mat >= w_inf) {
        stop("The maturity size of the smallest species w_mat must be ",
             "smaller than its maximum size w_inf")
    }
    if (!all(c(n, lambda, kappa, alpha, k_vb, beta, sigma, ks, f0) > 0)) {
        stop("The parameters n, lambda, kappa, alpha, k_vb, beta, sigma, ks and ",
             "f0, if supplied, need to be positive.")
    }
    
    ## Build Params Object ----
    erepro <- 0.1  # Will be changed later to achieve coexistence
    species_params <- data.frame(
        species = as.factor(1),
        w_min = w_min,
        w_inf = w_inf,
        w_mat = w_mat,
        w_min_idx = 1,
        k_vb =  k_vb,
        ks = ks,
        beta = beta,
        sigma = sigma,
        z0 = 0,
        alpha = alpha,
        erepro = erepro,
        stringsAsFactors = FALSE
    )
    params <-
        suppressMessages(newMultispeciesParams(
            species_params,
            min_w = w_min,
            no_w = no_w,
            max_w = w_inf,
            lambda = lambda,
            kappa = kappa,
            n = n,
            p = p,
            w_pp_cutoff = w_inf,
            plankton_dynamics = "plankton_constant"
        ))
    # No cannibalism
    params@interaction[] <- 0
    
    w <- params@w
    dw <- params@dw
    h <- params@species_params$h
    
    ## Construct steady state solution ----
    
    # Get constants for steady-state solution
    # Predation mortality rate coefficient
    mu0 <- get_power_law_mort(params)
    # Add backgound mortality rate
    mu0 <- mu0 / (1 - bmort_prop)
    hbar <- alpha * h * f0 - ks
    if (hbar < 0) {
        stop("The feeding level is not sufficient to maintain the fish.")
    }
    pow <- mu0 / hbar / (1 - n)
    if (pow < 1) {
        message("The ratio of death rate to growth rate is too small, leading to
                an accumulation of fish at their largest size.")
    }
    
    initial_n <- params@psi  # get array with correct dimensions and names
    initial_n[, ] <- 0
    mumu <- mu0 * w^(n - 1)  # Death rate
    params@mu_b[] <- mumu
    i_inf <- sum(params@w <= w_inf)  # index of asymptotic size
    idx <- 1:(i_inf - 1)
    idxs <- 1:i_inf
    gg <- hbar * w^n * (1 - params@psi[1, ])  # Growth rate
    # Steady state solution of the upwind-difference scheme used in project
    initial_n[1, idxs] <- c(1, cumprod(gg[idx] / ((gg + mumu * dw)[idx + 1])))
    
    # The plankton was already set up by newMultispeciesParams()
    initial_n_pp <- params@initial_n_pp
    
    # Normalise abundance so that the maximum ratio between the species
    # abundance and community abundance is 1/2
    fish <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)
    imax <- which.max(initial_n[1, ] / initial_n_pp[fish]) # index in fish spectrum
    pmax <- imax + length(params@w_full) - length(params@w) # corresponding plankton index
    initial_n <- initial_n / initial_n[1, imax] * initial_n_pp[pmax] / 2
    
    ## Set erepro to meet boundary condition ----
    rdi <- getRDI(params, initial_n, initial_n_pp)
    gg <- getEGrowth(params, initial_n, initial_n_pp)
    mumu <- getMort(params, initial_n, initial_n_pp)
    erepro_final <- 1:no_sp  # set up vector of right dimension
    for (i in (1:no_sp)) {
        gg0 <- gg[i, params@w_min_idx[i]]
        mumu0 <- mumu[i, params@w_min_idx[i]]
        DW <- params@dw[params@w_min_idx[i]]
        erepro_final[i] <- erepro * 
            (initial_n[i, params@w_min_idx[i]] *
                 (gg0 + DW * mumu0)) / rdi[i]
    }
    if (is.finite(rfac)) {
        # erepro has been multiplied by a factor of (rfac/(rfac-1)) to
        # compensate for using a stock recruitment relationship.
        erepro_final <- (rfac / (rfac - 1)) * erepro_final
    }
    params@species_params$erepro <- erepro_final
    # Record abundance of fish and plankton at steady state, as slots.
    params@initial_n <- initial_n
    params@initial_n_pp <- initial_n_pp
    # set rmax=fac*RDD
    # note that erepro has been multiplied by a factor of (rfac/(rfac-1)) to
    # compensate for using a stock recruitment relationship.
    params@species_params$R_max <-
        (rfac - 1) * getRDI(params, initial_n, initial_n_pp)
    
    return(params)
}


#' Retunes abundance of background species.
#' 
#' Rescales all background species in such a way that the total community
#' spectrum is as close to the Sheldon power law as possible. Background
#' species that are no longer needed are removed. The reproductive efficiencies
#' of all species are retuned.
#'
#' @param params A \linkS4class{MizerParams} object
#'   
#' @return An object of type \code{MizerParams}
#' @seealso markBackground
#' @export
retuneBackground <- function(params) {
    no_sp <- nrow(params@species_params)  # Number of species
    L <- is.na(params@A)
    
    # We find the abundance multipliers A_i so
    # that the integral of the square of the relative distance 
    # (sum_{i not in L} A_i*N_i(w) + sum_{i not in L} N_i(w) - sc(w))/sc(w) 
    # over w, between our limits, is minimized, where  L is the set of all
    # retuneable species.
    
    # ignore zero entries in params@sc and only use region above the smallest w_mat
    region <- params@sc > 0 & params@w > min(params@species_params$w_mat)
    sc <- params@sc[region]
    # rho is the total abundance of all the non-tunable species
    rho <- colSums(params@initial_n[!L, region, drop = FALSE])
    
    # Use Singular Value Decomposition to find optimal abundance multipliers.
    # See Numerical Recipes section 15.4.2
    #
    # Rescale by sc
    A <- t(sweep(params@initial_n[L, region, drop = FALSE], 2, sc, "/"))
    b <- (sc - rho) / sc
    
    sv <- svd(A)
    di <- 1/sv$d  # inverse of singular values
    di[di > 10^8] <- 0  # cut off
    x <- sweep(sv$v, 2, di, "*") %*% t(sv$u) %*% b
    A2 <- rep(1, no_sp) 
    A2[L] <- x
    
    # We may have to repeat this if any of the multipliers is negative or zero
    if (any(A2 <= 0)) {
        # Remove those species
        params <- removeSpecies(params, species = (A2 <= 0))
        # and try again retuning the remaining retunable species
        if (any(A2 > 0)) {
            params <- retuneBackground(params)
        } else {
            message("All background species have been removed.")
        }
    } else {
        # Use these abundance multipliers to rescale the abundance curves
        params@initial_n <- params@initial_n * A2
    }
    
    return(retuneReproductionEfficiency(params))
}

#' Removes species with abundance below a threshold
#' 
#' This species simply removes the low-abundance species from the params object. 
#' It does not recalculate the steady state for the remaining species or
#' retune their reproductive efficiencies.
#'
#' @param params A \linkS4class{MizerParams} object
#' @param cutoff Species with an abundance at maturity size that is less than 
#'               cutoff times community abundance will be removed. Default 1e-3.
#'   
#' @return An object of type \code{MizerParams}
#' @export
pruneSpecies <- function(params, cutoff = 1e-3) {
    no_sp <- nrow(params@species_params)  # Number of species
    # Determine which species need to be removed
    remove <- c()
    for (i in seq_along(params@species_params$species)) {
        # index of maturity size of this species
        w_mat_idx <- min(which(params@w > params@species_params$w_mat[i]))
        # If species abundance at maturity is less than cutoof * community 
        # abundance at that weight, then remove the species.
        if (params@initial_n[i, w_mat_idx] < params@sc[w_mat_idx] * cutoff) {
            remove <- c(remove, params@species_params$species[i])
        }
    }
    # Remove
    return(removeSpecies(params, remove))
}

#' Remove species from an ecosystem
#' 
#' This function simply removes all entries from the MizerParams object that
#' refer to the selected species. It does not recalculate the steady state for
#' the remaining species or retune their reproductive efficiency.
#' 
#' @param params A mizer params object for the original system.
#' @param species A vector of the names of the species to be deleted or a boolean
#'   vector indicating for each species whether it is to be removed (TRUE) or
#'   not.
#' 
#' @return An object of type \linkS4class{MizerParams}
#' @export
removeSpecies <- function(params, species) {
    no_sp <- length(params@w_min_idx)
    if (is.logical(species)) {
        if (length(species) != no_sp) {
            stop("The boolean species argument has the wrong length")
        }
    } else {
        species <- dimnames(params@initial_n)$sp %in% species
        if (length(species) == 0) {
            warning("The species argument matches none of the species in the params object")
            return(params)
        }
    }
    keep <- !species
    
    params@linecolour <- params@linecolour[!(names(params@linecolour) %in% 
                                                 params@species_params$species[species])]
    params@linetype <- params@linetype[!(names(params@linetype) %in% 
                                             params@species_params$species[species])]
    params@psi <- params@psi[keep, , drop = FALSE]
    params@maturity <- params@maturity[keep, , drop = FALSE]
    params@initial_n <- params@initial_n[keep, , drop = FALSE]
    params@intake_max <- params@intake_max[keep, , drop = FALSE]
    params@search_vol <- params@search_vol[keep, , drop = FALSE]
    params@metab <- params@metab[keep, , drop = FALSE]
    if (length(dim(params@ft_pred_kernel_e)) == 2) {
        params@ft_pred_kernel_e <- params@ft_pred_kernel_e[keep, , drop = FALSE]
    }
    if (length(dim(params@ft_pred_kernel_p)) == 2) {
        params@ft_pred_kernel_p <- params@ft_pred_kernel_p[keep, , drop = FALSE]
    }
    params@mu_b <- params@mu_b[keep, , drop = FALSE]
    params@species_params <- params@species_params[keep, , drop = FALSE]
    params@interaction <- params@interaction[keep, keep, drop = FALSE]
    params@selectivity <- params@selectivity[, keep, , drop = FALSE]
    params@catchability <- params@catchability[, keep, drop = FALSE]
    params@w_min_idx <- params@w_min_idx[keep]
    params@A <- params@A[keep]
    
    validObject(params)
    return(params)
}

#' Rescale Abundance
#' 
#' Multiplies the abundances of all or of selected species by given factors and
#' then retunes the reproductive efficiencies accordingly. 
#' 
#' Does not run the system to steady state. For that you should call 
#' \code{\link{steady}} explicitly afterwards.
#' 
#' @param params A mizer params object
#' @param factor The factor by which the abundance of each species is multiplied.
#'   This can be specified in two ways:
#'   \itemize{
#'   \item A named numeric vector where the name indicates the species and the
#'     value gives the factor for that species. Only the named species are 
#'     affected.
#'   \item  A number that gives the factor for all foreground species.
#'   }
#' 
#' @return An object of type \linkS4class{MizerParams}
#' @export
rescaleAbundance <- function(params, factor) {
    assert_that(is(params, "MizerParams"),
                is.numeric(factor),
                all(factor > 0))
    is_foreground <- !is.na(params@A)
    no_sp <- sum(is_foreground)
    if (length(factor) == 1 && length(names(factor)) == 0) {
        factor <- rep(factor, no_sp)
        names(factor) <- params@species_params$species[is_foreground]
    }
    to_rescale <- names(factor)
    wrong <- setdiff(to_rescale, params@species_params$species)
    if (length(wrong) > 0) {
        stop(paste(wrong, collapse = ", "),
             " do not exist.")
    }
    assert_that(length(to_rescale) == length(factor))

    params@initial_n[to_rescale, ] <- 
        params@initial_n[to_rescale, ] * factor
    
    return(retuneReproductionEfficiency(params))
}

#' Rescale System
#' 
#' The abundances in mizer and some rates depend on the size of the area to
#' which they refer. So they could be given per square meter or per square
#' kilometer or for an entire study area or any other choice of yours. This
#' function allows you to change the size by automatically changing the
#' abundances and rates accordingly.
#' 
#' If you rescale the system by a factor \eqn{c} then this function makes the
#' following rescalings in the params object:
#' \itemize{
#' \item The initial abundances \code{initial_n}, \code{initial_n_pp} and
#'   \code{initial_n_other} are rescaled by \eqn{c}.
#' \item The search volume is rescaled by \eqn{1/c}.
#' \item The plankton carrying capacity is rescaled by \eqn{c}
#' \item The maximum recruitment rate \eqn{R_{max}}, if used, is rescaled by 
#'   \eqn{c}.
#' }
#' The effect of this is that the dynamics of the rescaled system are identical
#' to those of the unscaled system, in the sense that it does not matter whether
#' one first calls \code{rescaleSystem} and then runs a simulation with
#' \code{project} or whether one first runs a simulation and then rescales the
#' resulting abundances.
#' 
#' Note that if you use non-standard plankton dynamics or other components then you
#' may need to rescale additional parameters that appear in those dynamics.
#' 
#' @param params A mizer params object
#' @param factor The factor by which the size is rescaled with respect to which 
#'   the abundances are given.
#' 
#' @return An object of type \linkS4class{MizerParams}
#' @export
rescaleSystem <- function(params, factor) {
    assert_that(is(params, "MizerParams"),
                is.number(factor),
                factor > 0)
    
    # Plankton replenishment rate
    params@cc_pp <- params@cc_pp * factor
    params@plankton_params$kappa <- params@plankton_params$kappa * factor
    
    # Rmax
    # r_max is a deprecated spelling of R_max. Get rid of it.
    if ("r_max" %in% names(params@species_params)) {
        params@species_params$R_max <- params@species_params$r_max
        params@species_params$r_max <- NULL
        message("The 'r_max' column has been renamed to 'R_max'.")
    }
    if ("R_max" %in% names(params@species_params)) {
        params@species_params$R_max <- params@species_params$R_max * factor
    }
    
    # Search volume
    params <- setSearchVolume(params, search_vol = params@search_vol / factor)
    if ("gamma" %in% names(params@species_params)) {
        params@species_params$gamma <- params@species_params$gamma / factor
    }
    
    # Initial values
    initial_n_other <- params@initial_n_other
    for (res in names(initial_n_other)) {
        initial_n_other[[res]] <- initial_n_other[[res]] * factor
    }
    params <- setInitial(params,
                         initial_n = params@initial_n * factor,
                         initial_n_pp = params@initial_n_pp * factor,
                         initial_n_other = initial_n_other)
    
    return(params)
}

#' Rename species
#' 
#' Changes the names of species in a MizerParams object
#' 
#' @param params A mizer params object
#' @param replace A named character vector, with new names as values, and old 
#'   names as names.
#' 
#' @return An object of type \linkS4class{MizerParams}
#' @export
#' @examples
#' \dontrun{
#' replace <- c(Cod = "Kabeljau", Haddock = "Schellfisch")
#' params <- renameSpecies(NS_params, replace)
#' params@species_params$species
#' }
renameSpecies <- function(params, replace) {
    replace[] <- as.character(replace)
    to_replace <- names(replace)
    species <- as.character(params@species_params$species)
    wrong <- setdiff(names(replace), species)
    if (length(wrong) > 0) {
        stop(paste(wrong, collapse = ", "),
             " do not exist.")
    }
    names(species) <- species
    species[to_replace] <- replace
    names(species) <- NULL
    rownames(params@species_params) <- species
    params@species_params$species <- species
    linenames <- names(params@linecolour)
    names(linenames) <- linenames
    linenames[to_replace] <- replace
    names(linenames) <- NULL
    names(params@linecolour) <- linenames
    names(params@linetype) <- linenames
    names(params@w_min_idx) <- species
    dimnames(params@maturity)$sp <- species
    dimnames(params@psi)$sp <- species
    dimnames(params@initial_n)$sp <- species
    dimnames(params@intake_max)$sp <- species
    dimnames(params@search_vol)$sp <- species
    dimnames(params@metab)$sp <- species
    if (length(dim(params@ft_pred_kernel_e)) == 2) {
        dimnames(params@ft_pred_kernel_e)$sp <- species
        dimnames(params@ft_pred_kernel_p)$sp <- species
    } else {
        dimnames(params@pred_kernel)$sp <- species
    }
    dimnames(params@mu_b)$sp <- species
    dimnames(params@interaction)$predator <- species
    dimnames(params@interaction)$prey <- species
    dimnames(params@selectivity)$sp <- species
    dimnames(params@catchability)$sp <- species
    
    validObject(params)
    return(params)
}


#' Add new species
#'
#' Takes a \linkS4class{MizerParams} object and adds additional species with
#' given parameters to the ecosystem. It sets the initial values for these new
#' species to its steady-state solution in the given initial state of the
#' existing ecosystem. This will be close to the true steady-state if the
#' abundances of the new species are sufficiently low. Hence the abundances of
#' the new species are set so that the maximal biomass density of each new
#' species lies at 1/100 of the community power law. The reproductive
#' efficiencies of the new species are set so as to keep them at that low level.
#' 
#' After adding the new species, the background species are not retuned and the
#' system is not run to steady state. You would have to call
#' \code{\link{retuneBackground}} and \code{\link{steady}} explicitly.
#' 
#' @param params A mizer params object for the original system. 
#' @param species_params The species parameters of the new species we
#'   want to add to the system.
#' @param interaction Interaction matrix. A square matrix giving either the
#'   interaction coefficients between all species or only those between the
#'   new species. In the latter case all interaction between an old and a new
#'   species are set to 1. If this argument is missing, all interactions 
#'   involving a new species are set to 1.
#' @inheritParams newTraitParams
#' 
#' @return An object of type \linkS4class{MizerParams}
#' @seealso \code{\link{removeSpecies}}
#' @export
#' @examples
#' \dontrun{
#' params <- newTraitParams(max_w_inf = 5000)
#' params <- markBackground(params)
#' a_m <- 0.0085
#' b_m <- 3.11
#' L_inf_m <- 24.3
#' L_mat <- 11.1
#' species_params <- data.frame(
#'     species = "mullet",
#'     w_min = 0.001, 
#'     w_inf = a_m*L_inf_m^b_m, 
#'     w_mat = a_m*L_mat^b_m, 
#'     beta = 283, 
#'     sigma = 1.8, 
#'     z0 = 0,
#'     alpha = 0.6,
#'     sel_func = "knife_edge", 
#'     knife_edge_size = 100, 
#'     gear = "knife_edge_gear",
#'     k = 0,
#'     k_vb = 0.6,
#'     a = a_m,
#'     b = b_m
#' )
#' params <- addSpecies(params, species_params)
#' plotSpectra(params)
#' sim <- project(params, t_max=50)
#' plotBiomass(sim)
#' }
addSpecies <- function(params, species_params, interaction,
                       initial_effort = NULL) {
    # check validity of parameters ----
    assert_that(is(params, "MizerParams"),
                is.data.frame(species_params))
    if (any(species_params$species %in% params@species_params$species)) {
        stop("You can not add species that are already there.")
    }
    no_old_sp <- nrow(params@species_params)
    old_sp <- 1:no_old_sp
    no_new_sp <- nrow(species_params)
    new_sp <- 1:no_new_sp + no_old_sp
    no_sp <- no_old_sp + no_new_sp
    if (missing(interaction)) {
        # keep existing interactions between old species and
        # set interactions involving new species to 1
        inter <- matrix(1, nrow = no_sp, ncol = no_sp)
        inter[old_sp, old_sp] <- params@interaction
    } else if (all(dim(interaction) == c(no_new_sp, no_new_sp))) {
        # keep existing interactions between old species,
        # set interactions involving an old and a new species to 1
        # and use supplied matrix for interaction among new species
        inter <- matrix(1, nrow = no_sp, ncol = no_sp)
        inter[old_sp, old_sp] <- params@interaction
        inter[new_sp, new_sp] <- interaction
    } else if (all(dim(interaction) != c(no_sp, no_sp))) {
        stop("interaction matrix has invalid dimensions.")
    }
    
    # combine species params ----

    # Move linecolour and linetype into species_params
    params@species_params$linetype <- 
        params@linetype[as.character(params@species_params$species)]
    params@species_params$linecolour <- 
        params@linecolour[as.character(params@species_params$species)]
    
    # Make sure that all columns exist in both data frames
    missing <- setdiff(names(params@species_params), names(species_params))
    species_params[missing] <- NA
    missing <- setdiff(names(species_params), names(params@species_params))
    params@species_params[missing] <- NA
    
    # add the new species (with parameters described by species_params), 
    # to make a larger species_params dataframe.
    combi_species_params <- rbind(params@species_params, species_params,
                                  stringsAsFactors = FALSE)
    # new params object ----
    # use dataframe and global settings from params to make a new MizerParams 
    # object.
    p <- newMultispeciesParams(
        combi_species_params,
        interaction = inter,
        min_w = min(params@w),
        max_w = max(params@w),
        min_w_pp = min(params@w_full),
        no_w = length(params@w),
        initial_effort = initial_effort
    )
    # Use the same plankton spectrum as params
    p@initial_n_pp <- params@initial_n_pp
    p@cc_pp <- params@cc_pp
    p@rr_pp <- params@rr_pp
    p@plankton_dynamics <- params@plankton_dynamics
    p@plankton_params <- params@plankton_params
    # Preserve comment
    comment(p) <- comment(params)
    
    # initial solution ----
    p@initial_n[old_sp, ] <- params@initial_n
    p@A[old_sp] <- params@A
    # Use the same psi and mu_b as before for old species
    p@psi[old_sp, ] <- params@psi
    p@sc <- params@sc
    p@mu_b[old_sp, ] <- params@mu_b
    # we assume same background death for all species
    p@mu_b[new_sp, ] <- rep(params@mu_b[1, ], each = no_new_sp)
    
    # Turn off self-interaction among the new species, so we can determine the
    # growth rates, and death rates induced upon them by the pre-existing species
    p@interaction[new_sp, new_sp] <- 0
    mumu <- getMort(p)
    gg <- getEGrowth(p)
    
    # Compute solution for new species
    for (i in new_sp) {
        g <- gg[i, ]
        mu <- mumu[i, ]
        w_inf_idx <- sum(p@w < p@species_params$w_inf[i])
        idx <- p@w_min_idx[i]:(w_inf_idx - 1)
        if (any(g[idx] == 0)) {
            stop("Can not compute steady state due to zero growth rates for ",
                 p@species_params$species[i])
        }
        p@initial_n[i, ] <- 0
        p@initial_n[i, p@w_min_idx[i]:w_inf_idx] <- 
            c(1, cumprod(g[idx] / ((g + mu * p@dw)[idx + 1])))
        
        # set low abundance ----
        # Normalise solution so that at its maximum it lies at 1/100 of the 
        # Sheldon spectrum.
        # We look at the maximum of abundance times w^lambda
        # because that is always an increasing function at small size.
        idx <- which.max(p@initial_n[i, ] * p@w^p@plankton_params$lambda)
        p@initial_n[i, ] <- p@initial_n[i, ] *
            p@plankton_params$kappa * p@w[idx]^(-p@plankton_params$lambda) / p@initial_n[i, idx] / 100
        p@A[i] <- sum(p@initial_n[i, ] * p@w * p@dw * p@maturity[i, ])
    }
    
    if (any(is.infinite(p@initial_n))) {
        stop("Candidate steady state holds infinities.")
    }
    if (any(is.na(p@initial_n) | is.nan(p@initial_n))) {
        stop("Candidate steady state holds non-numeric values.")
    }
    
    # Turn self interaction back on
    p@interaction[new_sp, new_sp] <- inter[new_sp, new_sp]
    
    # Retune reproductive efficiencies of new species
    p <- retuneReproductionEfficiency(p, p@species_params$species[new_sp])
    
    return(p)
}

#' Retune reproduction efficiency to maintain initial egg abundances
#' 
#' Sets the reproductive efficiency for all species so that the rate of egg
#' production exactly compensates for the loss from the first size class due
#' to growth and mortality. Sets the identical stock recruitment function.
#' 
#' @inheritParams steady
#' @param species A vector of the names of the species to be affected or a
#'   boolean vector indicating for each species whether it is to be affected
#'   (TRUE) or not. By default all species are affected
#' @return A MizerParams object
#' @export
retuneReproductionEfficiency <- function(params, 
                                         species = params@species_params$species) {
    assert_that(is(params, "MizerParams"))

    no_sp <- nrow(params@species_params)
    if (is.logical(species)) {
        if (length(species) != no_sp) {
            stop("The boolean species argument has the wrong length")
        }
    } else {
        species <- dimnames(params@initial_n)$sp %in% species
        if (length(species) == 0) {
            warning("The species argument matches none of the species in the params object")
            return(params)
        }
    }
    mumu <- getMort(params)
    gg <- getEGrowth(params)
    rdi <- getRDI(params)
    eff <- params@species_params$erepro
    for (i in (1:no_sp)[species]) {
        gg0 <- gg[i, params@w_min_idx[i]]
        mumu0 <- mumu[i, params@w_min_idx[i]]
        DW <- params@dw[params@w_min_idx[i]]
        if (!rdi[i] == 0) {
            eff[i] <- params@species_params$erepro[i] *
                (params@initial_n[i, params@w_min_idx[i]] *
                     (gg0 + DW * mumu0)) / rdi[i]
        }
        else {
            eff[i] <- 0.1
        }
    }
    params@species_params$erepro <- eff
    return(setReproduction(params, srr = "srrNone"))
}

#' Determine recruitment rate needed for initial egg abundance
#' 
#' @param params A MizerParams object
#' @return A vector of reproduction rates for all species
get_required_recruitment <- function(params) {
    assert_that(is(params, "MizerParams"))
    
    no_sp <- nrow(params@species_params)
    mumu <- getMort(params)
    gg <- getEGrowth(params)
    reproduction <- params@species_params$erepro # vector of correct length
    for (i in (1:no_sp)) {
        gg0 <- gg[i, params@w_min_idx[i]]
        mumu0 <- mumu[i, params@w_min_idx[i]]
        DW <- params@dw[params@w_min_idx[i]]
        reproduction[i] <- params@initial_n[i, params@w_min_idx[i]] *
            (gg0 + DW * mumu0)
    }
    return(reproduction)
}

#' Set maximum recruitment
#' 
#' Takes a MizerParams object with trivial stock recruitment function and sets 
#' Beverton-Holt stock recruitment with a maximum recruitment that is a chosen
#' factor \code{rfac} higher than the initial-state recruitment.
#' 
#' @param params A MizerParams object
#' @param rfac The factor by which the maximum recruitment should be higher than
#'   the initial-state recruitment
#' 
#' @return A MizerParams object
#' @export
setRmax <- function(params, rfac) {
    assert_that(is(params, "MizerParams"),
                is.numeric(rfac),
                length(rfac) %in% c(1, nrow(params@species_params)),
                all(rfac > 1))
    if (params@srr != "srrNone") {
        stop("setRmax can only be applied to params objects using the identity",
             " stock-recruitment function.")
    }
    # erepro needs to be divided by a factor of 1-1/rfac to
    # compensate for using a stock recruitment relationship
    # because RDD = (1-1/rfac) RDI
    params@species_params$erepro <- 
        params@species_params$erepro / (1 - 1/rfac)
    
    params@species_params$R_max <- params@species_params$w_inf
    params@species_params$R_max <- (rfac - 1) * getRDI(params)
    
    return(setReproduction(params, srr = "srrBevertonHolt"))
}


#' Designate species as background species
#'
#' Marks the specified set of species as background species. Background species
#' are handled differently in some plots and their abundance is automatically
#' adjusted in \code{\link{addSpecies}} to keep the community close to the
#' Sheldon spectrum.
#' 
#' @param object An object of class \linkS4class{MizerParams} or 
#'   \linkS4class{MizerSim}.
#' @param species Name or vector of names of the species to be designated as
#'   background species. By default this is set to all species.
#' @param ... Other arguments (unused)
#' 
#' @return An object of the same class as the \code{object} argument
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- newMultispeciesParams(NS_species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=20, t_save = 0.2)
#' sim <- markBackground(sim, species = c("Sprat", "Sandeel", 
#'                                        "N.pout", "Dab", "Saithe"))
#' plotSpectra(sim)
#' }
markBackground <- function(object, species) {
    if (is(object, "MizerSim")) {
        if (missing(species)) {
            species <- dimnames(object@params@initial_n)$sp
        }
        object@params@A[dimnames(object@params@initial_n)$sp %in% species] <- NA
    } else {
        if (missing(species)) {
            species <- dimnames(object@initial_n)$sp
        }
        object@A[dimnames(object@initial_n)$sp %in% species] <- NA
    }
    return(object)
}

#' Tune params object to be at steady state
#' 
#' This is done by running the dynamics for a specified number of years while
#' keeping recruitment and other components constant to reach steady state. Then
#' the reproductive efficiencies are retuned to achieve that level of
#' recruitment.
#' 
#' @param params A \linkS4class{MizerParams} object
#' @param t_max The maximum number of years to run the simulation. Default is 100.
#' @param t_per The simulation is broken up into shorter runs of `t_per` years,
#'   after each of which we check for convergence. Default value is 7.5. This
#'   should be chosen as an odd multiple of the timestep `dt` in order to be
#'   able to detect period 2 cycles.
#' @param tol The simulation stops when the relative change in the egg
#'   production RDI over `t_per years` is less than `tol` for every background
#'   species. Default value is 1/100.
#' @param dt The time step to use in `project()`.
#' @param return_sim If TRUE, the function returns the MizerSim object holding
#'   the result of the simulation run. If FALSE (default) the function returns
#'   a MizerParams object with the "initial" slots set to the steady state.
#' @param progress_bar A shiny progress object to implement
#'   a progress bar in a shiny app. Default FALSE.
#' @export
#' @md
#' @examples
#' \dontrun{
#' params <- newTraitParams()
#' params@species_params$gamma[5] <- 3000
#' params <- setSearchVolume(params)
#' params <- steady(params)
#' }
steady <- function(params, t_max = 100, t_per = 7.5, tol = 10^(-2),
                   dt = 0.1, return_sim = FALSE, progress_bar = TRUE) {
    assert_that(is(params, "MizerParams"),
                noNA(getRDD(params)))
    p <- params
    
    if (is(progress_bar, "Progress")) {
        # We have been passed a shiny progress object
        progress_bar$set(message = "Finding steady state", value = 0)
        proginc <- 1/ceiling(t_max/t_per)
    }
    
    # Force the recruitment to stay at the current level
    p@species_params$constant_recruitment <- getRDD(p)
    p@srr <- "srrConstant"
    old_rdi <- getRDI(p)
    rdi_limit <- old_rdi / 1e7
    # Force other componens to stay at current level
    old_other_dynamics <- p@other_dynamics
    for (res in names(p@other_dynamics)) {
        p@other_dynamics[[res]] <- constant_other(res)
    }
    
    n <- p@initial_n
    n_pp <- p@initial_n_pp
    n_other <- p@initial_n_other
    sim <- p
    for (ti in (1:ceiling(t_max/t_per))) {
        # advance shiny progress bar
        if (is(progress_bar, "Progress")) {
            progress_bar$inc(amount = proginc)
        }
        if (return_sim) {
            sim <- project(sim, dt = dt, t_max = t_per, t_save = t_per,
                           initial_n = n, initial_n_pp = n_pp, 
                           initial_n_other = n_other)
        } else {
            sim <- project(p, dt = dt, t_max = t_per, t_save = t_per,
                           initial_n = n, initial_n_pp = n_pp, 
                           initial_n_other = n_other)
        }
        no_t <- dim(sim@n)[1]
        n[] <- sim@n[no_t, , ]
        n_pp[] <- sim@n_pp[no_t, ]
        n_other <- sim@n_other[[no_t]]
        new_rdi <- getRDI(p, n, n_pp, n_other)
        deviation <- max(abs((new_rdi - old_rdi)/old_rdi))
        if (any(new_rdi < rdi_limit)) {
            if (return_sim) {
                message("One of the species is going extinct.")
                break
            }
            extinct <- p@species_params$species[new_rdi < rdi_limit]
            stop(paste(extinct, collapse = ", "),
                 " are going extinct.")
        }
        if (deviation < tol) {
            break
        }
        old_rdi <- new_rdi
    }
    if (deviation >= tol) {
        warning("Simulation run in steady() did not converge after ", 
                ti * t_per,
                " years. Residual relative rate of change = ", deviation)
    } else {
        message("Steady state was reached before ", ti * t_per, " years.")
    }
    
    # Restore original stock-recruitment relationship and other dynamics
    p@srr <- params@srr
    p@other_dynamics <- old_other_dynamics
    
    no_sp <- length(p@species_params$species)
    p@initial_n[] <- n
    p@initial_n_pp[] <- n_pp
    p@initial_n_other[] <- n_other

    # Retune the values of erepro so that we get the correct level of
    # recruitment
    p <- retuneReproductionEfficiency(p)
    
    if (return_sim) {
        sim@params <- p
        return(sim)
    } else {
        return(p)
    }
}

# Helper function to calculate the coefficient of the death rate created by
# a Sheldon spectrum of predators, assuming they have the same predation 
# parameters as the first species.
get_power_law_mort <- function(params) {
    params@interaction[] <- 0
    params@interaction[1, 1] <- 1
    params@initial_n[1, ] <- params@plankton_params$kappa * 
        params@w^(-params@plankton_params$lambda)
    return(getPredMort(params)[1, 1] / 
               params@w[[1]] ^ (1 + params@species_params$q[[1]] - 
                                    params@plankton_params$lambda))
}


# Helper function to keep other components constant
constant_other <- function(other_name) {
    force(other_name)
    function(params, n, n_pp, n_other, rates, t, dt, ...) n_other[other_name]
}
