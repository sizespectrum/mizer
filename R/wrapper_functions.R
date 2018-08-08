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
#'   \code{\link{set_community_model}} \tab Community model \tab 5 \cr
#'   \code{\link{set_trait_model}} \tab Trait-based model \tab 6 \cr
#'   \code{\link{set_scaling_model}} \tab Scale-invariant Trait-based model 
#'       \tab 7 \cr
#' }
#'
#'The file also contains a helper function \code{\link{retune_abundance}}.
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

#' Sets up parameters for a community-type model
#' 
#' This functions creates a \code{\linkS4class{MizerParams}} object so that
#' community-type models can be easily set up and run. A community model has
#' several features that distinguish it from the food-web type models. Only one
#' 'species' is resolved, i.e. one 'species' is used to represent the whole
#' community. The resource spectrum only extends to the start of the community
#' spectrum. Recruitment to the smallest size in the community spectrum is
#' constant and set by the user. As recruitment is constant, the proportion of
#' energy invested in reproduction (the slot \code{psi} of the returned 
#' \code{MizerParams} object) is set to 0. Standard metabolism has been turned 
#' off (the parameter \code{ks} is set to 0). Consequently, the growth rate is 
#' now determined solely by the assimilated food (see the package vignette for 
#' more details).
#' 
#' The function has many arguments, all of which have default values. The main 
#' arguments that the users should be concerned with are \code{z0}, 
#' \code{recruitment}, \code{alpha} and \code{f0} as these determine the average
#' growth rate of the community.
#' 
#' Fishing selectivity is modelled as a knife-edge function with one parameter, 
#' \code{knife_edge_size}, which determines the size at which species are 
#' selected.
#' 
#' The resulting \code{MizerParams} object can be projected forward using 
#' \code{project()} like any other \code{MizerParams} object. When projecting 
#' the community model it may be necessary to reduce \code{dt} to 0.1 to avoid 
#' any instabilities with the solver. You can check this by plotting the biomass
#' or abundance through time after the projection.
#' 
#' @param z0 The background mortality of the community. Default value is 0.1.
#' @param alpha The assimilation efficiency of the community. Default value 0.2
#' @param f0 The average feeding level of individuals who feed mainly on the 
#'   resource. This value is used to calculate the search rate parameter 
#'   \code{gamma} (see the package vignette). Default value is 0.7.
#' @param h The maximum food intake rate. Default value is 10.
#' @param beta The preferred predator prey mass ratio. Default value is 100.
#' @param sigma The width of the prey preference. Default value is 2.0.
#' @param q The search volume exponent. Default value is 0.8.
#' @param n The scaling of the intake. Default value is 2/3.
#' @param kappa The carrying capacity of the plankton spectrum. Default value
#'   is 1000.
#' @param lambda The exponent of the plankton spectrum. Default value is 2 + q
#'   - n.
#' @param r_pp Growth rate of the primary productivity. Default value is 10.
#' @param gamma Volumetric search rate. Estimated using \code{h}, \code{f0} and 
#'   \code{kappa} if not supplied.
#' @param recruitment The constant recruitment in the smallest size class of the
#'   community spectrum. This should be set so that the community spectrum 
#'   continues the plankton spectrum. Default value = \code{kappa} * 
#'   \code{min_w}^-\code{lambda}.
#' @param rec_mult Additional multiplier for the constant recruitment. Default 
#'   value is 1.
#' @param knife_edge_size The size at the edge of the knife-selectivity 
#'   function. Default value is 1000.
#' @param knife_is_min Is the knife-edge selectivity function selecting above 
#'   (TRUE) or below (FALSE) the edge. Default is TRUE.
#' @param max_w The maximum size of the community. The \code{w_inf} of the 
#'   species used to represent the community is set to 0.9 * this value. The 
#'   default value is 1e6.
#' @param min_w The minimum size of the community. Default value is 1e-3.
#' @param ... Other arguments to pass to the \code{MizerParams} constructor.
#' @export
#' @return An object of type \code{\linkS4class{MizerParams}}
#' @references K. H. Andersen,J. E. Beyer and P. Lundberg, 2009, Trophic and 
#'   individual efficiencies of size-structured communities, Proceedings of the 
#'   Royal Society, 276, 109-114
#' @examples
#' \dontrun{
#' params <- set_community_model(f0=0.7, z0=0.2, recruitment=3e7)
#' sim <- project(params, effort = 0, t_max = 100, dt=0.1)
#' plotBiomass(sim)
#' plotSpectra(sim)
#' }
set_community_model <- function(max_w = 1e6,
                                min_w = 1e-3,
                                z0 = 0.1,
                                alpha = 0.2,
                                h = 10,
                                beta = 100,
                                sigma = 2.0,
                                q = 0.8,
                                n = 2/3,
                                kappa = 1000,
                                lambda = 2 + q - n,
                                f0 = 0.7,
                                r_pp = 10,
                                gamma = NA,
                                knife_edge_size = 1000,
                                knife_is_min = TRUE,
                                recruitment = kappa * min_w^-lambda,
                                rec_mult = 1,
                                ...) {
    w_inf <- max_w * 0.9
    w_pp_cutoff <- min_w
    ks <- 0 # Turn off standard metabolism
    p <- n # But not used as ks = 0
    # Estimate gamma if not supplied
    if (is.na(gamma)) {
        gamma <- (f0 * h * beta^(2 - lambda)) /
            ((1 - f0) * sqrt(2 * pi) * kappa * sigma)
    }
    # Make the species data.frame
    com_params_df <- data.frame(
        species = "Community",
        w_inf = w_inf,
        w_mat = 1e12, # Has no affect as psi set to 0 but we set it to something 
                      # to help the constructor
        h = h, # max food intake
        gamma = gamma,# vol. search rate,
        ks = ks,# standard metabolism coefficient,
        beta = beta,
        sigma = sigma,
        z0 = z0, # background mortality
        alpha = alpha,
        erepro = 1, # not used
        sel_func = "knife_edge",
        knife_edge_size = knife_edge_size,
        knife_is_min = knife_is_min,
        constant_recruitment = recruitment * rec_mult # to be used in the SRR
    )
    # Set the recruitment function for constant recruitment
    constant_recruitment <- function(rdi, species_params){
        return(species_params$constant_recruitment)
    }
    com_params <- MizerParams(com_params_df, p = p, n = n, q = q, lambda = lambda, 
                              kappa = kappa, min_w = min_w, max_w = max_w, 
                              w_pp_cutoff = w_pp_cutoff, r_pp = r_pp, ...)
    com_params@srr <- constant_recruitment
    com_params@psi[] <- 0 # Need to force to be 0. Can try setting w_mat but 
                          # due to slope still not 0
    # Set w_mat to NA for clarity - it is not actually being used
    com_params@species_params$w_mat[] <- NA
    return(com_params)
}


#' Sets up parameters for a trait-based model
#' 
#' This functions creates a \code{MizerParams} object so that trait-based-type 
#' models can be easily set up and run. The trait-based size spectrum model can
#' be derived as a simplification of the general size-based model used in
#' \code{mizer}. All the species-specific parameters are the same, except for
#' the asymptotic size, which is considered the most important trait
#' characterizing a species. Other parameters are related to the asymptotic
#' size. For example, the size at maturity is given by w_inf * eta, where eta is
#' the same for all species. For the trait-based model the number of species is
#' not important. For applications of the trait-based model see Andersen &
#' Pedersen (2010). See the \code{mizer} vignette for more details and examples
#' of the trait-based model.
#' 
#' The function has many arguments, all of which have default values. Of
#' particular interest to the user are the number of species in the model and
#' the minimum and maximum asymptotic sizes. The asymptotic sizes of the species
#' are spread evenly on a logarithmic scale within this range.
#' 
#' The stock recruitment relationship is the default Beverton-Holt style. The
#' maximum recruitment is calculated using equilibrium theory (see Andersen &
#' Pedersen, 2010) and a multiplier, \code{k0}. Users should adjust \code{k0} to
#' get the spectra they want.
#' 
#' The factor for the search volume, \code{gamma}, is calculated using the
#' expected feeding level, \code{f0}.
#' 
#' Fishing selectivity is modelled as a knife-edge function with one parameter,
#' \code{knife_edge_size}, which is the size at which species are selected. Each
#' species can either be fished by the same gear (\code{knife_edge_size} has a
#' length of 1) or by a different gear (the length of \code{knife_edge_size} has
#' the same length as the number of species and the order of selectivity size is
#' that of the asymptotic size).
#' 
#' The resulting \code{MizerParams} object can be projected forward using
#' \code{project()} like any other \code{MizerParams} object. When projecting
#' the community model it may be necessary to reduce \code{dt} to 0.1 to avoid
#' any instabilities with the solver. You can check this by plotting the biomass
#' or abundance through time after the projection.
#' @param no_sp The number of species in the model. The default value is 10. The
#'   more species, the longer takes to run.
#' @param min_w_inf The asymptotic size of the smallest species in the
#'   community.
#' @param max_w_inf The asymptotic size of the largest species in the community.
#' @param no_w The number of size bins in the community spectrum.
#' @param min_w The smallest size of the community spectrum.
#' @param max_w The largest size of the community spectrum. Default value is the
#'   largest w_inf in the community x 1.1.
#' @param min_w_pp The smallest size of the plankton spectrum.
#' @param no_w_pp Obsolete argument that is no longer used because the number
#'    of plankton size bins is determined because all size bins have to
#'    be logarithmically equally spaced.
#' @param w_pp_cutoff The cut off size of the plankton spectrum. Default value
#'   is 1.
#' @param k0 Multiplier for the maximum recruitment. Default value is 50.
#' @param n Scaling of the intake. Default value is 2/3.
#' @param p Scaling of the standard metabolism. Default value is 0.75.
#' @param q Exponent of the search volume. Default value is 0.9.
#' @param eta Factor to calculate \code{w_mat} from asymptotic size.
#' @param r_pp Growth rate of the primary productivity. Default value is 4.
#' @param kappa Coefficient in abundance power law. Default value is
#'   0.005.
#' @param lambda Exponent of the abundance power law. Default value is (2+q-n).
#' @param alpha The assimilation efficiency of the community. The default value
#'   is 0.6
#' @param ks Standard metabolism coefficient. Default value is 4.
#' @param z0pre The coefficient of the background mortality of the community. z0
#'   = z0pre * w_inf ^ (n-1). The default value is 0.6.
#' @param h Maximum food intake rate. Default value is 30.
#' @param beta Preferred predator prey mass ratio. Default value is 100.
#' @param sigma Width of prey size preference. Default value is 1.3.
#' @param f0 Expected average feeding level. Used to set \code{gamma}, the
#'   factor for the search volume. The default value is 0.5.
#' @param gamma Volumetric search rate. Estimated using \code{h}, \code{f0} and
#'   \code{kappa} if not supplied.
#' @param knife_edge_size The minimum size at which the gear or gears select
#'   species. Must be of length 1 or no_sp.
#' @param gear_names The names of the fishing gears. A character vector, the
#'   same length as the number of species. Default is 1 - no_sp.
#' @param ... Other arguments to pass to the \code{MizerParams} constructor.
#' @export
#' @return An object of type \code{MizerParams}
#' @seealso \linkS4class{MizerParams}
#' @references K. H. Andersen and M. Pedersen, 2010, Damped trophic cascades
#'   driven by fishing in model marine ecosystems. Proceedings of the Royal
#'   Society V, Biological Sciences, 1682, 795-802.
#' @examples
#' \dontrun{
#' trait_params <- set_trait_model(no_sp = 15)
#' init_pop <- get_initial_n(trait_params, n0_mult = 0.001)
#' sim <- project(trait_params, effort = 0, t_max = 50, dt=0.2,
#'     initial_n = init_pop, t_save = 1)
#' plot(sim)
#' ## Set up industrial fishery that only fishes on species with w_inf <= 500 g
#' ## And where the selectivity of the industrial fishery = w_inf * 0.05
#' no_sp <- 10
#' min_w_inf <- 10
#' max_w_inf <- 1e5
#' w_inf <- 10^seq(from=log10(min_w_inf), to = log10(max_w_inf), length=no_sp)
#' knife_edges <- w_inf * 0.05
#' industrial_gears <- w_inf <= 500
#' other_gears <- w_inf > 500
#' gear_names <- rep("Industrial", no_sp)
#' gear_names[other_gears] <- "Other"
#' params_gear <- set_trait_model(no_sp = no_sp, min_w_inf = min_w_inf,
#'     max_w_inf = max_w_inf, knife_edge_size = knife_edges, 
#'     gear_names = gear_names)
#' ## Only turn on Industrial fishery. Set effort of the Other gear to 0
#' sim <- project(params_gear, t_max = 20, effort = c(Industrial = 1, Other = 0))
#' }
set_trait_model <- function(no_sp = 10,
                            min_w_inf = 10,
                            max_w_inf = 1e5,
                            no_w = 100,
                            min_w = 0.001,
                            max_w = max_w_inf * 1.1,
                            min_w_pp = 1e-10,
                            no_w_pp = NA,
                            w_pp_cutoff = 1,
                            k0 = 50, # recruitment adjustment parameter
                            n = 2/3,
                            p = 0.75,
                            q = 0.9, 
                            eta = 0.25,
                            r_pp = 4,
                            kappa = 0.005,
                            lambda = 2+q-n,
                            alpha = 0.6,
                            ks = 4,
                            z0pre = 0.6,
                            h = 30,
                            beta = 100,
                            sigma = 1.3,
                            f0 = 0.5,
                            gamma = NA,
                            knife_edge_size = 1000,
                            gear_names = "knife_edge_gear",
                            ...){
    if (!is.na(no_w_pp))
        warning("New mizer code does not support the parameter no_w_pp")
    
    w_inf <- 10^seq(from = log10(min_w_inf), to = log10(max_w_inf), length = no_sp)
    w_mat <- w_inf * eta
    
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
    
    # Make the species parameters data.frame
    trait_params_df <- data.frame(
        species = as.factor(1:no_sp),
        w_inf = w_inf,
        w_mat = w_mat,
        h = h, # max food intake
        gamma = gamma, # vol. search rate,
        ks = ks,# standard metabolism coefficient,
        beta = beta,
        sigma = sigma,
        z0 = z0pre * w_inf^(n - 1), # background mortality
        alpha = alpha,
        #r_max = r_max,
        sel_func = "knife_edge",
        knife_edge_size = knife_edge_size,
        gear = gear_names,
        erepro = 1
    )
    # Make the MizerParams
    trait_params <-
        MizerParams(
            trait_params_df,
            min_w = min_w,
            max_w = max_w,
            no_w = no_w,
            min_w_pp = min_w_pp,
            w_pp_cutoff = w_pp_cutoff,
            n = n,
            p = p,
            q = q,
            r_pp = r_pp,
            kappa = kappa,
            lambda = lambda
        ) 
    # Sort out maximum recruitment - see A&P 2009 Get max flux at recruitment
    # boundary, R_max R -> | -> g0 N0 R is egg flux, in numbers per time Actual
    # flux at recruitment boundary = RDD = NDD * g0 (where g0 is growth rate) So
    # in our BH SRR we need R_max comparable to RDI (to get RDD) R_max = N0_max
    # * g0 (g0 is the average growth rate of smallest size, i.e. at f0 = 0.5) N0
    # given by Appendix A of A&P 2010 - see Ken's email 12/08/13 Taken from
    # Ken's code 12/08/13 - equation in paper is wrong!
    alpha_p <- f0 * h * beta^(2 * n - q - 1) * 
        exp((2 * n * (q - 1) - q^2 + 1) * sigma^2 / 2)
    alpha_rec <- alpha_p / (alpha * h * f0 - ks)
    # Calculating dw using Ken's code - see Ken's email 12/08/13
    tmpA <- w_inf[1]
    tmpB <- (log10(w_inf[length(w_inf)]) - log10(w_inf[1])) / (no_sp - 1) # Difference between logged w_infs, fine
    dw_winf <- tmpB * tmpA * 10^(tmpB*( (1:no_sp) - 1)) # ?
    N0_max <- k0 * w_inf^(n*2-q-3+alpha_rec) * dw_winf  # Why * dw_winf, not / ? Ken confirms * in email
    # No need to include (1 - psi) in growth equation because allocation to reproduction at this size = 0, so 1 - psi = 1
    g0 <- (alpha * f0 * h * trait_params@w[1]^n - ks * trait_params@w[1]^p)
    r_max <- N0_max * g0
    
    trait_params@species_params$r_max <- r_max
    
    return(trait_params)
}


#' Sets up parameters for a scale free trait-based model
#'
#' This functions creates a \code{MizerParams} object so that scale free
#' trait-based-type models can be easily set up and run. The scale free
#' trait-based size spectrum model can be derived as a simplification of the
#' general size-based model used in \code{mizer}. All the species-specific
#' parameters are the same for all species, except for the egg size, maturity
#' size and asymptotic size. These differ over the species, but the ratio of egg
#' size to maturity size and the ratio of egg size to asymptotic size are the
#' same for each species. The asymptotic sizes of the species are spread evenly
#' on a logarithmic scale. See the \code{mizer} vignette and the Details section
#' below for more details and examples of the scale free trait-based model.
#' 
#' The scale free trait-based model is similar to the standard trait-based
#' model, with three main differences:
#' \enumerate{
#' \item We have an exact equation for a steady state of this system which is
#' often stable, even when we include no extra stabilization effects like
#' density dependence or stock recruitment relationships.
#' \item The egg size is proportional to the maturity size for each species
#' \item The parameters are chosen so that R_0 (the expected number of offspring
#' produced by an individual over a lifetime) is close to 1 for each species.
#' }
#'
#' The function has many arguments, all of which have default values. Of
#' particular interest to the user are the number of species in the model and
#' the minimum and maximum asymptotic sizes.
#'
#' The characteristic weights of the different species are defined by
#' min_egg, min_w_mat, min_w_inf, max_w_inf and no_sp, in the sense that the egg
#' weights of the no_sp species are logarithmically evenly spaced, ranging from
#' min_w=min_egg to max_w=max_w_inf. The maturity weights of the species can be
#' obtained by muliplying the egg_weights by min_w_mat/min_egg. The asymptotic
#' weights of the species can be obtained by multiplying the egg weights by
#' min_w_inf/min_egg.
#'
#' Although the scale free trait based model's default steady state is often
#' stable without imposing a stock recruitment relationship, the function can
#' set a Beverton-Holt type stock recruitment relationship that imposes a
#' maximal reproduction rate that is a multiple of the recruitment rate at
#' steady state. That multiple is set by the argument \code{rfac}.
#'
#' In addition to setting up the parameters, this function also evaluates the
#' analytic expression for a steady state of the scale free trait-based model
#' and sets it as the initial condition.
#'
#' The search rate coefficient \code{gamma} is calculated using the expected
#' feeding level, \code{f0}.
#'
#' The option of including fishing is given, but the steady state may lose its
#' natural stability if too much fishing is included. In such a case the user
#' may wish to include stablizing effects (like Rmax and chi) to ensure the
#' steady state is stable. Fishing selectivity is modelled as a knife-edge
#' function with one parameter, \code{knife_edge_size}, which is the size at
#' which species are selected. Each species can either be fished by the same
#' gear (\code{knife_edge_size} has a length of 1) or by a different gear (the
#' length of \code{knife_edge_size} has the same length as the number of species
#' and the order of selectivity size is that of the asymptotic size).
#'
#' The resulting \code{MizerParams} object can be projected forward using
#' \code{project()} like any other \code{MizerParams} object. When projecting
#' the model it may be necessary to reduce \code{dt} to 0.1 to avoid any
#' instabilities with the solver. You can check this by plotting the biomass or
#' abundance through time after the projection.
#'
#' @param no_sp The number of species in the model. The default value is 11.
#' @param min_w_inf The asymptotic size of the smallest species in the
#'   community. Default value is 10.
#' @param max_w_inf The asymptotic size of the largest species in the community.
#'   Default value is 1000.
#' @param min_egg The size of the the egg of the smallest species. Default value
#'   is 10^(-4).
#' @param min_w_mat The maturity size of the smallest species. Default value is
#'   10^(0.4),
#' @param no_w The number of size bins in the community spectrum. Default value
#'   is such that there are 100 bins for each factor of 10 in weight.
#' @param min_w_pp The smallest size of the plankton spectrum. Default value
#'   is min_egg/(beta*exp(5*sigma)) so that it covers the entire range of the
#'   feeding kernel of even the smallest fish larva.
#' @param w_pp_cutoff The largest size of the plankton spectrum. Default
#'   value is max_w_inf unless \code{perfect = TRUE} when it is Inf.
#' @param n Scaling of the intake. Default value is 2/3.
#' @param q Exponent of the search volume. Default value is 3/4 unless 
#'   \code{lambda} is provided, in which case this argument is ignored and
#'   q = lambda - 2 + n.
#' @param lambda Exponent of the abundance power law. If supplied, this 
#'   overrrules the \code{q} argument. Otherwise the default value is 2+q-n.
#' @param r_pp Growth rate of the primary productivity. Default value is 0.1.
#' @param kappa Coefficient in abundance power law. Default value is
#'   0.005.
#' @param alpha The assimilation efficiency of the community. The default value
#'   is 0.4.
#' @param ks Standard metabolism coefficient. Default value is 4.
#' @param h Maximum food intake rate. Default value is 30.
#' @param beta Preferred predator prey mass ratio. Default value is 100.
#' @param sigma Width of prey size preference. Default value is 1.3.
#' @param f0 Expected average feeding level. Used to set \code{gamma}, the
#'   coefficient in the search rate. The default value is 0.6.
#' @param knife_edge_size The minimum size at which the gear or gears select
#'   species. Must be of length 1 or no_sp. Default value is 100.
#' @param gear_names The names of the fishing gears. A character vector, the
#'   same length as the knife_edge_size parameter. Default value is
#'   "knife_edge_gear".
#' @param rfac The factor such that Rmax = rfac * R, where Rmax is the maximum
#'   recruitment allowed and R is the steady-state recruitment. Thus the larger
#'   \code{rfac} the less the impact of the non-linear stock-recruitment curve.
#'   The default is Inf.
#' @param perfect Boolean. Default FALSE. If TRUE then parameters are set so
#'   that the community abundance, growth before reproduction and death are
#'   perfect power laws.
#' @param ... Other arguments to pass to the \code{MizerParams} constructor.
#' @export
#' @return An object of type \code{MizerParams}
#' @seealso \linkS4class{MizerParams}
#' @examples
#' \dontrun{
#' s_params <- set_scaling_model()
#' sim <- project(s_params, t_max=5, effort = 0)
#' plotSpectra(sim)
#' }
set_scaling_model <- function(no_sp = 11,
                              min_w_inf = 10,
                              max_w_inf = 10 ^ 3,
                              min_egg = 10 ^ (-4),
                              min_w_mat = 10 ^ (0.4),
                              no_w = log10(max_w_inf / min_egg) * 100 + 1,
                              min_w_pp = min_egg / (beta * exp(5 * sigma)),
                              w_pp_cutoff = min_w_inf,
                              n = 2 / 3,
                              q = 3 / 4,
                              lambda = 2 + q - n,
                              r_pp = 0.1,
                              kappa = 0.005,
                              alpha = 0.4,
                              ks = 4,
                              h = 30,
                              beta = 100,
                              sigma = 1.3,
                              f0 = 0.6,
                              knife_edge_size = 100,
                              gear_names = "knife_edge_gear",
                              rfac = Inf,
                              perfect = FALSE,
                              ...) {
    if (hasArg(lambda)) {
        # The lambda argument overrules any q argument
        q <- lambda - 2 + n
    }
    # check validity of parameters
    if (rfac <= 1) {
        message("rfac needs to be larger than 1. Setting rfac=1.01")
        rfac <- 1.01
    }
    no_w <- round(no_w)
    if (no_w < 1) {
        stop("The number of size bins no_w must be a positive integer")
    }
    if (no_w < log10(max_w_inf/min_egg)*5) {
        no_w <- round(log10(max_w_inf / min_egg) * 5 + 1)
        message(paste("Increased no_w to", no_w, "so that there are 5 bins for an interval from w and 10w."))
    }
    if (no_w > 10000) {
        message("Running a simulation with", no_w, "size bins is going to be very slow.")
    }
    if (min_w_inf >= max_w_inf) {
        stop("The asymptotic size of the smallest species min_w_inf must be smaller than the asymptotic size of the largest species max_w_inf")
    }
    if (min_egg >= min_w_mat) {
        stop("The egg size of the smallest species min_egg must be smaller than its maturity size min_w_mat")
    }
    if (min_w_mat >= min_w_inf) {
        stop("The maturity size of the smallest species min_w_mat must be smaller than its maximum size min_w_inf")
    }
    no_sp <- as.integer(no_sp)
    if (no_sp < 2) {
        stop("The number of species must be at least 2.")
    }
    if (!all(c(n, q, r_pp, kappa, alpha, h, beta, sigma, ks, f0, knife_edge_size) > 0)) {
        stop("The parameters n, q, r_pp, kappa, alpha, h, beta, sigma, ks, f0 and knife_edge_size, if supplied, need to be positive.")
    }
    
    if (perfect) {
        w_pp_cutoff <- Inf
    }
    # Set exponents
    p <- n
    lambda <- 2 + q - n
    # Set grid points and characteristic sizes
    min_w <- min_egg
    max_w <- max_w_inf
    # min_egg and max_w already lie on grid points in w. 
    # Round min_w_mat up to the nearest grid point.
    delt <- (log10(max_w) - log10(min_w)) / (no_w - 1)
    v <- min_w_mat
    j <- 1 + ceiling((log10(v) - log10(min_w)) / delt)
    v <- 10 ^ (log10(min_w) + (j - 1) * delt)
    min_w_mat <- v
    # Round min_w_inf so that it is an integer multiple of the
    # species spacing away from max_w_inf
    j <- round((log10(max_w) - log10(min_w_inf)) / (delt * (no_sp - 1)))
    min_w_inf <- 10 ^ (log10(max_w) - j * (no_sp - 1) * delt)
    w_min_idx <- seq(1, by = j, length.out = no_sp)
    # Determine maximum egg size
    max_egg <- max_w * min_egg / min_w_inf
    log10_minimum_egg <- log10(min_egg)
    log10_maximum_egg <- log10(max_egg)
    # Determine logarithmic spacing of egg weights
    dist_sp <- (log10_maximum_egg - log10_minimum_egg) / (no_sp - 1)
    # Determine egg weights w_min for all species
    x_min <- seq(log10_minimum_egg, by = dist_sp, length.out = no_sp)
    w_min <- 10 ^ x_min
    # Use ratios to determine w_inf and w_mat from w_min
    w_inf <- w_min * min_w_inf / min_egg
    w_mat <- w_min * min_w_mat / min_egg
    # Build Params Object
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
        # not used but required
        knife_edge_size = knife_edge_size,
        gear = gear_names
    )
    params <-
        MizerParams(
            species_params,
            p = p,
            n = n,
            q = q,
            lambda = lambda,
            f0 = f0,
            kappa = kappa,
            min_w = min_w,
            max_w = max_w,
            no_w = no_w,
            min_w_pp = min_w_pp,
            w_pp_cutoff = max_w,
            r_pp = r_pp
        )
    # gamma is determined by MizerParams
    gamma <- params@species_params$gamma[1]
    w <- params@w
    dw <- params@dw
    
    # Get constants for steady-state solution
    mu0 <- (1 - f0) * sqrt(2 * pi) * kappa * gamma * sigma *
        (beta ^ (n - 1)) * exp(sigma ^ 2 * (n - 1) ^ 2 / 2)
    hbar <- alpha * h * f0 - ks
    if (hbar < 0) {
        stop("The feeding level is not sufficient to maintain the fish.")
    }
    pow <- mu0 / hbar / (1 - n)
    if (pow < 1) {
        message("The ratio of death rate to growth rate is too small, leading to
                an accumulation of fish at their largest size.")
    }
    
    # Create steady state solution n_exact for species 1
    # The following would calculate the analytic solution in the case
    # of the simplified discontinuous allocation to reproduction.
    # n_mult <- (1 - (w / w_inf[1]) ^ (1 - n)) ^ (pow - 1) * 
    #               (1 - (w_mat[1] / w_inf[1]) ^(1 - n)) ^ (-pow)
    # n_mult[w < w_mat[1]] <- 1
    # n_mult[w >= w_inf[1]] <- 0
    # n_exact <- ((w_min[1] / w) ^ (mu0 / hbar) / (hbar * w ^ n)) * n_mult
    # n_exact <- n_exact[w >= w_min[1] & w < w_inf[1]]
    # We instead evaluate the integral in the analytic solution numerically
    mumu <- mu0 * w^(n - 1)  # Death rate
    gg <- hbar * w^n * (1-params@psi[1, ])  # Growth rate
    
    w_inf_idx <- sum(w < w_inf[1])
    idx <- 1:(w_inf_idx - 1)
    # Compute integral in analytic solution to McKvF
    # We do not use this because it is not in agreement with scheme in project
    # w_inf_idx <- sum(w < w_inf[1])
    # idx <- 1:(w_inf_idx-1)
    # integrand <- dw[idx] * mumu[idx] / gg[idx]
    # n_exact <- exp(-cumsum(integrand)) / gg[idx]
    
    # Steady state solution of the upwind-difference scheme used in project
    n_exact <- c(1, cumprod(gg[idx] / ((gg + mumu * dw)[idx + 1])))
    
    # rescale fish abundance to line up with plankton spectrum
    mult <- kappa / 
        sum(n_exact * (w^(lambda - 1) * dw)[1:w_inf_idx])
    n_exact <- n_exact * mult * 
        (10^(dist_sp*(1-lambda)/2) - 10^(-dist_sp*(1-lambda)/2)) / (1-lambda)
    
    # Use n_exact as a template to create solution initial_n for all species
    initial_n <- params@psi  # get array with correct dimensions and names
    initial_n[, ] <- 0
    for (i in 1:no_sp) {
        # smallest index for species
        w_min_idx <- params@w_min_idx[i]
        # range of indices
        idxs <- w_min_idx:(w_min_idx + length(n_exact) - 1)  
        initial_n[i, idxs] <- n_exact * (w_min[1] / w_min[i]) ^ lambda
    }
    
    # Calculate the community spectrum
    sc <- colSums(initial_n)
    # The following was an attempt to calculate the community for a
    # finer-grained spectrum where species are spaced by just one
    # weight bracket. This level of detail is however not necessary.
    # Also in the below the start and end of the integration is not
    # quite right.
    # sc <- rep(0, no_w)
    # idxs <- seq_along(n_exact)
    # fac <- (w[1]/w[2])^lambda
    # for (i in seq_len(max(params@w_min_idx) - 1)) {
    #     sc[idxs] <- sc[idxs] + n_exact * fac^(i-1)
    #     idxs <- idxs+1
    # }
    # sc <- sc * 
    #     (10^(delt*(1-lambda)/2) - 10^(-delt*(1-lambda)/2)) /
    #     (10^(dist_sp*(1-lambda)/2) - 10^(-dist_sp*(1-lambda)/2))
    params@sc <- sc
    
    # Setup plankton
    plankton_vec <- (kappa * w ^ (-lambda)) - sc
    # Cut off plankton at w_pp_cutoff
    plankton_vec[w >= w_pp_cutoff] <- 0
    if (!perfect && any(plankton_vec < 0)) {
        # Do not allow negative plankton abundance
        message("Note: Negative plankton abundance values overwritten with zeros")
        plankton_vec[plankton_vec < 0] <- 0
    }
    # The cc_pp factor needs to be higher than the desired steady state in
    # order to compensate for predation mortality
    params@cc_pp[sum(params@w_full <= w[1]):length(params@cc_pp)] <-
        plankton_vec
    initial_n_pp <- params@cc_pp
    m2_background <- getPlanktonMort(params, initial_n, initial_n_pp)
    params@cc_pp <- (params@rr_pp + m2_background ) * initial_n_pp/params@rr_pp
    
    # Setup background death
    m2 <- getPredMort(params, initial_n, initial_n_pp)
    for (i in 1:no_sp) {
        # The steplike psi was only needed when we wanted to use the analytic
        # expression for the steady-state solution
        # params@psi[i,] <- (w / w_inf[i]) ^ (1 - n)
        # params@psi[i, w < (w_mat[i] - 1e-10)] <- 0
        # params@psi[i, w > (w_inf[i] - 1e-10)] <- 1
        params@mu_b[i,] <- mu0 * w ^ (n - 1) - m2[i, ]
        if (!perfect && any(params@mu_b[i,] < 0)) {
            message("Note: Negative background mortality rates overwritten with zeros")
            params@mu_b[i, params@mu_b[i,] < 0] <- 0
        }
    }
    # Set erepro to meet boundary condition
    rdi <- getRDI(params, initial_n, initial_n_pp)
    gg <- getEGrowth(params, initial_n, initial_n_pp)
    mumu <- getMort(params, initial_n, initial_n_pp, effort = 0)
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
    # Record abundance of fish and resources at steady state, as slots.
    params@initial_n <- initial_n
    params@initial_n_pp <- initial_n_pp
    # set rmax=fac*RDD
    # note that erepro has been multiplied by a factor of (rfac/(rfac-1)) to
    # compensate for using a stock recruitment relationship.
    params@species_params$r_max <-
        (rfac - 1) * getRDI(params, initial_n, initial_n_pp)
    return(params)
}


#' Retunes abundance of background species.
#' 
#' An unexported helper function.
#'
#' If N_i(w) is a steady state of the McKendrik-von Foerster (MVF) equation with
#' fixed growth and death rates, then A_i*N_i(w) is also a steady state, where
#' A_i is an abundance multiplier. When we add a foreground species to our
#' model, we want to choose new abundance multipliers of the background species
#' so that the community abundance after adding the new species is close to
#' the original community abundance, stored in \code{params@sc}.
#'
#' @param params A \linkS4class{MizerParams} object
#' @param retune A boolean vector that determines whether a species can be 
#'   retuned or not.
#' @param cutoff Species with an abundance at maturity size that is less than 
#'               cutoff times community abundance will be removed. Default 1e-3.
#'   
#' @return An object of type \code{MizerParams}
#' @seealso \linkS4class{MizerParams}
retune_abundance <- function(params, retune, cutoff = 1e-3) {
    no_sp <- length(params@species_params$species)  # Number of species
    if (length(retune) != no_sp) {
        stop("retune argument has the wrong length")
    }
    if (!any(retune)) {
        # nothing to retune
        return(params)
    }
    # We try to match the original abundance between the maturity size
    # of the smallest species and the maximum size of the largest species.
    # Determine the indices of these limits
    idx_start <- sum(params@w <= min(params@species_params$w_mat))
    idx_stop <- sum(params@w < max(params@species_params$w_inf))
    # More precisely, we find the abundance multipliers A_i so
    # that the integral of the square of the relative distance 
    # (sum_{i not in L} A_i*N_i(w) + sum_{i not in L} N_i(w) - sc(w))/sc(w) 
    # over w, between our limits, is minimized, where  L is the set of all
    # retuneable species.
    
    # deal with zero entries in params@sc
    nonzero <- params@sc > 0
    sc <- params@sc[nonzero]
    # rho is the total abundance of all the non-tunable species
    rho <- colSums(params@initial_n[!retune, nonzero, drop = FALSE])
    
    # Use Singular Value Decomposition to find optimal abundance multipliers.
    # See Numerical Recipes section 15.4.2
    #
    # Rescale by sc
    A <- t(sweep(params@initial_n[retune, nonzero, drop = FALSE], 2, sc, "/"))
    b <- (sc - rho) / sc
    
    sv <- svd(A)
    di <- 1/sv$d  # inverse of singular values
    di[di > 10^8] <- 0  # cut off
    x <- sweep(sv$v, 2, di, "*") %*% t(sv$u) %*% b
    A2 <- rep(1, no_sp) 
    A2[retune] <- x
    
    # We may have to repeat this if any of the multipliers is negative or zero
    if (any(A2 <= 0)) {
        # Remove those species
        params <- removeSpecies(params, A2 <= 0)
        # and try again retuning the remaining retunable species
        retune <- retune[A2 > 0]
        if (any(retune)) {
            params <- retune_abundance(params, retune)
        } else {
            message("All background species have been removed.")
        }
    } else {
        # Use these abundance multipliers to rescale the abundance curves
        params@initial_n <- params@initial_n * A2
        # update SSB
        params@A <- params@A * A2
    }
    # Remove low abundance species
    # TODO: this could be vectorised
    no_sp <- length(params@species_params$species)
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
    
    params <- removeSpecies(params, remove)
    
    return(params)
}

#' Remove species from an ecosystem
#' 
#' This method simply removes all entries from the MizerParams object that 
#' refer to the selected species. It does not recalculate the initial 
#' abundances.
#' 
#' @param params A mizer params object for the original system.
#' @param species A vector of the names of the species 
#'                to be deleted or a boolean vector indicating for each species 
#'                whether it is to be removed (TRUE) or not.
#' 
#' @return An object of type \linkS4class{MizerParams}
#' @export
removeSpecies <- function(params, species) {
    no_sp <- length(params@w_min_idx)
    if (is.logical(species)) {
        if (length(species) != no_sp) {
            stop("The boolean species argument has the wrong length")
        }
        remove <- species
    } else {
        remove <- dimnames(params@initial_n)$sp %in% species
        if (length(remove) == 0) {
            warning("The species argument matches none of the species in the params object")
            return(params)
        }
    }
    keep <- !remove
    p <- params
    p@psi <- p@psi[keep, , drop = FALSE]
    p@initial_n <- p@initial_n[keep, , drop = FALSE]
    p@intake_max <- p@intake_max[keep, , drop = FALSE]
    p@search_vol <- p@search_vol[keep, , drop = FALSE]
    p@metab <- p@metab[keep, , drop = FALSE]
    p@ft_pred_kernel_e <- p@ft_pred_kernel_e[keep, , drop = FALSE]
    p@ft_pred_kernel_p <- p@ft_pred_kernel_p[keep, , drop = FALSE]
    p@mu_b <- p@mu_b[keep, , drop = FALSE]
    p@species_params <- p@species_params[keep, , drop = FALSE]
    p@interaction <- p@interaction[keep, keep, drop = FALSE]
    p@selectivity <- p@selectivity[, keep, , drop = FALSE]
    p@catchability <- p@catchability[, keep, drop = FALSE]
    p@w_min_idx <- p@w_min_idx[keep]
    p@A <- p@A[keep]
    
    return(p)
}


#' Add more species into an ecosystem with background species.
#'
#' Takes a \linkS4class{MizerParams} object and adds an additional species with
#' given parameters to the ecosystem.
#'
#' @param params A mizer params object for the original system. 
#' @param species_params The species parameters of the new species we
#'   want to add to the system.
#' @param SSB The spawning stock biomass of the new species. If not provided, 
#'   the abundance of the new species will be chosen so that its maximal 
#'   biomass density lies at half the community power law.
#' @param rfac A number that determines the strength of the non-linearity in
#'   the Beverton-Holt stock-recruitment relationship. The maximal recruitment
#'   will be set to rfac times the normal steady-state recruitment.
#'   Default value is 10.
#' @param effort Default value is 0.
#' @param ... Other arguments (unused)
#' 
#' @return An object of type \linkS4class{MizerParams}
#' @export
#' @examples
#' \dontrun{
#' params <- set_scaling_model(max_w_inf = 5000)
#' params <- setBackground(params)
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
addSpecies <- function(params, species_params, SSB = NA,
                       rfac=10, effort = 0) {
    # The code adds a new species into the system, and sets its abundance to the
    # steady state in the system where the new species does not self interact. Then
    # the abundance multipliers of the background species are retuned to retain the
    # old aggregate abundance curve, using retune_abundance(). Then the values of
    # erepro are altered so that the resulting configuration satisfies the steady
    # state reproduction boundary condition. The idea is that if the params system
    # is at steady state, and if the death rates of pre-existing species are close
    # to what they where before the new species were added, and if the newly added
    # species is at a low enough abundance (i.e., if mult is low enough) that the
    # assumption of it being none self interacting is approximately valid, then the
    # abundance curves attached to the params object returned by addSpecies() will
    # be a steady state,
    #
    # Note that we are assuming that the first species is a background species, and
    # the last species is a foreground species, with abundance multiplier mult.
    
    # check validity of parameters
    if (rfac <= 1) {
        message("rfac needs to be larger than 1. Setting rfac=1.01")
        rfac <- 1.01
    }
    if (any(species_params$species %in% params@species_params$species)) {
        stop("You can not add species that are already there.")
    }
    
    # provide erepro column that is later overwritten
    species_params$erepro <- 0.1
    
    # Move linecolour and linetype into species_params
    params@species_params$linetype <- 
        params@linetype[params@species_params$species]
    params@species_params$linecolour <- 
        params@linecolour[params@species_params$species]
    # TODO: Check if we need to do this with selectivity as well
    
    # Make sure that all columns exist in both data frames
    missing <- setdiff(names(params@species_params), names(species_params))
    species_params[missing] <- NA
    missing <- setdiff(names(species_params), names(params@species_params))
    params@species_params[missing] <- NA
    
    # add the new species (with parameters described by species_params), 
    # to make a larger species_params dataframe.
    combi_species_params <- rbind(params@species_params, species_params)
    
    # use dataframe and global settings from params to make a new MizerParams 
    # object.
    p <- MizerParams(
        combi_species_params,
        p = params@p,
        n = params@n,
        q = params@q,
        lambda = params@lambda,
        f0 = params@f0,
        kappa = params@kappa,
        min_w = min(params@w),
        max_w = max(params@w),
        no_w = length(params@w),
        min_w_pp = min(params@w_full),
        w_pp_cutoff = max(params@w_full),
        r_pp = (params@rr_pp / (params@w_full ^ (params@p - 1)))[1]
    )
    # Use the same resource spectrum as params
    p@initial_n_pp <- params@initial_n_pp
    p@cc_pp <- params@cc_pp
    new_sp <- length(params@species_params$species) + 1
    no_sp <- new_sp
    # Initially use abundance curves for pre-existing species 
    # (we shall retune the abundance multipliers of such 
    # species from the background later)
    p@initial_n[1:(new_sp - 1), ] <- params@initial_n
    # Use the same psi and mu_b as before for old species
    p@psi[1:(new_sp - 1), ] <- params@psi
    p@sc <- params@sc
    p@mu_b[1:(new_sp - 1), ] <- params@mu_b
    p@mu_b[new_sp, ] <- params@mu_b[1, ]  # NOTE: we assume same
    # background death for all species
    p@srr <- params@srr
    
    # Turn off self-interaction of the new species, so we can determine the
    # growth rates, and death rates induced upon it by the pre-existing species
    p@interaction[new_sp, new_sp] <- 0
    # compute death rate for new species
    mumu <- getMort(p, p@initial_n, p@initial_n_pp, effort = effort)[new_sp, ]
    # compute growth rate for new species
    gg <- getEGrowth(p, p@initial_n, p@initial_n_pp)[new_sp, ]
    
    # Compute solution for new species
    w_inf_idx <- sum(p@w < p@species_params$w_inf[new_sp])
    idx <- p@w_min_idx[new_sp]:(w_inf_idx - 1)
    if (any(gg[idx] == 0)) {
        stop("Can not compute steady state due to zero growth rates")
    }
    p@initial_n[new_sp, ] <- 0
    p@initial_n[new_sp, p@w_min_idx[new_sp]:w_inf_idx] <- 
        c(1, cumprod(gg[idx] / ((gg + mumu * p@dw)[idx + 1])))
    if (any(is.infinite(p@initial_n))) {
        stop("Candidate steady state holds infinities")
    }
    if (any(is.na(p@initial_n) || is.nan(p@initial_n))) {
        stop("Candidate steady state holds none numeric values")
    }
    
    # Normalise solution
    if (is.na(SSB)) {
        # If spawning stock biomass of new species is not supplied, 
        # normalise solution so that at its maximum it lies at half the 
        # power law, and then calculate its SSB.
        # We choose the maximum of the biomass density in log space
        # because that is always an increasing function at small size.
        idx <- which.max(p@initial_n[new_sp, ] * p@w^p@lambda)
        p@initial_n[new_sp, ] <- p@initial_n[new_sp, ] *
            p@kappa * p@w[idx]^(-p@lambda) / p@initial_n[new_sp, idx] / 2
        SSB <- sum(p@initial_n[new_sp, ] * p@w * p@dw * p@psi[new_sp, ])
    } else {
        unnormalised_SSB <- sum(p@initial_n[new_sp,] * p@w * p@dw * 
                                    p@psi[new_sp, ])
        p@initial_n[new_sp, ] <- p@initial_n[new_sp, ] * SSB / unnormalised_SSB
    }
    p@A <- c(params@A, SSB)
    
    # Turn self interaction back on
    p@interaction[new_sp, new_sp] <- 1
    
    # Retune the abundance multipliers to recreate the aggregate abundance
    # spectrum of the old params object.
    # First identify the retunable species. These are all background
    # species except the largest one
    retune <- is.na(p@A)
    p <- retune_abundance(p, retune)
    no_sp <- length(p@species_params$species)
    
    
    # Retune the values of erepro, so that we are at steady state.
    # First get death, growth and reproduction rates
    mumu <- getMort(p, p@initial_n, p@initial_n_pp, effort = effort)
    gg <- getEGrowth(p, p@initial_n, p@initial_n_pp)
    rdi <- getRDI(p, p@initial_n, p@initial_n_pp)
    erepro_final <- 1:no_sp  # set up vector of right dimension
    for (i in (1:no_sp)) {
        gg0 <- gg[i, p@w_min_idx[i]]
        mumu0 <- mumu[i, p@w_min_idx[i]]
        DW <- p@dw[p@w_min_idx[i]]
        if (!rdi[i] == 0) {
            erepro_final[i] <- p@species_params$erepro[i] *
                (p@initial_n[i, p@w_min_idx[i]] *
                     (gg0 + DW * mumu0)) / rdi[i]
        }
        else {
            erepro_final[i] <- 0.1
        }
    }
    # erepro needs to be divided by a factor of 1-1/rfac to
    # compensate for using a stock recruitment relationship
    # because RDD = (1-1/rfac) RDI
    erepro_final <- erepro_final / (1 - 1/rfac)
    
    p@species_params$erepro <- erepro_final
    
    p@species_params$r_max <- p@species_params$w_inf
    # set rmax = rfac*RDD = (rfac - 1)*RDI
    p@species_params$r_max <-
        (rfac - 1) * getRDI(p, p@initial_n, p@initial_n_pp)
    return(p)
}


#' Designate species as background species
#'
#' Background species are handled differently in some plots and their
#' abundance is automatically adjusted in addSpecies() to keep the community
#' close to the Sheldon spectrum.
#' 
#' @param params An object of class \linkS4class{MizerParams} or 
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
#' params <- MizerParams(NS_species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=20, t_save = 0.2)
#' sim <- setBackground(sim, species = c("Sprat", "Sandeel", 
#'                                             "N.pout", "Dab", "Saithe"))
#' plotSpectra(sim)
#' }
setBackground <- function(params, species = dimnames(params@initial_n)$sp) {
    params@A[dimnames(params@initial_n)$sp %in% species] <- NA
    params@linecolour[species] <- "grey"
    return(params)
}

#' Tune params object to be at steady state
#' 
#' This is done by running the dynamics for a specified number of years
#' while keeping recruitment constant to reach steady state. Then
#' the reproductive efficiencies are retuned to achieve that level of
#' recruitment.
#' 
#' @param params A \linkS4class{MizerParams} object
#' @param effort The fishing effort. Default is 0
#' @param t_max The maximum number of years to run the simulation. Default is 50.
#' @param t_per The simulation is broken up into shorter runs of t_per years,
#'   after each of which we check for convergence. Default value is 2.
#' @param tol The simulation stops when the relative change in the egg
#'   production RDI over a t_per is less than rel_tol for every background
#'   species. Default value is 1/100.
#' @export
steady <- function(params, effort = 0, t_max = 50, t_per = 2, tol = 10^(-2),
                   shiny_progress = NULL) {
    p <- params
    
    if (hasArg(shiny_progress)) {
        # We have been passed a shiny progress object
        shiny_progress$set(message = "Finding steady state", value = 0)
        proginc <- 1/ceiling(t_max/t_per)
    }
    # Force the recruitment to stay at the current level
    rdd <- getRDD(p, p@initial_n, p@initial_n_pp)
    p@srr <- function(rdi, species_params) {rdd}
    
    n <- p@initial_n
    n_pp <- p@initial_n_pp
    old_rdi <- getRDI(p, n, n_pp)
    for (ti in (1:ceiling(t_max/t_per))){
        sim <- project(p, t_max = t_per, t_save = t_per, effort = effort, 
                       initial_n = n, initial_n_pp = n_pp)
        # advance shiny progress bar
        if (hasArg(shiny_progress)) {
            shiny_progress$inc(amount = proginc)
        }
        n <- sim@n[dim(sim@n)[1],,]
        n_pp <- sim@n_pp[dim(sim@n_pp)[1],]
        new_rdi <- getRDI(p, n, n_pp)
        deviation <- max(abs((new_rdi - old_rdi)/old_rdi)[!is.na(p@A)])
        if (deviation < tol) {
            break
        }
        old_rdi <- new_rdi
        
    }
    if (deviation >= tol) {
        warning(paste(
            "Simulation run in steady() did not converge after ", ti * t_per, 
            "years. Residual relative rate of change = ", deviation))
    } else {
        message(paste("Steady state was reached after ", ti * t_per, "years."))
    }
    
    # Restore original stock-recruitment relationship
    p@srr <- params@srr
    
    no_sp <- length(p@species_params$species)
    no_t <- dim(sim@n)[1]
    p@initial_n <- sim@n[no_t, , ]
    p@initial_n_pp <- sim@n_pp[no_t, ]
    
    # Retune the values of erepro so that we get the correct level of
    # recruitment
    mumu <- getMort(p, p@initial_n, p@initial_n_pp, effort = effort)
    gg <- getEGrowth(p, p@initial_n, p@initial_n_pp)
    rdd <- getRDD(p, p@initial_n, p@initial_n_pp)
    # TODO: vectorise this
    for (i in (1:no_sp)) {
        gg0 <- gg[i, p@w_min_idx[i]]
        mumu0 <- mumu[i, p@w_min_idx[i]]
        DW <- p@dw[p@w_min_idx[i]]
        p@species_params$erepro[i] <- p@species_params$erepro[i] *
            (p@initial_n[i, p@w_min_idx[i]] *
                 (gg0 + DW * mumu0)) / rdd[i]
    }
    
    return(p)
}