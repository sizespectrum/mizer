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
# Development has received funding from the European Commission’s Horizon 2020 
# Research and Innovation Programme under Grant Agreement No. 634495 
# for the project MINOUW (http://minouw-project.eu/).
# Distributed under the GPL 3 or later 
# Maintainer: Gustav Delius, University of York, <gustav.delius@york.ac.uk>

#### set_community_model ####
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
#' @param kappa The carrying capacity of the background spectrum. Default value
#'   is 1000.
#' @param lambda The exponent of the background spectrum. Default value is 2 + q
#'   - n.
#' @param r_pp Growth rate of the primary productivity. Default value is 10.
#' @param gamma Volumetric search rate. Estimated using \code{h}, \code{f0} and 
#'   \code{kappa} if not supplied.
#' @param recruitment The constant recruitment in the smallest size class of the
#'   community spectrum. This should be set so that the community spectrum 
#'   continues the background spectrum. Default value = \code{kappa} * 
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
                                lambda = 2+q-n,
                                f0 = 0.7,
                                r_pp = 10,
                                gamma = NA,
                                knife_edge_size = 1000,
                                knife_is_min = TRUE,
                                recruitment = kappa * min_w^-lambda,
                                rec_mult = 1,
                                ...
                                ){
    w_inf <- max_w * 0.9
    w_pp_cutoff <- min_w
    ks <- 0 # Turn off standard metabolism
    p <- n # But not used as ks = 0
    # Estimate gamma if not supplied
    if (is.na(gamma)){
        gamma <- (f0 * h * beta^(2-lambda)) / ((1-f0)*sqrt(2*pi)*kappa*sigma)
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
    com_params <- MizerParams(com_params_df, p=p, n=n,q=q, lambda = lambda, 
                              kappa = kappa, min_w = min_w, max_w = max_w, 
                              w_pp_cutoff = w_pp_cutoff, r_pp = r_pp, ...)
    com_params@srr <- constant_recruitment
    com_params@psi[] <- 0 # Need to force to be 0. Can try setting w_mat but 
                          # due to slope still not 0
    # Set w_mat to NA for clarity - it is not actually being used
    com_params@species_params$w_mat[] <- NA
    return(com_params)
}



#### set_trait_model ####
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
#' @param min_w_pp The smallest size of the background spectrum.
#' @param no_w_pp Obsolete argument that is no longer used because the number
#'    of plankton size bins is determined because all size bins have to
#'    be logarithmically equally spaced.
#' @param w_pp_cutoff The cut off size of the background spectrum. Default value
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
    # If not supplied, calculate gamma using equation 2.1 in A&P 2010
    # TODO: remove this here because it is already calculated in MizerParams()
    #       Having the same code in two locations is not a good idea
    if(is.na(gamma)){
        alpha_e <- sqrt(2*pi) * sigma * beta^(lambda-2) * 
            exp((lambda-2)^2 * sigma^2 / 2) # see A&P 2009
        gamma <- h * f0 / (alpha_e * kappa * (1-f0)) # see A&P 2009 
    }
    w_inf <- 10^seq(from=log10(min_w_inf), to = log10(max_w_inf), length=no_sp)
    w_mat <- w_inf * eta

    # Check gears
    if (length(knife_edge_size) > no_sp){
        stop("There cannot be more gears than species in the model")
    }
    if ((length(knife_edge_size) > 1) & (length(knife_edge_size) != no_sp)){
        warning("Number of gears is less than number of species so gear information is being recycled. Is this what you want?")
    }
    if ((length(gear_names) != 1) & (length(gear_names) != no_sp)){
        stop("Length of gear_names argument must equal the number of species.")
    }

    # Make the species parameters data.frame
    trait_params_df <- data.frame(
            species = 1:no_sp,
            w_inf = w_inf,
            w_mat = w_mat,
            h = h, # max food intake
            gamma = gamma, # vol. search rate,
            ks = ks,# standard metabolism coefficient,
            beta = beta,
            sigma = sigma,
            z0 = z0pre * w_inf^(n-1), # background mortality
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
    dw_winf <- tmpB * tmpA *10^(tmpB*((1:no_sp)-1)) # ?
    N0_max <- k0 * w_inf^(n*2-q-3+alpha_rec) * dw_winf  # Why * dw_winf, not / ? Ken confirms * in email
    # No need to include (1 - psi) in growth equation because allocation to reproduction at this size = 0, so 1 - psi = 1
    g0 <- (alpha * f0 * h * trait_params@w[1]^n - ks * trait_params@w[1]^p)
    r_max <- N0_max * g0

    trait_params@species_params$r_max <- r_max

    return(trait_params)
}


#### set_scaling_model ####
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
#' @param min_w_pp The smallest size of the background spectrum. Default value
#'   is min_egg/(beta*exp(5*sigma)) so that it covers the entire range of the
#'   feeding kernel of even the smallest fish larva.
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
                              ...) {
    if (hasArg(lambda)) {
        # The lambda argument overrules any q argument
        q <- lambda - 2 + n
    }
    # check validity of parameters
    if (rfac <= 1) {
        message("rfac can not be smaller than 1. Setting rfac=1.1")
        rfac <- 1.1
    }
    no_w <- round(no_w)
    if (no_w < 1) {
        stop("The number of size bins no_w must be a positive integer")
    }
    if (no_w < log10(max_w_inf/min_egg)*5) {
        no_w <- round(log10(max_w_inf/min_egg)*5+1)
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
        species = 1:no_sp,
        w_min = w_min,
        w_inf = w_inf,
        w_mat = w_mat,
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
    # Get constants for analytic solution
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
    n_mult <- (1 - (w / w_inf[1]) ^ (1 - n)) ^ (pow - 1) * 
                  (1 - (w_mat[1] / w_inf[1]) ^(1 - n)) ^ (-pow)
    n_mult[w < w_mat[1]] <- 1
    n_mult[w >= w_inf[1]] <- 0
    # Create steady state solution n_exact for species 1
    n_exact <- ((w_min[1] / w) ^ (mu0 / hbar) / (hbar * w ^ n)) * n_mult
    n_exact <- n_exact[w >= w_min[1] & w <= w_inf[1]]
    # rescale fish abundance to line up with background resource spectrum
    mult <- kappa / 
        sum(n_exact * (w^(lambda-1)*dw)[w >= w_min[1] & w <= w_inf[1]])
    n_exact <- n_exact * mult * 
        (10^(dist_sp*(1-lambda)/2) - 10^(-dist_sp*(1-lambda)/2)) / (1-lambda)
    # Use n_exact as a template to create solution initial_n for all species
    initial_n <- params@psi  # get array with correct dimensions and names
    initial_n[, ] <- 0
    for (i in 1:no_sp) {
        # smallest index for species
        w_min_idx <- params@species_params$w_min_idx[i]
        # range of indices
        idxs <- w_min_idx:(w_min_idx+length(n_exact)-1)  
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
    # for (i in seq_len(max(params@species_params$w_min_idx) - 1)) {
    #     sc[idxs] <- sc[idxs] + n_exact * fac^(i-1)
    #     idxs <- idxs+1
    # }
    # sc <- sc * 
    #     (10^(delt*(1-lambda)/2) - 10^(-delt*(1-lambda)/2)) /
    #     (10^(dist_sp*(1-lambda)/2) - 10^(-dist_sp*(1-lambda)/2))
    params@sc <- sc
    
    # Setup plankton
    plankton_vec <- (kappa * w ^ (-lambda)) - sc
    # Cut off plankton at maximum size of smallest species
    plankton_vec[w >= min_w_inf] <- 0
    # Do not allow negative plankton abundance
    if (sum(plankton_vec < 0) > 0) {
        message("Note: Negative abundance values in background resource overwritten with zeros")
        plankton_vec[plankton_vec < 0] <- 0
    }
    # The cc_pp factor needs to be higher than the desired steady state in
    # order to compensate for predation mortality
    params@cc_pp[sum(params@w_full <= w[1]):length(params@cc_pp)] <-
        plankton_vec
    initial_n_pp <- params@cc_pp
    m2_background <- getM2Background(params, initial_n, initial_n_pp)
    params@cc_pp <- (1 + m2_background / params@rr_pp) * initial_n_pp

    # Setup background death and steplike psi
    m2 <- getM2(params, initial_n, initial_n_pp)
    for (i in 1:no_sp) {
        params@psi[i,] <- (w / w_inf[i]) ^ (1 - n)
        params@psi[i, w < (w_mat[i] - 1e-10)] <- 0
        params@psi[i, w > (w_inf[i] - 1e-10)] <- 1
        params@mu_b[i,] <- mu0 * w ^ (n - 1) - m2[i, ]
        if (sum(params@mu_b[i,] < 0) > 0) {
            message("Note: Negative background mortality rates overwritten with zeros")
        }
        params@mu_b[i, params@mu_b[i,] < 0] <- 0
    }
    # Set erepro to meet boundary condition
    rdi <- getRDI(params, initial_n, initial_n_pp)
    gg <- getEGrowth(params, initial_n, initial_n_pp)
    mumu <- getZ(params, initial_n, initial_n_pp, effort = 0)
    erepro_final <- 1:no_sp  # set up vector of right dimension
    for (i in (1:no_sp)) {
        gg0 <- gg[i, params@species_params$w_min_idx[i]]
        mumu0 <- mumu[i, params@species_params$w_min_idx[i]]
        DW <- params@dw[params@species_params$w_min_idx[i]]
        erepro_final[i] <- erepro * 
            (initial_n[i, params@species_params$w_min_idx[i]] *
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
        (rfac - 1) * getRDI(params, initial_n, initial_n_pp)[,1]
    return(params)
}


#### retune_abundance ####
#' Retunes abundance multipliers of background species so aggregate abundance is
#' a power law.
#'
#' If N_i(w) is a steady state of the McKendrik-von Foerster (MVF) equation with
#' fixed growth and death rates, then A_i*N_i(w) is also a steady state, where
#' A_i is an abundance multiplier. When we add a foreground species to our
#' model, we want to choose new abundance multipliers of the background species
#' so that we the abundance, summed over all species and background resources,
#' is close to sc(w), which is the aggregate abundance of all but the last
#' species. We are assuming this last species is newly added, with A_i=1.
#'
#' retune_abundance operates of a params object, with a slot A. If i is a
#' background species, then A_i=NA, indicating we are allowed to retune the
#' abundance multiplier.
#'
#' @param params A mizer params object with an A slot with 1's for species we
#'   wish to hold fixed the abundance multiplier of, and NA's for species that
#'   we shall vary the abundance multiplier of.
#'   
#' @return An object of type \code{MizerParams}
#' @seealso \linkS4class{MizerParams}
#' @export
#' @examples
#' \dontrun{
#' params <- set_scaling_model()
#' params@A[] <- NA
#' params@A[length(params@A)] <- 1
#' retune_abundance(params)
#' }
retune_abundance <- function(params) {
    no_sp <- length(params@species_params$w_inf)
    # get a list of the species that we can tune abundance mulpliers of
    all_background <- is.na(params@A)
    largest_background <-
        which.max(params@species_params$w_inf[all_background])
    # we are assuming that the abundance multiplier of the largest backgroud
    # species should be held fixed at 1, if though it was initially an NA.
    # to represent this we make a new background indicator A2
    A2 <- params@A
    A2[largest_background] <- 1
    # we make a list L of species we will vary the abundance parameters of
    # (everything but largest background)
    L <- (1:no_sp)[is.na(A2)]
    # Determine the indices of the limits we shall integrate between
    idx_start <- sum(params@w <= min(params@species_params$w_mat))
    idx_stop <- sum(params@w < max(params@species_params$w_inf))
    # The problem is to vary the abundance multupliers of species in L so that,
    # between the limits, to sum of the abundances of the species is "close" to
    # the power law. More precisely, we find the abundance multipliers A_i so
    # that the integral of the square of the relative distance (sum_{i not in L}
    # A_i*N_i(w) + sum_{i not in L} N_i(w) - c(w))/c(w) over w, between our
    # limits, is minimized. Here c(w) is the sum of the abundances of all but
    # the last (newly added) species, and L is the set of all background
    # species, except the largest (that we keep the abundance multipliers of
    # fixed).
    # 
    #! how to define sc in general ? should it be smoothed ?
    # sc used to be defined as
    # sc <- colSums(params@initial_n[all_background, ])
    # but now we are assuming that the newly added species
    # is the last one, and we are retunning to
    # get similar to the aggregate abundance of the
    # species (1:(no_sp-1))
    sc <- params@sc
    # rho is the total abundance of all the species that have their abundance multipliers
    # held fixed.
    Lcomp <- (1:no_sp)[!is.na(A2)]
    rho <- colSums(params@initial_n[Lcomp, , drop=FALSE])
    # We solve a linear system to find the abundance multipliers, first we initialize
    # the matrix RR and vector QQ
    RR <- matrix(0, nrow = length(L), ncol = length(L))
    QQ <- (1:length(L))
    den <- sc ^ 2
    den[den == 0] <- 10 ^ (-50)
    # Next we fill out the values of QQ and RR
    for (i in (1:length(L))) {
        QQ[i] <-
            sum((params@initial_n[L[i],] * (sc - rho) * params@dw / (den))[idx_start:idx_stop])
        for (j in (1:length(L))) {
            RR[i, j] <-
                sum((
                    params@initial_n[L[i],] * params@initial_n[L[j],] * params@dw / (den)
                )[idx_start:idx_stop])
        }
    }
    # Now we solve the linear system to find the abundance multipliers that
    # yield our power law
    A2[L] <- solve(RR, QQ)
    print(A2)
    if (sum(A2 < 0) > 0) {
        stop('Abundance multipliers generated with negative entries')
        #! we should add an extra iteration to solve this issue of -ve
        # abundance multipliers by holding certain species off.
    }
    return(A2)
}


#### add_species ####
#' Add more species into an ecosystem with scaling background species.
#'
#' Takes a \linkS4class{MizerParams} object and adds additional species with
#' given parameters to the ecosystem.
#'
#' @param params A mizer params object for the original system. The A slot holds
#'   1's for foreground species we wish to hold fixed the abundance multiplier
#'   of, and NA's for for background species that we shall vary the abundance
#'   multiplier of.
#' @param species_params The species parameters of the foreground species we
#'   want to add to the system.
#' @param biomass The total biomass of members of the newly added species which
#'   have weights greater than min_w_observed. If not provided, the abundance
#'   of new species will be chosen so that the biomass of its mature fish is
#'   the same as that of the species with the next smaller maturity size.
#' @param min_w_observed The minimum weight of fish of the new species that
#'   contribute to the biomass which we setup for the new species. Default value
#'   is the maturity size of the new species.
#' @param rfac Default value is 10.
#' @param effort Default value is 0.
#' 
#' @export
#' @return An object of type \code{MizerParams}
#' @seealso \linkS4class{MizerParams}
#' @examples
#' \dontrun{
#' params <- set_scaling_model(max_w_inf = 5000)
#' params@A[] <- NA
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
#'     erepro = 0.1, 
#'     sel_func = "knife_edge", 
#'     knife_edge_size = 100, 
#'     gear = "knife_edge_gear",
#'     k = 0,
#'     r_max = Inf,
#'     k_vb = 0.6,
#'     a = a_m,
#'     b = b_m
#' )
#' params_out <- add_species(params, species_params, biomass = 3070953023, 
#'                           min_w_observed = 0)
#' sim <- project(params_out, t_max = 5, effort = 0)
#' plot(sim)
#' }
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
# abundance curves attached to the params object returned by add_species() will
# be a steady state,
#
# Note that we are assuming that the first species is a background species, and
# the last species is a foreground species, with abundance multiplier mult.
add_species <- function(params, species_params, biomass = NA, 
                        min_w_observed = species_params$w_mat, rfac=10, effort = 0) {
    
    # Set r_max to Inf if absent
    if (is.null(params@species_params$r_max)){
        params@species_params$r_max <- params@species_params$w_inf
        params@species_params$r_max[] <- Inf
    }
    if (is.null(species_params$r_max)){
        species_params$r_max <- Inf
    }
    
    # ensure dataframes to be merged have the same columns, 
    # regarding k_vb, a and b
    if (is.null(params@species_params$k_vb)&(!is.null(species_params$k_vb))){
        params@species_params$k_vb <- params@species_params$w_inf
        params@species_params$k_vb[] <- NA
    }
    if (!is.null(params@species_params$k_vb)&(is.null(species_params$k_vb))){
        species_params$k_vb <- NA
    }
    if ((!("a" %in% colnames(params@species_params)))&(("a" %in% colnames(species_params)))){
        params@species_params$a <- params@species_params$w_inf
        params@species_params$a[] <- NA
    }
    if ((("a" %in% colnames(params@species_params)))&(!("a" %in% colnames(species_params)))){
        species_params$a <- NA
    }
    if ((!("b" %in% colnames(params@species_params)))&(("b" %in% colnames(species_params)))){
        params@species_params$b <- params@species_params$w_inf
        params@species_params$b[] <- NA
    }
    if ((("b" %in% colnames(params@species_params)))&(!("b" %in% colnames(species_params)))){
        species_params$b <- NA
    }
    if ((!("l25" %in% colnames(params@species_params)))&(("l25" %in% colnames(species_params)))){
        params@species_params$l25 <- params@species_params$w_inf
        params@species_params$l25[] <- NA
    }
    if ((("l25" %in% colnames(params@species_params)))&(!("l25" %in% colnames(species_params)))){
        species_params$l25 <- NA
    }
    if ((!("l50" %in% colnames(params@species_params)))&(("l50" %in% colnames(species_params)))){
        params@species_params$l50 <- params@species_params$w_inf
        params@species_params$l50[] <- NA
    }
    if ((("l50" %in% colnames(params@species_params)))&(!("l50" %in% colnames(species_params)))){
        species_params$l50 <- NA
    }
    
    # calculate h if it is missing
    if (is.null(species_params$h) || is.na(species_params$h)){
        message("Note: \tNo h column in new species data frame so using f0 and k_vb to
                calculate it.")
        if(!("k_vb" %in% colnames(species_params))){
            stop("\t\tExcept I can't because there is no k_vb column in the new species data frame")
        }
        fc <- 0.2/species_params$alpha
        species_params$h <- 3*species_params$k_vb*(species_params$w_inf^(1/3))/
            (species_params$alpha*params@f0-0.2)
    }
    
    # calculate ks if it is missing
    if (is.null(species_params$ks) || is.na(species_params$ks)){
        message("Note: \tNo ks column in new species data frame. Setting ks = 0.2*h.")
        species_params$ks <- 0.2*species_params$h # mizer's default setting
    }
    
    # calculate gamma if it is missing
    if (is.null(species_params$gamma) || is.na(species_params$gamma)){
        message("Note: \tNo gamma column in new species data frame so using f0, h, beta, sigma, lambda and kappa to calculate it.")
        ae <- sqrt(2*pi) * species_params$sigma * 
            species_params$beta^(params@lambda-2) * 
            exp((params@lambda-2)^2 * species_params$sigma^2 / 2)
        species_params$gamma <- (species_params$h / (params@kappa * ae)) * 
            (params@f0 / (1 - params@f0))
    }
    
    # calculate w_min_idx if it is missing
    species_params$w_min_idx <- sum(params@w<=species_params$w_min)
    
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
    # Initially use abundance curves for pre-existing species 
    # (we shall retune the abundance multipliers of such 
    # species from the background later)
    p@initial_n[1:(new_sp - 1), ] <- params@initial_n
    # Use the same psi and mu_b as before for old species
    p@psi[1:(new_sp - 1), ] <- params@psi
    p@sc <- params@sc
    p@mu_b[1:(new_sp - 1), ] <- params@mu_b
    p@mu_b[new_sp, ] <- params@mu_b[1, ]
    p@srr <- params@srr
    
    # Turn off self-interaction of the new species, so we can determine the
    # growth rates, and death rates induced upon it by the pre-existing species
    p@interaction[new_sp, new_sp] <- 0
    # compute death rate for new species
    mumu <- getZ(p, p@initial_n, p@initial_n_pp, effort = effort)[new_sp, ]
    # compute growth rate for new species
    gg <- getEGrowth(p, p@initial_n, p@initial_n_pp)[new_sp, ]
    
    # Compute integral to solve MVF for new species
    w_inf_idx <- sum(p@w < p@species_params$w_inf[new_sp])
    if (any(gg[p@species_params$w_min_idx[new_sp]:w_inf_idx]==0)) {
        stop("Can not compute steady state due to zero growth rates")
    }
    integrand <-
        params@dw[p@species_params$w_min_idx[new_sp]:w_inf_idx] * 
        mumu[p@species_params$w_min_idx[new_sp]:w_inf_idx] /
        gg[p@species_params$w_min_idx[new_sp]:w_inf_idx]
    p@initial_n[new_sp, ] <- 0
    p@initial_n[new_sp, p@species_params$w_min_idx[new_sp]:w_inf_idx] <-
        exp(-cumsum(integrand)) / gg[p@species_params$w_min_idx[new_sp]:w_inf_idx]
    if (any(is.infinite(p@initial_n))) {
        stop("Candidate steady state holds infinities")
    }
    if (any(is.na(p@initial_n) || is.nan(p@initial_n))) {
        stop("Candidate steady state holds none numeric values")
    }
    
    # Normalise solution
    if (is.na(biomass)) {
        # If biomass of new species is not supplied, set it to the biomass of
        # the species with the next-smaller maturity size.
        idx <- which.min(p@species_params$w_mat < species_params$w_mat) - 1
        biomass <- sum(p@initial_n[idx, ] * p@w * p@dw * (p@w >= min_w_observed))
    }
    unnormalised_biomass <- sum(p@initial_n[new_sp,] * p@w * p@dw *
                                    (p@w >= min_w_observed))
    p@initial_n[new_sp, ] <- p@initial_n[new_sp, ] * biomass / 
        unnormalised_biomass
    
    # Turn self interaction back on
    p@interaction[new_sp, new_sp] <- 1
    # Arrange background inidicators A, so show the new species is in the
    # foreground
    p@A <- c(params@A, 1)
    # Retune the abundance multipliers to recreate the aggregate abundance
    # spectrum of the old params object.
    AA <- retune_abundance(p)
    # Use these abundance multipliers to rescale the abundance curves so that 
    # their aggregation is close to the power law
    p@initial_n <- p@initial_n * AA
    
    # Retune the values of erepro, so that we are at steady state.
    # First get death and growth rates
    mumu <- getZ(p, p@initial_n, p@initial_n_pp, effort = effort)
    gg <- getEGrowth(p, p@initial_n, p@initial_n_pp)
    rdi <- getRDI(p, p@initial_n, p@initial_n_pp)
    no_sp <- new_sp
    erepro_final <- 1:no_sp  # set up vector of right dimension
    for (i in (1:no_sp)) {
        gg0 <- gg[i, p@species_params$w_min_idx[i]]
        mumu0 <- mumu[i, p@species_params$w_min_idx[i]]
        DW <- params@dw[p@species_params$w_min_idx[i]]
        erepro_final[i] <- p@species_params$erepro[i] *
            (p@initial_n[i, p@species_params$w_min_idx[i]] *
                 (gg0 + DW * mumu0)) / rdi[i]
    }
    if (is.finite(rfac)) {
        # erepro has been multiplied by a factor of (rfac/(rfac-1)) to
        # compensate for using a stock recruitment relationship.
        erepro_final <- (rfac / (rfac - 1)) * erepro_final
    }
    p@species_params$erepro <- erepro_final
    
    p@species_params$r_max <- p@species_params$w_inf
    # set rmax=fac*RDD
    # note that erepro has been multiplied by a factor of (rfac/(rfac-1)) to
    # compensate for using a stock recruitment relationship.
    p@species_params$r_max <-
        (rfac - 1) * getRDI(p, p@initial_n, p@initial_n_pp)[,1]
    return(p)
}
