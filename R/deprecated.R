#' Deprecated obsolete function for setting up multispecies parameters
#' 
#' This function has been deprecated in favour of the function
#' [newMultispeciesParams()] that sets better default values.
#' 
#' @inheritParams newMultispeciesParams
#' @inheritParams newTraitParams
#' @param q Allometric exponent of search volume
#' @param ... Unused
#' @export
#' @family deprecated functions
set_multispecies_model <- 
    function(
        species_params,
        interaction = matrix(1,
                             nrow = nrow(species_params),
                             ncol = nrow(species_params)),
        min_w_pp = 1e-10,
        min_w = 0.001,
        max_w = max(species_params$w_inf)*1.1,
        no_w = 100,
        n = 2 / 3,
        q = 0.8,
        f0 = 0.6,
        kappa = 1e11,
        lambda = 2 + q - n,
        r_pp = 10,
        ...) {
    if (exists("no_w_pp")) {
        warning("New mizer code does not support the parameter no_w_pp")
    }
    
    # Need to correct for the fact that new mizer extends w_full to BELOW
    # min_w_pp
        dx <- log10(max_w / min_w) / (no_w - 1)
        min_w_pp <- min_w_pp * 10 ^ dx
        
    species_params[["q"]] <- q
    species_params[["f0"]] <- f0
    object <- species_params
    # old code from MizerParams() in version 1.0.1
    
    # Set default values for column values if missing
    # If no gear_name column in object, then named after species
    if (!("gear" %in% colnames(object)))
        object$gear <- object$species
    # If no k column (activity coefficient) in object, then set to 0
    if (!("k" %in% colnames(object)))
        object$k <- 0
    # If no alpha column in object, then set to 0.6
    # Should this be a column? Or just an argument?
    if (!("alpha" %in% colnames(object)))
        object$alpha <- 0.6
    # If no erepro column in object, then set to 1
    if (!("erepro" %in% colnames(object)))
        object$erepro <- 1
    # If no sel_func column in species_params, set to 'knife_edge'
    if (!("sel_func" %in% colnames(object))) {
        message("Note: No sel_func column in species data frame. Setting selectivity to be 'knife_edge' for all species.")
        object$sel_func <- 'knife_edge'
        # Set default selectivity size
        if (!("knife_edge_size" %in% colnames(object))) {
            message("Note: No knife_edge_size column in species data frame. Setting knife edge selectivity equal to w_mat.")
            object$knife_edge_size <- object$w_mat
        }
    }
    # If no catchability column in species_params, set to 1
    if (!("catchability" %in% colnames(object)))
        object$catchability <- 1
    # Sort out h column If not passed in directly, is calculated from f0 and
    # k_vb if they are also passed in
    if (!("h" %in% colnames(object))) {
        message("Note: No h column in species data frame so using f0 and k_vb to calculate it.")
        if (!("k_vb" %in% colnames(object))) {
            stop("Except I can't because there is no k_vb column in the species data frame")
        }
        object$h <- ((3 * object$k_vb) / (object$alpha * f0)) * (object$w_inf ^ (1/3))
    }
    # Sorting out gamma column
    if (!("gamma" %in% colnames(object))) {
        message("Note: No gamma column in species data frame so using f0, h, beta, sigma, lambda and kappa to calculate it.")
        ae <- sqrt(2*pi) * object$sigma * object$beta^(lambda-2) * exp((lambda-2)^2 * object$sigma^2 / 2)
        object$gamma <- (object$h / (kappa * ae)) * (f0 / (1 - f0))
    }
    # Sort out ks column
    if (!("ks" %in% colnames(object))) {
        message("Note: No ks column in species data frame so using ks = h * 0.2.")
        object$ks <- object$h * 0.2
    }
    
    # The m column did not exist in the old version, set it to 1
    if (!("m" %in% colnames(object))) {
        message("Note: No m column in species data frame so using m = 1.")
        object[["m"]] <- 1
    }
    
    return(newMultispeciesParams(object,
                                 interaction = interaction,
                                 n = n,
                                 kappa = kappa,
                                 lambda = lambda,
                                 r_pp = r_pp,
                                 max_w = max_w,
                                 min_w = min_w,
                                 min_w_pp = min_w_pp,
                                 no_w = no_w,
                                 ...))
}

#' Alias for set_multispecies_model
#' 
#' An alias provided for backward compatibility with mizer version <= 1.0
#' @inherit set_multispecies_model
#' @export
MizerParams <- set_multispecies_model


# Copied from version 1.0.1
#' Deprecated function for setting up parameters for a trait-based model
#' 
#' This function has been deprecated in favour of the function
#' [newTraitParams()] that sets better default values.
#' 
#' This functions creates a `MizerParams` object so that trait-based-type 
#' models can be easily set up and run. The trait-based size spectrum model can
#' be derived as a simplification of the general size-based model used in
#' `mizer`. The species-specific parameters are the same for all species, 
#' except for
#' the asymptotic size, which is considered the most important trait
#' characterizing a species. Other parameters are related to the asymptotic
#' size. For example, the size at maturity is given by \code{w_inf * eta}, 
#' where `eta` is
#' the same for all species. For the trait-based model the number of species is
#' not important. For applications of the trait-based model see Andersen &
#' Pedersen (2010). See the `mizer` vignette for more details and examples
#' of the trait-based model.
#' 
#' The function has many arguments, all of which have default values. Of
#' particular interest to the user are the number of species in the model and
#' the minimum and maximum asymptotic sizes. The asymptotic sizes of the species
#' are spread evenly on a logarithmic scale within this range.
#' 
#' The stock recruitment relationship is the default Beverton-Holt style. The
#' maximum recruitment is calculated using equilibrium theory (see Andersen &
#' Pedersen, 2010) and a multiplier, `k0`. Users should adjust `k0` to
#' get the spectra they want.
#' 
#' The factor for the search volume, `gamma`, is calculated using the
#' expected feeding level, `f0`.
#' 
#' Fishing selectivity is modelled as a knife-edge function with one parameter,
#' `knife_edge_size`, which is the size at which species are selected. Each
#' species can either be fished by the same gear (`knife_edge_size` has a
#' length of 1) or by a different gear (the length of `knife_edge_size` has
#' the same length as the number of species and the order of selectivity size is
#' that of the asymptotic size).
#' 
#' The resulting `MizerParams` object can be projected forward using
#' `project` like any other `MizerParams` object. When projecting
#' the community model it may be necessary to reduce `dt` to 0.1 to avoid
#' any instabilities with the solver. You can check this by plotting the biomass
#' or abundance through time after the projection.
#' 
#' @param no_sp The number of species in the model. The default value is 10. The
#'   more species, the longer takes to run.
#' @param min_w_inf The asymptotic size of the smallest species in the
#'   community.
#' @param max_w_inf The asymptotic size of the largest species in the community.
#' @param no_w The number of size bins in the community spectrum.
#' @param min_w The smallest size of the community spectrum.
#' @param max_w Obsolete argument because the maximum size of the consumer
#'   spectrum is set to max_w_inf.
#' @param min_w_pp Obsolete argument because the smallest resource size is set
#'   to the smallest size at which the consumers feed.
#' @param w_pp_cutoff The cut off size of the resource spectrum. Default value
#'    is 1.
#' @param k0 Multiplier for the maximum recruitment. Default value is 50.
#' @param n Scaling of the intake. Default value is 2/3.
#' @param p Scaling of the standard metabolism. Default value is 0.75.
#' @param q Exponent of the search volume. Default value is 0.9.
#' @param eta Factor to calculate `w_mat` from asymptotic size.
#' @param r_pp Growth rate parameter for the resource spectrum. Default value is 4.
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
#' @param f0 Expected average feeding level. Used to set `gamma`, the
#'   factor for the search volume. The default value is 0.5.
#' @param gamma Volumetric search rate. Estimated using `h`, `f0` and
#'   `kappa` if not supplied.
#' @param knife_edge_size The minimum size at which the gear or gears select
#'   species. Must be of length 1 or no_sp.
#' @param gear_names The names of the fishing gears. A character vector, the
#'   same length as the number of species. Default is 1 - no_sp.
#' @param ... Other arguments to pass to the `MizerParams` constructor.
#' 
#' @return An object of type `MizerParams`
#' @references K. H. Andersen and M. Pedersen, 2010, Damped trophic cascades
#'   driven by fishing in model marine ecosystems. Proceedings of the Royal
#'   Society V, Biological Sciences, 1682, 795-802.
#' @export
#' @family deprecated functions
set_trait_model <- function(no_sp = 10,
                            min_w_inf = 10,
                            max_w_inf = 1e5,
                            no_w = 100,
                            min_w = 0.001,
                            max_w = max_w_inf * 1.1,
                            min_w_pp = 1e-10,
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
    if (exists("no_w_pp")) {
        warning("New mizer code does not support the parameter no_w_pp")
    }
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
        species = as.factor(1:no_sp),
        w_inf = w_inf,
        w_mat = w_mat,
        w_min = min_w,
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
    
    trait_params@species_params$R_max <- r_max
    
    return(trait_params)
}


# Copied from version 1.0.1
#' Deprecated function for setting up parameters for a community-type model
#' 
#' This function has been deprecated in favour of the function
#' [newCommunityParams()] that sets better default values.
#' 
#' This functions creates a \code{\linkS4class{MizerParams}} object so that
#' community-type models can be easily set up and run. A community model has
#' several features that distinguish it from the food-web type models. Only one
#' 'species' is resolved, i.e. one 'species' is used to represent the whole
#' community. The resource spectrum only extends to the start of the community
#' spectrum. Recruitment to the smallest size in the community spectrum is
#' constant and set by the user. As recruitment is constant, the proportion of
#' energy invested in reproduction (the slot `psi` of the returned 
#' `MizerParams` object) is set to 0. Standard metabolism has been turned 
#' off (the parameter `ks` is set to 0). Consequently, the growth rate is 
#' now determined solely by the assimilated food (see the package vignette for 
#' more details).
#' 
#' The function has many arguments, all of which have default values. The main 
#' arguments that the users should be concerned with are `z0`, 
#' `recruitment`, `alpha` and `f0` as these determine the average
#' growth rate of the community.
#' 
#' Fishing selectivity is modelled as a knife-edge function with one parameter, 
#' `knife_edge_size`, which determines the size at which species are 
#' selected.
#' 
#' The resulting `MizerParams` object can be projected forward using 
#' \code{project()} like any other `MizerParams` object. When projecting 
#' the community model it may be necessary to keep a small time step size
#' `dt` of around 0.1 to avoid any instabilities with the solver. You can
#' check for these numerical instabilities by plotting the biomass or abundance
#' through time after the projection.
#' 
#' @param z0 The background mortality of the community. Default value is 0.1.
#' @param alpha The assimilation efficiency of the community. Default value 0.2
#' @param f0 The average feeding level of individuals who feed on a power-law 
#'   spectrum. This value is used to calculate the search rate parameter 
#'   `gamma` (see the package vignette). Default value is 0.7.
#' @param h The maximum food intake rate. Default value is 10.
#' @param beta The preferred predator prey mass ratio. Default value is 100.
#' @param sigma The width of the prey preference. Default value is 2.0.
#' @param q The search volume exponent. Default value is 0.8.
#' @param n The scaling of the intake. Default value is 2/3.
#' @param kappa The carrying capacity of the resource spectrum. Default value
#'   is 1000.
#' @param lambda The exponent of the resource spectrum. Default value is 2 + q
#'   - n.
#' @param r_pp Growth rate parameter for the resource spectrum. Default value is 10.
#' @param gamma Volumetric search rate. Estimated using `h`, `f0` and 
#'   `kappa` if not supplied.
#' @param recruitment The constant recruitment in the smallest size class of the
#'   community spectrum. This should be set so that the community spectrum 
#'   continues the resource spectrum. Default value = `kappa` * 
#'   `min_w`^-`lambda`.
#' @param rec_mult Additional multiplier for the constant recruitment. Default 
#'   value is 1.
#' @param knife_edge_size The size at the edge of the knife-selectivity 
#'   function. Default value is 1000.
#' @param knife_is_min Is the knife-edge selectivity function selecting above 
#'   (TRUE) or below (FALSE) the edge. Default is TRUE.
#' @param max_w The maximum size of the community. The `w_inf` of the 
#'   species used to represent the community is set to this value. The 
#'   default value is 1e6.
#' @param min_w The minimum size of the community. Default value is 1e-3.
#' @param min_w_pp The smallest size of the resource spectrum. 
#' @param ... Other arguments to pass to the `MizerParams` constructor.
#' @export
#' @return An object of type \code{\linkS4class{MizerParams}}
#' @references K. H. Andersen,J. E. Beyer and P. Lundberg, 2009, Trophic and 
#'   individual efficiencies of size-structured communities, Proceedings of the 
#'   Royal Society, 276, 109-114
#' @family deprecated functions
#' @examples
#' \dontrun{
#' params <- set_community_model(f0=0.7, z0=0.2, recruitment=3e7)
#' sim <- project(params, effort = 0, t_max = 100, dt=0.1)
#' plotBiomass(sim)
#' plotSpectra(sim)
#' }
set_community_model <- function(max_w = 1e6,
                                min_w = 1e-3,
                                min_w_pp = 1e-10,
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
        constant_reproduction = recruitment * rec_mult # to be used in the SRR
    )
    com_params <- MizerParams(com_params_df, p = p, n = n, q = q, lambda = lambda, 
                              kappa = kappa, min_w = min_w, max_w = max_w, 
                              min_w_pp = min_w_pp,
                              w_pp_cutoff = w_pp_cutoff, r_pp = r_pp, ...)
    com_params@rates_funcs$RDD <- "constantRDD"
    com_params@psi[] <- 0 # Need to force to be 0. Can try setting w_mat but 
    # due to slope still not 0
    # Set w_mat to NA for clarity - it is not actually being used
    com_params@species_params$w_mat[] <- NA
    return(com_params)
}

#### getPhiPrey ####
#' Get available energy
#' 
#' This is deprecated and is no longer used by the mizer project() method.
#' Calculates the amount \eqn{E_{a,i}(w)} of food exposed to each predator as
#' a function of predator size. 
#' 
#' @param object An \linkS4class{MizerParams} object
#' @param n A matrix of species abundances (species x size)
#' @param n_pp A vector of the background abundance by size
#' @param ... Other arguments (currently unused)
#'   
#' @return A two dimensional array (predator species x predator size)
#' @seealso [project()]
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' n <- sim@@n[21,,]
#' n_pp <- sim@@n_pp[21,]
#' getPhiPrey(params,n,n_pp)
#' }
getPhiPrey <- function(object, n, n_pp, ...) {
    phi_prey <- getEncounter(object, n, n_pp) / object@search_vol
}
