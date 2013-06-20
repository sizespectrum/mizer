#' Sets up parameters for a community-type model
#'
#' This functions creates a \code{MizerParams} object so that community-type models can be easily set up and run.
#' A community model has several features that distinguish it from the food-web type models.
#' Only one 'species' is resolved, i.e. one 'species' is used to represent the whole community.
#' The resource spectrum only extends to the start of the community spectrum.
#' Recruitment to the smallest size in the community spectrum is constant and set by the user.
#' As recruitment is constant, the proportion of energy invested in reproduction (the slot \code{psi} of the
#' returned \code{MizerParams} object) is set to 0.
#' Standard metabolism has been turned off (the parameter \code{ks} is set to 0).
#' Consequently, the growth rate is now determined solely by the assimilated food (see the package Vignette for more details).
#'
#' The function has many arguments, all of which have default values. The main arguments that the users should be concerned with are \code{z0}, \code{recruitment}, \code{alpha} and \code{f0} as these determine the average growth rate of the community.
#'
#' Fishing selectivity is modelled as a sigmoid function with two parameters: \code{l25} and \code{l50} which determine the lengths at which
#' selectivity is 0.25 and 0.5 respectively. Lengths are converted to weights using the default parameters a = 0.001 and b = 3.0.
#' 
#' The resulting \code{MizerParams} object can be projected forward using \code{project()} like any other \code{MizerParams} object.
#' When projecting the community model it may be necessary to reduce \code{dt} to 0.1 to avoid any instabilities with the solver. You can check this by plotting the biomass or abundance through time after the projection.
#' @param z0 The background mortality of the community. The default value is 0.1.
#' @param alpha The assimilation efficiency of the community. The default value is 0.2 (from Andersen et. al., 2009).
#' @param recruitment The constant recruitment in the smallest size class of the community spectrum. This should be set so that the community spectrum continues the background spectrum.
#' @param f0 The average feeding level of individuals who feed mainly on the resource. This value is to used to calculate the search rate parameter \code{ga,,a} (see the package Vignette). The default value is 0.7.
#' @param h The maximum food intake rate. The default value is 10.
#' @param beta The preferred predator prey mass ratio. The default value is 100.
#' @param sigma The width of the prey preference. The default value is 2.0.
#' @param q The search volume exponent. The default value is 0.8.
#' @param n The scaling of the intake. The default value is 2/3.
#' @param kappa The carrying capacity of the background spectrum. The default value is 1000.
#' @param lambda The exponent of the background spectrum. The default value is 2 + q - n.
#' @param a The length-weight coefficient. The default value is 0.01.
#' @param b The length-weight exponent. The default value is 3.0.
#' @param l25 The size at which fishing selectivity is 25%. The default value is 10.
#' @param l50 The size at which fishing selectivity is 50%. The default value is 40.
#' @param max_w The maximum size of the community. The \code{w_inf} of the species used to represent the community is set to 0.9 * this value. The default value is 1e6.
#' @param min_w The minimum size of the community. The default value is 1e-3.
#' @export
#' @return An object of type \code{MizerParams}
#' @seealso \code{\link{MizerParams}}
#' @references K. H. Andersen,J. E. Beyer and P. Lundberg, 2009, Trophic and individual efficiencies of size-structured communities, Proceedings of the Royal Society, 276, 109-114
#' @examples
#' params <- set_community_model(f0=0.7, z0=0.2, recruitment=3e7)
#' sim <- project(params, effort = 0, t_max = 100, dt=0.1)
#' plotBiomass(sim)
#' plotSpectra(sim)
set_community_model <- function(max_w = 1e6,
                                min_w = 1e-3,
                                z0 = 0.1,
                                recruitment = 4e7,
                                alpha = 0.2,
                                h = 10,
                                beta = 100,
                                sigma = 2.0,
                                q = 0.8,
                                n = 2/3,
                                kappa = 1000,
                                lambda = 2+q-n,
                                f0 = 0.7,
                                a = 0.01, b = 3.0, l25 = 10, l50 = 40,
                                ...
                                ){
    w_inf <- max_w * 0.9
    w_pp_cutoff <- min_w
    ks <- 0 # Turn off standard metabolism
    p <- n # But not used as ks = 0
    gamma <- (f0 * h * beta^(2-lambda)) / ((1-f0)*sqrt(2*pi)*kappa*sigma)
    # Make the species data.frame
    com_params_df <- data.frame(
        species = "Community",
        w_inf = w_inf,
        w_mat = 1e12, # Has no affect as psi set to 0 but we set it to something to help the constructor
        h = h, # max food intake
        gamma = gamma,# vol. search rate,
        ks = ks,# standard metabolism coefficient,
        beta = beta,
        sigma = sigma,
        z0 = z0, # background mortality
        alpha = alpha,
        l25 = l25, 
        l50 = l50, 
        a = a, 
        b = b, 
        erepro = 1, # not used
        constant_recruitment = recruitment # to be used in the SRR
    )
    # Set the recruitment function for constant recruitment
    constant_recruitment <- function(rdi, species_params){
        return(species_params$constant_recruitment)
    }
    com_params <- MizerParams(com_params_df, p=p, n=n,q=q, lambda = lambda, kappa = kappa, min_w = min_w, max_w = max_w, w_pp_cutoff = w_pp_cutoff, ...)
    com_params@srr <- constant_recruitment
    com_params@psi[] <- 0 # Need to force to be 0. Can try setting w_mat but due to slope still not 0
    # Set w_mat to NA for clarity - it is not actually being used
    com_params@species_params$w_mat[] <- NA
    return(com_params)
}





#' Sets up parameters for a trait-based model
#'
#' This functions creates a \code{MizerParams} object so that community-type models can be easily set up and run.
#' The trait-based size spectrum model can be derived as a simplification of the general size-based model used in \code{mizer}.
#' All the species-specific parameters are the same, except for the asymptotic size, which is considered the most important trait characterizing a species.
#' Other parameters are related to the asymptotic size.
#' For example, the size at maturity is given by w_inf * eta, where eta is the same for all species.
#' For the trait-based model the number of species is not important.
#' For applications of the trait-based model see Andersen & Pedersen (2010).
#' See the \code{mizer} vignette for more details and examples of the trait-based model.
#'
#' The function has many arguments, all of which have default values.
#' Of particular interest to the user are the number of species in the model and the minimum and maximum asymptotic sizes.
#' The asymptotic sizes of the species are spread evenly on a logarithmic scale within this range.
#'
#' The stock recruitment relationship in the trait-based model is a 'hockey-stick' (following Andersen and Pedersen, 2009).
#' The maximum value the recruitment can be controlled using the \code{k0} parameter.
#' Users should adjust this value to get the spectra they want.
#'
#' The factor for the search volume, \code{gamma}, is calculated using the expected feeding level, \code{f0}.
#' 
#' Fishing selectivity is modelled as a sigmoid function with two parameters: \code{l25} and \code{l50} which determine the lengths at which
#' selectivity is 0.25 and 0.5 respectively. Lengths are converted to weights using the default parameters a = 0.001 and b = 3.0.
#'
#' The resulting \code{MizerParams} object can be projected forward using \code{project()} like any other \code{MizerParams} object.
#' When projecting the community model it may be necessary to reduce \code{dt} to 0.1 to avoid any instabilities with the solver. You can check this by plotting the biomass or abundance through time after the projection.
#' @param no_sp The number of species in the model. The default value is 10. The more species, the longer takes to run.
#' @param min_w_inf The asymptotic size of the smallest species in the community.
#' @param max_w_inf The asymptotic size of the largest species in the community.
#' @param no_w The number of size bins in the community spectrum.
#' @param min_w The smallest size of the community spectrum.
#' @param max_w The largest size of the community spectrum. Default value is the largest w_inf in the community x 1.1.
#' @param min_w_pp The smallest size of the background spectrum.
#' @param no_w_pp The number of the extra size bins in the background spectrum (i.e. the difference between the number of sizes bins in the community spectrum and the full spectrum).
#' @param w_pp_cutoff The cut off size of the background spectrum. Default value is 1.
#' @param k0 Scaling parameter for the maximum recruitment.
#' @param n Scaling of the intake. Default value is 2/3. 
#' @param p Scaling of the standard metabolism. Default value is 0.7. 
#' @param q Exponent of the search volume. Default value is 0.8. 
#' @param eta Factor to calculate \code{w_mat} from asymptotic size.
#' @param r_pp Growth rate of the primary productivity. Default value is 10. 
#' @param kappa Carrying capacity of the resource spectrum. Default value is 1e11. 
#' @param lambda Exponent of the resource spectrum. Default value is (2+q-n). 
#' @param alpha The assimilation efficiency of the community. The default value is 0.2 (from Andersen et. al., 2010).
#' @param ks Standard metabolism coefficient. Default value is 4.
#' @param z0 The background mortality of the community. The default value is 0.6.
#' @param h Maximum food intake rate. Default value is 30.
#' @param beta Preferred predator prey mass ratio. Default value is 100.
#' @param sigma Width of prey size preference. Default value is 1.3.
#' @param f0 Expected average feeding level. Used to set \code{gamma}, the factor for the search volume.
#' @param l25 The size at which fishing selectivity is 25%. The default value is 10.
#' @param l50 The size at which fishing selectivity is 50%. The default value is 40.
#' @param a The length-weight coefficient. The default value is 0.01.
#' @param b The length-weight exponent. The default value is 3.0.
#' @export
#' @return An object of type \code{MizerParams}
#' @seealso \code{\link{MizerParams}}
#' @references K. H. Andersen and M. Pedersen, 2010, Damped trophic cascades driven by fishing in model marine ecosystems. Proceedings of the Royal Society V, Biological Sciences, 1682, 795-802.
#' @examples
#' trait_params <- set_trait_model(no_sp = 15)
#' init_pop <- get_initial_n(trait_params, n0_mult = 0.001)
#' sim <- project(trait_params, effort = 0, t_max = 50, dt=0.2, initial_n = init_pop, t_save = 1)
#' plot(sim)
set_trait_model <- function(no_sp = 10,
                            min_w_inf = 10,
                            max_w_inf = 1e5,
                            no_w = 100,
                            min_w = 0.001,
                            max_w = max_w_inf * 1.1,
                            min_w_pp = 1e-10,
                            no_w_pp = round(no_w)*0.3,
                            w_pp_cutoff = 1,
                            k0 = 100, # recruiment adjustment parameter
                            n = 2/3,
                            p = 0.75,
                            q = 0.9, 
                            eta = 0.25,
                            r_pp = 4,
                            kappa = 0.005,
                            lambda = 2+q-n,
                            alpha = 0.6,
                            ks = 4,
                            z0 = 0.6,
                            h = 30,
                            beta = 100,
                            sigma = 1.3,
                            f0 = 0.5,
                            l25 = 10,
                            l50 = 40,
                            a = 0.001, 
                            b = 3, 
                            ...){
    alpha_e <- sqrt(2*pi) * sigma * beta^(lambda-2) * exp((lambda-2)^2 * sigma^2 / 2) # see A&P 2009
    gamma <- h * f0 / (alpha_e * kappa * (1-f0)) # see A&P 2009 
    w_inf <- 10^seq(from=log10(min_w_inf), to = log10(max_w_inf), length=no_sp)
    w_mat <- w_inf * eta

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
            z0 = z0 * w_inf^(n-1), # background mortality
            alpha = alpha,
            l25 = l25,
            l50 = l50,
            a = a, 
            b = b, 
            erepro = 1 # not used but included out of necessity
    )
    # Make the MizerParams
    trait_params <- MizerParams(trait_params_df, min_w = min_w, max_w=max_w, no_w = no_w, min_w_pp = min_w_pp, w_pp_cutoff = w_pp_cutoff, n = n, p=p, q=q, r_pp=r_pp, kappa=kappa, lambda = lambda) 
    # Sort out maximum recruitment - see A&P 2009
    a_rec <- (f0* h) / (alpha * f0 * h - ks) * beta ^(2*n-q-1) * exp(2*n*(q-1)-q^2+1) * sigma^2 / 2
    dm_idx <- aaply(w_inf, 1, function(x) max(which(x >= trait_params@w)))
    dm <- trait_params@dw[dm_idx]
    # k0 is passed in by the user
    N0 <- k0 * w_inf^(2*n-q-2+a_rec)/ dm
    # Add to data frame so we can use it in the SRR
    trait_params@species_params$max_rec <- N0 #* 1e4
    hockey_stick_recruitment <- function(rdi, species_params){
        rdd <- pmin(rdi,species_params$max_rec) 
        return(rdd)
    }
    trait_params@srr <- hockey_stick_recruitment
    return(trait_params)
}

