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

