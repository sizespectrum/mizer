#' Get maximum predator/prey mass ratio for each species
#' 
#' Because of how we cut of the Gaussian feeding kernel, the maximum
#' predator/prey mass ratio is \eqn{\beta * \exp(3 \sigma)}.
#' 
#' This is used in \code{\link{emptyParams}} when setting up the size grid to
#' determine the smallest relevant plankton size. The reason we have a function
#' to calculate this is that it makes it easier for users to modify when they
#' choose their own feeding kernel.
#' 
#' @param species_params
#' 
#' @return Vector with the maximum predator/prey ratio for each species
#' @export
getMaxPPMR <- function(species_params) {
    #TODO: check validity of beta and sigma
    species_params$beta * exp(3 * species_params$sigma)
}


#' Set default feeding kernel
#' 
#' Still need to document
#' 
#' @param params MizerParams
#' @param store_kernel A boolean flag that determines whether the full
#'   feeding kernel is stored. If FALSE, only its Fourier transforms are stored.
#'   The default is TRUE if the number of size bins is no larger than 100 and
#'   FALSE otherwise.
#' 
#' @return MizerParams
#' @export
defaultPredKernel <- function(params, store_kernel) {
    species_params <- params@species_params
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
    Beta <- log(params@species_params$beta)
    sigma <- params@species_params$sigma
    # w_full has the weights from the smallest relevant plankton, to the largest fish
    x_full <- log(params@w_full)
    # We choose the origin of the x axis to be at the smallest plankton size
    x_full <- x_full - x_full[1]
    dx <- x_full[2] - x_full[1]
    # rr is the maximal log predator/prey mass ratio
    rr <- Beta + 3 * sigma
    ri <- floor(rr / dx)
    
    params@ft_pred_kernel_e <- matrix(0, nrow = no_sp, ncol = no_w_full)
    for (i in 1:no_sp) {
        # We compute the feeding kernel terms and their fft.
        phi <- exp(-(x_full - Beta[i])^2 / (2 * sigma[i]^2))
        phi[x_full > rr[i]] <- 0
        phi[1] <- 0
        # Fourier transform of feeding kernel for evaluating available energy
        params@ft_pred_kernel_e[i, ] <- fft(phi)
        # Fourier transform of feeding kernel for evaluating predation rate
        phi_p <- rep(0, no_w_full)
        phi_p[(no_w_full - ri[i] + 1):no_w_full] <- phi[(ri[i] + 1):2]
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


#' Set default resource encounter rate
#' 
#' @param params A MizerParams object
#' @param rho An array (species x resource) that gives the
#'   rate \eqn{\rho_{id}} that determines the rate at which species \eqn{i}
#'   encounters biomass of resource \eqn{d}. The rate is assumed to scale
#'   allometrically with size with exponent \code{n}. So the total contribution
#'   of the resources to the rate at which the individual encounters biomass is
#'   \deqn{\sum_d\rho_{id} w^n B_d,}
#'   where \eqn{B_d} is the biomass of the d-th unstructured resource component.
#'   See \code{\link{resource_dynamics}} help for more details.
#' 
#' @return A MizerParams object
#' @export
defaultResourceEncounter <- function(params, rho) {
    # Check validity of arguments
    assert_that(are_equal(length(dim(rho)), 2))
    no_res <- dim(rho)[2]
    assert_that(length(params@resource_dynamics) == no_res)
    # TODO: check that all parameters needed by resource dynamics functions
    # are included in resource_params
    if (dim(rho)[1] != no_sp) {
        stop("rho argument should have one row for each species.")
    }
    if (is.character(dimnames(rho)["res"])) {
        assert_that(are_equal(dimnames(rho)["res"],
                              names(resource_dynamics)))
    }
    params@rho[] <- outer(rho, params@w^params@n)
    params@initial_B[] <- rep(1, no_res)  # TODO: find better initial value
    
    return(params)
}


#' Set up plankton
#' 
#' @param params A MizerParams object
#' @param r_pp Growth rate of the primary productivity. Default is 10 g/year.
#' @param w_pp_cutoff The upper cut off size of the plankton spectrum. 
#'   Default is 10 g.
#' @param plankton_dynamics Function that determines plankton dynamics by
#'   calculating the plankton spectrum at the next time step from the current
#'   state. The default is \code{"\link{plankton_semichemostat}"}.
#' @param interaction_p Vector specifying for each species its interaction with
#'   plankton, similar to what the interaction matrix does for the interaction
#'   with other species. Entries should be numbers between 0 and 1. The 
#'   default is a vector of 1s.
#' 
#' @return A MizerParams object
#' @export
setPlankton <- function(params, r_pp = 10, w_pp_cutoff = 10,
                        plankton_dynamics = plankton_semichemostat,
                        interaction_p = rep(1, nrow(params@species_params))) {
    # TODO: check arguments
    # weight specific plankton growth rate
    params@rr_pp[] <- r_pp * params@w_full^(params@n - 1)
    # the plankton carrying capacity
    params@cc_pp[] <- params@kappa*params@w_full^(-params@lambda)
    params@cc_pp[params@w_full > w_pp_cutoff] <- 0
    params@plankton_dynamics <- plankton_dynamics
    params@interaction_p <- interaction_p
    
    return(params)
}


#' Set default reproduction proportion
#' 
#' Sets the proportion of the energy available for reproduction and growth that
#' is invested into reproduction as a function of the size of the individual.
#' 
#' The proportion is set to the product of a sigmoidal maturity ogive that 
#' gives the proportion of individuals of a given species and size that are
#' mature, and a factor that describes how investment into reproduction by mature
#' individuals scales with size. In formulas:
#' \deqn{\psi(w) = \left[1+\left(\frac{w}{w_{mat}}\right)^{-U}\right]^{-1}
#'   \left(\frac{w}{w_{inf}}\right)^{m-n}.}{[1+(w/w_mat)^(-U)]^(-1) * (w/w_inf)^(m - n)}
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
#' The result for \eqn{\psi(w)} for all species is stored in the \code{psi} slot
#' of the params object and this object is returned.
#' 
#' @param params A MizerParams object
#' 
#' @return A MizerParams object
#' @export
defaultReproProp <- function(params) {
    species_params <- params@species_params
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
    # Set w < 10% of w_mat to 0
    params@psi[unlist(
        tapply(params@w, 1:length(params@w),
               function(wx, w_mat) wx < (w_mat * 0.1),
               w_mat = species_params$w_mat))] <- 0
    # Set all w > w_inf to 1
    params@psi[unlist(
        tapply(params@w, 1:length(params@w),
               function(wx, w_inf) (wx/w_inf) > 1,
               w_inf = species_params$w_inf))] <- 1
    return(params)
}

#' Set default maximum intake rate
#' 
#' Still need to document
#' 
#' @param params MizerParams
#' 
#' @return MizerParams
#' @export
defaultIntakeMax <- function(params) {
    species_params <- params@species_params
    # If h column is not supplied, it is calculated from f0 and k_vb if they 
    # are supplied
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
    params@intake_max[] <- 
        unlist(tapply(params@w, 1:length(params@w), 
                      function(wx, h, n) h * wx^n,
                      h = species_params$h, n = params@n))
    return(params)
}

#' Set default search volume
#' 
#' Still need to document
#' 
#' @param params MizerParams
#' 
#' @return MizerParams
#' @export
defaultSearchVol <- function(params) {
    species_params <- params@species_params
    # Sorting out gamma column
    if (!("gamma" %in% colnames(species_params))) {
        species_params$gamma <- rep(NA, nrow(species_params))
    }
    if (any(is.na(species_params$gamma))) {
        message("Note: No gamma provided for some species, so using f0, h, beta, sigma, lambda and kappa to calculate it.")
        lm2 <- params@lambda - 2
        ae <- sqrt(2 * pi) * species_params$sigma * species_params$beta^lm2 *
            exp(lm2^2 * species_params$sigma^2 / 2) *
            # The factor on the following lines takes into account the cutoff
            # of the integral at 0 and at beta + 3 sigma
            (pnorm(3 - lm2 * species_params$sigma) + 
                 pnorm(log(species_params$beta)/species_params$sigma + 
                           lm2 * species_params$sigma) - 1)
        gamma <- (species_params$h / (params@kappa * ae)) * (params@f0 / (1 - params@f0))
        # Only overwrite missing gammas with calculated values
        missing <- is.na(species_params$gamma)
        if (any(is.na(gamma[missing]))) {
            stop("Could not calculate gamma.")
        }
        species_params$gamma[missing] <- gamma[missing]
        params@species_params$gamma <- species_params$gamma
    }
    params@search_vol[] <- 
        unlist(tapply(params@w, 1:length(params@w),
                      function(wx, gamma, q) gamma * wx^q,
                      gamma = species_params$gamma, q = params@q))
    return(params)
}

#' Set default metabolic rate
#' 
#' Still need to document
#' 
#' @param params MizerParams
#' 
#' @return MizerParams
#' @export
defaultMetab <- function(params) {
    species_params <- params@species_params
    
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
        unlist(tapply(params@w, 1:length(params@w),
                      function(wx, ks, k, p) ks * wx^p + k * wx,
                      ks = species_params$ks, k = species_params$k, p = params@p))
    return(params)
}

#' Set default background mortality rate
#' 
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
defaultBMort <- function(params, z0pre, z0exp) {
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


#' Construct MizerParams object for multispecies model
#'
#' Provides default functional forms for all slots in the
#' \linkS4class{MizerParams} object based on user-provided species parameters.
#' 
#' @inheritParams emptyParams
#' @param no_w_pp Obsolete argument that is no longer used because the number
#'    of plankton size bins is determined because all size bins have to
#'    be logarithmically equally spaced.
#' @param n Scaling of the intake. Default value is 2/3.
#' @param p Scaling of the standard metabolism. Default value is 0.7. 
#' @param q Exponent of the search volume. Default value is 0.8. 
#' @param kappa Carrying capacity of the plankton spectrum. Default value is
#'   1e11.
#' @param lambda Exponent of the plankton spectrum. Default value is (2+q-n).
#' @param f0 Average feeding level. Used to calculated \code{h} and \code{gamma}
#'   if those are not columns in the species data frame. Also requires
#'   \code{k_vb} (the von Bertalanffy K parameter) to be a column in the species
#'   data frame. If \code{h} and \code{gamma} are supplied then this argument is
#'   ignored. Default is 0.6.
#' @inheritParams defaultPredKernel
#' @inheritParams defaultResourceEncounter
#' @inheritParams setPlankton
#' @inheritParams defaultBMort
#' @param ... Additional arguments.
#'
#' @return An object of type \linkS4class{MizerParams}
#' 
#' @note The only essential argument is a data frame which contains the species
#'   data. The data frame is arranged species by parameter, so each column of
#'   the parameter data frame is a parameter and each row has the parameters for
#'   one of the species in the model.
#'   
#'   There are some essential columns that must be included in the species parameter
#'   data.frame and that do not have default values. Other columns do have
#'   default values, so that if they are not included in the species parameter
#'   data frame, they will be automatically added when the \code{MizerParams}
#'   object is created. See the accompanying vignette for details of these
#'   columns.
#'   
#' \if{html}{
#' The following table gives different components of a species parameters data.frame.
#' \tabular{lll}{
#'   Column name \tab Description \tab Default value \cr
#'   species \tab Name of the species \tab Compulsory (no default) \cr
#'   w_inf \tab The asymptotic mass of the species \tab Compulsory (no default) \cr
#'   w_mat \tab The maturation mass of the species \tab Compulsory (no default) \cr
#'   beta \tab Preferred predator prey mass ratio \tab Compulsory (no default) \cr
#'   sigma \tab Width of prey size preference \tab Compulsory (no default) \cr
#'   h \tab Maximum food intake rate. If this is not provided, it is calculated using the k_vb column. Therefore, either h or k_vb must be provided. \tab Optional (no default) \cr
#'   k_vb \tab The von Bertalanffy K parameter. Only used to calculate h if that column is not provided \tab Optional (no default) \cr
#'   gamma \tab Volumetric search rate. If this is not provided, it is calculated using the h column and other parameters. \tab Optional (no default) \cr
#'   ks \tab Standard metabolism coefficient \tab h*0.2 \cr
#'   z0 \tab Background mortality (constant for all sizes). If this is not provided then z0 is calculated as z0pre ∗ w_inf^z0exp. Also z0pre and z0exp have default values of 0.6 and -1/3 respectively. \tab Optional (no default) \cr
#'   k \tab Activity coefficient \tab 0 \cr
#'   alpha \tab Assimilation efficiency \tab 0.6 \cr
#'   erepro \tab Reproductive efficiency \tab 1 \cr
#'   w_min \tab The size class that recruits are placed in. \tab smallest size class of the species size spectrum
#' }}
#'   
#' @seealso \code{\link{project}} \linkS4class{MizerSim}
#' @export
#' @examples
#' data(NS_species_params_gears)
#' data(inter)
#' params <- set_multispecies_model(NS_species_params_gears, inter)
set_multispecies_model <- function(species_params,
                                   interaction = matrix(1,
                                                        nrow = nrow(species_params),
                                                        ncol = nrow(species_params)),
                                   no_w = 100,
                                   min_w = 0.001,
                                   max_w = NA,
                                   min_w_pp = NA,
                                   no_w_pp = NA,
                                   n = 2 / 3,
                                   p = 0.7,
                                   q = 0.8,
                                   r_pp = 10,
                                   kappa = 1e11,
                                   lambda = (2 + q - n),
                                   w_pp_cutoff = 10,
                                   f0 = 0.6,
                                   z0pre = 0.6,
                                   z0exp = n - 1,
                                   store_kernel = (no_w <= 100),
                                   plankton_dynamics = plankton_semichemostat,
                                   interaction_p = rep(1, nrow(species_params)),
                                   rho = NULL,
                                   resource_dynamics = list(),
                                   resource_params = list(),
                                   srr = srrBevertonHolt) {
    
    ## Make an empty MizerParams object of the right dimensions
    params <- emptyParams(species_params,
                          interaction = interaction,
                          no_w = no_w, 
                          min_w = min_w,  
                          max_w = max_w, 
                          min_w_pp = min_w_pp, 
                          no_w_pp = NA,
                          resource_dynamics = resource_dynamics,
                          resource_params = resource_params,
                          srr = srr)
    
    # Start filling the slots
    params@n <- n
    params@p <- p
    params@lambda <- lambda
    params@q <- q
    params@f0 <- f0
    params@kappa <- kappa
    
    params <- defaultPredKernel(params, store_kernel = store_kernel)
    params <- setPlankton(params, r_pp = r_pp, w_pp_cutoff = w_pp_cutoff,
                          plankton_dynamics = plankton_dynamics,
                          interaction_p = interaction_p)
    params <- defaultReproProp(params)
    params <- defaultIntakeMax(params)
    params <- defaultSearchVol(params)
    params <- defaultMetab(params)
    params <- defaultBMort(params, z0pre = z0pre, z0exp = z0exp)
    if (!is.null(rho)) {
        params <- defaultResourceEncounter(params, rho)
    }
    
    params@initial_n <- get_initial_n(params)
    params@initial_n_pp <- params@cc_pp
    params@A <- rep(1, nrow(species_params))
    
    return(params)
}

#' Alias for set_multispecies_model
#' 
#' An alias provided for backward compatibility with mizer version <= 1.0
#' @inherit set_multispecies_model
#' @export
MizerParams <- set_multispecies_model