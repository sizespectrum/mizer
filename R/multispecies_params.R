#' Construct MizerParams object for multispecies model
#'
#' Provides default functional forms for all slots in the
#' \linkS4class{MizerParams} object based on user-provided species parameters.
#' 
#' @inheritParams emptyParams
#' @param min_w_pp The smallest size of the plankton spectrum. By default this
#'   is set to the smallest value at which any of the consumers can feed.
#' @inheritParams setInteraction
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
#' @inheritParams setFishing
#' @inheritParams setSearchVolume
#' @inheritParams setIntakeMax
#' @inheritParams setMetab
#' @inheritParams setBMort
#' @inheritParams setReproduction
#' @inheritParams setPlankton
#' @inheritParams setResourceEncounter
#'
#' @return An object of type \linkS4class{MizerParams}
#' 
#' @section Species parameters:
#' The only essential argument is a data frame which contains the species data.
#' The data frame is arranged species by parameter, so each column of the
#' parameter data frame is a parameter and each row has the parameters for one
#' of the species in the model.
#'
#' There are some essential columns that must be included in the species
#' parameter data.frame and that do not have default values. Other columns do
#' have default values, so that if they are not included in the species
#' parameter data frame, they will be automatically added when the
#' \code{MizerParams} object is created. See the accompanying vignette for
#' details of these columns.
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
#' @inheritSection setInteraction Setting interactions
#' @inheritSection setPredKernel Setting predation kernel
#' @inheritSection setSearchVolume Setting search volume
#' @inheritSection setIntakeMax Setting maximum intake rate
#' @inheritSection setMetab Setting metabolic rate
#' @inheritSection setBMort Setting background mortality rate
#' @inheritSection setReproduction Setting reproduction
#' @inheritSection setPlankton Setting plankton dynamics
#' @inheritSection setResources Setting resource dynamics
#' @inheritSection setResourceEncounter Setting resource encounter rate
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
                                   interaction_p = rep(1, nrow(species_params)),
                                   no_w = 100,
                                   min_w = 0.001,
                                   max_w = NA,
                                   min_w_pp = NA,
                                   no_w_pp = NA,
                                   n = 2 / 3,
                                   p = 0.7,
                                   q = 0.8,
                                   kappa = 1e11,
                                   lambda = (2 + q - n),
                                   f0 = 0.6,
                                   # setPredKernel()
                                   pred_kernel = NULL,
                                   pred_kernel_type = "lognormal",
                                   store_kernel = FALSE,
                                   # setSearchVolume()
                                   gamma = NULL,
                                   # setIntakeMax()
                                   h = NULL,
                                   # setMetab()
                                   metab = NULL,
                                   # setBMort
                                   z0pre = 0.6,
                                   z0exp = n - 1,
                                   # setReproduction
                                   psi = NULL,
                                   # setPlankton
                                   r_pp = 10,
                                   w_pp_cutoff = 10,
                                   plankton_dynamics = plankton_semichemostat,
                                   # setResources
                                   resource_dynamics = list(),
                                   resource_params = list(),
                                   # setResourceEncounter
                                   rho = NULL,
                                   srr = srrBevertonHolt) {
    assert_that(is.data.frame(species_params))
    no_sp <- nrow(species_params)
    
    ## Determine min_w_pp ----
    # If not provided, set min_w_pp so that all fish have their full feeding 
    # kernel inside plankton spectrum
    getMaxPPMR <- get0(paste0(pred_kernel_type, "_max_ppmr"))
    if (is.function(getMaxPPMR)) {
        # First we need to set w_min if missing because we use it below
        if (!("w_min" %in% colnames(species_params))) {
            species_params$w_min <- rep(NA, no_sp)
        }
        missing <- is.na(species_params$w_min)
        if (any(missing)) {
            species_params$w_min[missing] <- min_w
        }
        min_w_feeding <- species_params$w_min / getMaxPPMR(species_params)
        if (is.na(min_w_pp)) {
            min_w_pp <- min(min_w_feeding)
        } else {
            hungry_sp <- species_params$species[min_w_feeding < min_w_pp]
            if (length(hungry_sp) > 0) {
                message(paste(
                    "Note: The following species have feeding kernels that extend",
                    "below the smallest plankton size specified by min_w_pp:",
                    toString(hungry_sp)))
            }
        }
    } else {
        if (is.na(min_w_pp)) {
            stop(paste0("You need to explicitly specify the minimum plankton ",
            "size via the min_w_pp argument or you need to define a function ",
            getMaxPPMR))
        }
    }
    
    ## Create MizerParams object ----
    params <- emptyParams(species_params,
                          no_w = no_w, 
                          min_w = min_w,  
                          max_w = max_w, 
                          min_w_pp = min_w_pp, 
                          no_w_pp = NA,
                          srr = srr)
    
    ## Fill the slots ----
    params@n <- n
    params@p <- p
    params@lambda <- lambda
    params@q <- q
    params@f0 <- f0
    params@kappa <- kappa
    
    params <- setInteraction(params, 
                             interaction = interaction, 
                             interaction_p = interaction_p)
    params <- setFishing(params)
    params <- setPredKernel(params, 
                            pred_kernel = pred_kernel,
                            pred_kernel_type = pred_kernel_type,
                            store_kernel = store_kernel)
    params <- setIntakeMax(params, 
                           h = h)
    params <- setMetab(params, 
                       metab = metab)
    params <- setBMort(params, 
                       z0pre = z0pre, 
                       z0exp = z0exp)
    params <- setSearchVolume(params, 
                              gamma = gamma)
    params <- setReproduction(params, 
                              psi = psi)
    params <- setPlankton(params, 
                          r_pp = r_pp, 
                          w_pp_cutoff = w_pp_cutoff,
                          plankton_dynamics = plankton_dynamics)
    params <- setResources(params,
                           resource_dynamics = resource_dynamics,
                           resource_params = resource_params)
    params <- setResourceEncounter(params, 
                                   rho = rho)
    
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