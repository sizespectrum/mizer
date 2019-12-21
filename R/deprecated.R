#' Deprecated obsolete function for setting up multispecies parameters
#' @inheritParams newMultispeciesParams
#' @export
#' @family deprecated functions
set_multispecies_model <- 
    function(
        species_params,
        interaction = matrix(1,
                             nrow = nrow(species_params),
                             ncol = nrow(species_params)),
        min_w_pp = 1e-10,
        max_w = max(species_params$w_inf)*1.1,
        n = 2 / 3,
        q = 0.8,
        f0 = 0.6,
        kappa = 1e11,
        lambda = 2 + q - n,
        ...) {
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
        message("\tNote: No sel_func column in species data frame. Setting selectivity to be 'knife_edge' for all species.")
        object$sel_func <- 'knife_edge'
        # Set default selectivity size
        if (!("knife_edge_size" %in% colnames(object))) {
            message("Note: \tNo knife_edge_size column in species data frame. Setting knife edge selectivity equal to w_mat.")
            object$knife_edge_size <- object$w_mat
        }
    }
    # If no catchability column in species_params, set to 1
    if (!("catchability" %in% colnames(object)))
        object$catchability <- 1
    # Sort out h column If not passed in directly, is calculated from f0 and
    # k_vb if they are also passed in
    if (!("h" %in% colnames(object))) {
        message("Note: \tNo h column in species data frame so using f0 and k_vb to calculate it.")
        if (!("k_vb" %in% colnames(object))) {
            stop("\t\tExcept I can't because there is no k_vb column in the species data frame")
        }
        object$h <- ((3 * object$k_vb) / (object$alpha * f0)) * (object$w_inf ^ (1/3))
    }
    # Sorting out gamma column
    if (!("gamma" %in% colnames(object))) {
        message("Note: \tNo gamma column in species data frame so using f0, h, beta, sigma, lambda and kappa to calculate it.")
        ae <- sqrt(2*pi) * object$sigma * object$beta^(lambda-2) * exp((lambda-2)^2 * object$sigma^2 / 2)
        object$gamma <- (object$h / (kappa * ae)) * (f0 / (1 - f0))
    }
    # Sort out ks column
    if (!("ks" %in% colnames(object))) {
        message("Note: \tNo ks column in species data frame so using ks = h * 0.2.")
        object$ks <- object$h * 0.2
    }
    
    # The m column did not exist in the old version, set it to 1
    if (!("m" %in% colnames(object))) {
        message("Note: \tNo m column in species data frame so using m = 1.")
        object$m <- 1
    }
    
    return(newMultispeciesParams(object,
                                 interaction = interaction,
                                 n = n,
                                 q = q,
                                 f0 = f0,
                                 kappa = kappa,
                                 lambda = lambda,
                                 max_w = max_w,
                                 min_w_pp = min_w_pp,
                                 ...))
}

#' Alias for set_multispecies_model
#' 
#' An alias provided for backward compatibility with mizer version <= 1.0
#' @inherit set_multispecies_model
#' @export
MizerParams <- set_multispecies_model


#' Deprecated obsolete function for setting up parameters for a trait-based
#' model
#' 
#' @inheritParams newTraitParams
#' @export
#' @family deprecated functions
# Copied from version 1.0.1
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
        species = as.factor(1:no_sp),
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

#' Deprecated obsolete function for setting up parameters for a community
#' model
#' 
#' @inheritParams newCommunityParams
#' @export
#' @family deprecated functions
# Copied from version 1.0.1
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
    com_params <- MizerParams(com_params_df, p = p, n = n, q = q, lambda = lambda, 
                              kappa = kappa, min_w = min_w, max_w = max_w, 
                              w_pp_cutoff = w_pp_cutoff, r_pp = r_pp, ...)
    com_params@srr <- "srrConstant"
    com_params@psi[] <- 0 # Need to force to be 0. Can try setting w_mat but 
    # due to slope still not 0
    # Set w_mat to NA for clarity - it is not actually being used
    com_params@species_params$w_mat[] <- NA
    return(com_params)
}