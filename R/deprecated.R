#' Deprecated obsolete function for setting up multispecies parameters
#' @inheritParams newMultispeciesParams
#' @export
#' @family functions for setting up models
set_multispecies_model <- 
    function(
        species_params,
        interaction = matrix(1,
                             nrow = nrow(species_params),
                             ncol = nrow(species_params)),
        n = 2 / 3,
        q = 0.8,
        f0 = 0.6,
        kappa = 1e11,
        lambda = 2 + q - n,
        ...) {
    object <- species_params
    # old code from MizerParams() in version 1.0.1
    
    # Set default values for column values if missing
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
                                 lambda = lambda, ...))
}


#' Alias for set_multispecies_model
#' 
#' An alias provided for backward compatibility with mizer version <= 1.0
#' @inherit set_multispecies_model
#' @export
MizerParams <- set_multispecies_model