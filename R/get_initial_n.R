#' Calculate initial population abundances
#' 
#' This function uses the model parameters and other parameters to calculate 
#' initial values for the species number densities. These initial 
#' abundances are currently quite arbitrary and not close to the steady state.
#' We intend to improve this in the future.
#' 
#' @param params The model parameters. An object of type \linkS4class{MizerParams}.
#' @param a A parameter with a default value of 0.35.
#' @param n0_mult Multiplier for the abundance at size 0. Default value is
#'   kappa/1000.
#' @export
#' @concept helper
#' @return A matrix (species x size) of population abundances.
#' @examples
#' \dontrun{
#' params <- newMultispeciesParams(NS_species_params_gears)
#' init_n <- get_initial_n(params)
#' }
get_initial_n <- function(params, n0_mult = NULL, a = 0.35) {
    if (!is(params,"MizerParams"))
        stop("params argument must of type MizerParams")
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)
    initial_n <- params@initial_n
    
    if (defaults_edition() < 2) {
        # N = N0 * Winf^(2*n-q-2+a) * w^(-n-a)
        # Reverse calc n and q from intake_max and search_vol slots (could add get_n function)
        n <- params@species_params[[1, "n"]]
        q <- params@species_params[[1, "q"]]
        # Guessing at a suitable n0 value based on kappa - this was figured out 
        # using trial and error and should be updated
        if (is.null(n0_mult)) {
            lambda <- 2 + q - n
            kappa <- params@cc_pp[1] / (params@w_full[1]^(-lambda))
            n0_mult <- kappa / 1000
        }
        initial_n[] <- unlist(tapply(params@w, 1:no_w, function(wx,n0_mult,w_max,a,n,q)
            n0_mult * w_max^(2 * n - q - 2 + a) * wx^(-n - a),
            n0_mult = n0_mult, w_max = params@species_params$w_max, a=a, n=n, q=q))
        #set densities at w > w_max to 0
        initial_n[unlist(tapply(params@w,1:no_w,function(wx,w_max) w_max<wx, w_max=params@species_params$w_max))] <- 0
        # Also any densities at w < w_min set to 0
        initial_n[unlist(tapply(params@w,1:no_w,function(wx,w_min)w_min>wx, w_min=params@species_params$w_min))] <- 0    
        return(initial_n)
    }
    
    p <- params
    p@initial_n[] <- 0
    p@initial_n_pp <- p@resource_params$kappa * 
        p@w_full ^ (-p@resource_params$lambda)
    p@interaction[] <- 0
    income <- getEReproAndGrowth(p) + p@metab
    for (i in seq_len(no_sp)) {
        # At small sizes the income should be A w^n. Determine A
        # Use w_min_idx + 1 in case user has implemented reduced growth
        # for the smallest size class (see e.g. #241)
        iw <- p@w_min_idx[i] + 1
        A <- income[i, iw] / (p@w[iw] ^ p@species_params[[i, "n"]])
        
        mort <- 0.4 * A * p@w ^ (p@species_params[[i, "n"]] - 1) + getFMort(p)
        growth <- getEGrowth(p)[i, ]
        
        idxs <- p@w_min_idx[i]:(min(which(c(growth, 0) <= 0)) - 1)
        idx <- idxs[1:(length(idxs) - 1)]
        # Steady state solution of the upwind-difference scheme used in project
        p@initial_n[i, idxs] <- 
            c(1, cumprod(growth[idx] / ((growth + mort * p@dw)[idx + 1])))
    }
    p <- matchBiomasses(p)
    return(p@initial_n)
}
