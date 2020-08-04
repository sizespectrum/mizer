#' Calculate initial population abundances for the community populations
#' 
#' This function uses the model parameters and other parameters to calculate 
#' initial population abundances for the community populations. These initial 
#' abundances should be reasonable guesses at the equilibrium values. The 
#' returned population can be passed to the `project` function.
#' 
#' @param params The model parameters. An object of type \linkS4class{MizerParams}.
#' @param a A parameter with a default value of 0.35.
#' @param n0_mult Multiplier for the abundance at size 0. Default value is
#'   kappa/1000.
#' @export
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
    initial_n <- array(NA, dim = c(no_sp, no_w))
    dimnames(initial_n) <- dimnames(params@intake_max)
    # N = N0 * Winf^(2*n-q-2+a) * w^(-n-a)
    # Reverse calc n and q from intake_max and search_vol slots (could add get_n function)
    n <- params@species_params$n[[1]]
    q <- params@species_params$q[[1]]
    # Guessing at a suitable n0 value based on kappa - this was figured out 
    # using trial and error and should be updated
    if (is.null(n0_mult)) {
        lambda <- 2 + q - n
        kappa <- params@cc_pp[1] / (params@w_full[1]^(-lambda))
        n0_mult <- kappa / 1000
    }
    initial_n[] <- unlist(tapply(params@w, 1:no_w, function(wx,n0_mult,w_inf,a,n,q)
        n0_mult * w_inf^(2 * n - q - 2 + a) * wx^(-n - a),
        n0_mult = n0_mult, w_inf = params@species_params$w_inf, a=a, n=n, q=q))
    #set densities at w > w_inf to 0
    initial_n[unlist(tapply(params@w,1:no_w,function(wx,w_inf) w_inf<wx, w_inf=params@species_params$w_inf))] <- 0
    # Also any densities at w < w_min set to 0
    initial_n[unlist(tapply(params@w,1:no_w,function(wx,w_min)w_min>wx, w_min=params@species_params$w_min))] <- 0    
    return(initial_n)
}
