

lognormal_pred_kernel <- function(ppmr, sp, params) {
    
    Beta <- log(params@species_params$beta[sp])
    sigma <- params@species_params$sigma[sp]
    phi <- exp(-(log(ppmr) - Beta)^2 / (2 * sigma^2))
    # rr is the maximal log predator/prey mass ratio
    rr <- exp(Beta + 3 * sigma)
    phi[ppmr > rr] <- 0
    # Do not allow feeding at own size
    phi[1] <- 0
    return(phi)
}