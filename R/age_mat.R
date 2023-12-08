#' Calculate age at maturity from von Bertalanffy growth parameters
#' 
#' This is not a good way to determine the age at maturity because the von
#' Bertalanffy growth curve is not reliable for larvae and juveniles. However
#' this was used in previous versions of mizer and is supplied for 
#' backwards compatibility.
#'
#' Uses the age at maturity that is implied by the von Bertalanffy growth curve
#' specified by the `w_inf`, `k_vb`, `t0`, `a` and `b` parameters in the
#' species_params data frame.
#'
#' If any of `k_vb` is missing for a species, the function returns NA for that
#' species. Default values of `b = 3` and `t0 = 0` are used if these are
#' missing. If `w_inf` is missing, `w_max` is used instead.
#'
#' @param object A MizerParams object or a species_params data frame
#' @return A named vector. The names are the species names and the values are
#'   the ages at maturity.
#' @concept helper
age_mat_vB <- function(object) {
    if (is(object, "MizerParams")) {
        sp <- object@species_params
    } else {
        if (!is.data.frame(object)) {
            stop("The first argument must be either a MizerParams object or a species_params data frame.")
        }
        sp <- completeSpeciesParams(object)
    }
    sp <- set_species_param_default(sp, "t0", 0)
    sp <- set_species_param_default(sp, "b", 3)
    sp <- set_species_param_default(sp, "k_vb", NA)
    sp <- set_species_param_default(sp, "w_inf", sp$w_max)
    
    a_mat <- -log(1 - (sp$w_mat / sp$w_inf) ^ (1/sp$b)) / sp$k_vb + sp$t0
    names(a_mat) <- sp$species
    a_mat
}

#' Calculate age at maturity
#' 
#' Uses the growth rate and the size at maturity to calculate the age at 
#' maturity
#' 
#' Using that by definition of the growth rate \eqn{g(w) = dw/dt} we have that
#' \deqn{\mathrm{age_{mat}} = \int_0^{w_{mat}.}\frac{dw}{g(w)}}{age_mat = \int_0^w_mat 1/g(w) dw.}
#' 
#' @param params A MizerParams object
#' @return A named vector. The names are the species names and the values are
#'   the ages at maturity.
#' @export
#' @concept helper
#' @examples 
#' age_mat(NS_params)
age_mat <- function(params) {
    assert_that(is(params, "MizerParams"))
    sp <- params@species_params
    no_sp <- nrow(sp)
    
    g <- getEGrowth(params)
    a_mat <- vector("double", no_sp)
    names(a_mat) <- sp$species
    for (i in seq_along(a_mat)) {
        sel <- params@w < sp$w_mat[i]
        a_mat[i] <- sum(params@dw[sel] / g[i, sel])
    }
    
    a_mat
}