#' Set reproduction parameters
#' 
#' Sets the proportion of the total energy available for reproduction and growth
#' that is invested into reproduction as a function of the size of the
#' individual and sets additional density dependence.
#' 
#' @section Setting reproduction:
#' 
#' For each species and at each size, the proportion \eqn{\psi}{psi} of the 
#' available energy 
#' that is invested into reproduction is the product of two factors: the
#' proportion `maturity` of individuals that are mature and the proportion
#' `repro_prop` of the energy available to a mature individual that is 
#' invested into reproduction.
#' 
#' \subsection{Maturity ogive}{
#' If the the proportion of individuals that are mature is not supplied via
#' the `maturity` argument , then it is set to a sigmoidal 
#' maturity ogive that changes from 0 to 1 at around the maturity size:
#' \deqn{{\tt maturity}(w) = \left[1+\left(\frac{w}{w_{mat}}\right)^{-U}\right]^{-1}.}{
#'   maturity(w) = [1+(w/w_mat)^(-U)]^(-1)}
#' (To avoid clutter, we are not showing the species index in the equations,
#' although each species has its own maturity ogive.)
#' The maturity weights are taken from the `w_mat` column of the 
#' species_params data frame. Any missing maturity weights are set to 1/4 of the
#' asymptotic weight in the `w_inf` column.
#' 
#' The exponent \eqn{U} determines the steepness of the maturity ogive. By
#' default it is chosen as \eqn{U = 10}, however this can be overridden by
#' including a column \code{w_mat25} in the species parameter dataframe that
#' specifies the weight at which 25% of individuals are mature, which sets
#' \eqn{U = \log(3) / \log(w_{mat} / w_{25}).}{U = log(3) / log(w_mat / w_25).}
#' 
#' The sigmoidal function given above would strictly reach 1 only asymptotically.
#' Mizer instead sets the function equal to 1 already at the species' 
#' maximum size, taken from the compulsory `w_inf` column in the
#' `species_params` data frame. Also, for computational simplicity, any 
#' proportion smaller than `1e-8` is set to `0`.
#' }
#' 
#' \subsection{Investment into reproduction}{
#' If the the energy available to a mature individual that is 
#' invested into reproduction is not supplied via the `repro_prop` argument,
#' it is set to the allometric form
#' \deqn{{\tt repro\_prop}(w) = \left(\frac{w}{w_{inf}}\right)^{m-n}.}{
#'   repro_prop(w) = (w/w_inf)^(m - n).}
#' Here \eqn{n} is the scaling exponent of the energy income rate. Hence
#' the exponent \eqn{m} determines the scaling of the investment into
#' reproduction for mature individuals. By default it is chosen to be 
#' \eqn{m = 1} so that the rate at which energy is invested into reproduction 
#' scales linearly with the size. This default can be overridden by including a 
#' column `m` in the species parameter dataframe. The asymptotic sizes
#' are taken from the compulsory `w_inf` column in the species_params
#' data frame.
#' 
#' So finally we have
#' \deqn{\psi(w) = {\tt maturity}(w){\tt repro\_prop}(w)}{psi(w) = maturity(w) * repro_prop(w)}
#' }
#' 
#' \subsection{Reproductive efficiency}{
#' The reproductive efficiency \eqn{\epsilon}, i.e., the proportion of energy allocated to
#' reproduction that results in egg biomass, is set through the `erepro`
#' column in the species_params data frame. If that is not provided, the default
#' is set to 1 (which you will want to override). The offspring biomass divided
#' by the egg biomass gives the rate of egg production, returned by
#' [getRDI()]:
#' \deqn{R_{di} = \frac{\epsilon}{2 w_{min}} \int N(w)  E_r(w) \psi(w) \, dw}{R_di = (\epsilon/(2 w_min)) \int N(w)  E_r(w) \psi(w) dw}
#' }
#' 
#' \subsection{Density dependence}{
#' The stock-recruitment relationship is an emergent phenomenon in mizer, with
#' several sources of density dependence. Firstly, the amount of energy invested
#' into reproduction depends on the energy income of the spawners, which is
#' density-dependent due to competition for prey. Secondly, the proportion of
#' larvae that grow up to recruitment size depends on the larval mortality,
#' which depends on the density of predators, and on larval growth rate, which
#' depends on density of prey.
#' 
#' Finally, to encode all the density dependence in the stock-recruitment
#' relationship that is not already included in the other two sources of density
#' dependence, mizer puts the the density-independent rate of egg production
#' through a density-dependence function. The result is returned by
#' [getRDD()]. The name of the density-dependence function is
#' specified by the `RDD` argument. The default is the Beverton-Holt
#' function [BevertonHoltRDD()], which requires an `R_max` column
#' in the species_params data frame giving the maximum egg production rate. If
#' this column does not exist, it is initialised to `Inf`, leading to no
#' density-dependence. Other functions provided by mizer are
#' [RickerRDD()] and [SheperdRDD()] and you can easily use
#' these as models for writing your own functions.
#' }
#' @param params A MizerParams object
#' @param maturity Optional. An array (species x size) that holds the proportion
#'   of individuals of each species at size that are mature. If not supplied, a
#'   default is set as described in the section "Setting reproduction".
#' @param repro_prop Optional. An array (species x size) that holds the
#'   proportion of consumed energy that a mature individual allocates to
#'   reproduction for each species at size. If not supplied, a default is set as
#'   described in the section "Setting reproduction".
#' @param RDD The name of the function calculating the density-dependent 
#'   reproduction rate from the density-independent rate. Defaults to 
#'   "[BevertonHoltRDD()]".
#' @param ... Unused
#' 
#' @return For `setReproduction()`:
#'   The updated MizerParams object. Because of the way the R language
#'   works, `setReproduction()` does not make the changes to the params object
#'   that you pass to it but instead returns a new params object. So to affect
#'   the change you call the function in the form
#'   `params <- setReproduction(params, ...)`.
#' @export
#' @family functions for setting parameters
#' @examples
#' \donttest{
#' # Plot maturity and reproduction ogives for Cod in North Sea model
#' maturity <- getMaturityProportion(NS_params)["Cod", ]
#' repro_prop <- getReproductionProportion(NS_params)["Cod", ]
#' df <- data.frame(Size = w(NS_params), 
#'                  Reproduction = repro_prop, 
#'                  Maturity = maturity, 
#'                  Total = maturity * repro_prop)
#' dff <- melt(df, id.vars = "Size", 
#'             variable.name = "Type", 
#'             value.name = "Proportion")
#' library(ggplot2)
#' ggplot(dff) + geom_line(aes(x = Size, y = Proportion, colour = Type))
#' }
setReproduction <- function(params, maturity = NULL, repro_prop = NULL,
                            RDD = NULL, ...) {
    # check arguments ----
    assert_that(is(params, "MizerParams"))
    if (is.null(RDD)) RDD <- params@rates_funcs[["RDD"]]
    assert_that(is.string(RDD),
                exists(RDD),
                is.function(get(RDD)))
    species_params <- params@species_params
    
    # Check maximum sizes
    if (!("w_inf" %in% colnames(species_params))) {
        stop("The maximum sizes of the species must be specified in the w_inf ",
             "column of the species parameter data frame.")
    }
    missing <- is.na(species_params$w_inf)
    if (any(missing)) {
        stop("The following species are missing data for their maximum size w_inf: ",
             toString(species_params$species[missing]))
    }
    if (any(species_params$w_inf <= species_params$w_min)) {
        stop("Some of the asymptotic sizes are smaller than the egg sizes.")
    }
    # # Round maximum sizes to nearest grid point
    # for (i in seq_along(species_params$w_inf)) {
    #     idx <- which.min(abs(species_params$w_inf[i] - params@w))
    #     params@species_params$w_inf[i] < params@w[idx]
    # }
    
    # set maturity proportion ----
    if (!is.null(maturity)) {
        assert_that(is.array(maturity),
                    identical(dim(maturity), dim(params@psi)))
        if (!is.null(dimnames(maturity)) && 
            !all(dimnames(maturity)[[1]] == species_params$species)) {
            stop("You need to use the same ordering of species as in the ",
                 "params object: ", toString(species_params$species))
        }
    } else {
        # Check maturity sizes
        if (!("w_mat" %in% colnames(species_params))) {
            species_params$w_mat <- rep(NA, nrow(species_params))
        }
        missing <- is.na(species_params$w_mat)
        if (any(missing)) {
            message("Note: The following species were missing data for ",
                    "their maturity size w_mat: ",
                    toString(species_params$species[missing]),
                    ". These have been set to 1/4 w_inf.")
            species_params$w_mat[missing] <- species_params$w_inf[missing] / 4
        }
        assert_that(all(species_params$w_mat > species_params$w_min))
        
        # Set defaults for w_mat25
        species_params <- set_species_param_default(
            species_params, "w_mat25",       
            species_params$w_mat / (3 ^ (1 / 10)))
        # Check w_mat25
        assert_that(all(species_params$w_mat25 > species_params$w_min))
        assert_that(all(species_params$w_mat25 < species_params$w_mat))
        params@species_params$w_mat25 <- species_params$w_mat25
        
        maturity <- params@maturity  # To get the right dimensions
        maturity[] <- 
            unlist(
                tapply(params@w, seq_along(params@w),
                       function(wx, w_inf, w_mat, w_mat25) {
                           U <- log(3) / log(w_mat / w_mat25)
                           return((1 + (wx / w_mat)^-U)^-1)
                       },
                       w_inf = species_params$w_inf,
                       w_mat = species_params$w_mat,
                       w_mat25 = species_params$w_mat25
                )
            )
        
        # For reasons of efficiency we next set all very small values to 0 
        maturity[maturity < 1e-8] <- 0
        
        # If maturity is protected by a comment, keep the old value
        if (!is.null(comment(params@maturity))) {
            if (!isTRUE(all.equal(params@maturity, maturity,
                                  check.attributes = FALSE))) {
                message("The maturity ogive has been commented and therefore will ",
                        "not be recalculated from the species parameters.")
            }
            maturity <- params@maturity
        }
    }
    assert_that(all(maturity >= 0 & maturity <= 1))
    
    # Need to update psi because it contains maturity as a factor
    if (!isTRUE(all.equal(params@maturity, maturity,
                          check.attributes = FALSE))) {
        params@psi[] <- params@psi / params@maturity * maturity
        params@psi[is.nan(params@psi)] <- 0
    }
    
    params@maturity[] <- maturity
    comment(params@maturity) <- comment(maturity)
    
    # set reproduction proportion ----
    if (!is.null(repro_prop)) {
        assert_that(is.array(repro_prop),
                    identical(dim(repro_prop), dim(params@psi)))
        if (!is.null(dimnames(repro_prop)) && 
            !all(dimnames(repro_prop)[[1]] == species_params$species)) {
            stop("You need to use the same ordering of species as in the ",
                 "params object: ", toString(species_params$species))
        }
    } else {
        # Set defaults for m
        params <- set_species_param_default(params, "m", 1)
        
        repro_prop <- array(
            unlist(
                tapply(params@w, seq_along(params@w),
                       function(wx, w_inf, mn) (wx / w_inf)^(mn),
                       w_inf = params@species_params$w_inf,
                       mn = params@species_params$m - params@species_params$n
                )
            ), dim = c(nrow(species_params), length(params@w)))
    }
    
    psi <- params@maturity * repro_prop
    # psi should never be larger than 1
    psi[params@psi > 1] <- 1
    # Set psi for all w > w_inf to 1
    psi[outer(species_params$w_inf, params@w, "<")] <- 1
    assert_that(all(psi >= 0 & psi <= 1))
    
    if (!is.null(comment(params@psi))) {
        if (!isTRUE(all.equal(params@psi, psi,
                              check.attributes = FALSE))) {
            message("The reproductive proportion has been commented and therefore ",
                    "will not be recalculated from the species parameters.")
        }
    } else {
        params@psi[] <- psi
        comment(params@psi) <- comment(repro_prop)
    }
    
    # If no erepro (reproductive efficiency), then set to 1
    params <- set_species_param_default(params, "erepro", 1)
    assert_that(all(params@species_params$erepro > 0))
    
    # RDD function is currently called only with three arguments
    if (!all(names(formals(RDD)) %in%  c("rdi", "species_params", "t", "..."))) {
        stop("Arguments of RDD function can only contain 'rdi', 'species_params' and `t`.")
    }
    if (!all(c("rdi", "...") %in% names(formals(RDD)))) {
        stop("The RDD function needs to have at least arguments `rdi` and `...`.")
    }
    params@rates_funcs$RDD <- RDD
    if (identical(params@rates_funcs$RDD, "BevertonHoltRDD")) {
        
        # for legacy reasons (R_max used to be called r_max):
        if ("r_max" %in% names(params@species_params)) {
            params@species_params$R_max <- params@species_params$r_max
            params@species_params$r_max <- NULL
            message("The 'r_max' column has been renamed to 'R_max'.")
        }
        
        params <- set_species_param_default(params, "R_max", Inf)
    }
    
    return(params)
}

#' @rdname setReproduction
#' @return For `getMaturityProportion()`: 
#'   An array (species x size) that holds the proportion
#'   of individuals of each species at size that are mature.
#' @export
getMaturityProportion <- function(params) {
    assert_that(is(params, "MizerParams"))
    params@maturity
}

#' @rdname setReproduction
#' @return For `getReproductionProportion`:
#'   An array (species x size) that holds the
#'   proportion of consumed energy that a mature individual allocates to
#'   reproduction for each species at size. For sizes where the maturity
#'   proportion is zero, also the reproduction proportion is returned as zero.
#' @export
getReproductionProportion <- function(params) {
    assert_that(is(params, "MizerParams"))
    repro_prop <- params@psi / params@maturity
    repro_prop[is.nan(repro_prop)] <- 0
    comment(repro_prop) <- comment(params@psi)
    repro_prop
}
