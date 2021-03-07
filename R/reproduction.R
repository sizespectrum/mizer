#' Get density-independent rate of reproduction needed to project standard
#' mizer model
#'
#' Calculates the density-independent rate of total egg production 
#' \eqn{R_{di}}{R_di} (units 1/year) before density dependence, by species. 
#' You would not usually call this
#' function directly but instead use [getRDI()], which then calls this
#' function unless an alternative function has been registered, see below.
#' 
#' This rate is obtained by taking the per capita rate \eqn{E_r(w)\psi(w)} at
#' which energy is invested in reproduction, as calculated by [getERepro()],
#' multiplying it by the number of individuals\eqn{N(w)} and integrating over all sizes
#' \eqn{w} and then multiplying by the reproductive efficiency \eqn{\epsilon}
#' and dividing by the egg size `w_min`, and by a factor of two to account for
#' the two sexes:
#' \deqn{R_{di} = \frac{\epsilon}{2 w_{min}} \int N(w)  E_r(w) \psi(w) \, dw}{R_di = (\epsilon/(2 w_min)) \int N(w)  E_r(w) \psi(w) dw}
#' 
#' Used by [getRDD()] to calculate the actual, density dependent rate.
#' See [setReproduction()] for more details.
#' 
#' 
#' @section Your own reproduction function:
#' By default [getRDI()] calls [mizerRDI()]. However you can
#' replace this with your own alternative reproduction function. If 
#' your function is called `"myRDI"` then you register it in a MizerParams
#' object `params` with
#' ```
#' params <- setRateFunction(params, "RDI", "myRDI")
#' ```
#' Your function will then be called instead of [mizerRDI()], with the
#' same arguments. For an example of an alternative reproduction function
#' see `constantEggRDI()`.
#'
#' @inheritParams mizerRates
#' @param e_repro An array (species x size) with the energy available for
#'   reproduction as calculated by [getERepro()].
#' @param e_growth An array (species x size) with the energy available for
#'   growth as calculated by [getEGrowth()]. Unused.
#' @param mort An array (species x size) with the mortality rate as calculated
#'   by [getMort()]. Unused.
#'
#' @return A numeric vector with the rate of egg production for each species.
#' @export
#' @family mizer rate functions
mizerRDI <- function(params, n, n_pp, n_other, t,
                     e_growth, mort, e_repro, ...) {
    # Calculate total energy from per capita energy
    e_repro_pop <- drop((e_repro * n) %*% params@dw)
    # Assume sex_ratio = 0.5
    rdi <- 0.5 * (e_repro_pop * params@species_params$erepro) /
        params@w[params@w_min_idx]
    return(rdi)
}

#' Choose egg production to keep egg density constant
#'
#' `r lifecycle::badge("experimental")`
#' The new egg production is set to compensate for the loss of individuals from
#' the smallest size class through growth and mortality. The result should not
#' be modified by density dependence, so this should be used together with
#' the `noRDD()` function, see example.
#'
#' @param params A MizerParams object
#' @param n A matrix of species abundances (species x size).
#' @param e_growth A two dimensional array (species x size) holding the energy
#'   available for growth as calculated by [mizerEGrowth()].
#' @param mort A two dimensional array (species x size) holding the mortality
#'   rate as calculated by [mizerMort()].
#' @param ... Unused
#'
#' @export
#' @family functions calculating density-dependent reproduction rate
#' @examples
#' \dontrun{
#' # choose an example params object
#' params <- NS_params
#' # We set the reproduction rate functions
#' params <- setRateFunction(params, "RDI", "constantEggRDI")
#' params <- setRateFunction(params, "RDD", "noRDD")
#' # Now the egg density should stay fixed no matter how we fish
#' sim <- project(params, effort = 10, progress_bar = FALSE)
#' # To check that indeed the egg densities have not changed, we first construct
#' # the indices for addressing the egg densities
#' no_sp <- nrow(params@@species_params)
#' idx <- (params@w_min_idx - 1) * no_sp + (1:no_sp)
#' # Now we can check equality between egg densities at the start and the end
#' all.equal(finalN(sim)[idx], initialN(params)[idx])
#' }
constantEggRDI <- function(params, n, e_growth, mort, ...) {
    no_sp <- nrow(params@species_params) # number of species
    # Hacky shortcut to access the correct element of a 2D array using 1D notation
    idx <- (params@w_min_idx - 1) * no_sp + (1:no_sp)
    rdi <- n[idx] * (e_growth[idx] + mort[idx] * params@dw[params@w_min_idx])
    rdi
}


#' Beverton Holt function to calculate density-dependent reproduction rate
#'
#' Takes the density-independent rates \eqn{R_{di}}{R_di} of egg production (as
#' calculated by [getRDI()]) and returns
#' reduced, density-dependent reproduction rates \eqn{R_{dd}}{R_dd} given as 
#' \deqn{R_{dd} = R_{di}
#' \frac{R_{max}}{R_{di} + R_{max}}}{R_dd = R_di R_max/(R_di + R_max)} where
#' \eqn{R_{max}}{R_max} are the maximum possible reproduction rates that must be
#' specified in a column in the species parameter dataframe.
#' (All quantities in the above equation are species-specific but we dropped
#' the species index for simplicity.)
#'
#' This is only one example of a density-dependence. You can write your own
#' function based on this example, returning different density-dependent
#' reproduction rates. Three other examples provided are [RickerRDD()],
#' [SheperdRDD()], [noRDD()] and [constantRDD()]. For more explanation see
#' [setReproduction()].
#'
#' @param rdi Vector of density-independent reproduction rates 
#'   \eqn{R_{di}}{R_di} for all species.
#' @param species_params A species parameter dataframe. Must contain a column
#'   R_max holding the maximum reproduction rate \eqn{R_{max}}{R_max} for each species.
#' @param ... Unused
#'
#' @return Vector of density-dependent reproduction rates.
#' @export
#' @family functions calculating density-dependent reproduction rate
BevertonHoltRDD <- function(rdi, species_params, ...) {
    if (!("R_max" %in% names(species_params))) {
        stop("The R_max column is missing in species_params.")
    }
    return(rdi / (1 + rdi/species_params$R_max))
}

#' Ricker function to calculate density-dependent reproduction rate
#' 
#' `r lifecycle::badge("experimental")`
#' Takes the density-independent rates \eqn{R_{di}}{R_di} of egg production and 
#' returns reduced, density-dependent rates \eqn{R_{dd}}{R_dd} given as
#' \deqn{R_{dd} = R_{di} \exp(- b R_{di})}{R_dd = R_di exp(- b R_di)}
#' 
#' @param rdi Vector of density-independent reproduction rates 
#'   \eqn{R_{di}}{R_di} for all species.
#' @param species_params A species parameter dataframe. Must contain a column
#'   `ricker_b` holding the coefficient b.
#' @param ... Unused
#' 
#' @return Vector of density-dependent reproduction rates.
#' @export
#' @family functions calculating density-dependent reproduction rate
RickerRDD <- function(rdi, species_params, ...) {
    if (!("ricker_b" %in% names(species_params))) {
        stop("The ricker_b column is missing in species_params")
    }
    return(rdi * exp(-species_params$ricker_b * rdi))
}

#' Sheperd function to calculate density-dependent reproduction rate
#' 
#' `r lifecycle::badge("experimental")`
#' Takes the density-independent rates \eqn{R_{di}}{R_di} of egg production and returns
#' reduced, density-dependent rates \eqn{R_{dd}}{R_dd} given as
#' \deqn{R_{dd} = \frac{R_{di}}{1+(b\ R_{di})^c}}{R_dd = R_di / (1 + (b R_di)^c)}
#' 
#' @param rdi Vector of density-independent reproduction rates 
#'   \eqn{R_{di}}{R_di} for all species.
#' @param species_params A species parameter dataframe. Must contain columns
#'   `sheperd_b` and `sheperd_c` with the parameters b and c.
#' @param ... Unused
#' 
#' @return Vector of density-dependent reproduction rates.
#' @export
#' @family functions calculating density-dependent reproduction rate
SheperdRDD <- function(rdi, species_params, ...) {
    if (!all(c("sheperd_b", "sheperd_c") %in% names(species_params))) {
        stop("The species_params dataframe must contain columns sheperd_b and sheperd_c.")
    }
    return(rdi / (1 + (species_params$sheperd_b * rdi)^species_params$sheperd_c))
}

#' Give density-independent reproduction rate
#' 
#' Simply returns its `rdi` argument.
#' 
#' @param rdi Vector of density-independent reproduction rates 
#'   \eqn{R_{di}}{R_di} for all species.
#' @param ... Not used.
#' 
#' @return Vector of density-dependent reproduction rates.
#' @export
#' @family functions calculating density-dependent reproduction rate
noRDD <- function(rdi, ...) {
    return(rdi)
}

#' Give constant reproduction rate
#'
#' `r lifecycle::badge("experimental")`
#' Simply returns the value from `species_params$constant_reproduction`.
#'
#' @param rdi Vector of density-independent reproduction rates 
#'   \eqn{R_{di}}{R_di} for all species.
#' @param species_params A species parameter dataframe. Must contain a column
#'   `constant_reproduction`.
#' @param ... Unused
#'
#' @return Vector `species_params$constant_reproduction`
#' @export
#' @family functions calculating density-dependent reproduction rate
constantRDD <- function(rdi, species_params, ...){
    return(species_params$constant_reproduction)
}
