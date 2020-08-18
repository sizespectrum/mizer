# Copyright 2012 Finlay Scott and Julia Blanchard.
# Copyright 2018 Gustav Delius and Richard Southwell.
# Development has received funding from the European Commission's Horizon 2020 
# Research and Innovation Programme under Grant Agreement No. 634495 
# for the project MINOUW (http://minouw-project.eu/).
# Distributed under the GPL 3 or later 
# Maintainer: Gustav Delius, University of York, <gustav.delius@york.ac.uk>


#' Get all rates needed to project standard mizer model
#' 
#' Calls other rate functions in sequence and collects the results in a list.
#' 
#' By default this function returns a list with the following components:
#'   \itemize{
#'     \item encounter from [mizerEncounter()]
#'     \item feeding_level from [mizerFeedingLevel()]
#'     \item e from [mizerEReproAndGrowth()]
#'     \item e_repro from [mizerERepro()]
#'     \item e_growth from [mizerEGrowth()]
#'     \item pred_rate from [mizerPredRate()]
#'     \item pred_mort from [mizerPredMort()]
#'     \item f_mort from [mizerFMort()]
#'     \item mort from [mizerMort()]
#'     \item rdi from [mizerRDI()]
#'     \item rdd from [BevertonHoltRDD()]
#'     \item resource_mort from [mizerResourceMort()]
#'   }
#' However you can replace any of these rate functions by your own rate
#' function if you wish, see [setRateFunction()] for details.
#' 
#' @param params A \linkS4class{MizerParams} object
#' @param n A matrix of species abundances (species x size).
#' @param n_pp A vector of the resource abundance by size
#' @param n_other A list of abundances for other dynamical components of the
#'   ecosystem
#' @param t The time for which to do the calculation (Not used by standard
#'   mizer rate functions but useful for extensions with time-dependent
#'   parameters.)
#' @param effort The effort for each fishing gear
#' @param rates_fns Named list of the functions to call to calculate the rates.
#'   Note that this list holds the functions themselves, not their names.
#' @param ... Unused
#' @export
#' @family mizer rate functions
mizerRates <- function(params, n, n_pp, n_other,
                       t = 0, effort, rates_fns, ...) {
    r <- list()
    
    ## Growth ----
    # Calculate rate E_{e,i}(w) of encountered food
    r$encounter <- rates_fns$Encounter(
        params, n = n, n_pp = n_pp, n_other = n_other, t = t, ...)
    # Calculate feeding level f_i(w)
    r$feeding_level <- rates_fns$FeedingLevel(
        params, n = n, n_pp = n_pp, n_other = n_other, 
        encounter = r$encounter, t = t, ...)
    # Calculate the energy available for reproduction and growth
    r$e <- rates_fns$EReproAndGrowth(
        params, n = n, n_pp = n_pp, n_other = n_other,
        encounter = r$encounter, feeding_level = r$feeding_level, t = t, ...)
    # Calculate the energy for reproduction
    r$e_repro <- rates_fns$ERepro(
        params, n = n, n_pp = n_pp, n_other = n_other, 
        e = r$e, t = t, ...)
    # Calculate the growth rate g_i(w)
    r$e_growth <- rates_fns$EGrowth(
        params, n = n, n_pp = n_pp, n_other = n_other,
        e_repro = r$e_repro, e = r$e, t = t, ...)
    
    ## Mortality ----
    # Calculate the predation rate
    r$pred_rate <- rates_fns$PredRate(
        params, n = n, n_pp = n_pp, n_other = n_other,
        feeding_level = r$feeding_level, t = t, ...)
    # Calculate predation mortality on fish \mu_{p,i}(w)
    r$pred_mort <- rates_fns$PredMort(
        params, n = n, n_pp = n_pp, n_other = n_other,
        pred_rate = r$pred_rate, t = t, ...)
    # Calculate fishing mortality
    r$f_mort <- rates_fns$FMort(
        params, n = n, n_pp = n_pp, n_other = n_other, 
        effort = effort, t = t, 
        e_growth = r$e_growth, pred_mort = r$pred_mort, ...)
    # Calculate total mortality \mu_i(w)
    r$mort <- rates_fns$Mort(
        params, n = n, n_pp = n_pp, n_other = n_other,
        f_mort = r$f_mort, pred_mort = r$pred_mort, t = t, ...)
    
    ## Reproduction ----
    # R_di
    r$rdi <- rates_fns$RDI(
        params, n = n, n_pp = n_pp, n_other = n_other,
        e_growth = r$e_growth,
        mort = r$mort,
        e_repro = r$e_repro, t = t, ...)
    # R_dd
    r$rdd <- rates_fns$RDD(
        rdi = r$rdi, species_params = params@species_params, ...)
    
    ## Resource ----
    # Calculate mortality on the resource spectrum
    r$resource_mort <- rates_fns$ResourceMort(
        params, n = n, n_pp = n_pp, n_other = n_other,
        pred_rate = r$pred_rate, t = t, ...)
    
    return(r)
}

#' Get encounter rate needed to project standard mizer model
#' 
#' Calculates the rate \eqn{E_i(w)} at which a predator of species \eqn{i} and
#' weight \eqn{w} encounters food (grams/year). You would not usually call this
#' function directly but instead use [getEncounter()], which then calls this
#' function unless an alternative function has been registered, see below.
#' 
#' @section Predation encounter:
#' The encounter rate \eqn{E_i(w)} at which a predator of species \eqn{i}
#' and weight \eqn{w} encounters food has contributions from the encounter of
#' fish prey and of resource. This is determined by summing over all prey
#' species and the resource spectrum and then integrating over all prey sizes
#' \eqn{w_p}, weighted by predation kernel \eqn{\phi(w,w_p)}:
#' \deqn{
#' E_i(w) = \gamma_i(w) \int 
#' \left( \theta_{ip} N_R(w_p) + \sum_{j} \theta_{ij} N_j(w_p) \right) 
#' \phi_i(w,w_p) w_p \, dw_p.
#' }{\gamma_i(w) \int 
#' ( \theta_{ip} N_R(w_p) + \sum_{j} \theta_{ij} N_j(w_p) ) 
#' \phi_i(w,w_p) w_p dw_p.}
#' Here \eqn{N_j(w)} is the abundance density of species \eqn{j} and
#' \eqn{N_R(w)} is the abundance density of resource.
#' The overall prefactor \eqn{\gamma_i(w)} determines the predation power of the
#' predator. It could be interpreted as a search volume and is set with the
#' [setSearchVolume()] function. The predation kernel
#' \eqn{\phi(w,w_p)} is set with the [setPredKernel()] function. The
#' species interaction matrix \eqn{\theta_{ij}} is set with [setInteraction()]
#' and the resource interaction vector \eqn{\theta_{ip}} is taken from the
#' `interaction_resource` column in `params@species_params`.
#' 
#' @section Details:
#' The encounter rate is multiplied by \eqn{1-f_0} to obtain the consumption rate,
#' where \eqn{f_0} is the feeding level calculated with [getFeedingLevel()].
#' This is used by the [project()] function for performing simulations.
#' 
#' The function returns values also for sizes outside the size-range of the
#' species. These values should not be used, as they are meaningless.
#' 
#' If your model contains additional components that you added with 
#' [setComponent()] and for which you specified an `encounter_fun` function then
#' the encounters of these components will be included in the returned value.
#' 
#' @section Your own encounter function:
#' By default [getEncounter()] calls [mizerEncounter()]. However you can
#' replace this with your own alternative encounter function. If 
#' your function is called `"myEncounter"` then you register it in a MizerParams
#' object `params` with
#' ```
#' params <- setRateFunction(params, "Encounter", "myEncounter")
#' ```
#' Your function will then be called instead of [mizerEncounter()], with the
#' same arguments.
#' 
#' @inheritParams mizerRates
#' @param ... Unused
#'   
#' @return A named two dimensional array (predator species x predator size) with
#'   the encounter rates.
#' @export
#' @family mizer rate functions
mizerEncounter <- function(params, n, n_pp, n_other, t, ...) {

    # idx_sp are the index values of params@w_full such that
    # params@w_full[idx_sp] = params@w
    idx_sp <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)
    
    # If the feeding kernel does not have a fixed predator/prey mass ratio
    # then the integral is not a convolution integral and we can not use fft.
    # In this case we use the code from mizer version 0.3
    if (length(params@ft_pred_kernel_e) == 1) {
        # n_eff_prey is the total prey abundance by size exposed to each
        # predator (prey not broken into species - here we are just working out
        # how much a predator eats - not which species are being eaten - that is
        # in the mortality calculation
        # \sum_j \theta_{ij} N_j(w_p) w_p dw_p
        n_eff_prey <- sweep(params@interaction %*% n, 2, 
                            params@w * params@dw, "*", check.margin = FALSE) 
        # pred_kernel is predator species x predator size x prey size
        # So multiply 3rd dimension of pred_kernel by the prey biomass density
        # Then sum over 3rd dimension to get consumption rate of each predator by 
        # predator size
        # This line is a bottle neck
        phi_prey_species <- rowSums(sweep(
            params@pred_kernel[, , idx_sp, drop = FALSE],
            c(1, 3), n_eff_prey, "*", check.margin = FALSE), dims = 2)
        # Eating the background
        # This line is a bottle neck
        phi_prey_background <- params@species_params$interaction_resource *
            rowSums(sweep(
            params@pred_kernel, 3, params@dw_full * params@w_full * n_pp,
            "*", check.margin = FALSE), dims = 2)
        encounter <- params@search_vol * (phi_prey_species + phi_prey_background)
    } else {
        prey <- outer(params@species_params$interaction_resource, n_pp)
        prey[, idx_sp] <- prey[, idx_sp] + params@interaction %*% n
        # The vector prey equals everything inside integral (3.4) except the feeding
        # kernel phi_i(w_p/w).
        prey <- sweep(prey, 2, params@w_full * params@dw_full, "*")
        # Eq (3.4) is then a convolution integral in terms of prey[w_p] and phi[w_p/w].
        # We approximate the integral by the trapezoidal method. Using the
        # convolution theorem we can evaluate the resulting sum via fast fourier
        # transform.
        # mvfft() does a Fourier transform of each column of its argument, but
        # we need the Fourier transforms of each row, so we need to apply mvfft()
        # to the transposed matrices and then transpose again at the end.
        avail_energy <- Re(base::t(mvfft(base::t(params@ft_pred_kernel_e) * 
                                             mvfft(base::t(prey)),
                                   inverse = TRUE))) / length(params@w_full)
        # Only keep the bit for fish sizes
        avail_energy <- avail_energy[, idx_sp, drop = FALSE]
        # Due to numerical errors we might get negative or very small entries that
        # should be 0
        avail_energy[avail_energy < 1e-18] <- 0
        
        encounter <- params@search_vol * avail_energy
    }
    
    # Add contributions from other components
    for (i in seq_along(params@other_encounter)) {
        encounter <- encounter + 
            do.call(params@other_encounter[[i]], 
                    list(params = params,
                         n = n, n_pp = n_pp, n_other = n_other,
                         component = names(params@other_encounter)[[i]], ...))
    }
    return(encounter)
}

#' Get feeding level needed to project standard mizer model
#' 
#' You would not usually call this function directly but instead use
#' [getFeedingLevel()], which then calls this function unless an alternative
#' function has been registered, see below.
#' 
#' @section Feeding level:
#' The feeding level \eqn{f_i(w)} is the
#' proportion of its maximum intake rate at which the predator is actually
#' taking in fish. It is calculated from the encounter rate \eqn{E_i} and the
#' maximum intake rate \eqn{h_i(w)} as
#' \deqn{f_i(w) = \frac{E_i(w)}{E_i(w)+h_i(w)}.}{E_i(w)/(E_i(w)+h_i(w)).}
#' The encounter rate \eqn{E_i} is passed as an argument or calculated with
#' [getEncounter()]. The maximum intake rate \eqn{h_i(w)} is
#' taken from the `params` object, and is set with 
#' [setMaxIntakeRate()].
#' As a consequence of the above expression for the feeding level,
#' \eqn{1-f_i(w)} is the proportion of the food available to it that the
#' predator actually consumes.
#' 
#' @section Your own feeding level function:
#' By default [getFeedingLevel()] calls [mizerFeedingLevel()]. However you can
#' replace this with your own alternative feeding level function. If 
#' your function is called `"myFeedingLevel"` then you register it in a MizerParams
#' object `params` with
#' ```
#' params <- setRateFunction(params, "FeedingLevel", "myFeedingLevel")
#' ```
#' Your function will then be called instead of [mizerFeedingLevel()], with the
#' same arguments.
#'
#' @seealso The feeding level is used in [mizerEReproAndGrowth()] and in
#' [mizerPredRate()].
#' 
#' @inheritParams mizerEncounter
#' @param encounter A two dimensional array (predator species x predator size) 
#'   with the encounter rate.
#'
#' @return A two dimensional array (predator species x predator size) with the
#'   feeding level.
#' 
#' @export
#' @family mizer rate functions
mizerFeedingLevel <- function(params, n, n_pp, n_other, t, encounter, ...) {
    return(encounter / (encounter + params@intake_max))
}

#' Get energy rate available for reproduction and growth  needed to project 
#' standard mizer model
#'
#' Calculates the energy rate
#' \eqn{E_{r.i}(w)} (grams/year) available to an
#' individual of species i and size w for reproduction and growth after
#' metabolism and movement have been accounted for.
#' You would not usually call this function directly but instead use
#' [getEReproAndGrowth()], which then calls this function unless an alternative
#' function has been registered, see below. 
#' 
#' @inheritParams mizerRates
#' @param encounter An array (species x size) with the encounter rate as
#'   calculated by [getEncounter()].
#' @param feeding_level An array (species x size) with the feeding level as
#'   calculated by [getFeedingLevel()].
#'
#' @return A two dimensional array (species x size) holding
#' \deqn{E_{r.i}(w) = \max(0, \alpha_i\, (1 - {\tt feeding\_level}_i(w))\, 
#'                            {\tt encounter}_i(w) - {\tt metab}_i(w)).}{
#'   E_{r.i}(w) = max(0, alpha_i * (1 - feeding_level_i(w)) * 
#'                       encounter_i(w) - metab_i(w)).}
#' Due to the form of the feeding level, calculated by
#' [getFeedingLevel()], this can also be expressed as
#' \deqn{E_{r.i}(w) = \max(0, \alpha_i\, {\tt feeding\_level}_i(w)\, 
#'                            h_i(w) - {\tt metab}_i(w))}{
#'   E_{r.i}(w) = max(0, alpha_i * feeding_level_i(w) * 
#'                       h_i(w) - metab_i(w))}
#' where \eqn{h_i} is the maximum intake rate, set with 
#' [setMaxIntakeRate()].
#' The assimilation rate \eqn{\alpha_i} is taken from the species parameter
#' data frame in `params`. The metabolic rate `metab` is taken from 
#' `params` and set with [setMetabolicRate()].
#' 
#' The return value can be negative, which means that the energy intake does not
#' cover the cost of metabolism and movement.
#' 
#' @section Your own energy rate function:
#' By default [getEReproAndGrowth()] calls [mizerEReproAndGrowth()]. However you
#' can replace this with your own alternative energy rate function. If 
#' your function is called `"myEReproAndGrowth"` then you register it in a
#' MizerParams object `params` with
#' ```
#' params <- setRateFunction(params, "EReproAndGrowth", "myEReproAndGrowth")
#' ```
#' Your function will then be called instead of [mizerEReproAndGrowth()], with
#' the same arguments.
#' 
#' @export
#' @family mizer rate functions
mizerEReproAndGrowth <- function(params, n, n_pp, n_other, t, encounter,
                                 feeding_level, ...) {
    
    sweep((1 - feeding_level) * encounter, 1,
          params@species_params$alpha, "*", check.margin = FALSE) - 
        params@metab
}

#' Get energy rate available for reproduction needed to project standard mizer 
#' model
#'
#' Calculates the energy rate (grams/year) available for reproduction after
#' growth and metabolism have been accounted for.
#' You would not usually call this
#' function directly but instead use [getERepro()], which then calls this
#' function unless an alternative function has been registered, see below.
#' 
#' @section Your own reproduction rate function:
#' By default [getERepro()] calls [mizerERepro()]. However you can
#' replace this with your own alternative reproduction rate function. If 
#' your function is called `"myERepro"` then you register it in a MizerParams
#' object `params` with
#' ```
#' params <- setRateFunction(params, "ERepro", "myERepro")
#' ```
#' Your function will then be called instead of [mizerERepro()], with the
#' same arguments.
#' 
#' @inheritParams mizerRates
#' @param e A two dimensional array (species x size) holding the energy available
#'   for reproduction and growth as calculated by [mizerEReproAndGrowth()].
#'
#' @return A two dimensional array (species x size) holding
#' \deqn{\psi_i(w)E_{r.i}(w)}
#' where \eqn{E_{r.i}(w)} is the rate at which energy becomes available for
#' growth and reproduction, calculated with [mizerEReproAndGrowth()],
#' and \eqn{\psi_i(w)} is the proportion of this energy that is used for
#' reproduction. This proportion is taken from the `params` object and is
#' set with [setReproduction()].
#' @export
#' @family mizer rate functions
mizerERepro <- function(params, n, n_pp, n_other, t, e, ...) {
    # Because getEReproAndGrowth can return negative values, 
    # we add an extra line here 
    e[e < 0] <- 0 # Do not allow negative growth
    
    params@psi * e
}

#' Get energy rate available for growth needed to project standard mizer model
#'
#' Calculates the energy rate \eqn{g_i(w)} (grams/year) available by species and
#' size for growth after metabolism, movement and reproduction have been
#' accounted for. Used by [project()] for performing simulations.
#' You would not usually call this
#' function directly but instead use [getEGrowth()], which then calls this
#' function unless an alternative function has been registered, see below.
#' 
#' @section Your own growth rate function:
#' By default [getEGrowth()] calls [mizerEGrowth()]. However you can
#' replace this with your own alternative growth rate function. If 
#' your function is called `"myEGrowth"` then you register it in a MizerParams
#' object `params` with
#' ```
#' params <- setRateFunction(params, "EGrowth", "myEGrowth")
#' ```
#' Your function will then be called instead of [mizerEGrowth()], with the
#' same arguments.
#' 
#' @inheritParams mizerRates
#' @param e The energy available for reproduction and growth as calculated by
#'   [getEReproAndGrowth()].
#' @param e_repro The energy available for reproduction as calculated by
#'   [getERepro()].
#'   
#' @return A two dimensional array (species x size) with the growth rates.
#' @export
#' @family mizer rate functions
mizerEGrowth <- function(params, n, n_pp, n_other, t, e_repro, e, ...) {
    # Because getEReproAndGrowth can return negative values, we add an 
    # extra line here 
    e[e < 0] <- 0 # Do not allow negative growth
    
    # energy for growth is intake - energy for reproduction
    e - e_repro
}


#' Get predation rate needed to project standard mizer model
#' 
#' Calculates the potential rate (in units 1/year) at which a prey individual of
#' a given size \eqn{w} is killed by predators from species \eqn{j}. In formulas
#' \deqn{{\tt pred\_rate}_j(w_p) = \int \phi_j(w,w_p) (1-f_j(w)) 
#'   \gamma_j(w) N_j(w) \, dw.}{pred_rate_j(w_p) = \int\phi_i(w,w_p) (1-f_i(w)) 
#'   \gamma_i(w) N_i(w) dw.}
#' This potential rate is used in the function [mizerPredMort()] to
#' calculate the realised predation mortality rate on the prey individual.
#' You would not usually call this
#' function directly but instead use [getPredRate()], which then calls this
#' function unless an alternative function has been registered, see below.
#' 
#' @section Your own predation rate function:
#' By default [getPredRate()] calls [mizerPredRate()]. However you can
#' replace this with your own alternative predation rate function. If 
#' your function is called `"myPredRate"` then you register it in a MizerParams
#' object `params` with
#' ```
#' params <- setRateFunction(params, "PredRate", "myPredRate")
#' ```
#' Your function will then be called instead of [mizerPredRate()], with
#' the same arguments.
#'
#' @inheritParams mizerEReproAndGrowth
#'   
#' @return A named two dimensional array (predator species x prey size) with the
#'   predation rate, where the prey size runs over fish community plus resource
#'   spectrum.
#' @export
#' @family mizer rate functions
mizerPredRate <- function(params, n, n_pp, n_other, t, feeding_level, ...) {
    no_sp <- dim(params@interaction)[1]
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)
    
    # If the feeding kernel does not have a fixed predator/prey mass ratio
    # then the integral is not a convolution integral and we can not use fft.
    # In this case we use the code from mizer version 0.3
    if (length(params@ft_pred_kernel_p) == 1) {
        n_total_in_size_bins <- sweep(n, 2, params@dw, '*', check.margin = FALSE)
        # The next line is a bottle neck
        pred_rate <- sweep(params@pred_kernel, c(1,2),
                           (1 - feeding_level) * params@search_vol * 
                               n_total_in_size_bins,
                           "*", check.margin = FALSE)
        # integrate over all predator sizes
        pred_rate <- colSums(aperm(pred_rate, c(2, 1, 3)), dims = 1)
        return(pred_rate)
    }

    # Get indices of w_full that give w
    idx_sp <- (no_w_full - no_w + 1):no_w_full
    # We express the result as a a convolution  involving
    # two objects: Q[i,] and ft_pred_kernel_p[i,].
    # Here Q[i,] is all the integrand of (3.12) except the feeding kernel
    # and theta
    Q <- matrix(0, nrow = no_sp, ncol = no_w_full)
    # We fill the end of each row of Q with the proper values
    Q[, idx_sp] <- sweep( (1 - feeding_level) * params@search_vol * n, 2,
                         params@dw, "*")

    # We do our spectral integration in parallel over the different species
    pred_rate <- Re(t(mvfft(t(params@ft_pred_kernel_p) *
                                 mvfft(t(Q)), inverse = TRUE))) / no_w_full
    # Due to numerical errors we might get negative or very small entries that
    # should be 0
    pred_rate[pred_rate < 1e-18] <- 0
    
    return(pred_rate * params@ft_mask)
}


#' Get total predation mortality rate needed to project standard mizer model
#'
#' Calculates the total predation mortality rate \eqn{\mu_{p,i}(w_p)} (in units
#' of 1/year) on each prey species by prey size:
#' \deqn{\mu_{p.i}(w_p) = \sum_j {\tt pred\_rate}_j(w_p)\, \theta_{ji}.}{
#'   \mu_{p.i}(w_p) = \sum_j pred_rate_j(w_p) \theta_{ji}.}
#' You would not usually call this
#' function directly but instead use [getPredMort()], which then calls this
#' function unless an alternative function has been registered, see below.
#' 
#' @section Your own predation mortality function:
#' By default [getPredMort()] calls [mizerPredMort()]. However you can
#' replace this with your own alternative predation mortality function. If 
#' your function is called `"myPredMort"` then you register it in a MizerParams
#' object `params` with
#' ```
#' params <- setRateFunction(params, "PredMort", "myPredMort")
#' ```
#' Your function will then be called instead of [mizerPredMort()], with the
#' same arguments.
#' 
#' @inheritParams mizerRates
#' @param pred_rate A two dimensional array (predator species x predator size)
#'   with the feeding level.
#'
#' @return A two dimensional array (prey species x prey size) with the predation
#'   mortality
#' @family mizer rate functions
#' @export
mizerPredMort <- function(params, n, n_pp, n_other, t, pred_rate, ...) {
    idx_sp <- (length(params@w_full) - 
                   length(params@w) + 1):length(params@w_full)
    return((base::t(params@interaction) %*% pred_rate)[, idx_sp, drop = FALSE])
}

#' Get the fishing mortality by time, gear, species and size needed to project 
#' standard mizer model
#'
#' Calculates the fishing mortality rate \eqn{F_{g,i,w}} by gear, species and
#' size at each time step in the `effort` argument (in units 1/year).
#' This is a helper function for [mizerFMort()].
#' 
#' @inheritParams mizerRates
#' @param effort A vector with the effort for each fishing gear.
#'   
#' @return An three dimensional array (gear x species x size) with the
#'    fishing mortality 
#' @note Here: fishing mortality = catchability x selectivity x effort.
#' 
#' @export
#' @family mizer rate functions
mizerFMortGear <- function(params, effort) {
    # Streamlined for speed increase - note use of recycling
    out <- params@selectivity
    out[] <- effort * c(params@catchability) * c(params@selectivity)
    return(out)
}


#' Get the total fishing mortality rate from all fishing gears by time, species
#' and size needed to project standard mizer model
#' 
#' Calculates the total fishing mortality  (in units 1/year) from all gears by
#' species and size at each time step in the `effort` argument.
#' The total fishing mortality is just the sum of the fishing mortalities
#' imposed by each gear, \eqn{\mu_{f.i}(w)=\sum_g F_{g,i,w}}.
#' You would not usually call this
#' function directly but instead use [getFMort()], which then calls this
#' function unless an alternative function has been registered, see below.
#' 
#' @section Your own fishing mortality function:
#' By default [getFMort()] calls [mizerFMort()]. However you can
#' replace this with your own alternative fishing mortality function. If 
#' your function is called `"myFMort"` then you register it in a MizerParams
#' object `params` with
#' ```
#' params <- setRateFunction(params, "FMort", "myFMort")
#' ```
#' Your function will then be called instead of [mizerFMort()], with the
#' same arguments.
#' 
#' @inheritParams mizerRates
#' @param effort A vector with the effort for each fishing gear.
#' @param e_growth An array (species x size) with the energy available for
#'   growth as calculated by [getEGrowth()]. Unused.
#' @param pred_mort A two dimensional array (species x size) with the predation
#'   mortality as calculated by [getPredMort()]. Unused.
#'
#' @return An array (species x size) with the fishing mortality.
#' @note Here: fishing mortality = catchability x selectivity x effort.
#' @export
#' @family mizer rate functions
mizerFMort <- function(params, n, n_pp, n_other, t, effort,
                       e_growth, pred_mort, ...) {
    colSums(mizerFMortGear(params, effort))
}

#' Get total mortality rate needed to project standard mizer model
#'
#' Calculates the total mortality rate \eqn{\mu_i(w)}  (in units 1/year) on each
#' species by size from predation mortality, background mortality and fishing
#' mortality.
#' You would not usually call this
#' function directly but instead use [getMort()], which then calls this
#' function unless an alternative function has been registered, see below.
#' 
#' If your model contains additional components that you added with 
#' [setComponent()] and for which you specified a `mort_fun` function then
#' the mortality inflicted by these components will be included in the returned
#' value.
#' 
#' @section Your own mortality function:
#' By default [getMort()] calls [mizerMort()]. However you can
#' replace this with your own alternative mortality function. If 
#' your function is called `"myMort"` then you register it in a MizerParams
#' object `params` with
#' ```
#' params <- setRateFunction(params, "Mort", "myMort")
#' ```
#' Your function will then be called instead of [mizerMort()], with the
#' same arguments.
#'
#' @inheritParams mizerRates
#' @param f_mort A two dimensional array (species x size) with the fishing
#'   mortality
#' @param pred_mort A two dimensional array (species x size) with the predation
#'   mortality
#' @param ... Unused
#'
#' @return A named two dimensional array (species x size) with the total
#'   mortality rates.
#' @export
#' @family mizer rate functions
mizerMort <- function(params, n, n_pp, n_other, t, f_mort, pred_mort, ...){
    mort <- pred_mort + params@mu_b + f_mort
    # Add contributions from other components
    for (i in seq_along(params@other_mort)) {
        mort <- mort + 
            do.call(params@other_mort[[i]], 
                    list(params = params,
                         n = n, n_pp = n_pp, n_other = n_other, t = t,
                         component = names(params@other_mort)[[i]], ...))
    }
    return(mort)
}


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
#' same arguments.
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
#' reproduction rates. Two other examples provided are [RickerRDD()]
#' and [SheperdRDD()]. For more explanation see
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
#' Takes the density-independent rates \eqn{R_{di}}{R_di} of egg production and returns
#' reduced, density-dependent rates \eqn{R_{dd}}{R_dd} given as
#' \deqn{R_{dd} = R_{di} \exp{- b R_{di}}}{R_dd = R_di \exp{- b R_di}}
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


#' Get predation mortality rate for resource needed to project standard mizer 
#' model
#' 
#' Calculates the predation mortality rate \eqn{\mu_p(w)} on the resource
#' spectrum by resource size (in units 1/year).
#' You would not usually call this
#' function directly but instead use [getResourceMort()], which then calls this
#' function unless an alternative function has been registered, see below.
#' 
#' @section Your own resource mortality function:
#' By default [getResourceMort()] calls [mizerResourceMort()]. However you can
#' replace this with your own alternative resource mortality function. If 
#' your function is called `"myResourceMort"` then you register it in a MizerParams
#' object `params` with
#' ```
#' params <- setRateFunction(params, "ResourceMort", "myResourceMort")
#' ```
#' Your function will then be called instead of [mizerResourceMort()], with the
#' same arguments.
#' 
#' @inheritParams mizerRates
#' @param pred_rate A two dimensional array (predator species x prey size) with
#'   the predation rate, where the prey size runs over fish community plus
#'   resource spectrum.
#'
#' @return A vector of mortality rate by resource size.
#' @family mizer rate functions
#' @export
mizerResourceMort <- 
    function(params, n, n_pp, n_other, t, pred_rate, ...) {
        
        return(as.vector(params@species_params$interaction_resource %*% pred_rate))
    }