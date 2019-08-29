#' Functions for calculating rates used for projecting
#'
#' The functions defined in the file project_methods calculate the various
#' quantities needed to project the size-spectra forward in time, using the
#' model described in section 3 of the mizer vignette.
#'
#' @section List of functions:
#' In this list we relate the functions in this file to the quantities named in
#' the mizer vignette.
#' \tabular{llll}{
#'   Function \tab Expression \tab Description \tab Section in vignette\cr
#'   \code{\link{getRates}} \tab \tab All of the below \cr
#'   \code{\link{getEncounter}} \tab \eqn{E_{e.i}(w)} \tab Encounter rate \tab 3.2 \cr
#'   \code{\link{getFeedingLevel}} \tab \eqn{f_i(w)} \tab Feeding level \tab 3.3 \cr
#'   \code{\link{getPredRate}} \tab \eqn{\phi_i(w_p/w) (1-f_i(w)) \gamma_i w^q N_i(w) dw} \tab Predation \tab 3.7 \cr
#'   \code{\link{getPredMort}} \tab \eqn{\mu_{p.i}(w)} \tab Predation mortality \tab 3.7 \cr
#'   \code{\link{getPlanktonMort}} \tab \eqn{\mu_{p}(w)} \tab Mortality on plankton \tab 3.8 \cr
#'   \code{\link{getFMortGear}} \tab \eqn{F_{g,i}(w)} \tab Fishing mortality by gear \tab 8.3 \cr
#'   \code{\link{getFMort}} \tab \eqn{\mu_{f.i}(w)} \tab Total fishing mortality \tab 8.3 \cr
#'   \code{\link{getMort}} \tab \eqn{\mu_{i}(w)} \tab Total mortality \tab 3.7 \cr
#'   \code{\link{getEReproAndGrowth}} \tab \eqn{E_{r.i}(w)} \tab Energy put into growth and reproduction \tab 3.4 \cr
#'   \code{\link{getERepro}} \tab \eqn{\psi_i(w)E_{r.i}(w)} \tab Energy put reproduction\tab 3.5 \cr
#'   \code{\link{getEGrowth}} \tab \eqn{g_i(w)} \tab Energy put growth \tab 3.4 \cr
#'   \code{\link{getRDI}} \tab \eqn{R_{p.i}} \tab Egg production \tab 3.5 \cr
#'   \code{\link{getRDD}} \tab \eqn{R_i} \tab Recruitment \tab 3.6 \cr
#' }
#'
#' @name project_methods
NULL

# Copyright 2012 Finlay Scott and Julia Blanchard.
# Copyright 2018 Gustav Delius and Richard Southwell.
# Development has received funding from the European Commission's Horizon 2020 
# Research and Innovation Programme under Grant Agreement No. 634495 
# for the project MINOUW (http://minouw-project.eu/).
# Distributed under the GPL 3 or later 
# Maintainer: Gustav Delius, University of York, <gustav.delius@york.ac.uk>


#' Get all rates
#' 
#' Calls all the other rate functions in sequence and collects the results in a
#' list.
#' 
#' @param params A \linkS4class{MizerParams} object
#' @param n A matrix of species abundances (species x size)
#' @param n_pp A vector of the plankton abundance by size
#' @param B A vector of biomasses of unstructured resource components
#' @param effort The effort for each fishing gear. Default 1.
#' @param sex_ratio Proportion of the population that is female. Default value
#'   is 0.5.
#' 
#' @return A list with the following components:
#'   \itemize{
#'     \item encounter from \code{\link{getEncounter}}
#'     \item feeding_level from \code{\link{getFeedingLevel}}
#'     \item pred_rate from \code{\link{getPredRate}}
#'     \item pred_mort from \code{\link{getPredMort}}
#'     \item mort from \code{\link{getMort}}
#'     \item plankton_mort from \code{\link{getPlanktonMort}}
#'     \item e from \code{\link{getEReproAndGrowth}}
#'     \item e_repro from \code{\link{getERepro}}
#'     \item e_growth from \code{\link{getEGrowth}}
#'     \item rdi from \code{\link{getRDI}}
#'     \item rdd from \code{\link{getRDD}}
#'   }
#' @export
#' @family rate functions
getRates <- function(params, n = params@initial_n, 
                     n_pp = params@initial_n_pp,
                     B = params@initial_B,
                     effort = 1, sex_ratio = 0.5) {
    r <- list()
    # Calculate rate E_{e,i}(w) of encountered food
    r$encounter <- getEncounter(params, n = n, n_pp = n_pp, B = B)
    # Calculate feeding level f_i(w)
    r$feeding_level <- getFeedingLevel(params, n = n, n_pp = n_pp, B = B,
                                       encounter = r$encounter)
    # Calculate the predation rate
    r$pred_rate <- getPredRate(params, n = n, n_pp = n_pp, B = B,
                               feeding_level = r$feeding_level)
    # Calculate predation mortality on fish \mu_{p,i}(w)
    r$pred_mort <- getPredMort(params, pred_rate = r$pred_rate)
    # Calculate total mortality \mu_i(w)
    r$mort <- getMort(params, n = n, n_pp = n_pp, B = B, 
                      effort = effort, m2 = r$pred_mort)
    # Calculate mortality on the plankton spectrum
    r$plankton_mort <- getPlanktonMort(params, n = n, n_pp = n_pp, B = B,
                                       pred_rate = r$pred_rate)
    # Calculate the energy available for reproduction and growth
    r$e <- getEReproAndGrowth(params, n = n, n_pp = n_pp, B = B,
                              encounter = r$encounter,
                              feeding_level = r$feeding_level)
    # Calculate the energy for reproduction
    r$e_repro <- getERepro(params, n = n, n_pp = n_pp, B = B, e = r$e)
    # Calculate the growth rate g_i(w)
    r$e_growth <- getEGrowth(params, n = n, n_pp = n_pp, B = B, 
                             e_repro = r$e_repro, e = r$e)
    # R_{p,i}
    r$rdi <- getRDI(params, n = n, n_pp = n_pp, B = B, 
                    e_repro = r$e_repro, sex_ratio = sex_ratio)
    # R_i
    r$rdd <- params@srr(rdi = r$rdi, species_params = params@species_params)
    
    return(r)
}

#' Get encounter rate
#' 
#' Calculates the rate \eqn{E_i(w)} at which a predator of species \eqn{i} and
#' weight \eqn{w} encounters food (grams/year).
#' 
#' @section Predation encounter:
#' The encounter rate has contributions from the encounter of fish prey, of
#' plankton, and of other resources. The contribution from fish and plankton
#' is determined by summing over all prey species and the resource spectrum and
#' then integrating over all prey sizes \eqn{w_p}, weighted by predation kernel 
#' \eqn{\phi(w,w_p)}:
#' \deqn{
#' E_{e.i}(w) = \gamma_i(w) \int 
#' \left( \theta_{ip} N_R(w_p) + \sum_{j} \theta_{ij} N_j(w_p) \right) 
#' \phi_i(w,w_p) w_p \, dw_p.
#' }{\gamma_i(w) \int 
#' ( \theta_{ip} N_R(w_p) + \sum_{j} \theta_{ij} N_j(w_p) ) 
#' \phi_i(w,w_p) w_p dw_p.}
#' Here \eqn{N_j(w)} is the abundance density of species \eqn{j} and
#' \eqn{N_R(w)} is the abundance density of plankton.
#' The overall prefactor \eqn{\gamma_i(w)} determines the predation power of the
#' predator. It could be interpreted as a search volume and is changed with the
#' \code{\link{setSearchVolume}} function. The predation kernel
#' \eqn{\phi(w,w_p)} is changed with the \code{\link{setPredKernel}} function. The
#' species interaction matrix \eqn{\theta_{ij}} and the plankton interaction
#' vector \eqn{\theta_{ip}} are changed with \code{\link{setInteraction}}.
#' 
#' @section Resource encounter:
#' In addition to the contribution from predation on fish prey and plankton,
#' the food encounter rate may have a contribution from unstructured resource
#' components. This takes the form
#' \deqn{E_{u.i} = \sum_d \rho_{id}(w) B_d.}
#' where \eqn{B_d} is the biomass of the d-th unstructured resource component
#' and \eqn{\rho_{id}(w)} is a parameter that therefore determines the rate at
#' which a predator of species \eqn{i} and size \eqn{w} encounters biomass from
#' the d-th unstructured resource component. This is changed with
#' \code{\link{setResourceEncounter}}.
#' 
#' @section Details:
#' The total encounter rate is the sum of the contribution from fish and
#' plankton and the contribution from unstructured resources, if any:
#' \deqn{E_i(w)=E_{e.i}(w)+E_{u.i}(w).}
#'
#' The encounter rate is multiplied by \eqn{1-f_0} to obtain the consumption rate,
#' where \eqn{f_0} is the feeding level calculated with \code{\link{getFeedingLevel}}.
#' This is used by the \code{\link{project}} function for performing simulations.
#' 
#' The function returns values also for sizes outside the size-range of the
#' species. These values should not be used, as they are meaningless.
#' 
#' @param params An \linkS4class{MizerParams} object
#' @param n A matrix of species abundances (species x size)
#' @param n_pp A vector of the plankton abundance by size
#' @param B A vector of biomasses of unstructured resource components
#'   
#' @return A two dimensional array (predator species x predator size)
#' @export
#' @family rate functions
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- set_multispecies_model(NS_species_params_gears, inter)
#' # Run simulation with constant fishing effort for all gears for 20 years
#' sim <- project(params, t_max = 20, effort = 0.5)
#' n <- sim@@n[21, , ]
#' n_pp <- sim@@n_pp[21, ]
#' getEncounter(params, n, n_pp)
#' }
getEncounter <- function(params, n = params@initial_n, 
                         n_pp = params@initial_n_pp,
                         B = params@initial_B) {

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
        phi_prey_background <- params@species_params$interaction_p *
            rowSums(sweep(
            params@pred_kernel, 3, params@dw_full * params@w_full * n_pp,
            "*", check.margin = FALSE), dims = 2)
        encounter <- params@search_vol * (phi_prey_species + phi_prey_background)
    } else {
        prey <- outer(params@species_params$interaction_p, n_pp)
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
        avail_energy <- Re(t(mvfft(t(params@ft_pred_kernel_e) * mvfft(t(prey)),
                                   inverse = TRUE))) / length(params@w_full)
        # Only keep the bit for fish sizes
        avail_energy <- avail_energy[, idx_sp, drop = FALSE]
        # Due to numerical errors we might get negative or very small entries that
        # should be 0
        avail_energy[avail_energy < 1e-18] <- 0
        
        encounter <- params@search_vol * avail_energy
    }
    dimnames(encounter) <- dimnames(params@metab)
    
    # Add contribution from unstructured resources
    # Can't use rowSums or colSums unfortunately because
    # the resource index that we want to sum over is the middle index.
    for (u in seq_along(B)) {
        encounter[] <- encounter + params@rho[, u, ] * B[u]
    }
    return(encounter)
}


#' Get feeding level
#'
#' Calculates the feeding level \eqn{f_i(w)} by predator size based on food
#' availability, search volume and maximum intake. The feeding level is the
#' proportion of the encountered food that is actually consumed. It is 
#' defined in terms of the encounter rate \eqn{E_i} and the maximum intake 
#' rate \eqn{h_i(w)} as
#' \deqn{f_i(w) = \frac{E_i(w)}{E_i(w)+h_i(w)}}{E_i(w)/(E_i(w)+h_i(w))}
#' The feeding rate is used in \code{\link{getEReproAndGrowth}} and in
#' \code{\link{getPredRate}}.
#' 
#' @param object A \code{MizerParams} or \code{MizerSim} object
#' @param n A matrix of species abundance (species x size). Only used if
#'   \code{object} argument is of type \code{MizerParams}.
#' @param n_pp A vector of the plankton abundance by size. Only used if
#'   \code{object} argument is of type \code{MizerParams}.
#' @param B A vector of biomasses of unstructured resource components
#' @param encounter The encounter rate matrix (optional) of dimension no.
#'   species x no. size bins. If not passed in, it is calculated internally
#'   using \code{\link{getEncounter}}. Only used if \code{object}
#'   argument is of type \code{MizerParams}.
#' @param time_range Subset the returned fishing mortalities by time. The time
#'   range is either a vector of values, a vector of min and max time, or a
#'   single value. Default is the whole time range. Only used if the
#'   \code{object} argument is of type \code{MizerSim}.
#' @param drop should extra dimensions of length 1 in the output be dropped,
#'   simplifying the output. Defaults to FALSE.
#'
#' @return If a \code{MizerParams} object is passed in, the function returns a two
#'   dimensional array (predator species x predator size) based on the
#'   abundances also passed in.
#'   If a \code{MizerSim} object is passed in, the function returns a three
#'   dimensional array (time step x predator species x predator size) with the
#'   feeding level calculated at every time step in the simulation.
#'   If `drop = TRUE` then the dimension of length 1 will be removed from the
#'   returned array.
#'   
#' @seealso \code{\link{getEncounter}}
#' @export
#' @family rate functions
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- set_multispecies_model(NS_species_params_gears, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the feeding level at one time step
#' n <- sim@@n[21,,]
#' n_pp <- sim@@n_pp[21,]
#' fl <- getFeedingLevel(params,n,n_pp)
#' # Get the feeding level at all saved time steps
#' fl <- getFeedingLevel(sim)
#' # Get the feeding level for time 15 - 20
#' fl <- getFeedingLevel(sim, time_range = c(15,20))
#' }
getFeedingLevel <- function(object, n, n_pp, B, encounter,
                            time_range, drop = FALSE) {
    if (is(object, "MizerParams")) {
        params <- object
        if (missing(encounter)) {
            if (missing(n)) n <- params@initial_n
            if (missing(n_pp)) n_pp <- params@initial_n_pp
            if (missing(B)) B <- params@initial_B
            encounter <- getEncounter(params, n, n_pp, B)
        }
        # Check dims of encounter
        if (!all(dim(encounter) == c(nrow(params@species_params),
                                  length(params@w)))) {
            stop("encounter argument must have dimensions: no. species (",
                 nrow(params@species_params), ") x no. size bins (",
                 length(params@w), ")")
        }
        # calculate feeding level
        f <- encounter / (encounter + params@intake_max)
        return(f)
    } else {
        sim <- object
        if (missing(time_range)) {
            time_range <- dimnames(sim@n)$time
        }
        time_elements <- get_time_elements(sim, time_range)
        feed_time <- plyr::aaply(which(time_elements), 1, function(x) {
            # Necessary as we only want single time step but may only have 1
            # species which makes using drop impossible
            n <- array(sim@n[x, , ], dim = dim(sim@n)[2:3])
            dimnames(n) <- dimnames(sim@n)[2:3]
            B <- sim@params@initial_B
            B[] <- sim@B[x, ]
            feed <- getFeedingLevel(sim@params, n = n,
                                    n_pp = sim@n_pp[x, ],
                                    B = B)
            return(feed)
            }, .drop = drop)
        return(feed_time)
    }
}



#' Get predation rate
#' 
#' Calculates the potential rate (in units 1/year) at which a prey individual of
#' a given size \eqn{w} is killed by predators from species \eqn{j}. In formulas
#' \deqn{{\tt pred\_rate}_j(w_p) = \int \phi_j(w,w_p) (1-f_j(w)) 
#'   \gamma_j(w) N_j(w) \, dw.}{pred_rate_j(w_p) = \int\phi_i(w,w_p) (1-f_i(w)) 
#'   \gamma_i(w) N_i(w) dw.}
#' This potential rate is used in the function \code{\link{getPredMort}} to
#' calculate the realised predation mortality rate on the prey individual.
#'
#' @param params An \linkS4class{MizerParams} object
#' @param n A matrix of species abundances (species x size)
#' @param n_pp A vector of the plankton abundance by size
#' @param B A vector of biomasses of unstructured resource components
#' @param feeding_level The current feeding level (optional). A matrix of size
#'   no. species x no. size bins. If not supplied, is calculated internally
#'   using the \code{\link{getFeedingLevel}} function.
#'   
#' @return A two dimensional array (predator species x prey size), 
#'   where the prey size runs over fish community plus plankton spectrum.
#' @export
#' @family rate functions
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- set_multispecies_model(NS_species_params_gears, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the feeding level at one time step
#' n <- sim@@n[21,,]
#' n_pp <- sim@@n_pp[21,]
#' getPredRate(params,n,n_pp)
#' }

getPredRate <- function(params, n = params@initial_n, 
                        n_pp = params@initial_n_pp,
                        B = params@initial_B,
                        feeding_level = getFeedingLevel(params, n = n,
                                                        n_pp = n_pp, B = B)
                        ) {
    no_sp <- dim(params@interaction)[1]
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)
    if (!all(dim(feeding_level) == c(no_sp, no_w))) {
        stop("feeding_level argument must have dimensions: no. species (",
             no_sp, ") x no. size bins (", no_w, ")")
    }
    
    # If the feeding kernel does not have a fixed predator/prey mass ratio
    # then the integral is not a convolution integral and we can not use fft.
    # In this case we use the code from mizer version 0.3
    if (length(params@ft_pred_kernel_p) == 1) {
        n_total_in_size_bins <- sweep(n, 2, params@dw, '*', check.margin = FALSE)
        # The next line is a bottle neck
        pred_rate <- sweep(params@pred_kernel, c(1,2),
                           (1-feeding_level) * params@search_vol * 
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
    
    dimnames(pred_rate) <- list(sp = params@species_params$species,
                                w_prey = names(n_pp))
    return(pred_rate)
}


#' Get total predation mortality rate
#'
#' Calculates the total predation mortality rate \eqn{\mu_{p,i}(w_p)} (in units
#' of 1/year) on each prey species by prey size:
#' \deqn{\mu_{p.i}(w_p) = \sum_j {\tt pred\_rate}_j(w_p)\, \theta_{ji}.}{
#'   \mu_{p.i}(w_p) = \sum_j pred_rate_j(w_p) \theta_{ji}.}
#' 
#' @param object A \code{MizerParams} or \code{MizerSim} object.
#' @param n A matrix of species abundance (species x size). Only used if
#'   \code{object} argument is of type \code{MizerParams}.
#' @param n_pp A vector of the plankton abundance by size. Only used if
#'   \code{object} argument is of type \code{MizerParams}.
#' @param B A vector of biomasses of unstructured resource components. Only used
#'   if \code{object} argument is of type \code{MizerParams}.
#' @param pred_rate An array of predation rates of dimension no. sp x no.
#'   community size bins x no. of size bins in whole spectra (i.e. community +
#'   plankton, the w_full slot). The array is optional. If it is not provided
#'   it is calculated by the \code{\link{getPredRate}} function.
#' @param time_range Subset the returned fishing mortalities by time. The time
#'   range is either a vector of values, a vector of min and max time, or a
#'   single value. Default is the whole time range. Only used if the
#'   \code{object} argument is of type \code{MizerSim}.
#' @param drop Only used when object is of type \code{MizerSim}. Should
#'   dimensions of length 1 in the output be dropped, simplifying the output.
#'   Defaults to TRUE
#'
#' @return
#'   If a \code{MizerParams} object is passed in, the function returns a two
#'   dimensional array (prey species x prey size) based on the abundances also
#'   passed in. If a \code{MizerSim} object is passed in, the function returns a
#'   three dimensional array (time step x prey species x prey size) with the
#'   predation mortality calculated at every time step in the simulation.
#'   Dimensions may be dropped if they have length 1 unless `drop = FALSE`.
#' @family rate functions
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- set_multispecies_model(NS_species_params_gears, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get predation mortality at one time step
#' n <- sim@@n[21,,]
#' n_pp <- sim@@n_pp[21,]
#' getPredMort(params,n,n_pp)
#' # Get predation mortality at all saved time steps
#' getPredMort(sim)
#' # Get predation mortality over the time 15 - 20
#' getPredMort(sim, time_range = c(15,20))
#' }
getPredMort <- function(object, n, n_pp, B, 
                        pred_rate, time_range, drop = TRUE) {
    if (is(object, "MizerParams")) {
        params <- object
        if (missing(n)) n <- params@initial_n
        if (missing(n_pp)) n_pp <- params@initial_n_pp
        if (missing(B)) B <- params@initial_B
        if (missing(pred_rate)) {
            feeding_level <- getFeedingLevel(params, n = n, n_pp = n_pp, B = B)

            pred_rate <- getPredRate(params, n = n, n_pp = n_pp, 
                                     B = B, feeding_level = feeding_level)
        }
        idx_sp <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)

        m2 <- (t(params@interaction) %*% pred_rate)[, idx_sp, drop = FALSE]
        return(m2)
    } else {
        sim <- object
        if (missing(time_range)) {
            time_range <- dimnames(sim@n)$time
        }
        time_elements <- get_time_elements(sim, time_range)
        m2_time <- plyr::aaply(which(time_elements), 1, function(x) {
            n <- array(sim@n[x, , ], dim = dim(sim@n)[2:3])
            dimnames(n) <- dimnames(sim@n)[2:3]
            B <- sim@params@initial_B
            B[] <- sim@B[x, ]
            m2 <- getPredMort(sim@params, n = n, 
                              n_pp = sim@n_pp[x, ], B = B)
            return(m2)
        }, .drop = drop)
        return(m2_time)
    }
}

#' Alias for getPredMort
#' 
#' An alias provided for backward compatibility with mizer version <= 1.0
#' @inherit getPredMort
#' @export
getM2 <- getPredMort


#' Get predation mortality rate for plankton
#' 
#' Calculates the predation mortality rate \eqn{\mu_p(w)} on the plankton
#' spectrum by plankton size (in units 1/year).
#' 
#' Used by the \code{project} function for running size based simulations.
#' 
#' @param params An \linkS4class{MizerParams} object
#' @param n A matrix of species abundances (species x size)
#' @param n_pp A vector of the plankton abundance by size
#' @param B A vector of biomasses of unstructured resource components
#' @param pred_rate An array of predation rates of dimension no. sp x no.
#'   community size bins x no. of size bins in whole spectra (i.e. community +
#'   plankton, the w_full slot). The array is optional. If it is not provided
#'   it is calculated by the \code{\link{getPredRate}} function.
#'
#' @return A vector of mortality rate by plankton size.
#' @family rate functions
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- set_multispecies_model(NS_species_params_gears, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get plankton mortality at one time step
#' n <- sim@@n[21,,]
#' n_pp <- sim@@n_pp[21,]
#' getPlanktonMort(params,n,n_pp)
#' }
getPlanktonMort <- 
    function(params, n = params@initial_n, 
             n_pp = params@initial_n_pp,
             B = params@initial_B,
             pred_rate = getPredRate(params, n = n, n_pp = n_pp, B = B)) {

    if ( (!all(dim(pred_rate) ==
               c(nrow(params@species_params), length(params@w_full)))) |
         (length(dim(pred_rate)) != 2)) {
        stop("pred_rate argument must have 2 dimensions: no. species (",
             nrow(params@species_params),
             ") x no. size bins in community + plankton (",
             length(params@w_full), ")")
    }
    return(as.vector(params@species_params$interaction_p %*% pred_rate))
}

#' Alias for getPlanktonMort
#' 
#' An alias provided for backward compatibility with mizer version <= 1.0
#' @inherit getPlanktonMort
#' @export
getM2Background <- getPlanktonMort

#' Get the fishing mortality by time, gear, species and size
#'
#' Calculates the fishing mortality rate \eqn{F_{g,i,w}} by gear, species and
#' size at each time step in the \code{effort} argument (in units 1/year).
#' Used by the \code{project} function to perform simulations.
#' 
#' @param object A \code{MizerParams} object or a \code{MizerSim} object.
#' @param effort The effort of each fishing gear. Only needed if the object
#'   argument is of class \code{MizerParams}. See notes below.
#' @param time_range Subset the returned fishing mortalities by time. The time
#'   range is either a vector of values, a vector of min and max time, or a
#'   single value. Default is the whole time range. Only used if the
#'   \code{object} argument is of type \code{MizerSim}.
#'   
#' @return An array. If the effort argument has a time dimension, or a
#'   \code{MizerSim} is passed in, the output array has four dimensions (time x
#'   gear x species x size). If the effort argument does not have a time
#'   dimension (i.e. it is a vector or a single numeric), the output array has
#'   three dimensions (gear x species x size).
#' @note Here: fishing mortality = catchability x selectivity x effort.
#' 
#' The \code{effort} argument is only used if a \code{MizerParams} object is
#' passed in. The \code{effort} argument can be a two dimensional array (time x
#' gear), a vector of length equal to the number of gears (each gear has a
#' different effort that is constant in time), or a single numeric value (each
#' gear has the same effort that is constant in time). The order of gears in the
#' \code{effort} argument must be the same the same as in the \code{MizerParams}
#' object.
#' 
#' If the object argument is of class \code{MizerSim} then the effort slot of
#' the \code{MizerSim} object is used and the \code{effort} argument is not
#' used.
#' @export
#' @family rate functions
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- set_multispecies_model(NS_species_params_gears, inter)
#' # Get the fishing mortality when effort is constant
#' # for all gears and time:
#' getFMortGear(params, effort = 1)
#' # Get the fishing mortality when effort is different
#' # between the four gears but constant in time:
#' getFMortGear(params, effort = c(0.5,1,1.5,0.75))
#' # Get the fishing mortality when effort is different
#' # between the four gears and changes with time:
#' effort <- array(NA, dim = c(20,4))
#' effort[,1] <- seq(from=0, to = 1, length=20)
#' effort[,2] <- seq(from=1, to = 0.5, length=20)
#' effort[,3] <- seq(from=1, to = 2, length=20)
#' effort[,4] <- seq(from=2, to = 1, length=20)
#' getFMortGear(params, effort=effort)
#' # Get the fishing mortality using the effort already held in a MizerSim object.
#' sim <- project(params, t_max = 20, effort = 0.5)
#' getFMortGear(sim)
#' getFMortGear(sim, time_range=c(10,20))
#' }
#' 
getFMortGear <- function(object, effort = 1, time_range) {
    if (is(object, "MizerSim")) {
        sim <- object
        if (missing(time_range)) {
            time_range <- dimnames(sim@effort)$time
        }
        time_elements <- get_time_elements(sim, time_range, slot_name = "effort")
        f_mort_gear <- getFMortGear(sim@params, sim@effort)
        return(f_mort_gear[time_elements, , , , drop = FALSE])
    } else {
        params <- object
        if (is(effort, "numeric")) {
            no_gear <- dim(params@catchability)[1]
            # If a single value, just repeat it for all gears
            if (length(effort) == 1) {
                effort <- rep(effort, no_gear)
            }
            if (length(effort) != no_gear) {
                stop("Effort must be a single value or a vector as long as the number of gears\n")
            }
            # Streamlined for speed increase - note use of recycling
            out <- params@selectivity
            out[] <- effort * c(params@catchability) * c(params@selectivity)
            return(out)
        } else {
            # assuming effort is a matrix, and object is of MizerParams class
            no_gear <- dim(params@catchability)[1]
            if (dim(effort)[2] != no_gear)
                stop("Effort array must have a single value or a vector as long as the number of gears for each time step\n")
            # Make the output array - note that we put time as last dimension
            # and then aperm before returning. This is because of the order of
            # the values when we call the other getFMortGear function.
            # Fill it up by calling the other function and passing in each line
            # of the effort matrix
            out <- array(NA, dim = c(dim(params@selectivity), dim(effort)[1]),
                         dimnames = c(dimnames(params@selectivity),
                                      list(time = dimnames(effort)[[1]])))
            out[] <- apply(effort, 1, function(x) getFMortGear(params, x))
            out <- aperm(out, c(4, 1, 2, 3))
            return(out)
        }
    }
}


#' Get the total fishing mortality rate from all fishing gears by time, species
#' and size.
#' 
#' Calculates the total fishing mortality  (in units 1/year) from all gears by
#' species and size at each time step in the \code{effort} argument.
#' The total fishing mortality is just the sum of the fishing mortalities
#' imposed by each gear, \eqn{\mu_{f.i}(w)=\sum_g F_{g,i,w}}.
#' 
#' @param object A \code{MizerParams} object or a \code{MizerSim} object
#' @param effort The effort of each fishing gear. Only needed if the object
#'   argument is of class \code{MizerParams}. See notes below.
#' @param time_range Subset the returned fishing mortalities by time. The time
#'   range is either a vector of values, a vector of min and max time, or a
#'   single value. Default is the whole time range. Only used if the
#'   \code{object} argument is of type \code{MizerSim}.
#' @param drop Only used when object is of type \code{MizerSim}. Should
#'   dimensions of length 1 be dropped, e.g. if your community only has one
#'   species it might make presentation of results easier. Default is TRUE
#'
#' @return An array. If the effort argument has a time dimension, or object is
#'   of class \code{MizerSim}, the output array has three dimensions (time x
#'   species x size). If the effort argument does not have a time dimension, the
#'   output array has two dimensions (species x size).
#' @note Here: fishing mortality = catchability x selectivity x effort.
#'
#' The \code{effort} argument is only used if a \code{MizerParams} object is
#' passed in. The \code{effort} argument can be a two dimensional array (time x
#' gear), a vector of length equal to the number of gears (each gear has a
#' different effort that is constant in time), or a single numeric value (each
#' gear has the same effort that is constant in time). The order of gears in the
#' \code{effort} argument must be the same the same as in the \code{MizerParams}
#' object.
#'
#' If the object argument is of class \code{MizerSim} then the effort slot of
#' the \code{MizerSim} object is used and the \code{effort} argument is not
#' used.
#' @export
#' @family rate functions
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- set_multispecies_model(NS_species_params_gears, inter)
#' # Get the total fishing mortality when effort is constant for all 
#' # gears and time:
#' getFMort(params, effort = 1)
#' # Get the total fishing mortality when effort is different
#' # between the four gears but constant in time:
#' getFMort(params, effort = c(0.5,1,1.5,0.75))
#' # Get the total fishing mortality when effort is different
#' # between the four gears and changes with time:
#' effort <- array(NA, dim = c(20,4))
#' effort[,1] <- seq(from=0, to = 1, length=20)
#' effort[,2] <- seq(from=1, to = 0.5, length=20)
#' effort[,3] <- seq(from=1, to = 2, length=20)
#' effort[,4] <- seq(from=2, to = 1, length=20)
#' getFMort(params, effort=effort)
#' # Get the total fishing mortality using the effort already held in a 
#' # MizerSim object.
#' sim <- project(params, t_max = 20, effort = 0.5)
#' getFMort(sim)
#' getFMort(sim, time_range = c(10,20))
#' }
getFMort <- function(object, effort, time_range, drop=TRUE){
    if (is(object, "MizerParams")) {
        params <- object
        if (is(effort, "numeric")) {
            f_mort_gear <- getFMortGear(params, effort)
            f_mort <- colSums(f_mort_gear)
            return(f_mort)
        } else {
            #assuming effort is a matrix
            f_mort_gear <- getFMortGear(params, effort)
            f_mort <- apply(f_mort_gear, c(1, 3, 4), sum)
            return(f_mort)
        }
    } else {
        #case where object is MizerSim, and we use effort from there
        sim <- object
        if (missing(time_range)) {
            time_range <- dimnames(sim@effort)$time
        }
        time_elements <- get_time_elements(sim, time_range, slot_name = "effort")
        f_mort <- getFMort(sim@params, sim@effort)
        return(f_mort[time_elements, , , drop = drop])
    }}


#' Get total mortality rate
#'
#' Calculates the total mortality rate \eqn{\mu_i(w)}  (in units 1/year) on each
#' species by size from predation mortality, background mortality and fishing
#' mortality for a single time step.
#' @param params An \linkS4class{MizerParams} object
#' @param n A matrix of species abundances (species x size)
#' @param n_pp A vector of the plankton abundance by size
#' @param B A vector of biomasses of unstructured resource components
#' @param effort A numeric vector of the effort by gear or a single numeric
#'   effort value which is used for all gears. Default 1.
#' @param m2 A two dimensional array of predation mortality (optional). Has
#'   dimensions no. sp x no. size bins in the community. If not supplied is
#'   calculated using the \code{\link{getPredMort}} function.
#'
#' @return A two dimensional array (prey species x prey size). 
#'
#' @export
#' @seealso \code{\link{getPredMort}}, \code{\link{getFMort}}
#' @family rate functions
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- set_multispecies_model(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the total mortality at a particular time step
#' getMort(params,sim@@n[21,,],sim@@n_pp[21,],effort=0.5)
#' }
getMort <- function(params, n = params@initial_n, 
                    n_pp = params@initial_n_pp,
                    B = params@initial_B,
                    effort = 1, 
                    m2 = getPredMort(params, n = n, n_pp = n_pp, B = B)) {
    if (!all(dim(m2) == c(nrow(params@species_params), length(params@w)))) {
        stop("m2 argument must have dimensions: no. species (",
             nrow(params@species_params), ") x no. size bins (",
             length(params@w), ")")
    }
    return(m2 + params@mu_b + getFMort(params, effort = effort))
}

#' Alias for getMort
#' 
#' An alias provided for backward compatibility with mizer version <= 1.0
#' @inherit getMort
#' @export
getZ <- getMort


#' Get energy rate available for reproduction and growth
#'
#' Calculates the energy rate (grams/year) available by species and size for
#' reproduction and growth after metabolism and movement have been accounted
#' for. 
#' 
#' @param params An \linkS4class{MizerParams} object
#' @param n A matrix of species abundances (species x size)
#' @param n_pp A vector of the plankton abundance by size
#' @param B A vector of biomasses of unstructured resource components
#' @param encounter The encounter rate matrix (optional) of dimension no.
#'   species x no. size bins. If not passed in, it is calculated internally
#'   using the \code{\link{getEncounter}} function.
#' @param feeding_level The current feeding level (optional). A matrix of size
#'   no. species x no. size bins. If not supplied, is calculated internally
#'   using the \code{\link{getFeedingLevel}} function.
#'
#' @return A two dimensional array (species x size) holding
#' \deqn{E_{r.i}(w) = \max(0, \alpha_i\, (1 - {\tt feeding\_level}_i(w))\, 
#'                            {\tt encounter}_i(w) - {\tt metab}_i(w)).}{
#'   E_{r.i}(w) = max(0, \alpha_i * (1 - feeding_level_i(w)) * 
#'                       encounter_i(w) - metab_i(w)).}
#' The assimilation rate \eqn{\alpha_i} is taken from the species parameter
#' data frame in \code{params}. The metabolic rate \code{metab} is taken from 
#' \code{params}. 
#' @export
#' @seealso \code{\link{getERepro}} and \code{\link{getEGrowth}}
#' @family rate functions
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- set_multispecies_model(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the energy at a particular time step
#' getEReproAndGrowth(params,sim@@n[21,,],sim@@n_pp[21,])
#' }
getEReproAndGrowth <- function(params, n = params@initial_n, 
                               n_pp = params@initial_n_pp,
                               B = params@initial_B,
                               encounter = getEncounter(params, n = n,
                                                        n_pp = n_pp, B = B),
                               feeding_level = getFeedingLevel(params, n = n,
                                                               n_pp = n_pp, B = B,
                                                               encounter = encounter)) {
    if (!all(dim(feeding_level) == c(nrow(params@species_params), length(params@w)))) {
        stop("feeding_level argument must have dimensions: no. species (",
             nrow(params@species_params), ") x no. size bins (",
             length(params@w), ")")
    }
    # assimilated intake
    e <- sweep((1 - feeding_level) * encounter, 1,
               params@species_params$alpha, "*", check.margin = FALSE)
    # Subtract metabolism
    e <- e - params@metab
    e[e < 0] <- 0 # Do not allow negative growth
    return(e)
}


#' Get energy rate available for reproduction
#'
#' Calculates the energy rate (grams/year) available by species and size for
#' reproduction after metabolism and movement have been accounted for:
#' \eqn{\psi_i(w)E_{r.i}(w)}. Used by the \code{project} function for performing
#' simulations.
#' 
#' @param params An \linkS4class{MizerParams} object
#' @param n A matrix of species abundances (species x size)
#' @param n_pp A vector of the plankton abundance by size
#' @param B A vector of biomasses of unstructured resource components
#' @param e The energy available for reproduction and growth (optional). A
#'   matrix of size no. species x no. size bins. If not supplied, is calculated
#'   internally using \code{\link{getEReproAndGrowth}}.
#'
#' @return A two dimensional array (prey species x prey size) 
#' @export
#' @seealso \code{\link{getERepro}} and \code{\link{getEReproAndGrowth}}.
#' @family rate functions
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- set_multispecies_model(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the energy at a particular time step
#' getERepro(params,sim@@n[21,,],sim@@n_pp[21,])
#' }
getERepro <- function(params, n = params@initial_n, 
                      n_pp = params@initial_n_pp,
                      B = params@initial_B,
                      e = getEReproAndGrowth(params, n = n, n_pp = n_pp, B = B)) {
    if (!all(dim(e) == c(nrow(params@species_params), length(params@w)))) {
        stop("e argument must have dimensions: no. species (",
             nrow(params@species_params), ") x no. size bins (",
             length(params@w), ")")
    }
    e_repro <- params@psi * e
    return(e_repro)
}

#' Alias for getERepro
#' 
#' An alias provided for backward compatibility with mizer version <= 1.0
#' @inherit getERepro
#' @export
getESpawning <- getERepro


#' Get energy rate available for growth
#'
#' Calculates the energy rate \eqn{g_i(w)} (grams/year) available by species and
#' size for growth after metabolism, movement and reproduction have been
#' accounted for. Used by \code{\link{project}} for performing simulations.
#' @param params A \linkS4class{MizerParams} object
#' @param n A matrix of species abundances (species x size)
#' @param n_pp A vector of the plankton abundance by size
#' @param B A vector of biomasses of unstructured resource components
#' @param e The energy available for reproduction and growth (optional, although
#'   if specified, e_repro must also be specified). A matrix of size no.
#'   species x no. size bins. If not supplied, is calculated internally using
#'   \code{\link{getEReproAndGrowth}}.
#' @param e_repro The energy available for reproduction (optional, although if
#'   specified, e must also be specified). A matrix of size no. species x no.
#'   size bins. If not supplied, is calculated internally using
#'   \code{\link{getERepro}}.
#'   
#' @return A two dimensional array (prey species x prey size) 
#' @export
#' @seealso \code{\link{getERepro}}, \code{\link{getEReproAndGrowth}}
#' @family rate functions
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- set_multispecies_model(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the energy at a particular time step
#' getEGrowth(params,sim@@n[21,,],sim@@n_pp[21,])
#' }
getEGrowth <- function(params, n = params@initial_n, 
                       n_pp = params@initial_n_pp,
                       B = params@initial_B,
                       e_repro = getERepro(params, n = n, n_pp = n_pp, B = B),
                       e=getEReproAndGrowth(params, n = n, n_pp = n_pp, B = B)) {
    if (!all(dim(e_repro) == c(nrow(params@species_params), length(params@w)))) {
        stop("e_repro argument must have dimensions: no. species (",
             nrow(params@species_params), ") x no. size bins (",
             length(params@w), ")")
    }
    if (!all(dim(e) == c(nrow(params@species_params), length(params@w)))) {
        stop("e argument must have dimensions: no. species (",
             nrow(params@species_params), ") x no. size bins (",
             length(params@w), ")")
    }
    # Assimilated intake less activity and metabolism
    # energy for growth is intake - energy for reproduction
    e_growth <- e - e_repro
    return(e_growth)
}


#' Get density independent rate of egg production
#'
#' Calculates the density independent rate of egg production \eqn{R_{p.i}}
#' (units 1/year) before density dependence, by species. Used by
#' \code{\link{getRDD}} to calculate the actual density dependent rate.
#' See \code{\link{setRecruitment}} for more details.
#' 
#' @param params An \linkS4class{MizerParams} object
#' @param n A matrix of species abundances (species x size)
#' @param n_pp A vector of the plankton abundance by size
#' @param B A vector of biomasses of unstructured resource components
#' @param e_repro The energy available for reproduction (optional). A matrix of
#'   size no. species x no. size bins. If not supplied, is calculated internally
#'   using \code{\link{getERepro}}.
#' @param sex_ratio Proportion of the population that is female. Default value
#'   is 0.5.
#'   
#' @return A numeric vector the length of the number of species 
#' @export
#' @seealso \code{\link{getRDD}}
#' @family rate functions
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- set_multispecies_model(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the recruitment at a particular time step
#' getRDI(params,sim@@n[21,,],sim@@n_pp[21,])
#' }
getRDI <- function(params, n = params@initial_n, 
                   n_pp = params@initial_n_pp,
                   B = params@initial_B,
                   e_repro = getERepro(params, n = n, n_pp = n_pp, B = B),
                   sex_ratio = 0.5) {
    if (!all(dim(e_repro) == c(nrow(params@species_params), length(params@w)))) {
        stop("e_repro argument must have dimensions: no. species (",
             nrow(params@species_params), ") x no. size bins (",
             length(params@w), ")")
    }
    e_repro_pop <- drop( (e_repro * n) %*% params@dw)
    rdi <- sex_ratio * (e_repro_pop * params@species_params$erepro) /
        params@w[params@w_min_idx]
    return(rdi)
}


#' Get density dependent rate of larvae production
#'
#' Calculates the density dependent rate of larvae production \eqn{R_i} (units
#' 1/year) for each species. This is the flux entering the smallest size class
#' of each species. The density dependent rate is the density independent
#' rate obtained with \code{\link{getRDI}} after it has been put through the 
#' density dependent "stock-recruitment" relationship function. See
#' \code{\link{setRecruitment}} for more details.
#' 
#' @param params An \linkS4class{MizerParams} object
#' @param n A matrix of species abundances (species x size)
#' @param n_pp A vector of the plankton abundance by size
#' @param B A vector of biomasses of unstructured resource components
#' @param sex_ratio Proportion of the population that is female. Default value
#'   is 0.5.
#' @param rdi A vector of density independent recruitment for each species. 
#'   If not specified rdi is calculated internally using
#'   \code{\link{getRDI}}.
#'   
#' @return A numeric vector the length of the number of species. 
#' @export
#' @seealso \code{\link{getRDI}}
#' @family rate functions
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- set_multispecies_model(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the rate at a particular time step
#' getRDD(params,sim@@n[21,,],sim@@n_pp[21,])
#' }
getRDD <- function(params, n = params@initial_n, 
                   n_pp = params@initial_n_pp,
                   B = params@initial_B, sex_ratio = 0.5,
                   rdi = getRDI(params, n = n, n_pp = n_pp, B = B, sex_ratio = sex_ratio)) {
    rdd <- params@srr(rdi = rdi, species_params = params@species_params)
    return(rdd)
}

#' Get_time_elements
#'
#' Internal function to get the array element references of the time dimension
#' for the time based slots of a MizerSim object.
#' @param sim A MizerSim object
#' @param time_range The time_range can be character or numeric.
#' @param slot_name Necessary to include a slot_name argument because the effort
#'   and abundance slots have different time dimensions
#' @export
#' @concept helper
#' @keywords internal
get_time_elements <- function(sim, time_range, slot_name = "n"){
    assert_that(is(sim, "MizerSim"))
    time_range <- range(as.numeric(time_range))
    # Check that time range is even in object
    sim_times <- as.numeric(dimnames(sim@effort)$time)
    sim_time_range <- range(sim_times)
    if ( (time_range[1] < sim_time_range[1]) |
         (time_range[2] > sim_time_range[2]))
        stop("Time range is outside the time range of the model")
    time_elements <- (sim_times >= time_range[1]) & (sim_times <= time_range[2])
    names(time_elements) <- dimnames(sim@effort)$time
    return(time_elements)
}
