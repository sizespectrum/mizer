#' Methods used for projecting
#'
#' The functions defined in the file project_methods calculate the various
#' quantities needed to project the size-spectra forward in time, using the
#' model described in section 3 of the mizer vignette.
#'
#' @section List of functions:
#' In this list we relate the functions in this file to the quantities named in
#' the mizer vignette.
#' \tabular{llll}{
#'   Method name \tab Expression \tab Description \tab Section in vignette\cr
#'   \code{\link{getAvailEnergy}} \tab \eqn{E_{a.i}(w)} \tab Available energy \tab 3.2 \cr
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

#' Get available energy
#' 
#' Calculates the amount \eqn{E_{a,i}(w)} of food exposed to each predator as
#' a function of predator size. 
#' 
#' This function is used by the \code{\link{project}} method for
#' performing simulations.
#' 
#' The function returns values also for sizes outside the size-range of the
#' species. These values should not be trusted, as they are meaningless.
#' 
#' @param object An \linkS4class{MizerParams} object
#' @param n A matrix of species abundances (species x size)
#' @param n_pp A vector of the plankton abundance by size
#'   
#' @return A two dimensional array (predator species x predator size)
#' @seealso \code{\link{project}}
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' n <- sim@@n[21,,]
#' n_pp <- sim@@n_pp[21,]
#' getAvailEnergy(params,n,n_pp)
#' }
getAvailEnergy <- function(object, n, n_pp) {

    # idx_sp are the index values of object@w_full such that
    # object@w_full[idx_sp] = object@w
    idx_sp <- (length(object@w_full) - length(object@w) + 1):length(object@w_full)
    
    # If the feeding kernel does not have a fixed predator/prey mass ratio
    # then the integral is not a convolution integral and we can not use fft.
    # In this case we use the code from mizer version 0.3
    if (length(object@ft_pred_kernel_e) == 1) {
        # n_eff_prey is the total prey abundance by size exposed to each
        # predator (prey not broken into species - here we are just working out
        # how much a predator eats - not which species are being eaten - that is
        # in the mortality calculation
        n_eff_prey <- sweep(object@interaction %*% n, 2, 
                            object@w * object@dw, "*", check.margin = FALSE) 
        # pred_kernel is predator species x predator size x prey size
        # So multiply 3rd dimension of pred_kernel by the prey abundance
        # Then sum over 3rd dimension to get total eaten by each predator by 
        # predator size
        # This line is a bottle neck
        phi_prey_species <- rowSums(sweep(
            object@pred_kernel[, , idx_sp, drop = FALSE],
            c(1, 3), n_eff_prey, "*", check.margin = FALSE), dims = 2)
        # Eating the background
        # This line is a bottle neck
        phi_prey_background <- rowSums(sweep(
            object@pred_kernel, 3, object@dw_full * object@w_full * n_pp,
            "*", check.margin = FALSE), dims = 2)
        return(phi_prey_species + phi_prey_background)
    }

    prey <- matrix(0, nrow = dim(n)[1], ncol = length(object@w_full))
    # Looking at Equation (3.4), for available energy in the mizer vignette,
    # we have, for our predator species i, that prey[k] equals
    # the sum over all species j of fish, of theta_{i,j}*N_j(wFull[k])
    prey[, idx_sp] <- object@interaction %*% n
    # The vector f2 equals everything inside integral (3.4) except the feeding
    # kernel phi_i(w_p/w).
    prey <- sweep(sweep(prey, 2, n_pp, "+"), 2, 
                object@w_full * object@dw_full, "*")
    # Eq (3.4) is then a convolution integral in terms of prey[w_p] and phi[w_p/w].
    # We approximate the integral by the trapezoidal method. Using the
    # convolution theorem we can evaluate the resulting sum via fast fourier
    # transform.
    # mvfft() does a Fourier transform of each column of its argument, but
    # we need the Fourier transforms of each row, so we need to apply mvfft()
    # to the transposed matrices and then transpose again at the end.
    avail_energy <- Re(t(mvfft(t(object@ft_pred_kernel_e) * mvfft(t(prey)),
                               inverse = TRUE))) / length(object@w_full)
    # Only keep the bit for fish sizes
    avail_energy <- avail_energy[, idx_sp, drop = FALSE]
    # Due to numerical errors we might get negative or very small entries that
    # should be 0
    avail_energy[avail_energy < 1e-18] <- 0
    
    dimnames(avail_energy) <- dimnames(object@metab)
    return(avail_energy)
}

#' Alias for getAvailEnergy
#' 
#' An alias provided for backward compatibility with mizer version <= 1.0
#' @inherit getAvailEnergy
#' @export
getPhiPrey <- getAvailEnergy


#' Get feeding level
#'
#' Calculates the feeding level \eqn{f_i(w)} by predator size based on food
#' availability, search volume and maximum intake. The feeding level is the
#' proportion of the encountered food that is actually consumed. This method is
#' used by the \code{\link{project}} method for performing simulations.
#' @param object A \code{MizerParams} or \code{MizerSim} object
#' @param n A matrix of species abundance (species x size). Only used if
#'   \code{object} argument is of type \code{MizerParams}.
#' @param n_pp A vector of the plankton abundance by size. Only used if
#'   \code{object} argument is of type \code{MizerParams}.
#' @param avail_energy The available energy matrix (optional) of dimension no.
#'   species x no. size bins. If not passed in, it is calculated internally
#'   using the \code{\link{getAvailEnergy}} method. Only used if \code{object}
#'   argument is of type \code{MizerParams}.
#' @param time_range Subset the returned fishing mortalities by time. The time
#'   range is either a vector of values, a vector of min and max time, or a
#'   single value. Default is the whole time range. Only used if the
#'   \code{object} argument is of type \code{MizerSim}.
#' @param drop should extra dimensions of length 1 in the output be dropped,
#'   simplifying the output. Defaults to FALSE.
#'
#' @note If a \code{MizerParams} object is passed in, the method returns a two
#'   dimensional array (predator species x predator size) based on the
#'   abundances also passed in.
#'
#'   If a \code{MizerSim} object is passed in, the method returns a three
#'   dimensional array (time step x predator species x predator size) with the
#'   feeding level calculated at every time step in the simulation.
#' @seealso \code{\link{getAvailEnergy}}
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
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
getFeedingLevel <- function(object, n, n_pp, avail_energy,
                            time_range, drop = FALSE) {
    if (is(object, "MizerParams")) {
        if (missing(avail_energy)) {
            avail_energy <- getAvailEnergy(object, n, n_pp)
        }
        # Check dims of avail_energy

        if (!all(dim(avail_energy) == c(nrow(object@species_params),
                                        length(object@w)))) {
            stop("avail_energy argument must have dimensions: no. species (",
                 nrow(object@species_params), ") x no. size bins (",
                 length(object@w), ")")
        }
        # encountered food = available food * search volume
        encount <- object@search_vol * avail_energy
        # calculate feeding level
        f <- encount / (encount + object@intake_max)
        return(f)
    } else {
        if (missing(time_range)) {
            time_range <- dimnames(object@n)$time
        }
        time_elements <- get_time_elements(object, time_range)
        feed_time <- aaply(which(time_elements), 1, function(x) {
            # Necessary as we only want single time step but may only have 1
            # species which makes using drop impossible
            n <- array(object@n[x, , ], dim = dim(object@n)[2:3])
            dimnames(n) <- dimnames(object@n)[2:3]
            feed <- getFeedingLevel(object@params, n = n, n_pp = object@n_pp[x, ])
            return(feed)
            }, .drop = drop)
        return(feed_time)
    }
}



#' Get predation rate
#' 
#' Calculates the potential rate at which a prey individual of a given size 
#' \eqn{w} is killed by predators from species \eqn{i}. 
#' In formulas \deqn{\int\phi_i(w_p/w) (1-f_i(w)) \gamma_i w^q N_i(w) dw}
#' This potential rate is used in the function \code{\link{getPredMort}} to
#' calculate the realised predation mortality rate on the prey individual.
#' @param object A \code{MizerParams} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the plankton abundance by size.
#' @param feeding_level The current feeding level (optional). A matrix of size
#'   no. species x no. size bins. If not supplied, is calculated internally
#'   using the \code{getFeedingLevel()} method.
#'   
#' @return A two dimensional array (predator species x prey size), 
#'   where the prey size runs over fish community plus plankton spectrum.
#' @export
#' @seealso \code{\link{getPredMort}} and
#'   \code{\link{getFeedingLevel}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the feeding level at one time step
#' n <- sim@@n[21,,]
#' n_pp <- sim@@n_pp[21,]
#' getPredRate(params,n,n_pp)
#' }

getPredRate <- function(object, n,  n_pp,
                        feeding_level = getFeedingLevel(object, n = n, n_pp = n_pp)
                        ) {
    no_sp <- dim(object@interaction)[1]
    no_w <- length(object@w)
    no_w_full <- length(object@w_full)
    if (!all(dim(feeding_level) == c(no_sp, no_w))) {
        stop("feeding_level argument must have dimensions: no. species (",
             no_sp, ") x no. size bins (", no_w, ")")
    }
    
    # If the feeding kernel does not have a fixed predator/prey mass ratio
    # then the integral is not a convolution integral and we can not use fft.
    # In this case we use the code from mizer version 0.3
    if (length(object@ft_pred_kernel_p) == 1) {
        n_total_in_size_bins <- sweep(n, 2, object@dw, '*', check.margin = FALSE)
        # The next line is a bottle neck
        pred_rate <- sweep(object@pred_kernel, c(1,2),
                           (1-feeding_level) * object@search_vol * 
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
    Q[, idx_sp] <- sweep( (1 - feeding_level) * object@search_vol * n, 2,
                         object@dw, "*")

    # We do our spectral integration in parallel over the different species
    pred_rate <- Re(t(mvfft(t(object@ft_pred_kernel_p) *
                                 mvfft(t(Q)), inverse = TRUE))) / no_w_full
    # Due to numerical errors we might get negative or very small entries that
    # should be 0
    pred_rate[pred_rate < 1e-18] <- 0
    
    dimnames(pred_rate) <- list(sp = object@species_params$species,
                                w_prey = names(n_pp))
    return(pred_rate)
}


#' get predation mortality rate
#'
#' Calculates the total predation mortality rate \eqn{\mu_{p,i}(w_p)} on each
#' prey species by prey size. This method is used by the \code{\link{project}}
#' method for performing simulations.
#' @param object A \code{MizerParams} or \code{MizerSim} object.
#' @param n A matrix of species abundance (species x size). Only used if
#'   \code{object} argument is of type \code{MizerParams}.
#' @param n_pp A vector of the plankton abundance by size. Only used if
#'   \code{object} argument is of type \code{MizerParams}.
#' @param pred_rate An array of predation rates of dimension no. sp x no.
#'   community size bins x no. of size bins in whole spectra (i.e. community +
#'   plankton, the w_full slot). The array is optional. If it is not provided
#'   it is calculated by the \code{getPredRate()} method.
#' @param time_range Subset the returned fishing mortalities by time. The time
#'   range is either a vector of values, a vector of min and max time, or a
#'   single value. Default is the whole time range. Only used if the
#'   \code{object} argument is of type \code{MizerSim}.
#' @param drop Only used when object is of type \code{MizerSim}. Should
#'   dimensions of length 1 in the output be dropped, simplifying the output.
#'   Defaults to TRUE
#'
#' @return
#'   If a \code{MizerParams} object is passed in, the method returns a two
#'   dimensional array (prey species x prey size) based on the abundances also
#'   passed in. If a \code{MizerSim} object is passed in, the method returns a
#'   three dimensional array (time step x prey species x prey size) with the
#'   predation mortality calculated at every time step in the simulation.
#' @seealso \code{\link{getPredRate}} and \code{\link{project}}.
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
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
getPredMort <- function(object, n, n_pp, pred_rate, time_range, drop = TRUE) {
    if (is(object, "MizerParams")) {
        if (missing(pred_rate)) {
            feeding_level <- getFeedingLevel(object, n = n, n_pp = n_pp)

            pred_rate <- getPredRate(object = object, n = n,
                                     n_pp = n_pp, feeding_level = feeding_level)
        }
        idx_sp <- (length(object@w_full) - length(object@w) + 1):length(object@w_full)

        m2 <- (t(object@interaction) %*% pred_rate)[, idx_sp, drop = FALSE]
        return(m2)
    } else {
        if (missing(time_range)) {
            time_range <- dimnames(object@n)$time
        }
        time_elements <- get_time_elements(object, time_range)
        m2_time <- aaply(which(time_elements), 1, function(x) {
            n <- array(object@n[x, , ], dim = dim(object@n)[2:3])
            dimnames(n) <- dimnames(object@n)[2:3]
            m2 <- getPredMort(object@params, n = n, n_pp = object@n_pp[x, ])
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
#' spectrum by plankton size. Used by the \code{project} method for running size
#' based simulations.
#' @param object A \code{MizerParams} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the plankton abundance by size.
#' @param pred_rate An array of predation rates of dimension no. sp x no.
#'   community size bins x no. of size bins in whole spectra (i.e. community +
#'   plankton, the w_full slot). The array is optional. If it is not provided
#'   it is calculated by the \code{getPredRate()} method.
#'
#' @return A vector of mortality rate by plankton size.
#' @seealso \code{\link{project}} and \code{\link{getPredMort}}.
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get plankton mortality at one time step
#' n <- sim@@n[21,,]
#' n_pp <- sim@@n_pp[21,]
#' getPlanktonMort(params,n,n_pp)
#' }
getPlanktonMort <- 
    function(object, n, n_pp,
             pred_rate = getPredRate(object, n = n, n_pp = n_pp)) {

    if ( (!all(dim(pred_rate) ==
               c(nrow(object@species_params), length(object@w_full)))) |
         (length(dim(pred_rate)) != 2)) {
        stop("pred_rate argument must have 2 dimensions: no. species (",
             nrow(object@species_params),
             ") x no. size bins in community + plankton (",
             length(object@w_full), ")")
    }
    return(colSums(pred_rate))
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
#' size at each time step in the \code{effort} argument. 
#' Used by the \code{project} method to perform simulations.
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
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
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
getFMortGear <- function(object, effort, time_range) {
    if (is(object, "MizerSim")) {
        if (missing(time_range)) {
            time_range <- dimnames(object@effort)$time
        }
        time_elements <- get_time_elements(object, time_range, slot_name = "effort")
        f_mort_gear <- getFMortGear(object@params, object@effort)
        return(f_mort_gear[time_elements, , , , drop = FALSE])
    } else {
        if (is(effort, "numeric")) {
            no_gear <- dim(object@catchability)[1]
            # If a single value, just repeat it for all gears
            if (length(effort) == 1) {
                effort <- rep(effort, no_gear)
            }
            if (length(effort) != no_gear) {
                stop("Effort must be a single value or a vector as long as the number of gears\n")
            }
            # Streamlined for speed increase - note use of recycling
            out <- object@selectivity
            out[] <- effort * c(object@catchability) * c(object@selectivity)
            return(out)
        } else {
            # assuming effort is a matrix, and object is of MizerParams class
            no_gear <- dim(object@catchability)[1]
            if (dim(effort)[2] != no_gear)
                stop("Effort array must have a single value or a vector as long as the number of gears for each time step\n")
            # Make the output array - note that we put time as last dimension
            # and then aperm before returning. This is because of the order of
            # the values when we call the other getFMortGear method.
            # Fill it up by calling the other method and passing in each line
            # of the effort matrix
            out <- array(NA, dim = c(dim(object@selectivity), dim(effort)[1]),
                         dimnames = c(dimnames(object@selectivity),
                                      list(time = dimnames(effort)[[1]])))
            out[] <- apply(effort, 1, function(x) getFMortGear(object, x))
            out <- aperm(out, c(4, 1, 2, 3))
            return(out)
        }
    }
}


#' Get the total fishing mortality rate from all fishing gears by time, species
#' and size.
#' 
#' Calculates the total fishing mortality from all gears by species and size at 
#' each time step in the \code{effort} argument.
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
#' @seealso \code{getFMortGear}, \code{project}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
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
        if (is(effort, "numeric")) {
            f_mort_gear <- getFMortGear(object, effort)
            f_mort <- colSums(f_mort_gear)
            return(f_mort)
        } else {
            #assuming effort is a matrix
            f_mort_gear <- getFMortGear(object, effort)
            f_mort <- apply(f_mort_gear, c(1, 3, 4), sum)
            return(f_mort)
        }
    } else {
        #case where object is mizersim, and we use effort from there
        if (missing(time_range)) {
            time_range <- dimnames(object@effort)$time
        }
        time_elements <- get_time_elements(object, time_range, slot_name = "effort")
        f_mort <- getFMort(object@params, object@effort)
        return(f_mort[time_elements, , , drop = drop])
    }}


#' Get total mortality rate
#'
#' Calculates the total mortality rate \eqn{\mu_i(w)} on each species by size
#' from predation mortality, background mortality and fishing mortality
#' for a single time step.
#' @param object A \code{MizerParams} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the plankton abundance by size.
#' @param effort A numeric vector of the effort by gear or a single numeric
#'   effort value which is used for all gears.
#' @param m2 A two dimensional array of predation mortality (optional). Has
#'   dimensions no. sp x no. size bins in the community. If not supplied is
#'   calculated using the \code{getPredMort()} method.
#'
#' @return A two dimensional array (prey species x prey size). 
#'
#' @export
#' @seealso \code{\link{getPredMort}}, \code{\link{getFMort}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the total mortality at a particular time step
#' getMort(params,sim@@n[21,,],sim@@n_pp[21,],effort=0.5)
#' }
getMort <- function(object, n, n_pp, effort,
                 m2 = getPredMort(object, n = n, n_pp = n_pp)){
    if (!all(dim(m2) == c(nrow(object@species_params), length(object@w)))) {
        stop("m2 argument must have dimensions: no. species (",
             nrow(object@species_params), ") x no. size bins (",
             length(object@w), ")")
    }
    return(m2 + object@mu_b + getFMort(object, effort = effort))
}

#' Alias for getMort
#' 
#' An alias provided for backward compatibility with mizer version <= 1.0
#' @inherit getMort
#' @export
getZ <- getMort


#' Get energy after metabolism and movement
#'
#' Calculates the energy rate available by species and size for reproduction and
#' growth after metabolism and movement have been accounted for: \eqn{E_{r.i}(w)}.
#' Used by the \code{project} method for performing simulations.
#' @param object A \code{MizerParams} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the plankton abundance by size.
#' @param feeding_level The current feeding level (optional). A matrix of size
#'   no. species x no. size bins. If not supplied, is calculated internally
#'   using the \code{getFeedingLevel()} method.
#'
#' @return A two dimensional array (species x size) 
#' @export
#' @seealso \code{\link{project}} and \code{\link{getFeedingLevel}}.
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the energy at a particular time step
#' getEReproAndGrowth(params,sim@@n[21,,],sim@@n_pp[21,])
#' }
getEReproAndGrowth <- function(object, n, n_pp,
                               feeding_level = getFeedingLevel(object, n = n,
                                                               n_pp = n_pp)) {
    if (!all(dim(feeding_level) == c(nrow(object@species_params), length(object@w)))) {
        stop("feeding_level argument must have dimensions: no. species (",
             nrow(object@species_params), ") x no. size bins (",
             length(object@w), ")")
    }
    # assimilated intake
    e <- sweep(feeding_level * object@intake_max, 1,
               object@species_params$alpha, "*", check.margin = FALSE)
    # Subtract metabolism
    e <- e - object@metab
    e[e < 0] <- 0 # Do not allow negative growth
    return(e)
}


#' Get energy rate available for reproduction
#'
#' Calculates the energy rate available by species and size for reproduction
#' after metabolism and movement have been accounted for:
#' \eqn{\psi_i(w)E_{r.i}(w)}. Used by the \code{project} method for performing
#' simulations.
#' @param object A \code{MizerParams} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the plankton abundance by size.
#' @param e The energy available for reproduction and growth (optional). A
#'   matrix of size no. species x no. size bins. If not supplied, is calculated
#'   internally using the \code{getEReproAndGrowth()} method.
#'
#' @return A two dimensional array (prey species x prey size) 
#' @export
#' @seealso \code{\link{project}} and \code{\link{getEReproAndGrowth}}.
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the energy at a particular time step
#' getERepro(params,sim@@n[21,,],sim@@n_pp[21,])
#' }
getERepro <- function(object, n, n_pp,
                         e = getEReproAndGrowth(object, n = n, n_pp = n_pp)) {
    if (!all(dim(e) == c(nrow(object@species_params), length(object@w)))) {
        stop("e argument must have dimensions: no. species (",
             nrow(object@species_params), ") x no. size bins (",
             length(object@w), ")")
    }
    e_repro <- object@psi * e
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
#' Calculates the energy rate \eqn{g_i(w)} available by species and size for
#' growth after metabolism, movement and reproduction have been accounted for.
#' Used by the \code{\link{project}} method for performing simulations.
#' @param object A \linkS4class{MizerParams} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the plankton abundance by size.
#' @param e The energy available for reproduction and growth (optional, although
#'   if specified, e_repro must also be specified). A matrix of size no.
#'   species x no. size bins. If not supplied, is calculated internally using
#'   the \code{\link{getEReproAndGrowth}} method.
#' @param e_repro The energy available for reproduction (optional, although if
#'   specified, e must also be specified). A matrix of size no. species x no.
#'   size bins. If not supplied, is calculated internally using the
#'   \code{\link{getERepro}} method.
#'   
#' @return A two dimensional array (prey species x prey size) 
#' @export
#' @seealso \code{\link{project}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the energy at a particular time step
#' getEGrowth(params,sim@@n[21,,],sim@@n_pp[21,])
#' }
getEGrowth <- function(object, n, n_pp,
                       e_repro = getERepro(object, n = n, n_pp = n_pp),
                       e=getEReproAndGrowth(object, n = n, n_pp = n_pp)) {
    if (!all(dim(e_repro) == c(nrow(object@species_params), length(object@w)))) {
        stop("e_repro argument must have dimensions: no. species (",
             nrow(object@species_params), ") x no. size bins (",
             length(object@w), ")")
    }
    if (!all(dim(e) == c(nrow(object@species_params), length(object@w)))) {
        stop("e argument must have dimensions: no. species (",
             nrow(object@species_params), ") x no. size bins (",
             length(object@w), ")")
    }
    # Assimilated intake less activity and metabolism
    # energy for growth is intake - energy for reproduction
    e_growth <- e - e_repro
    return(e_growth)
}


#' Get density independent recruitment
#'
#' Calculates the density independent recruitment (total egg production)
#' \eqn{R_{p.i}} before density dependence, by species. Used by the
#' \code{project} method for performing simulations.
#' @param object A \code{MizerParams} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the plankton abundance by size.
#' @param e_repro The energy available for reproduction (optional). A matrix of
#'   size no. species x no. size bins. If not supplied, is calculated internally
#'   using the \code{\link{getERepro}} method.
#' @param sex_ratio Proportion of the population that is female. Default value
#'   is 0.5.
#'   
#' @return A numeric vector the length of the number of species 
#' @export
#' @seealso \code{\link{project}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the recruitment at a particular time step
#' getRDI(params,sim@@n[21,,],sim@@n_pp[21,])
#' }
getRDI <- function(object, n, n_pp,
                   e_repro = getERepro(object, n = n, n_pp = n_pp),
                   sex_ratio = 0.5) {
    if (!all(dim(e_repro) == c(nrow(object@species_params), length(object@w)))) {
        stop("e_repro argument must have dimensions: no. species (",
             nrow(object@species_params), ") x no. size bins (",
             length(object@w), ")")
    }
    e_repro_pop <- drop( (e_repro * n) %*% object@dw)
    rdi <- sex_ratio * (e_repro_pop * object@species_params$erepro) /
        object@w[object@w_min_idx]
    return(rdi)
}


#' Get density dependent recruitment
#'
#' Calculates the density dependent recruitment (total egg production) \eqn{R_i}
#' for each species. This is the flux entering the smallest size class of each
#' species. The density dependent recruitment is the density independent
#' recruitment after it has been put through the density dependent
#' stock-recruitment relationship function. This method is used by the
#' \code{project} method for performing simulations.
#' @param object An \code{MizerParams} object
#' @param n A matrix of species abundance (species x size)
#' @param n_pp A vector of the plankton abundance by size
#' @param rdi A vector of density independent recruitment for each species. 
#'   If not specified rdi is calculated internally using
#'   the \code{\link{getRDI}} method.
#'   
#' @return A numeric vector the length of the number of species. 
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the energy at a particular time step
#' getRDD(params,sim@@n[21,,],sim@@n_pp[21,])
#' }
getRDD <- function(object, n, n_pp, sex_ratio = 0.5,
                   rdi = getRDI(object, n = n, n_pp = n_pp, sex_ratio = sex_ratio)) {
    rdd <- object@srr(rdi = rdi, species_params = object@species_params)
    return(rdd)
}

# get_time_elements
# internal function to get the array element references of the time dimension
# for the time based slots of a MizerSim object
# time_range can be character or numeric
# Necessary to include a slot_name argument because the effort and abundance
# slots have different time dimensions
get_time_elements <- function(sim, time_range, slot_name = "n"){
    if (!is(sim, "MizerSim"))
        stop("First argument to get_time_elements function must be of class MizerSim")
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
