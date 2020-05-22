# Copyright 2012 Finlay Scott and Julia Blanchard.
# Copyright 2018 Gustav Delius and Richard Southwell.
# Development has received funding from the European Commission's Horizon 2020 
# Research and Innovation Programme under Grant Agreement No. 634495 
# for the project MINOUW (http://minouw-project.eu/).
# Distributed under the GPL 3 or later 
# Maintainer: Gustav Delius, University of York, <gustav.delius@york.ac.uk>

#' Get encounter rate
#' 
#' Returns the rate at which a predator of species \eqn{i} and
#' weight \eqn{w} encounters food (grams/year).
#' By default this function uses the function [mizerEncounter()] to calculate
#' the encounter rate, but this can be overruled via [setRateFunction()].
#' 
#' @inherit mizerEncounter
#' 
#' @export
#' @family rate functions
#' @examples
#' \dontrun{
#' params <- newMultispeciesParams(NS_species_params_gears, inter)
#' # Run simulation with constant fishing effort for all gears for 20 years
#' sim <- project(params, t_max = 20, effort = 0.5)
#' n <- sim@@n[21, , ]
#' n_pp <- sim@@n_pp[21, ]
#' getEncounter(params, n, n_pp)
#' }
getEncounter <- function(params, n = initialN(params), 
                         n_pp = initialNResource(params),
                         n_other = initialNOther(params)) {

    f <- get(params@rates_funcs$Encounter)
    f(params, n = n, n_pp = n_pp, n_other = n_other)
}


#' Get feeding level
#'
#' Returns the feeding level.
#' By default this function uses the function [mizerFeedingLevel()] to calculate
#' the feeding level, but this can be overruled via [setRateFunction()].
#' 
#' @inherit mizerFeedingLevel
#' 
#' @param object A `MizerParams` object or a `MizerSim` object
#' @inheritParams mizerRates
#' @inheritParams get_time_elements
#' @param encounter The encounter rate matrix (optional) of dimension no.
#'   species x no. size bins. If not passed in, it is calculated internally
#'   using [getEncounter()]. Only used if `object`
#'   argument is of type `MizerParams`.
#' @param drop If `TRUE` then any dimension of length 1 will be removed
#'   from the returned array.
#'
#' @return If a `MizerParams` object is passed in, the function returns a two
#'   dimensional array (predator species x predator size) based on the
#'   abundances also passed in.
#'   If a `MizerSim` object is passed in, the function returns a three
#'   dimensional array (time step x predator species x predator size) with the
#'   feeding level calculated at every time step in the simulation.
#'   If \code{drop = TRUE} then the dimension of length 1 will be removed from
#'   the returned array.
#' 
#' @export
#' @family rate functions
#' @examples
#' \dontrun{
#' params <- newMultispeciesParams(NS_species_params_gears, inter)
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
getFeedingLevel <- function(object, n, n_pp, n_other, encounter,
                            time_range, drop = FALSE) {
    if (is(object, "MizerParams")) {
        params <- object
        if (missing(encounter)) {
            if (missing(n)) n <- params@initial_n
            if (missing(n_pp)) n_pp <- params@initial_n_pp
            if (missing(n_other)) n_other <- params@initial_n_other
            encounter <- getEncounter(params, n, n_pp, n_other)
        }
        # Check dims of encounter
        if (!all(dim(encounter) == c(nrow(params@species_params),
                                  length(params@w)))) {
            stop("encounter argument must have dimensions: no. species (",
                 nrow(params@species_params), ") x no. size bins (",
                 length(params@w), ")")
        }
        # calculate feeding level
        f <- get(params@rates_funcs$FeedingLevel)
        return(f(params, encounter = encounter))
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
            n_other <- sim@n_other[x, ]
            names(n_other) <- dimnames(sim@n_other)$component
            feed <- getFeedingLevel(sim@params, n = n,
                                    n_pp = sim@n_pp[x, ],
                                    n_other = n_other)
            return(feed)
            }, .drop = drop)
        return(feed_time)
    }
}


#' Get critical feeding level
#' 
#' The critical feeding level is the feeding level at which the food intake is
#' just high enough to cover the metabolic costs, with nothing left over for
#' growth or reproduction. 
#' 
#' @param params A MizerParams object
#' @return A matrix (species x size) with the critical feeding level
#' @export
getCriticalFeedingLevel <- function(params) {
  validObject(params)
  params@metab/params@intake_max/params@species_params$alpha
}

#' Get predation rate
#' 
#' Calculates the potential rate (in units 1/year) at which a prey individual of
#' a given size \eqn{w} is killed by predators from species \eqn{j}. In formulas
#' \deqn{{\tt pred\_rate}_j(w_p) = \int \phi_j(w,w_p) (1-f_j(w)) 
#'   \gamma_j(w) N_j(w) \, dw.}{pred_rate_j(w_p) = \int\phi_i(w,w_p) (1-f_i(w)) 
#'   \gamma_i(w) N_i(w) dw.}
#' This potential rate is used in the function [getPredMort()] to
#' calculate the realised predation mortality rate on the prey individual.
#'
#' @inheritParams mizerRates
#' @param feeding_level The current feeding level (optional). A matrix of size
#'   no. species x no. size bins. If not supplied, is calculated internally
#'   using the [getFeedingLevel()] function.
#'   
#' @return A two dimensional array (predator species x prey size), 
#'   where the prey size runs over fish community plus resource spectrum.
#' @export
#' @family rate functions
#' @examples
#' \dontrun{
#' params <- newMultispeciesParams(NS_species_params_gears, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the feeding level at one time step
#' n <- sim@@n[21,,]
#' n_pp <- sim@@n_pp[21,]
#' getPredRate(params,n,n_pp)
#' }
getPredRate <- function(params, n = initialN(params), 
                        n_pp = initialNResource(params),
                        n_other = initialNOther(params),
                        feeding_level = getFeedingLevel(params, n = n,
                                                        n_pp = n_pp, 
                                                        n_other = n_other)
                        ) {
    no_sp <- dim(params@interaction)[1]
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)
    if (!all(dim(feeding_level) == c(no_sp, no_w))) {
        stop("feeding_level argument must have dimensions: no. species (",
             no_sp, ") x no. size bins (", no_w, ")")
    }
    
    f <- get(params@rates_funcs$PredRate)
    pred_rate <- f(params, n = n, n_pp = n_pp, n_other = n_other, 
                   feeding_level = feeding_level)
    dimnames(pred_rate) <- list(sp = dimnames(params@initial_n)$sp,
                                w_prey = as.character(signif(params@w_full, 3)))
    pred_rate
}

#' Get total predation mortality rate
#'
#' Calculates the total predation mortality rate \eqn{\mu_{p,i}(w_p)} (in units
#' of 1/year) on each prey species by prey size:
#' \deqn{\mu_{p.i}(w_p) = \sum_j {\tt pred\_rate}_j(w_p)\, \theta_{ji}.}{
#'   \mu_{p.i}(w_p) = \sum_j pred_rate_j(w_p) \theta_{ji}.}
#' 
#' @inheritParams getFeedingLevel
#' @param pred_rate An array of predation rates of dimension no. sp x no.
#'   community size bins x no. of size bins in whole spectra (i.e. community +
#'   resource, the w_full slot). The array is optional. If it is not provided
#'   it is calculated by the [getPredRate()] function.
#'
#' @return
#'   If a `MizerParams` object is passed in, the function returns a two
#'   dimensional array (prey species x prey size) based on the abundances also
#'   passed in. If a `MizerSim` object is passed in, the function returns a
#'   three dimensional array (time step x prey species x prey size) with the
#'   predation mortality calculated at every time step in the simulation.
#'   Dimensions may be dropped if they have length 1 unless `drop = FALSE`.
#' @family rate functions
#' @export
#' @examples
#' \dontrun{
#' params <- newMultispeciesParams(NS_species_params_gears, inter)
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
getPredMort <- function(object, n, n_pp, n_other,
                        pred_rate, time_range, drop = TRUE) {
    if (is(object, "MizerParams")) {
        params <- object
        if (missing(n)) n <- params@initial_n
        if (missing(n_pp)) n_pp <- params@initial_n_pp
        if (missing(n_other)) n_other <- params@initial_n_other
        if (missing(pred_rate)) {
            feeding_level <- getFeedingLevel(params, n = n, n_pp = n_pp, 
                                             n_other = n_other)

            pred_rate <- getPredRate(params, n = n, n_pp = n_pp, 
                                     n_other = n_other,
                                     feeding_level = feeding_level)
        }
        
        f <- get(params@rates_funcs$PredMort)
        pred_mort <- f(params, n = n, n_pp = n_pp, n_other = n_other, 
                       pred_rate = pred_rate)
        dimnames(pred_mort) <- list(prey = dimnames(params@initial_n)$sp,
                                    w_prey = dimnames(params@initial_n)$w)
        pred_mort
    } else {
        sim <- object
        if (missing(time_range)) {
            time_range <- dimnames(sim@n)$time
        }
        time_elements <- get_time_elements(sim, time_range)
        pred_mort_time <- plyr::aaply(which(time_elements), 1, function(x) {
            n <- array(sim@n[x, , ], dim = dim(sim@n)[2:3])
            dimnames(n) <- dimnames(sim@n)[2:3]
            n_other <- sim@n_other[x, ]
            names(n_other) <- dimnames(sim@n_other)$component
            n_pp <- sim@n_pp[x, ]
            return(getPredMort(sim@params, n = n, 
                               n_pp = n_pp, n_other = n_other))
        }, .drop = drop)
        return(pred_mort_time)
    }
}

#' Alias for getPredMort
#' 
#' An alias provided for backward compatibility with mizer version <= 1.0
#' @inherit getPredMort
#' @export
getM2 <- getPredMort


#' Get predation mortality rate for resource
#' 
#' Calculates the predation mortality rate \eqn{\mu_p(w)} on the resource
#' spectrum by resource size (in units 1/year).
#' 
#' Used by the `project` function for running size based simulations.
#' 
#' @inheritParams mizerRates
#' @param pred_rate An array of predation rates of dimension no. sp x no.
#'   community size bins x no. of size bins in whole spectra (i.e. community +
#'   resource, the w_full slot). The array is optional. If it is not provided
#'   it is calculated by the [getPredRate()] function.
#'
#' @return A vector of mortality rate by resource size.
#' @family rate functions
#' @export
#' @examples
#' \dontrun{
#' params <- newMultispeciesParams(NS_species_params_gears, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get resource mortality at one time step
#' n <- sim@@n[21,,]
#' n_pp <- sim@@n_pp[21,]
#' getResourceMort(params,n,n_pp)
#' }
getResourceMort <- 
    function(params, n = initialN(params), 
             n_pp = initialNResource(params),
             n_other = initialNOther(params),
             pred_rate = getPredRate(params, n = n, n_pp = n_pp, 
                                     n_other = n_other)) {

    if ( (!all(dim(pred_rate) ==
               c(nrow(params@species_params), length(params@w_full)))) |
         (length(dim(pred_rate)) != 2)) {
        stop("pred_rate argument must have 2 dimensions: no. species (",
             nrow(params@species_params),
             ") x no. size bins in community + resource (",
             length(params@w_full), ")")
    }
    
    f <- get(params@rates_funcs$ResourceMort)
    f(params, n = n, n_pp = n_pp, n_other = n_other, pred_rate = pred_rate)
}

#' Alias for getResourceMort
#' 
#' An alias provided for backward compatibility with mizer version <= 1.0
#' @inherit getResourceMort
#' @export
getM2Background <- getResourceMort


#' Get the fishing mortality by time, gear, species and size
#'
#' Calculates the fishing mortality rate \eqn{F_{g,i,w}} by gear, species and
#' size at each time step in the `effort` argument (in units 1/year).
#' Used by the `project` function to perform simulations.
#' 
#' @param object A `MizerParams` object or a `MizerSim` object.
#' @param effort The effort for each fishing gear. See notes below.
#' @param time_range Subset the returned fishing mortalities by time. The time
#'   range is either a vector of values, a vector of min and max time, or a
#'   single value. Default is the whole time range. Only used if the
#'   `object` argument is of type `MizerSim`.
#'   
#' @return An array. If the effort argument has a time dimension, or a
#'   `MizerSim` is passed in, the output array has four dimensions (time x
#'   gear x species x size). If the effort argument does not have a time
#'   dimension (i.e. it is a vector or a single numeric), the output array has
#'   three dimensions (gear x species x size).
#' @note Here: fishing mortality = catchability x selectivity x effort.
#' 
#' The `effort` argument is only used if a `MizerParams` object is
#' passed in. The `effort` argument can be a two dimensional array (time x
#' gear), a vector of length equal to the number of gears (each gear has a
#' different effort that is constant in time), or a single numeric value (each
#' gear has the same effort that is constant in time). The order of gears in the
#' `effort` argument must be the same the same as in the `MizerParams`
#' object. If the `effort` argument is not supplied, its value is taken
#' from the `@initial_effort` slot in the params object.
#' 
#' If the object argument is of class `MizerSim` then the effort slot of
#' the `MizerSim` object is used and the `effort` argument is not
#' used.
#' @export
#' @family rate functions
#' @examples
#' \dontrun{
#' params <- newMultispeciesParams(NS_species_params_gears, inter)
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
        sim <- object
        if (missing(time_range)) {
            time_range <- dimnames(sim@effort)$time
        }
        time_elements <- get_time_elements(sim, time_range, slot_name = "effort")
        f_mort_gear <- getFMortGear(sim@params, sim@effort)
        return(f_mort_gear[time_elements, , , , drop = FALSE])
    } else {
        params <- object
        if (missing(effort)) {
            effort <- params@initial_effort
        }
        if (is(effort, "numeric")) {
            no_gear <- dim(params@catchability)[1]
            # If a single value, just repeat it for all gears
            if (length(effort) == 1) {
                effort <- rep(effort, no_gear)
            }
            if (length(effort) != no_gear) {
                stop("Effort must be a single value or a vector as long as the number of gears\n")
            }
            
            mizerFMortGear(params, effort = effort)
        
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
            out[] <- apply(effort, 1, function(x) mizerFMortGear(params, x))
            out <- aperm(out, c(4, 1, 2, 3))
            return(out)
        }
    }
}


#' Get the total fishing mortality rate from all fishing gears by time, species
#' and size.
#' 
#' Calculates the total fishing mortality  (in units 1/year) from all gears by
#' species and size at each time step in the `effort` argument.
#' The total fishing mortality is just the sum of the fishing mortalities
#' imposed by each gear, \eqn{\mu_{f.i}(w)=\sum_g F_{g,i,w}}.
#' 
#' @param object A `MizerParams` object or a `MizerSim` object
#' @param effort The effort of each fishing gear. Only used if the object
#'   argument is of class `MizerParams`. See notes below.
#' @param time_range Subset the returned fishing mortalities by time. The time
#'   range is either a vector of values, a vector of min and max time, or a
#'   single value. Default is the whole time range. Only used if the
#'   `object` argument is of type `MizerSim`.
#' @param drop Only used when object is of type `MizerSim`. Should
#'   dimensions of length 1 be dropped, e.g. if your community only has one
#'   species it might make presentation of results easier. Default is TRUE
#'
#' @return An array. If the effort argument has a time dimension, or object is
#'   of class `MizerSim`, the output array has three dimensions (time x
#'   species x size). If the effort argument does not have a time dimension, the
#'   output array has two dimensions (species x size).
#' @note Here: fishing mortality = catchability x selectivity x effort.
#'
#' The `effort` argument is only used if a `MizerParams` object is
#' passed in. The `effort` argument can be a two dimensional array (time x
#' gear), a vector of length equal to the number of gears (each gear has a
#' different effort that is constant in time), or a single numeric value (each
#' gear has the same effort that is constant in time). The order of gears in the
#' `effort` argument must be the same the same as in the `MizerParams`
#' object.
#'
#' If the object argument is of class `MizerSim` then the effort slot of
#' the `MizerSim` object is used and the `effort` argument is not
#' used.
#' @export
#' @family rate functions
#' @examples
#' \dontrun{
#' params <- newMultispeciesParams(NS_species_params_gears, inter)
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
getFMort <- function(object, effort, time_range, drop = TRUE){
    if (is(object, "MizerParams")) {
        params <- object
        if (missing(effort)) {
          effort <- params@initial_effort
        }
        no_gears <- dim(params@catchability)[[1]]
        f <- get(params@rates_funcs$FMort)
        if (length(dim(effort)) == 2) {
            f_mort <- t(apply(effort, 1, function(x) f(params, x)))
            dim(f_mort) <- c(dim(effort)[[1]], dim(params@initial_n))
            dimnames(f_mort) <- c(list(time = dimnames(effort)[[1]]),
                              dimnames(params@initial_n))
            return(f_mort)
        } else if (length(effort) == 1) {
            return(f(params, rep(effort, no_gears)))
        } else if (length(effort) == no_gears) {
            return(f(params, effort))
        } else {
            stop("Invalid effort argument")
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
#' @inheritParams mizerRates
#' @param effort A numeric vector of the effort by gear or a single numeric
#'   effort value which is used for all gears.
#' @param f_mort A two dimensional array (species x size) of fishing mortality
#'   (optional).If not supplied is calculated using the [getFMort()]
#'   function.
#' @param pred_mort A two dimensional array of predation mortality (optional). Has
#'   dimensions no. sp x no. size bins in the community. If not supplied is
#'   calculated using the [getPredMort()] function.
#'
#' @return A two dimensional array (prey species x prey size). 
#'
#' @export
#' @seealso [getPredMort()], [getFMort()]
#' @family rate functions
#' @examples
#' \dontrun{
#' params <- newMultispeciesParams(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the total mortality at a particular time step
#' getMort(params,sim@@n[21,,],sim@@n_pp[21,],effort=0.5)
#' }
getMort <- function(params, 
                    n = initialN(params), 
                    n_pp = initialNResource(params),
                    n_other = initialNOther(params),
                    effort = getInitialEffort(params),
                    f_mort = getFMort(params, effort),
                    pred_mort = getPredMort(params, n = n, n_pp = n_pp, 
                                     n_other = n_other)){
    if (!all(dim(f_mort) == c(nrow(params@species_params), length(params@w)))) {
        stop("f_mort argument must have dimensions: no. species (",
             nrow(params@species_params), ") x no. size bins (",
             length(params@w), ")")
    }
    if (!all(dim(pred_mort) == c(nrow(params@species_params), length(params@w)))) {
        stop("pred_mort argument must have dimensions: no. species (",
             nrow(params@species_params), ") x no. size bins (",
             length(params@w), ")")
    }
  
    f <- get(params@rates_funcs$Mort)
    z <- f(params, n = n, n_pp = n_pp, n_other = n_other, 
           f_mort = f_mort, pred_mort = pred_mort)
    dimnames(z) <- list(prey = dimnames(params@initial_n)$sp,
                        w_prey = dimnames(params@initial_n)$w)
    return(z)
}

#' Alias for getMort
#' 
#' An alias provided for backward compatibility with mizer version <= 1.0
#' @inherit getMort
#' @export
getZ <- getMort


#' Get energy rate available for reproduction and growth
#'
#' Calculates the energy rate \eqn{E_{r.i}(w)} (grams/year) available for
#' reproduction and growth after metabolism and movement have been accounted
#' for. 
#' 
#' @inheritParams mizerRates
#' @param encounter The encounter rate matrix (optional) of dimension no.
#'   species x no. size bins. If not passed in, it is calculated internally
#'   using the [getEncounter()] function.
#' @param feeding_level The current feeding level (optional). A matrix of size
#'   no. species x no. size bins. If not supplied, is calculated internally
#'   using the [getFeedingLevel()] function.
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
#' @export
#' @seealso The parts of this energy rate that is invested into growth is
#'   calculated with [getERepro()] and the part that is invested into
#'   reproduction is calculated with [getEGrowth()].
#' @family rate functions
#' @examples
#' \dontrun{
#' params <- newMultispeciesParams(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the energy at a particular time step
#' getEReproAndGrowth(params,sim@@n[21,,],sim@@n_pp[21,])
#' }
getEReproAndGrowth <- function(params, n = initialN(params), 
                               n_pp = initialNResource(params),
                               n_other = initialNOther(params),
                               encounter = getEncounter(params, n = n,
                                                        n_pp = n_pp, 
                                                        n_other = n_other),
                               feeding_level = getFeedingLevel(params, n = n,
                                                               n_pp = n_pp, 
                                                               n_other = n_other,
                                                               encounter = encounter)) {
    if (!all(dim(encounter) == c(nrow(params@species_params), length(params@w)))) {
        stop("encounter argument must have dimensions: no. species (",
             nrow(params@species_params), ") x no. size bins (",
             length(params@w), ")")
    }
    if (!all(dim(feeding_level) == c(nrow(params@species_params), length(params@w)))) {
        stop("feeding_level argument must have dimensions: no. species (",
             nrow(params@species_params), ") x no. size bins (",
             length(params@w), ")")
    }
    f <- get(params@rates_funcs$EReproAndGrowth)
    f(params, n = n, n_pp = n_pp, n_other = n_other, 
      encounter = encounter, feeding_level = feeding_level)
}


#' Get energy rate available for reproduction
#'
#' Calculates the energy rate (grams/year) available for reproduction after
#' growth and metabolism have been accounted for.
#' 
#' @inheritParams mizerRates
#' @param e The energy available for reproduction and growth (optional). A
#'   matrix of size no. species x no. size bins. If not supplied, is calculated
#'   internally using [getEReproAndGrowth()].
#'
#' @return A two dimensional array (prey species x prey size) holding
#' \deqn{\psi_i(w)E_{r.i}(w)}
#' where \eqn{E_{r.i}(w)} is the rate at which energy becomes available for
#' growth and reproduction, calculated with [getEReproAndGrowth()],
#' and \eqn{\psi_i(w)} is the proportion of this energy that is used for
#' reproduction. This proportion is taken from the `params` object and is
#' set with [setReproduction()].
#' @export
#' @family rate functions
#' @examples
#' \dontrun{
#' params <- newMultispeciesParams(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the energy at a particular time step
#' getERepro(params,sim@@n[21,,],sim@@n_pp[21,])
#' }
getERepro <- function(params, n = initialN(params), 
                      n_pp = initialNResource(params),
                      n_other = initialNOther(params),
                      e = getEReproAndGrowth(params, n = n, n_pp = n_pp,
                                             n_other = n_other)) {
    if (!all(dim(e) == c(nrow(params@species_params), length(params@w)))) {
        stop("e argument must have dimensions: no. species (",
             nrow(params@species_params), ") x no. size bins (",
             length(params@w), ")")
    }
    f <- get(params@rates_funcs$ERepro)
    f(params, n = n, n_pp = n_pp, n_other = n_other, e = e)
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
#' accounted for. Used by [project()] for performing simulations.
#' 
#' @inheritParams mizerRates
#' @param e The energy available for reproduction and growth (optional, although
#'   if specified, e_repro must also be specified). A matrix of size no.
#'   species x no. size bins. If not supplied, is calculated internally using
#'   [getEReproAndGrowth()].
#' @param e_repro The energy available for reproduction (optional, although if
#'   specified, e must also be specified). A matrix of size no. species x no.
#'   size bins. If not supplied, is calculated internally using
#'   [getERepro()].
#'   
#' @return A two dimensional array (prey species x prey size) 
#' @export
#' @seealso [getERepro()], [getEReproAndGrowth()]
#' @family rate functions
#' @examples
#' \dontrun{
#' params <- newMultispeciesParams(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the energy at a particular time step
#' getEGrowth(params,sim@@n[21,,],sim@@n_pp[21,])
#' }
getEGrowth <- function(params, n = initialN(params), 
                       n_pp = initialNResource(params),
                       n_other = initialNOther(params),
                       e_repro = getERepro(params, n = n, n_pp = n_pp, 
                                           n_other = n_other),
                       e = getEReproAndGrowth(params, n = n, n_pp = n_pp, 
                                              n_other = n_other)) {
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
    f <- get(params@rates_funcs$EGrowth)
    f(params, n = n, n_pp = n_pp, n_other = n_other, e_repro = e_repro, e = e)
}


#' Get density independent rate of egg production
#'
#' Calculates the density independent rate of egg production \eqn{R_{p.i}}
#' (units 1/year) before density dependence, by species. Used by
#' [getRDD()] to calculate the actual density dependent rate.
#' See [setReproduction()] for more details.
#' 
#' @inheritParams mizerRates
#' @param e_repro The energy available for reproduction (optional). A matrix of
#'   size no. species x no. size bins. If not supplied, is calculated internally
#'   using [getERepro()].
#'   
#' @return A numeric vector the length of the number of species 
#' @export
#' @seealso [getRDD()]
#' @family rate functions
#' @examples
#' \dontrun{
#' params <- newMultispeciesParams(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the density-independent reproduction rate at a particular time step
#' getRDI(params,sim@@n[21,,],sim@@n_pp[21,])
#' }
getRDI <- function(params, n = initialN(params), 
                   n_pp = initialNResource(params),
                   n_other = initialNOther(params),
                   e_repro = getERepro(params, n = n, n_pp = n_pp, 
                                       n_other = n_other)) {
    if (!all(dim(e_repro) == c(nrow(params@species_params), length(params@w)))) {
        stop("e_repro argument must have dimensions: no. species (",
             nrow(params@species_params), ") x no. size bins (",
             length(params@w), ")")
    }
    f <- get(params@rates_funcs$RDI)
    f(params, n = n, n_pp = n_pp, n_other = n_other, e_repro = e_repro)
}


#' Get density dependent reproduction rate
#'
#' Calculates the density dependent rate of egg production \eqn{R_i} (units
#' 1/year) for each species. This is the flux entering the smallest size class
#' of each species. The density dependent rate is the density independent
#' rate obtained with [getRDI()] after it has been put through the 
#' density dependence function. This is the Beverton-Holt function
#' [BevertonHoltRDD()] by default, but this can be changed. See
#' [setReproduction()] for more details.
#' 
#' @inheritParams mizerRates
#' @param rdi A vector of density-independent reproduction rates for each
#'   species. If not specified, rdi is calculated internally using
#'   [getRDI()].
#'   
#' @return A numeric vector the length of the number of species. 
#' @export
#' @seealso [getRDI()]
#' @family rate functions
#' @examples
#' \dontrun{
#' params <- newMultispeciesParams(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the rate at a particular time step
#' getRDD(params,sim@@n[21,,],sim@@n_pp[21,])
#' }
getRDD <- function(params, n = initialN(params), 
                   n_pp = initialNResource(params),
                   n_other = initialNOther(params),
                   rdi = getRDI(params, n = n, n_pp = n_pp, 
                                n_other = n_other)) {
    # Avoid getting into infinite loops
    if (params@rates_funcs$RDD == "getRDD") {
        stop('"getRDD" is not a valid name for the function giving the density',
             'dependent reproductive rate.')
    }
    f <- get(params@rates_funcs$RDD)
    f(rdi = rdi, species_params = params@species_params)
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
