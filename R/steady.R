#' Measure distance between current and previous state in terms of RDI
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' This function can be used in [projectToSteady()] to decide when sufficient
#' convergence to steady state has been achieved.
#'
#' @param params MizerParams
#' @param current A named list with entries `n`, `n_pp` and `n_other`
#'   describing the current state
#' @param previous A named list with entries `n`, `n_pp` and `n_other`
#'   describing the previous state
#' @return The largest absolute relative change in rdi:
#'   `max(abs((current_rdi - previous_rdi) / previous_rdi))`
#' @family distance functions
#' @concept helper
#' @export
distanceMaxRelRDI <- function(params, current, previous) {
    UseMethod("distanceMaxRelRDI")
}
#' @export
distanceMaxRelRDI.MizerParams <- function(params, current, previous) {
    current_rdi <- getRDI(params, n = current$n, n_pp = current$n_pp,
                          n_other = current$n_other)
    previous_rdi <- getRDI(params, n = previous$n, n_pp = previous$n_pp,
                           n_other = previous$n_other)
    max(abs((current_rdi - previous_rdi) / previous_rdi))
}

#' Measure distance between current and previous state in terms of fish abundances
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Calculates the sum squared difference between log(N) in current and previous
#' state. This function can be used in [projectToSteady()] to decide when
#' sufficient convergence to steady state has been achieved.
#'
#' @param params MizerParams
#' @param current A named list with entries `n`, `n_pp` and `n_other`
#'   describing the current state
#' @param previous A named list with entries `n`, `n_pp` and `n_other`
#'   describing the previous state
#' @return The sum of squares of the difference in the logs of the (nonzero)
#'   fish abundances n:
#'   `sum((log(current$n) - log(previous$n))^2)`
#' @family distance functions
#' @concept helper
#' @export
distanceSSLogN <- function(params, current, previous) {
    UseMethod("distanceSSLogN")
}
#' @export
distanceSSLogN.MizerParams <- function(params, current, previous) {
    sel <- current$n > 0 & previous$n > 0
    sum((log(current$n[sel]) - log(previous$n[sel]))^2)
}

#' Project to steady state
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Run the full dynamics, as in [project()], but stop once the change has slowed
#' down sufficiently, in the sense that the distance between states at
#' successive time steps is less than `tol`. You determine how the distance is
#' calculated.
#'
#' @inheritParams steady
#' @param effort The fishing effort to be used throughout the simulation.
#'   This must be a vector or list with one named entry per fishing gear.
#' @param distance_func A function that will be called after every `t_per` years
#'   with both the previous and the new state and that should return a number
#'   that in some sense measures the distance between the states. By default
#'   this uses the function [distanceSSLogN()] that you can use as a model for your
#'   own distance function.
#' @param ... Further arguments will be passed on to your distance function.
#' 
#' @return A MizerParams or a MizerSim object
#' @seealso [distanceSSLogN()], [distanceMaxRelRDI()]
#' @export
projectToSteady <- function(params,
                            effort = params@initial_effort,
                            distance_func = distanceSSLogN,
                            t_per = 1.5,
                            t_max = 100,
                            dt = 0.1,
                            tol = 0.1 * t_per,
                            return_sim = FALSE,
                            progress_bar = TRUE, ...) {
    UseMethod("projectToSteady")
}
#' @export
projectToSteady.MizerParams <- function(params,
                            effort = params@initial_effort,
                            distance_func = distanceSSLogN,
                            t_per = 1.5,
                            t_max = 100,
                            dt = 0.1,
                            tol = 0.1 * t_per,
                            return_sim = FALSE,
                            progress_bar = TRUE, ...) {
    params <- validParams(params)
    effort <- validEffortVector(effort, params = params)
    params@initial_effort <- effort
    assert_that(t_max >= t_per,
                tol > 0)
    if ((t_per < dt) || !isTRUE(all.equal((t_per - round(t_per / dt) * dt), 0))) {
        stop("t_per must be a positive multiple of dt")
    }
    t_dimnames <-  seq(0, t_max, by = t_per)
    
    if (is(progress_bar, "Progress")) {
        # We have been passed a shiny progress object
        progress_bar$set(message = "Finding steady state", value = 0)
        proginc <- 1 / ceiling(t_max/t_per)
    }
    
    if (return_sim) {
        # create MizerSim object
        sim <- MizerSim(params, t_dimnames =  t_dimnames)
        sim@n[1, , ] <- params@initial_n
        sim@n_pp[1, ] <- params@initial_n_pp
        sim@n_other[1, ] <- params@initial_n_other
        sim@effort[1, ] <- params@initial_effort
    }
    
    # get functions
    resource_dynamics_fn <- get(params@resource_dynamics)
    other_dynamics_fns <- lapply(params@other_dynamics, get)
    rates_fns <- lapply(params@rates_funcs, get)
    r <- rates_fns$Rates(
        params, n = params@initial_n,
        n_pp = params@initial_n_pp,
        n_other = params@initial_n_other,
        t = 0, 
        effort = effort, rates_fns = rates_fns, ...)
    
    previous <- list(n = params@initial_n,
                     n_pp = params@initial_n_pp,
                     n_other = params@initial_n_other,
                     rates = r)
    
    for (i in 2:length(t_dimnames)) {
        # advance shiny progress bar
        if (is(progress_bar, "Progress")) {
            progress_bar$inc(amount = proginc)
        }
        current <- project_simple(params, n = previous$n, n_pp = previous$n_pp,
                                  n_other = previous$n_other, t = 0,
                                  dt = dt, steps = round(t_per / dt),
                                  effort = params@initial_effort,
                                  resource_dynamics_fn = resource_dynamics_fn,
                                  other_dynamics_fns = other_dynamics_fns,
                                  rates_fns = rates_fns)
        if (return_sim) {
            # Store result
            sim@n[i, , ] <- current$n
            sim@n_pp[i, ] <- current$n_pp
            sim@n_other[i, ] <- unserialize(serialize(current$n_other, NULL))
            sim@effort[i, ] <- params@initial_effort
        }
        
        # Species with no reproduction are going extinct, so stop.
        extinct <- is.na(current$rates$rdd) | current$rates$rdd <= 1e-20
        if (any(extinct)) {
            warning(paste(params@species_params$species[extinct], collapse = ", "),
                    " are going extinct.")
            success <- FALSE
            distance <- NA
            break
        }
        
        distance <- distance_func(params,
                                  current = current,
                                  previous = previous, ...)
        success <- distance < tol
        if (success == TRUE) {
            break
        }
        previous <- current
    }
    if (!success) {
        message("Simulation run did not converge after ",
                (i - 1) * t_per,
                " years. Value returned by the distance function was: ",
                distance)
    } else {
        message("Convergence was achieved in ", (i - 1) * t_per, " years.")
    }
    
    params@initial_n[] <- current$n
    params@initial_n_pp[] <- current$n_pp
    params@initial_n_other[] <- current$n_other
    
    if (return_sim) {
        sim@params <- params
        sel <- 1:i
        sim@n <- sim@n[sel, , , drop = FALSE]
        sim@n_pp <- sim@n_pp[sel, , drop = FALSE]
        sim@n_other <- sim@n_other[sel, , drop = FALSE]
        sim@effort <- sim@effort[sel, , drop = FALSE]
        return(sim)
    } else {
        params@time_modified <- lubridate::now()
        return(params)
    }
}

#' Set initial values to a steady state for the model
#'
#' The steady state is found by running the dynamics while keeping reproduction,
#' resource and other components constant until the size spectra no longer
#' change much (or until time `t_max` is reached, if earlier).
#' 
#' If the model use Beverton-Holt reproduction then the reproduction parameters
#' are set to values that give the level of reproduction observed in that
#' steady state. The `preserve` argument can be used to specify which of the 
#' reproduction parameters should be preserved.
#' 
#' @param params A \linkS4class{MizerParams} object
#' @param t_max The maximum number of years to run the simulation. Default is 100.
#' @param t_per The simulation is broken up into shorter runs of `t_per` years,
#'   after each of which we check for convergence. Default value is 1.5. This
#'   should be chosen as an odd multiple of the timestep `dt` in order to be
#'   able to detect period 2 cycles.
#' @param dt The time step to use in `project()`.
#' @param tol The simulation stops when the relative change in the egg
#'   production RDI over `t_per` years is less than `tol` for every species.
#' @param return_sim If TRUE, the function returns the MizerSim object holding
#'   the result of the simulation run, saved at intervals of `t_per`. If FALSE (default) the function returns
#'   a MizerParams object with the "initial" slots set to the steady state.
#' @param preserve `r lifecycle::badge("experimental")`
#'   Specifies whether the `reproduction_level` should be preserved (default)
#'   or the maximum reproduction rate `R_max` or the reproductive
#'   efficiency `erepro`. See [setBevertonHolt()] for an explanation
#'   of the `reproduction_level`.
#' @param progress_bar A shiny progress object to implement a progress bar in a
#'   shiny app. Default FALSE.
#' @return A MizerParams or a MizerSim object
#' @export
#' @examples
#' \donttest{
#' params <- newTraitParams()
#' species_params(params)$gamma[5] <- 3000
#' params <- steady(params)
#' plotSpectra(params)
#' }
steady <- function(params, t_max = 100, t_per = 1.5, dt = 0.1,
                   tol = 0.1 * dt, return_sim = FALSE, 
                   preserve = c("reproduction_level", "erepro", "R_max"),
                   progress_bar = TRUE) {
    UseMethod("steady")
}

#' @export
steady.MizerParams <- function(params, t_max = 100, t_per = 1.5, dt = 0.1,
                   tol = 0.1 * dt, return_sim = FALSE, 
                   preserve = c("reproduction_level", "erepro", "R_max"),
                   progress_bar = TRUE) {
    
    if (params@rates_funcs$RDD == "BevertonHoltRDD") {
        preserve <- match.arg(preserve)
        old_reproduction_level <- getReproductionLevel(params)
        old_R_max <- params@species_params$R_max
        old_erepro <- params@species_params$erepro
    }
    
    # Force the reproduction to stay at the current level
    params@species_params$constant_reproduction <- getRDD(params)
    old_rdd_fun <- params@rates_funcs$RDD
    params@rates_funcs$RDD <- "constantRDD"
    
    # Force resource to stay at current level
    old_resource_dynamics <- params@resource_dynamics
    params@resource_dynamics <- "resource_constant"
    
    # Force other components to stay at current level
    old_other_dynamics <- params@other_dynamics
    for (res in names(params@other_dynamics)) {
        params@other_dynamics[[res]] <- "constant_other"
    }
    
    object <- projectToSteady(params,
                              distance_func = distanceMaxRelRDI,
                              t_per = t_per,
                              t_max = t_max,
                              dt = dt,
                              tol = tol,
                              return_sim = return_sim,
                              progress_bar = progress_bar)
    if (return_sim) {
        params <- object@params
    } else {
        params <- object
    }
    # Restore original RDD and dynamics
    params@rates_funcs$RDD <- old_rdd_fun
    params@other_dynamics <- old_other_dynamics
    params@species_params$constant_reproduction <- NULL
    
    # Set resource dynamics
    params <- setResource(params, resource_dynamics = old_resource_dynamics)
    
    if (params@rates_funcs$RDD == "BevertonHoltRDD") {
        if (preserve == "reproduction_level") {
            params <- setBevertonHolt(params, 
                                      reproduction_level = old_reproduction_level)
        } else if (preserve == "R_max") {
            params <- setBevertonHolt(params, 
                                      R_max = old_R_max)
        } else {
            params <- setBevertonHolt(params, erepro = old_erepro)
        }
    }
    
    if (return_sim) {
        object@params <- params
        return(object)
    } else {
        params@time_modified <- lubridate::now()
        return(params)
    }
}


#' Helper function to keep other components constant
#' 
#' @param params MizerParams object
#' @param n_other Abundances of other components
#' @param component Name of the component that is being updated
#' @param ... Unused
#' @return The current value of the component
#' @export
#' @concept helper
constant_other <- function(params, n_other, component, ...) {
    n_other[[component]]
}

#' Helper function to assure validity of species argument
#' 
#' If the species argument contains invalid species, then these are
#' ignored but a warning is issued.
#' 
#' @param object A MizerSim or MizerParams object from which the species
#'   should be selected.
#' @param species The species to be selected. Optional. By default all target
#'   species are selected. A vector of species names, or a
#'   numeric vector with the species indices, or a logical vector indicating for
#'   each species whether it is to be selected (TRUE) or not. 
#' @param return.logical Whether the return value should be a logical vector.
#'   Default FALSE.
#' @param error_on_empty Whether to throw an error if there are zero valid
#'   species. Default FALSE.
#'   
#' @return A vector of species names, in the same order as specified in the
#'   'species' argument. If 'return.logical = TRUE' then a logical vector is
#'   returned instead, with length equal to the number of species, with
#'   TRUE entry for each selected species.
#' @export
#' @concept helper
valid_species_arg <- function(object, species = NULL, return.logical = FALSE,
                              error_on_empty = FALSE) {
    if (is(object, "MizerSim")) {
        params <- object@params
    } else if (is(object, "MizerParams")) {
        params <- object
    } else {
        stop("The first argument must be a MizerSim or MizerParams object.")
    }
    assert_that(is.flag(return.logical),
                is.flag(error_on_empty))
    all_species <- dimnames(params@initial_n)$sp
    no_sp <- nrow(params@species_params)
    # Set species if missing to list of all non-background species
    if (is.null(species)) {
        species <- dimnames(params@initial_n)$sp[!is.na(params@A)]
        if (length(species) == 0) {  # There are no non-background species.
            if (error_on_empty) {
                stop("No species have been selected.")
            }
            if (return.logical) {
                return(rep(FALSE, no_sp))
            } else {
                vector("character")
            }
        }
    }
    if (is.logical(species)) {
        if (length(species) != no_sp) {
            stop("The boolean `species` argument has the wrong length.")
        }
        if (!any(species) && error_on_empty) {
            stop("No species have been selected.")
        }
        if (return.logical) {
            return(species)
        }
        return(all_species[species])
    }
    if (is.numeric(species)) {
        if (!all(species %in% (1:no_sp))) {
            warning("A numeric 'species' argument should only contain the ",
                    "integers 1 to ", no_sp, ".")
        }
        species.logical <- 1:no_sp %in% species
        if (!any(species.logical) && error_on_empty) {
            stop("No species have been selected.")
        }
        if (return.logical) {
            return(species.logical)
        }
        return(all_species[species[species %in% (1:no_sp)]])
    }
    invalid <- setdiff(species, all_species)
    if (length(invalid) > 0) {
        warning("The following species do not exist: ", 
                toString(invalid), ".")
    }
    species <- intersect(species, all_species)
    if (length(species) == 0 && error_on_empty) {
        stop("No species have been selected.")
    }
    if (return.logical) {
        return(all_species %in% species)
    }
    species
}

#' Helper function to assure validity of gears argument
#' 
#' If the gears argument contains invalid gears, then these are
#' ignored but a warning is issued.
#' 
#' @param object A MizerSim or MizerParams object from which the gears
#'   should be selected.
#' @param gears The gears to be selected. Optional. By default all gears
#'   are selected. A vector of gear names.
#' @param error_on_empty Whether to throw an error if there are zero valid
#'   gears. Default FALSE.
#'   
#' @return A vector of gear names.
#' @export
#' @concept helper
valid_gears_arg <- function(object, gears = NULL,
                            error_on_empty = FALSE) {
    if (is(object, "MizerSim")) {
        params <- object@params
    } else if (is(object, "MizerParams")) {
        params <- object
    } else {
        stop("The first argument must be a MizerSim or MizerParams object.")
    }
    assert_that(is.flag(error_on_empty))
    all_gears <- unique(params@gear_params$gear)
    no_gear <- length(all_gears)
    # Set gears if missing to list of all gears
    if (is.null(gears)) {
        gears <- all_gears
    }
    invalid <- setdiff(gears, all_gears)
    if (length(invalid) > 0) {
        warning("The following gears do not exist: ", 
                toString(invalid), ".")
    }
    gears <- intersect(gears, all_gears)
    if (length(gears) == 0 && error_on_empty) {
        stop("No gears have been selected.")
    }
    gears
}
