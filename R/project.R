# Project function for the size based modelling package mizer

# Copyright 2012 Finlay Scott and Julia Blanchard.
# Copyright 2018 Gustav Delius and Richard Southwell.
# Development has received funding from the European Commission's Horizon 2020 
# Research and Innovation Programme under Grant Agreement No. 634495 
# for the project MINOUW (http://minouw-project.eu/).
# Distributed under the GPL 3 or later 
# Maintainer: Gustav Delius, University of York, <gustav.delius@york.ac.uk>

#' @useDynLib mizer
#' @importFrom Rcpp sourceCpp
NULL


#' Project size spectrum forward in time
#' 
#' Runs the size spectrum model simulation.
#' The function returns an object of type
#' \linkS4class{MizerSim} that can then be explored with a range of
#' [summary_functions], [indicator_functions] and 
#' [plotting_functions].
#' 
#' @param object Either a \linkS4class{MizerParams} object or a 
#'   \linkS4class{MizerSim} object (which contains a `MizerParams` object).
#' @param effort The effort of each fishing gear through time. See notes below.
#' @param t_max The number of years the projection runs for. The default value is
#'   100. However, this argument is ignored if an array is used for the
#'   `effort` argument. See notes below.
#' @param dt Time step of the solver. The default value is 0.1.
#' @param t_save The frequency with which the output is stored. The default
#'   value is 1. Must be an integer multiple of dt.
#' @param t_start The the year of the start of the simulation. The simulation
#'   will cover the period from `t_start` to \code{t_start + t_max}.
#'   Defaults to 0. Ignored if an array is used for the `effort`
#'   argument or a `MizerSim` for the `object` argument.
#' @param initial_n Deprecated. The initial abundances of species. Instead of
#'   using this argument you should set `initialN(params)` to the desired value.
#' @param initial_n_pp Deprecated. The initial abundances of resource. Instead
#'   of using this argument you should set `initialNResource(params)` to the
#'   desired value.
#' @param append A boolean that determines whether the new simulation results
#'   are appended to the previous ones. Only relevant if `object` is a
#'   `MizerSim` object. Default = TRUE.
#' @param progress_bar Either a boolean value to determine whether a progress
#'   bar should be shown in the console, or a shiny Progress object to implement 
#'   a progress bar in a shiny app.
#' @param ... Other arguments will be passed to rate functions.
#' 
#' @note The `effort` argument specifies the level of fishing effort during
#' the simulation. If it is not supplied, the initial effort stored in the params
#' object is used. The effort can be specified in three different ways: 
#' \itemize{ 
#' \item A single numeric value. This specifies the effort of all fishing gears
#' which is constant through time (i.e. all the gears have the same constant
#' effort).
#' \item A numerical vector which has the same length as the number of fishing
#' gears. The vector must be named and the names must correspond to the gear
#' names in the `MizerParams` object. The values in the vector specify the
#' constant fishing effort of each of the fishing gears, i.e. the effort is
#' constant through time but each gear may have a different fishing effort. 
#' \item A numerical array with dimensions time x gear. This specifies the
#' fishing effort of each gear at each time step.  The first dimension, time,
#' must be named numerically and increasing. The second dimension of the array
#' must be named and the names must correspond to the gear names in the
#' `MizerParams` argument. The value for the effort for a particular time
#' is used during the interval from that time to the next time in the array.
#' }
#' 
#' If effort is specified as an array then the smallest time in the array is 
#' used as the initial time for the simulation. Otherwise the initial time is
#' set to the final time of the previous simulation if `object` is a 
#' `MizerSim` object or to `t_start` otherwise. Also, if the effort is
#' an array then the `t_max` argument is ignored and the maximum simulation
#' time is the largest time of the effort array.
#' 
#' If the `object` argument is of class `MizerSim` then the initial
#' values for the simulation are taken from the final values in the 
#' `MizerSim` object and the corresponding arguments to this function will
#' be ignored.
#' 
#' @return An object of class \linkS4class{MizerSim}.
#' 
#' @export
#' @examples
#' \dontrun{
#' # Data set with different fishing gears
#' params <- newMultispeciesParams(NS_species_params_gears, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # With constant fishing effort which is different for each gear
#' effort <- c(Industrial = 0, Pelagic = 1, Beam = 0.5, Otter = 0.5)
#' sim <- project(params, t_max = 20, effort = effort)
#' # With fishing effort that varies through time for each gear
#' gear_names <- c("Industrial","Pelagic","Beam","Otter")
#' times <- seq(from = 1, to = 10, by = 1)
#' effort_array <- array(NA, dim = c(length(times), length(gear_names)),
#'     dimnames = list(time = times, gear = gear_names))
#' effort_array[,"Industrial"] <- 0.5
#' effort_array[,"Pelagic"] <- seq(from = 1, to = 2, length = length(times))
#' effort_array[,"Beam"] <- seq(from = 1, to = 0, length = length(times))
#' effort_array[,"Otter"] <- seq(from = 1, to = 0.5, length = length(times))
#' sim <- project(params, effort = effort_array)
#' }
project <- function(object, effort,
                    t_max = 100, dt = 0.1, t_save = 1, t_start = 0,
                    initial_n, initial_n_pp,
                    append = TRUE,
                    progress_bar = TRUE, ...) {
    
    # Set and check initial values ----
    if (is(object, "MizerSim")) {
        validObject(object)
        params <- object@params
        no_t <- dim(object@n)[[1]]
        initial_n <- params@initial_n # Needed to get the right dimensions
        initial_n[] <- object@n[no_t, , ]
        initial_n_pp <- params@initial_n_pp # Needed to get the right dimensions
        initial_n_pp[] <- object@n_pp[no_t, ]
        initial_n_other <- object@n_other[no_t, ]
        t_start <- as.numeric(dimnames(object@n)[[1]][[no_t]])
    } else if (is(object, "MizerParams")) {
        params <- validParams(object)
        if (missing(initial_n))       initial_n <- params@initial_n
        if (missing(initial_n_pp)) initial_n_pp <- params@initial_n_pp
        initial_n_other <- params@initial_n_other
    } else {
        stop("The `object` argument must be either a MizerParams or a MizerSim object.")
    }
    params@initial_n[] <- initial_n
    params@initial_n_pp[] <- initial_n_pp
    params@initial_n_other <- initial_n_other
    
    no_sp <- length(params@w_min_idx)
    assert_that(is.array(initial_n),
                is.numeric(initial_n),
                are_equal(dim(initial_n), c(no_sp, length(params@w))))
    assert_that(is.numeric(initial_n_pp),
                length(initial_n_pp) == length(params@w_full))
    
    assert_that(is.null(initial_n_other) || is.list(initial_n_other))
    other_names <- names(params@other_dynamics)
    if (length(other_names) > 0) {
        if (is.null(names(initial_n_other))) {
            stop("The initial_n_other needs to be a named list")
        }
        if (!setequal(names(initial_n_other), other_names)) {
            stop("The names of the entries in initial_n_other do not match ",
                 "the names of the other components of the model.")
        }
    }
    
    # Set effort array ----
    if (missing(effort)) effort <- params@initial_effort
    # Do we need to create an effort array?
    if (is.null(dim(effort))) { # effort is a vector
        no_gears <- dim(params@catchability)[1]
        if ((length(effort) > 1) & (length(effort) != no_gears)) {
            stop("Effort vector must be the same length as the number of fishing gears\n")
        }
        # If more than 1 gear need to check that gear names match
        gear_names <- dimnames(params@catchability)[[1]]
        effort_gear_names <- names(effort)
        if (length(effort) == 1 && is.null(effort_gear_names)) {
            effort_gear_names <- gear_names
        }
        if (!all(gear_names %in% effort_gear_names)) {
            stop("Gear names in the MizerParams object (", 
                 paste(gear_names, collapse = ", "), 
                 ") do not match those in the effort vector.")
        }
        # Set up the effort array transposed so we can use the recycling rules
        time_dimnames <- c(t_start, t_start + t_max)
        effort <- t(array(effort, dim = c(no_gears, 2), 
                          dimnames = list(gear = effort_gear_names, 
                                          time = time_dimnames)))
    }
    
    # Check that number and names of gears in effort array is same as in 
    # MizerParams object
    no_gears <- dim(params@catchability)[1]
    if (dim(effort)[2] != no_gears) {
        stop("The number of gears in the effort array (length of the second dimension = ", 
             dim(effort)[2], 
             ") does not equal the number of gears in the MizerParams object (", 
             no_gears, ").")
    }
    gear_names <- dimnames(params@catchability)[[1]]
    if (!all(gear_names %in% dimnames(effort)[[2]])) {
        stop("Gear names in the MizerParams object (", 
             paste(gear_names, collapse = ", "), 
             ") do not match those in the effort array.")
    }
    # Sort effort array to match order in MizerParams
    effort <- effort[, gear_names, drop = FALSE]
    
    if (is.null(dimnames(effort)[[1]])) {
        stop("The time dimname of the effort argument must be numeric.")
    }
    time_effort <- as.numeric(dimnames(effort)[[1]])
    if (any(is.na(time_effort))) {
        stop("The time dimname of the effort argument must be numeric.")
    }
    if (is.unsorted(time_effort)) {
        stop("The time dimname of the effort argument should be increasing.")
    }
    t_start <- time_effort[[1]]
    t_end <- time_effort[[length(time_effort)]]
    t_max <- t_end - t_start
    
    # Make the MizerSim object with the right size ----
    # We only save every t_save years
    # no point running a simulation with no saved results
    if (t_max < t_save) {
        t_save <- t_max
    }
    # Divisibility test needs to be careful about machine rounding errors,
    # see https://github.com/sizespectrum/mizer/pull/2
    if ((t_save < dt) || !isTRUE(all.equal((t_save - round(t_save / dt) * dt), 0)))
        stop("t_save must be a positive multiple of dt")
    skip <- round(t_save/dt)
    t_dimnames <- seq(t_start, t_end, by = t_save)
    sim <- MizerSim(params, t_dimnames = t_dimnames)
    
    ## Initialise ----
    # get functions
    resource_dynamics_fn <- get(params@resource_dynamics)
    other_dynamics_fns <- lapply(params@other_dynamics, get)
    rates_fns <- lapply(params@rates_funcs, get)
    
    # Set up progress bar
    if (is(progress_bar, "Progress")) {
        # We have been passed a shiny progress object
        progress_bar$set(message = "Running simulation", value = 0)
        proginc <- 1/length(t_dimnames)
    } else if (progress_bar == TRUE) {
        pb <- progress::progress_bar$new(
            format = "[:bar] :percent ETA: :eta",
            total = length(t_dimnames), width = 60)
        pb$tick(0)
    }
    
    t <- t_start  # keep track of time
    i_save_time <- 2 # index into which to save the results
    
    # set up variables for keeping track of effort changes
    i_effort <- 1
    current_effort <- effort[1, ]
    t_next_effort <- time_effort[[2]] - 1e-8  # to avoid rounding problems
    
    # Set initial population and effort
    n_list <- list(n = initial_n, n_pp = initial_n_pp, n_other = initial_n_other)
    sim@n[1, , ] <- initial_n 
    sim@n_pp[1, ] <- initial_n_pp
    sim@n_other[1, ] <- initial_n_other
    sim@effort[1, ] <- current_effort
    
    ## Loop over time ----
    for (i in 2:length(t_dimnames)) {
        
        n_list <- project_simple(
            params, n = n_list$n, n_pp = n_list$n_pp, n_other = n_list$n_other,
            t = t, dt = dt, steps = skip, effort = current_effort,
            resource_dynamics_fn = resource_dynamics_fn,
            other_dynamics_fns = other_dynamics_fns,
            rates_fns = rates_fns, ...)
        t <- t + t_save
        
        # Advance progress bar
        if (is(progress_bar, "Progress")) {
            progress_bar$inc(amount = proginc)
        } else if (progress_bar == TRUE) {
            pb$tick()
        }
        
        # Do we need new effort?
        if (t >= t_next_effort) {
            i_effort <- i_effort + 1
            current_effort <- effort[i_effort, ]
            # We have to be careful when we reach the final time
            t_next_effort <- ifelse(i_effort < length(time_effort),
                                    time_effort[[i_effort + 1]] - 1e-8,
                                    t_max + 100)
        }
        
        # Store result
        sim@n[i, , ] <- n_list$n
        sim@n_pp[i, ] <- n_list$n_pp
        sim@n_other[i, ] <- n_list$n_other
        sim@effort[i, ] <- current_effort
    }
    
    # append to previous simulation ----
    if (is(object, "MizerSim") && append) {
        no_t_old <- dim(object@n)[1]
        no_t <- length(t_dimnames)
        new_t_dimnames <- c(as.numeric(dimnames(object@n)[[1]]),
                            t_dimnames[2:no_t])
        new_sim <- MizerSim(params, t_dimnames = new_t_dimnames)
        old_indices <- 1:no_t_old
        new_indices <- seq(from = no_t_old + 1, length.out = no_t - 1)
        new_sim@n[old_indices, , ]  <- object@n
        new_sim@n[new_indices, , ]  <- sim@n[2:no_t, , ]
        new_sim@n_pp[old_indices, ] <- object@n_pp
        new_sim@n_pp[new_indices, ] <- sim@n_pp[2:no_t, ]
        new_sim@n_other[old_indices, ]  <- object@n_other
        new_sim@n_other[new_indices, ]  <- sim@n_other[2:no_t, ]
        new_sim@effort[old_indices, ] <- object@effort
        new_sim@effort[new_indices, ] <- sim@effort[2:no_t, ]
        return(new_sim)
    }
    return(sim)
}

#' Project abundances by a given number of time steps into the future
#' 
#' This is an internal function used by the user-facing `project()` function.
#' It is of potential interest only to mizer extension authors.
#' 
#' The function does not check its arguments because it is meant to be as fast
#' as possible to allow it to be used in a loop. For example, it is called in
#' `project()` once for every saved value. The function also does not save its
#' intermediate results but only returns the result at time `t + dt * steps`.
#' During this time it uses the constant fishing effort `effort`.
#' 
#' The functional arguments can be calculated from slots in the `params` object
#' with
#' ```
#' resource_dynamics_fn <- get(params@resource_dynamics)
#' other_dynamics_fns <- lapply(params@other_dynamics, get)
#' rates_fns <- lapply(params@rates_funcs, get)
#' ```
#' The reason the function does not do that itself is to shave 20 microseconds
#' of its running time, which pays when the function is called hundreds of
#' times in a row.
#' 
#' This function is also used in `steady()`. In between calls to 
#' `project_simple()` the `steady()` function checks whether the values are
#' still changing significantly, so that it can stop when a steady state has
#' been approached. Mizer extension packages might have a similar need to run
#' a simulation repeatedly for short periods to run some other code in
#' between. Because this code may want to use the values of the rates at the
#' final time step, these too are included in the returned list.
#' 
#' @param params A MizerParams object.
#' @param n An array (species x size) with the number density at start of
#'   simulation.
#' @param n_pp A vector (size) with the resource number density at start of
#'   simulation.
#' @param n_other A named list with the abundances of other components at start
#'   of simulation.
#' @param t Time at the start of the simulation.
#' @param dt Size of time step.
#' @param steps The number of time steps by which to project.
#' @param effort The fishing effort to be used throughout the simulation. This
#'   must be a vector or list with one named entry per fishing gear.
#' @param resource_dynamics_fn The function for the resource
#'   dynamics. See Details.
#' @param other_dynamics_fns List with the functions for the
#'   dynamics of the other components. See Details.
#' @param rates_fns List with the functions for calculating
#'   the rates. See Details.
#' @param ... Other arguments that are passed on to the rate functions.
#' @return List with the final values of `n`, `n_pp` and `n_other`, `rates`.
#' 
#' @export
project_simple <- function(params, n, n_pp, n_other, t, dt, steps, 
                           effort, resource_dynamics_fn, other_dynamics_fns,
                           rates_fns, ...) {    
    # Handy things ----
    no_sp <- nrow(params@species_params) # number of species
    no_w <- length(params@w) # number of fish size bins
    idx <- 2:no_w
    # Hacky shortcut to access the correct element of a 2D array using 1D notation
    # This references the egg size bracket for all species, so for example
    # n[w_min_idx_array_ref] = n[,w_min_idx]
    w_min_idx_array_ref <- (params@w_min_idx - 1) * no_sp + (1:no_sp)
    # Matrices for solver
    a <- matrix(0, nrow = no_sp, ncol = no_w)
    b <- matrix(0, nrow = no_sp, ncol = no_w)
    S <- matrix(0, nrow = no_sp, ncol = no_w)

    # Loop over time steps ----
    for (i_time in 1:steps) {
        r <- rates_fns$Rates(
            params, n = n, n_pp = n_pp, n_other = n_other,
            t = t, effort = effort, rates_fns = rates_fns, ...)
        
        # Update time
        t <- t + dt
        
        # * Update other components ----
        n_other_current <- n_other  # So that the resource dynamics can still 
        # use the current value
        for (component in names(params@other_dynamics)) {
            n_other[[component]] <-
                other_dynamics_fns[[component]](
                    params,
                    n = n,
                    n_pp = n_pp,
                    n_other = n_other_current,
                    rates = r,
                    t = t,
                    dt = dt,
                    component = component,
                    ...
                )
        }
        
        # * Update resource ----
        n_pp <- resource_dynamics_fn(params, n = n, n_pp = n_pp,
                                     n_other = n_other_current, rates = r,
                                     t = t, dt = dt, ...)
        
        # * Update species ----
        # a_{ij} = - g_i(w_{j-1}) / dw_j dt
        a[, idx] <- sweep(-r$e_growth[, idx - 1, drop = FALSE] * dt, 2,
                          params@dw[idx], "/")
        # b_{ij} = 1 + g_i(w_j) / dw_j dt + \mu_i(w_j) dt
        b[] <- 1 + sweep(r$e_growth * dt, 2, params@dw, "/") + r$mort * dt
        # S_{ij} <- N_i(w_j)
        S[,idx] <- n[, idx, drop = FALSE]
        # Update first size group of n
        n[w_min_idx_array_ref] <-
            (n[w_min_idx_array_ref] + r$rdd * dt / 
                 params@dw[params@w_min_idx]) /
            b[w_min_idx_array_ref]
        # Update n
        # for (i in 1:no_sp) # number of species assumed small, so no need to 
        #                      vectorize this loop over species
        #     for (j in (params@w_min_idx[i]+1):no_w)
        #         n[i,j] <- (S[i,j] - A[i,j]*n[i,j-1]) / B[i,j]
        # This is implemented via Rcpp
        n <- inner_project_loop(no_sp = no_sp, no_w = no_w, n = n,
                                A = a, B = b, S = S,
                                w_min_idx = params@w_min_idx)
    }
    return(list(n = n, n_pp = n_pp, n_other = n_other, rates = r))
}
