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
#' @param t_max The number of years the projection runs for. The default value
#'   is 100. When an effort array is supplied, this argument can be used to
#'   extend the simulation beyond the times specified in the effort array.
#'   See notes below.
#' @param dt Time step of the solver. The default value is 0.1.
#' @param t_save The frequency with which the output is stored. The default
#'   value is 1. When an effort array is supplied, this argument can be used
#'   to control the times at which the simulation results are saved.
#'   See notes below.
#' @param t_start The the year of the start of the simulation. The simulation
#'   will cover the period from `t_start` to \code{t_start + t_max}.
#'   Defaults to 0. Ignored if an array is used for the `effort`
#'   argument or a `MizerSim` for the `object` argument.
#' @param initial_n `r lifecycle::badge("deprecated")` The initial abundances of
#'   species. Instead of using this argument you should set `initialN(params)`
#'   to the desired value.
#' @param initial_n_pp `r lifecycle::badge("deprecated")` The initial abundances
#'   of resource. Instead of using this argument you should set
#'   `initialNResource(params)` to the desired value.
#' @param append A boolean that determines whether the new simulation results
#'   are appended to the previous ones. Only relevant if `object` is a
#'   `MizerSim` object. Default = TRUE.
#' @param progress_bar Either a boolean value to determine whether a progress
#'   bar should be shown in the console, or a shiny Progress object to implement
#'   a progress bar in a shiny app.
#' @param ... Other arguments will be passed to rate functions.
#'
#' @note The `effort` argument specifies the level of fishing effort during the
#'   simulation. If it is not supplied, the initial effort stored in the params
#'   object is used. The effort can be specified in four different ways:
#' \itemize{
#' \item A single numeric value. This specifies the effort of all fishing gears
#' which is constant through time (i.e. all the gears have the same constant
#' effort).
#' \item A named vector whose names match with existing gear names.
#'  The values in the vector specify the constant fishing effort for those
#'  fishing gears, i.e. the effort is constant through time. The
#'   effort for gears that are not included in the effort vector is set to 0.
#' \item A numerical vector which has the same length as the number of fishing
#' gears. The values in the vector specify the
#' constant fishing effort of each of the fishing gears, with the ordering
#' assumed to be the same as in the MizerParams object.
#' \item A numerical array with dimensions time x gear. This specifies the
#' fishing effort of each gear at each time step.  The first dimension, time,
#' must be named numerically and increasing. The second dimension of the array
#' must be named and the names must correspond to the gear names in the
#' `MizerParams` object. The value for the effort for a particular time
#' is used during the interval from that time to the next time in the array.
#' }
#'
#' If effort is specified as an array then the smallest time in the array is
#' used as the initial time for the simulation. Otherwise the initial time is
#' set to the final time of the previous simulation if `object` is a
#' `MizerSim` object or to `t_start` otherwise.
#'
#' When an effort array is provided, the `t_max` argument can be used to
#' extend the simulation beyond the last time specified in the effort array.
#' In this case, the effort values from the last time in the array will be
#' used for the extended period. The `t_save` argument can be used to specify
#' the frequency at which simulation results are saved. If `t_save` is not
#' supplied, the results will be saved at the times specified in the effort
#' array. If both `t_max` and `t_save` are provided with an effort array,
#' effort values will be interpolated (using step function) or extrapolated
#' (using the last known value) as needed for the new time points. The
#' `t_start` argument continues to be ignored when an effort array is supplied.
#'
#' Note that if `t_max` or `t_save` are specified, the time grid for the
#' simulation is resampled based on `t_save`. This means that if the time
#' points in the `effort` array are irregular and do not align with the new
#' grid, those specific time points may be lost and the effort values at the
#' new grid points will be calculated via interpolation.
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
#' \donttest{
#' params <- NS_params
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # With constant fishing effort which is different for each gear
#' effort <- c(Industrial = 0, Pelagic = 1, Beam = 0.5, Otter = 0.5)
#' sim <- project(params, t_max = 20, effort = effort)
#' # With fishing effort that varies through time for each gear
#' gear_names <- c("Industrial", "Pelagic", "Beam", "Otter")
#' times <- seq(from = 1, to = 10, by = 1)
#' effort_array <- array(NA,
#'     dim = c(length(times), length(gear_names)),
#'     dimnames = list(time = times, gear = gear_names)
#' )
#' effort_array[, "Industrial"] <- 0.5
#' effort_array[, "Pelagic"] <- seq(from = 1, to = 2, length = length(times))
#' effort_array[, "Beam"] <- seq(from = 1, to = 0, length = length(times))
#' effort_array[, "Otter"] <- seq(from = 1, to = 0.5, length = length(times))
#' sim <- project(params, effort = effort_array)
#' # Extend a simulation beyond the effort array times
#' # Effort values from the final time are used for the extension
#' sim <- project(params, effort = effort_array, t_max = 15)
#' # Control save times with an effort array using t_save
#' sim <- project(params, effort = effort_array, t_save = 2)
#' }
#' @rdname project
#' @export
project <- function(object, effort,
                    t_max = 100, dt = 0.1, t_save = 1, t_start = 0,
                    initial_n, initial_n_pp,
                    append = TRUE,
                    progress_bar = TRUE, ...) {
    UseMethod("project")
}

#' @rdname project
#' @export
project.MizerParams <- function(object, effort,
                                t_max = 100, dt = 0.1, t_save = 1, t_start = 0,
                                initial_n, initial_n_pp,
                                append = TRUE,
                                progress_bar = TRUE, ...) {
    params <- validParams(object)
    # Set and check initial values ----
    assert_that(t_max > 0)
    if (!missing(initial_n)) params@initial_n[] <- initial_n
    if (!missing(initial_n_pp)) params@initial_n_pp[] <- initial_n_pp

    initial_n <- params@initial_n
    initial_n_pp <- params@initial_n_pp
    initial_n_other <- params@initial_n_other

    no_sp <- length(params@w_min_idx)
    assert_that(
        is.array(initial_n),
        is.numeric(initial_n),
        are_equal(dim(initial_n), c(no_sp, length(params@w)))
    )
    assert_that(
        is.numeric(initial_n_pp),
        length(initial_n_pp) == length(params@w_full)
    )

    assert_that(is.null(initial_n_other) || is.list(initial_n_other))
    other_names <- names(params@other_dynamics)
    if (length(other_names) > 0) {
        if (is.null(names(initial_n_other))) {
            stop("The initial_n_other needs to be a named list")
        }
        if (!setequal(names(initial_n_other), other_names)) {
            stop(
                "The names of the entries in initial_n_other do not match ",
                "the names of the other components of the model."
            )
        }
    }

    # Set effort array ----
    if (missing(effort)) effort <- params@initial_effort

    if (is.null(dim(effort))) { # effort is a vector or scalar
        # Set up the effort array transposed so we can use the recycling rules
        # no point running a simulation with no saved results
        if (t_max < t_save) {
            t_save <- t_max
        }
        times <- seq(t_start, t_start + t_max, by = t_save)
        effort <- validEffortVector(effort, params)
        effort <- t(array(effort,
            dim = c(length(effort), length(times)),
            dimnames = list(
                gear = names(effort),
                time = times
            )
        ))
    } else {
        effort <- validEffortArray(effort, params)

        # When an effort array is supplied, we now respect t_max and t_save
        # if they are explicitly provided by the user
        effort_times <- as.numeric(dimnames(effort)[[1]])
        t_start_effort <- effort_times[1]
        t_end_effort <- effort_times[length(effort_times)]

        # Check if t_max or t_save were explicitly provided by the user
        t_max_provided <- !missing(t_max)
        t_save_provided <- !missing(t_save)

        if (t_max_provided || t_save_provided) {
            # Determine the end time
            if (t_max_provided) {
                t_end <- t_start_effort + t_max
            } else {
                # Keep the end time from the effort array
                t_end <- t_end_effort
            }

            # Determine the save frequency
            if (t_save_provided) {
                # Use the provided t_save
                save_freq <- t_save
            } else {
                # Infer from effort array times or use default
                if (length(effort_times) > 1) {
                    # Use the first interval as the save frequency
                    save_freq <- effort_times[2] - effort_times[1]
                } else {
                    save_freq <- 1 # default
                }
            }

            # Ensure we don't use a t_save larger than the duration
            duration <- t_end - t_start_effort
            if (duration < save_freq) {
                save_freq <- duration
            }

            # Generate new time points
            times <- seq(t_start_effort, t_end, by = save_freq)

            # Interpolate effort values for the new time points
            gear_names <- dimnames(effort)[[2]]
            new_effort <- array(NA,
                dim = c(length(times), length(gear_names)),
                dimnames = list(time = times, gear = gear_names)
            )

            for (g in seq_along(gear_names)) {
                # Use approx with rule=2 for constant extrapolation beyond range
                # method="constant" means the effort value applies from that time
                # until the next time point (step function)
                new_effort[, g] <- stats::approx(
                    x = effort_times,
                    y = effort[, g],
                    xout = times,
                    method = "constant",
                    rule = 2,
                    f = 0
                )$y
            }
            effort <- new_effort
        }
    }

    times <- as.numeric(dimnames(effort)[[1]])

    # Make the MizerSim object with the right size ----
    # We only save every t_save years
    sim <- MizerSim(params, t_dimnames = times)
    # Set initial population and effort
    sim@n[1, , ] <- initial_n
    sim@n_pp[1, ] <- initial_n_pp
    sim@n_other[1, ] <- initial_n_other
    sim@effort <- effort

    ## Initialise ----
    # get functions
    resource_dynamics_fn <- get(params@resource_dynamics)
    other_dynamics_fns <- lapply(params@other_dynamics, get)
    rates_fns <- lapply(params@rates_funcs, get)

    # Set up progress bar
    if (is(progress_bar, "Progress")) {
        # We have been passed a shiny progress object
        progress_bar$set(message = "Running simulation", value = 0)
        proginc <- 1 / length(times)
    } else if (progress_bar == TRUE) {
        pb <- progress::progress_bar$new(
            format = "[:bar] :percent ETA: :eta",
            total = length(times), width = 60
        )
        pb$tick(0)
    }

    n_list <- list(
        n = initial_n, n_pp = initial_n_pp,
        n_other = unserialize(serialize(initial_n_other, NULL))
    )
    t <- times[[1]]

    ## Loop over time ----
    for (i in 2:length(times)) {
        # number of time steps between saved times
        steps <- round((times[[i]] - t) / dt)
        # advance to next saved time
        n_list <- project_simple(
            params,
            n = n_list$n, n_pp = n_list$n_pp, n_other = n_list$n_other,
            t = t, dt = dt, steps = steps,
            effort = effort[i - 1, ],
            resource_dynamics_fn = resource_dynamics_fn,
            other_dynamics_fns = other_dynamics_fns,
            rates_fns = rates_fns, ...
        )
        # Calculate start time for next iteration
        # The reason we don't simply use the next entry in `times` is that
        # those entries may not be separated by exact multiples of dt.
        t <- t + steps * dt
        # Advance progress bar
        if (is(progress_bar, "Progress")) {
            progress_bar$inc(amount = proginc)
        } else if (progress_bar == TRUE) {
            pb$tick()
        }

        # Store result
        sim@n[i, , ] <- n_list$n
        sim@n_pp[i, ] <- n_list$n_pp
        sim@n_other[i, ] <- unserialize(serialize(n_list$n_other, NULL))
    }

    return(sim)
}

#' @rdname project
#' @export
project.MizerSim <- function(object, effort,
                             t_max = 100, dt = 0.1, t_save = 1, t_start = 0,
                             initial_n, initial_n_pp,
                             append = TRUE,
                             progress_bar = TRUE, ...) {
    validObject(object)
    params <- setInitialValues(object@params, object)
    t_start <- getTimes(object)[idxFinalT(object)]

    sim <- project(params,
        effort = effort, t_max = t_max, dt = dt,
        t_save = t_save, t_start = t_start,
        initial_n = initial_n, initial_n_pp = initial_n_pp,
        progress_bar = progress_bar, ...
    )

    if (append) {
        times <- as.numeric(dimnames(sim@n)[[1]])
        no_t_old <- dim(object@n)[1]
        no_t <- length(times)
        new_t_dimnames <- c(
            as.numeric(dimnames(object@n)[[1]]),
            times[2:no_t]
        )
        new_sim <- MizerSim(params, t_dimnames = new_t_dimnames)
        old_indices <- 1:no_t_old
        new_indices <- seq(from = no_t_old + 1, length.out = no_t - 1)
        new_sim@n[old_indices, , ] <- object@n
        new_sim@n[new_indices, , ] <- sim@n[2:no_t, , ]
        new_sim@n_pp[old_indices, ] <- object@n_pp
        new_sim@n_pp[new_indices, ] <- sim@n_pp[2:no_t, ]
        new_sim@n_other[old_indices, ] <- object@n_other
        new_sim@n_other[new_indices, ] <- sim@n_other[2:no_t, ]
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
#' @concept helper
project_simple <- function(params, n, n_pp, n_other, effort, t, dt, steps,
                           resource_dynamics_fn, other_dynamics_fns,
                           rates_fns, ...) {
    UseMethod("project_simple")
}

#' @export
project_simple.MizerParams <-
    function(params,
             n = params@initial_n,
             n_pp = params@initial_n_pp,
             n_other = params@initial_n_other,
             effort = params@initial_effort,
             t = 0, dt = 0.1, steps,
             resource_dynamics_fn = get(params@resource_dynamics),
             other_dynamics_fns = lapply(params@other_dynamics, get),
             rates_fns = lapply(params@rates_funcs, get), ...) {
        # Handy things ----
        no_sp <- nrow(params@species_params) # number of species
        no_w <- length(params@w) # number of fish size bins
        idx <- 2:no_w
        # Hacky shortcut to access the correct element of a 2D array using 1D
        # notation
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
                params,
                n = n, n_pp = n_pp, n_other = n_other,
                t = t, effort = effort, rates_fns = rates_fns, ...
            )

            # * Update other components ----
            n_other_new <- list() # So that the resource dynamics can still
            # use the current value
            for (component in names(params@other_dynamics)) {
                n_other_new[[component]] <-
                    other_dynamics_fns[[component]](
                        params,
                        n = n,
                        n_pp = n_pp,
                        n_other = n_other,
                        rates = r,
                        t = t,
                        dt = dt,
                        component = component,
                        ...
                    )
            }

            # * Update resource ----
            n_pp <- resource_dynamics_fn(params,
                n = n, n_pp = n_pp,
                n_other = n_other, rates = r,
                t = t, dt = dt,
                resource_rate = params@rr_pp,
                resource_capacity = params@cc_pp, ...
            )

            # * Update species ----
            n <- project_n(params, r, n, dt, a, b, S, idx, w_min_idx_array_ref,
                           no_sp, no_w)

            # * Update time ----
            t <- t + dt
        }

        return(list(n = n, n_pp = n_pp, n_other = n_other_new, rates = r))
    }

validEffortArray <- function(effort, params) {
    # Check that number and names of gears in effort array is same as in
    # MizerParams object
    no_gears <- dim(params@catchability)[1]
    if (dim(effort)[2] != no_gears) {
        stop(
            "The number of gears in the effort array (length of the second dimension = ",
            dim(effort)[2],
            ") does not equal the number of gears in the MizerParams object (",
            no_gears, ")."
        )
    }
    gear_names <- dimnames(params@catchability)[[1]]
    if (!all(gear_names %in% dimnames(effort)[[2]])) {
        stop(
            "Gear names in the MizerParams object (",
            paste(gear_names, collapse = ", "),
            ") do not match those in the effort array."
        )
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

    # Replace any NA's with default value
    effort_default <- ifelse(defaults_edition() < 2, 0, 1)
    effort[is.na(effort)] <- effort_default

    names(dimnames(effort)) <- c("time", "gear")
    effort
}
