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
#' \code{\link{summary_functions}} and \code{\link{plotting_functions}}.
#' 
#' @param object Either a \linkS4class{MizerParams} object or a 
#'   \linkS4class{MizerSim} object (which contains a \code{MizerParams} object).
#' @param effort The effort of each fishing gear through time. See notes below.
#' @param t_max The number of years the projection runs for. The default value is
#'   100. However, this argument is ignored if an array is used for the
#'   \code{effort} argument. See notes below.
#' @param dt Time step of the solver. The default value is 0.1.
#' @param t_save The frequency with which the output is stored. The default
#'   value is 1. Must be an integer multiple of dt.
#' @param t_start The the year of the start of the simulation. The simulation
#'   will cover the period from \code{t_start} to \code{t_start + t_max}.
#'   Defaults to 0. Ignored if an array is used for the \code{effort}
#'   argument or a \code{MizerSim} for the \code{object} argument.
#' @param initial_n The initial abundances of species. A matrix with dimensions 
#'   species x size. The order of species must be the same as in the 
#'   \code{MizerParams} argument. Ignored if the \code{object} argument is a 
#'   MizerSim object, but overrules the \code{initial_n} slot if \code{object} 
#'   is a \code{MizerParams} object.
#' @param initial_n_pp The initial abundances of plankton. A numeric vector.
#'   Ignored if the \code{object} argument is a MizerSim object, but overrules
#'   the \code{initial_n_pp} slot if \code{object} is a \code{MizerParams}
#'   object.
#' @param initial_B The initial biomasses of the unstructured resources. It
#'   should be a named vector with one entry for each resource. Ignored if the
#'   \code{object} argument is a MizerSim object, but overrules the
#'   \code{initial_B} slot if \code{object} is a \code{MizerParams} object.
#' @param append A boolean that determines whether the new simulation results
#'   are appended to the previous ones. Only relevant if \code{object} is a
#'   \code{MizerSim} object. Default = TRUE.
#' @param progress_bar Either a boolean value to determine whether a progress
#'   bar should be shown in the console of a shiny progress object to implement 
#'   a progress bar in a shiny app
#' @param ... Currently unused.
#' 
#' @note The \code{effort} argument specifies the level of fishing effort during
#' the simulation. It can be specified in three different ways: \itemize{ \item
#' A single numeric value. This specifies the effort of all fishing gears which
#' is constant through time (i.e. all the gears have the same constant effort). 
#' \item A numerical vector which has the same length as the number of fishing
#' gears. The vector must be named and the names must correspond to the gear
#' names in the \code{MizerParams} object. The values in the vector specify the
#' constant fishing effort of each of the fishing gears, i.e. the effort is
#' constant through time but each gear may have a different fishing effort. 
#' \item A numerical array with dimensions time x gear. This specifies the
#' fishing effort of each gear at each time step.  The first dimension, time,
#' must be named numerically and contiguously. The second dimension of the array
#' must be named and the names must correspond to the gear names in the
#' \code{MizerParams} argument. The value for the effort for a particular time
#' is used during the interval from that time to the next time in the array.}
#' 
#' If effort is specified as an array then the smallest time in the array is 
#' used as the initial time for the simulation. Otherwise the initial time is
#' set to the final time of the previous simulation if \code{object} is a 
#' \code{MizerSim} object or to \code{t_start} otherwise. Also, if the effort is
#' an array then the \code{t_max} argument is ignored and the maximum simulation
#' time is the largest time of the effort array.
#' 
#' If the \code{object} argument is of class \code{MizerSim} then the initial
#' values for the simulation are taken from the final values in the 
#' \code{MizerSim} object and the corresponding arguments to this function will
#' be ignored.
#' 
#' @return An object of class \linkS4class{MizerSim}.
#' 
#' @export
#' @seealso \code{\link{MizerParams}}, \code{\link{summary_functions}} and 
#'   \code{\link{plotting_functions}}
#' @examples
#' \dontrun{
#' # Data set with different fishing gears
#' data(NS_species_params_gears)
#' data(inter)
#' params <- set_multispecies_model(NS_species_params_gears, inter)
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
project <- function(object, effort = 1,
                    t_max = 100, dt = 0.1, t_save = 1, t_start = 0,
                    initial_n = params@initial_n,
                    initial_n_pp = params@initial_n_pp,
                    initial_B = params@initial_B,
                    append = TRUE,
                    progress_bar = TRUE, ...) {
    validObject(object)
    if (is(object, "MizerSim")) {
        params <- object@params
        no_t <- dim(object@B)[1]
        initial_n <- params@initial_n # Needed to get the right dimensions
        initial_n[] <- object@n[no_t, , ]
        initial_n_pp <- params@initial_n_pp # Needed to get the right dimensions
        initial_n_pp[] <- object@n_pp[no_t, ]
        initial_B <- params@initial_B # Needed to get the right dimensions
        initial_B[] <- object@B[no_t, ]
        t_start <- as.numeric(dimnames(object@n)[[1]][no_t])
    } else {
        params <- object
        if (missing(initial_n))       initial_n <- params@initial_n
        if (missing(initial_n_pp)) initial_n_pp <- params@initial_n_pp
        if (missing(initial_B))       initial_B <- params@initial_B
        t_start <- 0
    }
    params@initial_n[] <- initial_n
    params@initial_n_pp[] <- initial_n_pp
    params@initial_B[] <- initial_B
    
    no_sp <- length(params@w_min_idx)
    assert_that(is.array(initial_n),
                are_equal(dim(initial_n), c(no_sp, length(params@w))))
    assert_that(is.vector(initial_n_pp),
                length(initial_n_pp) == length(params@w_full))
    if (length(params@resource_dynamics) > 0) {
        if (!is.character(names(initial_B))) {
            stop("The initial_B needs to be a named vector")
        }
        if (!setequal(names(initial_B), names(params@resource_dynamics))) {
            stop("The names of the entries in initial_B do not match the names of the unstructured resource components of the model.")
        }
    }
    
    # Create effort array ----
    # Do we need to create an effort array?
    if (is.vector(effort)) {
        no_gears <- dim(params@catchability)[1]
        if ((length(effort) > 1) & (length(effort) != no_gears)) {
            stop("Effort vector must be the same length as the number of fishing gears\n")
        }
        # If more than 1 gear need to check that gear names match
        gear_names <- dimnames(params@catchability)[[1]]
        effort_gear_names <- names(effort)
        if (length(effort) == 1 & is.null(effort_gear_names)) {
            effort_gear_names <- gear_names
        }
        if (!all(gear_names %in% effort_gear_names)) {
            stop(paste0("Gear names in the MizerParams object (", 
                        paste(gear_names, collapse = ", "), 
                        ") do not match those in the effort vector."))
        }
        # Set up the effort array transposed so we can use the recycling rules
        time_dimnames <- signif(seq(from = t_start, 
                                    to = t_start + t_max, 
                                    by = dt), 3)
        effort <- t(array(effort, dim = c(no_gears, length(time_dimnames)), 
                          dimnames = list(gear = effort_gear_names, 
                                          time = time_dimnames)))
    }
    
    # Check that number and names of gears in effort array is same as in 
    # MizerParams object
    no_gears <- dim(params@catchability)[1]
    if (dim(effort)[2] != no_gears) {
        stop(paste0("The number of gears in the effort array (length of the second dimension = ", 
                   dim(effort)[2], 
                   ") does not equal the number of gears in the MizerParams object (", 
                   no_gears, ")."))
    }
    gear_names <- dimnames(params@catchability)[[1]]
    if (!all(gear_names %in% dimnames(effort)[[2]])) {
        stop(paste0("Gear names in the MizerParams object (", 
                    paste(gear_names, collapse = ", "), 
                    ") do not match those in the effort array."))
    }
    # Sort effort array to match order in MizerParams
    effort <- effort[, gear_names, drop = FALSE]
    
    # Blow up time dimension of effort array
    # i.e. effort might have been passed in using time steps of 1, but actual 
    # dt = 0.1, so need to blow up
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
    t_end <- time_effort[length(time_effort)]
    # Blow up effort so that rows are dt spaced
    time_effort_dt <- seq(from = time_effort[1], to = t_end, by = dt)
    effort_dt <- t(array(NA, dim = c(length(time_effort_dt), dim(effort)[2]), 
                         dimnames = list(time = time_effort_dt,
                                         dimnames(effort)[[2]])))
    for (i in 1:(length(time_effort) - 1)) {
        effort_dt[,time_effort_dt >= time_effort[i]] <- effort[i,]
    }
    effort_dt <- t(effort_dt)
    
    # Make the MizerSim object with the right size ----
    # We only save every t_save steps
    # Divisibility test needs to be careful about machine rounding errors,
    # see https://github.com/sizespectrum/mizer/pull/2
    if ((t_save < dt) || !isTRUE(all.equal((t_save - round(t_save / dt) * dt), 0)))
        stop("t_save must be a positive multiple of dt")
    t_skip <- round(t_save/dt)
    t_dimnames_index <- seq(1, to = length(time_effort_dt), by = t_skip)
    t_dimnames <- time_effort_dt[t_dimnames_index]
    sim <- MizerSim(params, t_dimnames = t_dimnames) 
    # Fill up the effort array
    sim@effort[] <- effort_dt[t_dimnames_index,]
    
    ## Initialise ----
    # Set initial population
    sim@n[1, , ] <- initial_n 
    sim@n_pp[1, ] <- initial_n_pp
    sim@B[1, ] <- initial_B
    
    # Handy things
    no_sp <- nrow(sim@params@species_params) # number of species
    no_w <- length(sim@params@w) # number of fish size bins
    idx <- 2:no_w
    # Hacky shortcut to access the correct element of a 2D array using 1D notation
    # This references the egg size bracket for all species, so for example
    # n[w_minidx_array_ref] = n[,w_min_idx]
    w_min_idx_array_ref <- (sim@params@w_min_idx - 1) * no_sp + (1:no_sp)
    
    # sex ratio - DO SOMETHING LATER WITH THIS
    sex_ratio <- 0.5
    
    # Matrices for solver
    a <- matrix(0, nrow = no_sp, ncol = no_w)
    b <- matrix(0, nrow = no_sp, ncol = no_w)
    S <- matrix(0, nrow = no_sp, ncol = no_w)
    
    # initialise n n_pp and B
    # We want the first time step only but cannot use drop as there may only 
    # be a single species
    n <- array(sim@n[1, , ], dim = dim(sim@n)[2:3])
    dimnames(n) <- dimnames(sim@n)[2:3]
    n_pp <- sim@n_pp[1, ]
    B <- initial_B
    
    # Set up progress bar
    if (progress_bar == TRUE) {
        pb <- progress::progress_bar$new(
            format = "[:bar] :percent ETA: :eta",
            total = length(t_dimnames_index), width = 60)
    }
    if (is(progress_bar, "Progress")) {
        # We have been passed a shiny progress object
        progress_bar$set(message = "Running simulation", value = 0)
        proginc <- 1/length(t_dimnames_index)
    }
    
    ## Loop over time ----
    t <- 0  # keep track of time
    t_steps <- dim(effort_dt)[1] - 1
    for (i_time in 1:t_steps) {
        r <- getRates(params, n = n, n_pp = n_pp, B = B,
                      effort = effort_dt[i_time,])
        
        # Update unstructured resource biomasses
        B_current <- B  # So that the plankton dynamics can still use the 
                        # current value
        if (length(params@resource_dynamics) > 0) {
            for (res in names(params@resource_dynamics)) {
                B[res] <-
                    params@resource_dynamics[[res]](
                        params,
                        n = n,
                        n_pp = n_pp,
                        B = B_current,
                        rates = r,
                        dt = dt
                    )
            }
        }
        
        # Update plankton
        n_pp <- params@plankton_dynamics(params, n = n, n_pp = n_pp, 
                                         B = B_current, rates = r, dt = dt)
        
        # Iterate species one time step forward:
        # See Ken's PDF
        # a_{ij} = - g_i(w_{j-1}) / dw_j dt
        a[, idx] <- sweep(-r$e_growth[, idx - 1, drop = FALSE] * dt, 2,
                          sim@params@dw[idx], "/")
        # b_{ij} = 1 + g_i(w_j) / dw_j dt + \mu_i(w_j) dt
        b[, idx] <- 1 + sweep(r$e_growth[, idx, drop = FALSE] * dt, 2, 
                              sim@params@dw[idx], "/") +
                        r$mort[, idx, drop = FALSE] * dt
        # S_{ij} <- N_i(w_j)
        S[,idx] <- n[, idx, drop = FALSE]
        # Boundary condition upstream end (recruitment)
        b[w_min_idx_array_ref] <- 1 + r$e_growth[w_min_idx_array_ref] * dt /
                                        sim@params@dw[sim@params@w_min_idx] +
                                    r$mort[w_min_idx_array_ref] * dt
        # Update first size group of n
        n[w_min_idx_array_ref] <-
            (n[w_min_idx_array_ref] + r$rdd * dt / 
                 sim@params@dw[sim@params@w_min_idx]) /
            b[w_min_idx_array_ref]
        # Update n
        # for (i in 1:no_sp) # number of species assumed small, so no need to 
        #                      vectorize this loop over species
        #     for (j in (sim@params@w_min_idx[i]+1):no_w)
        #         n[i,j] <- (S[i,j] - A[i,j]*n[i,j-1]) / B[i,j]
        # This is implemented via Rcpp
        n <- inner_project_loop(no_sp = no_sp, no_w = no_w, n = n,
                                A = a, B = b, S = S,
                                w_min_idx = sim@params@w_min_idx)
        
        # Update time
        t <- t + dt
        # Store results only every t_step steps.
        store <- t_dimnames_index %in% (i_time + 1)
        if (any(store)) {
            # Advance progress bar
            if (is(progress_bar, "Progress")) {
                progress_bar$inc(amount = proginc)
            }
            if (progress_bar == TRUE) {
                pb$tick()
            }
            # Store result
            sim@n[which(store), , ] <- n
            sim@n_pp[which(store), ] <- n_pp
            sim@B[which(store), ] <- B
        }
    }
    if (is(object, "MizerSim") && append) {
        # append to previous simulation ----
        no_t_old <- dim(object@n)[1]
        no_t <- length(t_dimnames)
        new_t_dimnames <- c(as.numeric(dimnames(object@n)[[1]]),
                            t_dimnames[2:length(t_dimnames)])
        new_sim <- MizerSim(params, t_dimnames = new_t_dimnames)
        old_indices <- 1:no_t_old
        new_indices <- seq(from = no_t_old + 1, length.out = no_t - 1)
        new_sim@n[old_indices, , ]  <- object@n
        new_sim@n[new_indices, , ]  <- sim@n[2:no_t, , ]
        new_sim@n_pp[old_indices, ] <- object@n_pp
        new_sim@n_pp[new_indices, ] <- sim@n_pp[2:no_t, ]
        new_sim@B[old_indices, ]    <- object@B
        new_sim@B[new_indices, ]    <- sim@B[2:no_t, ]
        new_sim@effort[old_indices, ] <- object@effort
        new_sim@effort[new_indices, ] <- sim@effort[2:no_t, ]
        return(new_sim)
    }
    return(sim)
}


#' Calculate initial population abundances for the community populations
#' 
#' This function uses the model parameters and other parameters to calculate 
#' initial population abundances for the community populations. These initial 
#' abundances should be reasonable guesses at the equilibrium values. The 
#' returned population can be passed to the \code{project} function.
#' 
#' @param params The model parameters. An object of type \linkS4class{MizerParams}.
#' @param a A parameter with a default value of 0.35.
#' @param n0_mult Multiplier for the abundance at size 0. Default value is
#'   kappa/1000.
#' @export
#' @return A matrix (species x size) of population abundances.
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' params <- set_multispecies_model(NS_species_params_gears)
#' init_n <- get_initial_n(params)
#' }
get_initial_n <- function(params, n0_mult = NULL, a = 0.35) {
    if (!is(params,"MizerParams"))
        stop("params argument must of type MizerParams")
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)
    initial_n <- array(NA, dim = c(no_sp, no_w))
    dimnames(initial_n) <- dimnames(params@intake_max)
    # N = N0 * Winf^(2*n-q-2+a) * w^(-n-a)
    # Reverse calc n and q from intake_max and search_vol slots (could add get_n function)
    n <- (log(params@intake_max[,1] / params@species_params$h) / log(params@w[1]))[1]
    q <- (log(params@search_vol[,1] / params@species_params$gamma) / log(params@w[1]))[1]
    # Guessing at a suitable n0 value based on kappa - this was figured out using trial and error and should be updated
    if (is.null(n0_mult)) {
        lambda <- 2 + q - n
        kappa <- params@cc_pp[1] / (params@w_full[1]^(-lambda))
        n0_mult <- kappa / 1000
    }
    initial_n[] <- unlist(tapply(params@w, 1:no_w, function(wx,n0_mult,w_inf,a,n,q)
        n0_mult * w_inf^(2 * n - q - 2 + a) * wx^(-n - a),
        n0_mult = n0_mult, w_inf = params@species_params$w_inf, a=a, n=n, q=q))
    #set densities at w > w_inf to 0
    initial_n[unlist(tapply(params@w,1:no_w,function(wx,w_inf) w_inf<wx, w_inf=params@species_params$w_inf))] <- 0
    # Also any densities at w < w_min set to 0
    initial_n[unlist(tapply(params@w,1:no_w,function(wx,w_min)w_min>wx, w_min=params@species_params$w_min))] <- 0    
    return(initial_n)
}
