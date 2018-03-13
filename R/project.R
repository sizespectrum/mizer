# Project method for the size based modelling package mizer

# Copyright 2012 Finlay Scott and Julia Blanchard. 
# Distributed under the GPL 3 or later 
# Maintainer: Gustav Delius, University of York, <gustav.delius@york.ac.uk>

# project can dispatch with effort being different classes (missing, numeric,
# array).

#' @useDynLib mizer
#' @importFrom Rcpp sourceCpp
NULL

#' project method for the size based modelling
#' 
#' Runs the size-based model simulation and projects the size based model
#' through time. \code{project} is called using an object of type
#' \linkS4class{MizerParams} and an object that contains the effort of the
#' fishing gears through time. The method returns an object of type
#' \linkS4class{MizerSim} which can then be explored with a range of summary and
#' plotting methods.
#' 
#' @param object A \code{MizerParams} object
#' @param effort The effort of each fishing gear through time. See notes below.
#' @param t_max The maximum time the projection runs for. The default value is
#'   100. However, this argument is not needed if an array is used for the
#'   \code{effort} argument, in which case this argument is ignored. See notes
#'   below.
#' @param dt Time step of the solver. The default value is 0.1.
#' @param t_save The frequency with which the output is stored. The default
#'   value is 1. Must be an integer multiple of dt.
#' @param initial_n The initial populations of the species. By default the 
#'   \code{initial_n} slot of the \code{\link{MizerParams}} argument is used.
#'   See the notes below.
#' @param initial_n_pp The initial population of the background spectrum. It
#'   should be a numeric vector of the same length as the \code{w_full} slot of
#'   the \code{MizerParams} argument. By default the \code{initial_n_pp} slot of the
#'   \code{\link{MizerParams}} argument is used.
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
#' \item A numerical array with dimensions time step x gear. This specifies the
#' fishing effort of each gear at each time step.  The first dimension, time,
#' must be named numerically and contiguously. The second dimension of the array
#' must be named and the names must correspond to the gear names in the
#' \code{MizerParams} argument. The value for the effort for a particular time
#' is used during the interval from that time to the next time in the array.}
#' 
#' If effort is specified as an array then the smallest time in the array is 
#' used as the initial time for the simulation. Otherwise the initial time is
#' set to 0. Also, if the effort is an array then the \code{t_max} argument is 
#' ignored and the maximum simulation time is the largest time of the effort
#' array.
#' 
#' The \code{initial_n} argument is a matrix with dimensions species x size. 
#' It specifies the abundances of the species at the initial time. The
#' order of species must be the same as in the \code{MizerParams} argument. If
#' the initial population is not specified, the argument is set by default by
#' the \code{get_initial_n} function which is set up for a North Sea model.
#' @return An object of type \code{MizerSim}.
#' @export
#' @seealso \code{\link{MizerParams}}
#' @examples
#' \dontrun{
#' # Data set with different fishing gears
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
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
setGeneric('project', function(object, effort, ...)
    standardGeneric('project'))

# No effort is specified - default is to set an effort of 0
# All other arguments passed as ...

#' Project without an effort argument.
#' @rdname project
setMethod('project', signature(object='MizerParams', effort='missing'),
    function(object, ...){
	res <- project(object, effort=0, ...)
	return(res)
})

#' Project with a constant effort.
#' @rdname project
setMethod('project', signature(object='MizerParams', effort='numeric'),
    function(object, effort,  t_max = 100, dt = 0.1, ...){
	no_gears <- dim(object@catchability)[1]
	if ((length(effort)>1) & (length(effort) != no_gears))
	    stop("Effort vector must be the same length as the number of fishing gears\n")
    # If more than 1 gear need to check that gear names match
	gear_names <- dimnames(object@catchability)[[1]]
    effort_gear_names <- names(effort)
    if (length(effort) == 1 & is.null(effort_gear_names)){
        effort_gear_names <- gear_names
    }
    if(!all(gear_names %in% effort_gear_names)){
        gear_names_error_message <- paste("Gear names in the MizerParams object (", paste(gear_names, collapse=", "), ") do not match those in the effort vector.", sep="")
        stop(gear_names_error_message)
    }
	# Set up the effort array transposed so we can use the recycling rules
    time_dimnames <- signif(seq(from=0, to=t_max, by=dt), 3)
	effort_array <- t(array(effort, dim=c(no_gears,length(time_dimnames)), dimnames=list(gear=effort_gear_names,time=time_dimnames)))
	res <- project(object, effort_array, dt=dt, ...)
	return(res)
})

#' Project with time varying effort
#' @rdname project
setMethod('project', signature(object='MizerParams', effort='array'),
    function(object, effort, t_save=1, dt=0.1, initial_n=object@initial_n, initial_n_pp=object@initial_n_pp,  ...){
        validObject(object)
        # Check that number and names of gears in effort array is same as in MizerParams object
        no_gears <- dim(object@catchability)[1]
        if(dim(effort)[2] != no_gears){
            no_gears_error_message <- paste("The number of gears in the effort array (length of the second dimension = ", dim(effort)[2], ") does not equal the number of gears in the MizerParams object (", no_gears, ").", sep="")
            stop(no_gears_error_message)
        }
        gear_names <- dimnames(object@catchability)[[1]]
        if(!all(gear_names %in% dimnames(effort)[[2]])){
            gear_names_error_message <- paste("Gear names in the MizerParams object (", paste(gear_names, collapse=", "), ") do not match those in the effort array.", sep="")
            stop(gear_names_error_message)
        }
        # Sort effort array to match order in MizerParams
        effort <- effort[,gear_names, drop=FALSE]

        # Blow up time dimension of effort array
        # i.e. effort might have been passed in using time steps of 1, but actual dt = 0.1, so need to blow up
        if (is.null(dimnames(effort)[[1]])){
            stop("The time dimname of the effort argument must be numeric.")
        }
        if (any(is.na(as.numeric(dimnames(effort)[[1]])))){
            stop("The time dimname of the effort argument must be numeric.")
        }
        time_effort <- as.numeric(dimnames(effort)[[1]])
        if (is.unsorted(time_effort)) {
            stop("The time dimname of the effort argument should be increasing.")
        }
        t_max <- time_effort[length(time_effort)]
        # Blow up effort so that rows are dt spaced
        time_effort_dt <- seq(from = time_effort[1], to = t_max, by = dt)
        effort_dt <- t(array(NA, dim = c(length(time_effort_dt), dim(effort)[2]), dimnames=list(time = time_effort_dt, dimnames(effort)[[2]])))
        for (i in 1:(length(time_effort)-1)){
            effort_dt[,time_effort_dt >= time_effort[i]] <- effort[i,]
        }
        effort_dt <- t(effort_dt)

        # Make the MizerSim object with the right size
        # We only save every t_save steps
        # Divisibility test needs to be careful about machine rounding errors,
        # see https://github.com/sizespectrum/mizer/pull/2
        if((t_save < dt) || !isTRUE(all.equal((t_save - round(t_save / dt) * dt), 0)))
            stop("t_save must be a positive multiple of dt")
        t_skip <- round(t_save/dt)
        t_dimnames_index <- seq(1, to = length(time_effort_dt), by = t_skip)
        t_dimnames <- time_effort_dt[t_dimnames_index]
        sim <- MizerSim(object, t_dimnames = t_dimnames) 
        # Fill up the effort array
        sim@effort[] <- effort_dt[t_dimnames_index,]

        # Set initial population
        sim@n[1,,] <- initial_n 
        sim@n_pp[1,] <- initial_n_pp

        # Handy things
        no_sp <- nrow(sim@params@species_params) # number of species
        no_w <- length(sim@params@w) # number of fish size bins
        idx <- 2:no_w
        # If no w_min_idx column in species_params, issue error
        if (!("w_min_idx" %in% names(sim@params@species_params))) {
            stop("w_min_idx column missing in species params")
        }
        # Hacky shortcut to access the correct element of a 2D array using 1D notation
        # This references the egg size bracket for all species, so for example
        # n[w_minidx_array_ref] = n[,w_min_idx]
        w_min_idx_array_ref <- (sim@params@species_params$w_min_idx-1) * no_sp + (1:no_sp)

        # sex ratio - DO SOMETHING LATER WITH THIS
        sex_ratio <- 0.5

        # Matrices for solver
        A <- matrix(0,nrow=no_sp,ncol=no_w)
        B <- matrix(0,nrow=no_sp,ncol=no_w)
        S <- matrix(0,nrow=no_sp,ncol=no_w)

        # initialise n and nPP
        # We want the first time step only but cannot use drop as there may only be a single species
        n <- array(sim@n[1,,],dim=dim(sim@n)[2:3])
        dimnames(n) <- dimnames(sim@n)[2:3]
        n_pp <- sim@n_pp[1,]
        t_steps <- dim(effort_dt)[1]-1
        # Set up progress bar
        pb <- progress::progress_bar$new(
            format = "[:bar] :percent ETA: :eta",
            total = t_steps, width= 60)
        for (i_time in 1:t_steps){
            # Advance progress bar
            pb$tick()
            # Do it piece by piece to save repeatedly calling methods
            # Calculate amount E_{a,i}(w) of available food
            phi_prey <- getPhiPrey(sim@params, n=n, n_pp=n_pp)
            # Calculate amount f_i(w) of food consumed
            feeding_level <- getFeedingLevel(sim@params, n=n, n_pp=n_pp, phi_prey=phi_prey)
            # Calculate the predation rate
            pred_rate <- getPredRate(sim@params, n=n, n_pp=n_pp, feeding_level=feeding_level)
            # Calculate predation mortality on fish \mu_{p,i}(w)
            #m2 <- getM2(sim@params, n=n, n_pp=n_pp, pred_rate=pred_rate)
            m2 <- getM2(sim@params, pred_rate=pred_rate)
            # Calculate total mortality \mu_i(w)
            z <- getZ(sim@params, n=n, n_pp=n_pp, effort=effort_dt[i_time,], m2=m2)
            # Calculate predation mortality on the background spectrum
            m2_background <- getM2Background(sim@params, n=n, n_pp=n_pp, pred_rate=pred_rate)
            # Calculate the resources available for reproduction and growth
            e <- getEReproAndGrowth(sim@params, n=n, n_pp=n_pp, feeding_level=feeding_level)
            # Calculate the resources for reproduction
            e_spawning <- getESpawning(sim@params, n=n, n_pp=n_pp, e=e)
            # Calculate the growth rate g_i(w)
            e_growth <- getEGrowth(sim@params, n=n, n_pp=n_pp, e_spawning=e_spawning, e=e)
            # R_{p,i}
            rdi <- getRDI(sim@params, n=n, n_pp=n_pp, e_spawning=e_spawning, sex_ratio=sex_ratio)
            # R_i
            rdd <- getRDD(sim@params, n=n, n_pp=n_pp, rdi=rdi, sex_ratio=sex_ratio)

            # Iterate species one time step forward:
            # See Ken's PDF
            # A_{ij} = - g_i(w_{j-1}) / dw_j dt
            A[,idx] <- sweep(-e_growth[,idx-1,drop=FALSE]*dt, 2, sim@params@dw[idx], "/")
            # B_{ij} = 1 + g_i(w_j) / dw_j dt + \mu_i(w_j) dt
            B[,idx] <- 1 + sweep(e_growth[,idx,drop=FALSE]*dt,2,sim@params@dw[idx],"/") + z[,idx,drop=FALSE]*dt
            # S_{ij} <- N_i(w_j)
            S[,idx] <- n[,idx,drop=FALSE]
            # Boundary condition upstream end (recruitment)
            B[w_min_idx_array_ref] <- 1+e_growth[w_min_idx_array_ref]*dt/sim@params@dw[sim@params@species_params$w_min_idx]+z[w_min_idx_array_ref]*dt
            # Update first size group of n
            n[w_min_idx_array_ref] <- (n[w_min_idx_array_ref] + rdd*dt/sim@params@dw[sim@params@species_params$w_min_idx]) / B[w_min_idx_array_ref]
            # Update n
            # for (i in 1:no_sp) # number of species assumed small, so no need to vectorize this loop over species
            #     for (j in (sim@params@species_params$w_min_idx[i]+1):no_w)
            #         n[i,j] <- (S[i,j] - A[i,j]*n[i,j-1]) / B[i,j]
            
            n <- inner_project_loop(no_sp=no_sp, no_w=no_w, n=n, A=A, B=B, S=S,
                                    w_min_idx=sim@params@species_params$w_min_idx)

            # Dynamics of background spectrum uses a semi-chemostat model (de Roos - ask Ken)
            # We use the exact solution under the assumption of constant mortality during timestep
            tmp <- (sim@params@rr_pp * sim@params@cc_pp / (sim@params@rr_pp + m2_background))
            n_pp <- tmp - (tmp - n_pp) * exp(-(sim@params@rr_pp+m2_background)*dt)

            # Store results only every t_step steps.
            store <- t_dimnames_index %in% (i_time+1)
            if (any(store)){
                sim@n[which(store),,] <- n 
                sim@n_pp[which(store),] <- n_pp
            }
        }
        # and end
        return(sim)
    }
)

#' Calculate initial population abundances for the community populations
#' 
#' This function uses the model parameters and other parameters to calculate 
#' initial population abundances for the community populations. These initial 
#' abundances should be reasonable guesses at the equilibrium values. The 
#' returned population can be passed to the \code{project} method.
#' 
#' @param params The model parameters. An object of type \code{MizerParams}.
#' @param a A parameter with a default value of 0.35.
#' @param n0_mult Multiplier for the abundance at size 0. Default value is
#'   kappa/1000.
#' @export
#' @return A matrix (species x size) of population abundances.
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' params <- MizerParams(NS_species_params_gears)
#' init_n <- get_initial_n(params)
#' }
get_initial_n<- function(params, n0_mult = NULL, a = 0.35){
    if (!is(params,"MizerParams"))
        stop("params argument must of type MizerParams")
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)
    initial_n <- array(NA, dim=c(no_sp,no_w))
    dimnames(initial_n) <- dimnames(params@intake_max)
    # N = N0 * Winf^(2*n-q-2+a) * w^(-n-a)
    # Reverse calc n and q from intake_max and search_vol slots (could add get_n as method)
    n <- (log(params@intake_max[,1] / params@species_params$h) / log(params@w[1]))[1]
    q <- (log(params@search_vol[,1] / params@species_params$gamma) / log(params@w[1]))[1]
    # Guessing at a suitable n0 value based on kappa - this was figured out using trial and error and should be updated
    if (is.null(n0_mult)){
        lambda <- 2+q-n
        kappa <- params@cc_pp[1] / (params@w_full[1]^(-lambda))
        n0_mult <- kappa / 1000
    }
    initial_n[] <- unlist(tapply(params@w,1:no_w,function(wx,n0_mult,w_inf,a,n,q)
        n0_mult * w_inf^(2*n-q-2+a) * wx^(-n-a),
        n0_mult=n0_mult, w_inf=params@species_params$w_inf, a=a, n=n, q=q))
     #set densities at w > w_inf to 0
    initial_n[unlist(tapply(params@w,1:no_w,function(wx,w_inf) w_inf<wx, w_inf=params@species_params$w_inf))] <- 0
    # Also any densities at w < w_min set to 0
    initial_n[unlist(tapply(params@w,1:no_w,function(wx,w_min)w_min>wx, w_min=params@species_params$w_min))] <- 0    
    return(initial_n)
}
