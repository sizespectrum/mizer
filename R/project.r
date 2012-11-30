# Project method for the size based modelling package mizer

# Copyright 2012 Finlay Scott and Julia Blanchard. Distributed under the GPL 2 or later
# Maintainer: Finlay Scott, CEFAS

# project can dispatch with effort being different classes
# But meaty project dispatches with effort being an array (or data.frame?)

#' project method for the size based modelling
#'
#' Projects the size based model through time.
#' \code{project()} is called using an object of type \code{MizerParams} and an object that stores the effort of the fishing gears through time. 
#'
#' @param object A \code{MizerParams} object
#' @param effort The effort of each fishing gear through time. Can be an array (time x gear), a vector (of length ngears where each gear has a constant effort through time), or a single numeric value (each gear has the same effort, constant through time).
#' @param t_max The maximum time the projection runs for params. Not needed if an array is used for the \code{effort} argument.
#' @param dt Time step of the solver
#' @param t_save Store the output at every \code{t_save} 
#' @param initial_n the initial populations of the species
#' @param initial_n_pp the initial population of the background spectrum 
#'
#' @return An object of type of \code{MizerSim}
#' @note details and whatnot about setting up t from effort
#' @export
#' @docType methods
#' @rdname project-methods
setGeneric('project', function(object, effort, ...)
    standardGeneric('project'))

# No effort is specified - default is to set an effort of 1
# All other arguments passed as ...

#' @rdname project-methods
#' @aliases project,MizerParams,missing-method
setMethod('project', signature(object='MizerParams', effort='missing'),
    function(object, ...){
	res <- project(object, effort=1, ...)
	return(res)
})

#' @rdname project-methods
#' @aliases project,MizerParams,numeric-method
setMethod('project', signature(object='MizerParams', effort='numeric'),
    function(object, effort,  t_max = 100, dt = 1, ...){
	if (!all.equal(t_max %% dt, 0))
	    stop("t_max must be divisible by dt with no remainder")
	no_gears <- dim(object@catchability)[1]
	if ((length(effort)>1) & (length(effort) != no_gears))
	    stop("Effort vector must be the same length as the number of fishing gears\n")
	# Set it up transposed so we can use the recycling rules
	effort_array <- t(array(effort, dim=c(no_gears,t_max / dt), dimnames=list(gear=dimnames(object@catchability)$gear,t_step=signif(seq(from=dt,to=t_max,by=dt),3))))
	res <- project(object,effort_array, dt=dt, ...)
	return(res)
})

# n is included an argument to help set initial abundances
# Makes sense for user to use the same one used in sbmParams but as it's only for setting initial pop, it doesn't matter

#' @rdname project-methods
#' @aliases project,MizerParams,array-method
setMethod('project', signature(object='MizerParams', effort='array'),
    function(object, effort, t_save=1, dt=1, initial_n=getInitialN(object), initial_n_pp=object@cc_pp,  ...){
	validObject(object)
	#args <- list(...)
	# Use the effort array to pull out time info
	t_steps <- dim(effort)[1]
	t_max <- t_steps * dt
	if (!all.equal(t_save %% dt, 0))
	    stop("t_save must be divisible by dt with no remainder")
	# Make the sbmSim object
	sim <- MizerSim(object, t_max = t_max, t_save=t_save)
	# Set the effort slot - the effort slot only stores effort every tSave
	effort_ref <- seq(from = t_save / dt, to = t_max/dt, by = t_save/dt)
	sim@effort[] <- effort[effort_ref,]

	# Set initial population
	sim@n[1,,] <- initial_n 
	sim@n_pp[1,] <- initial_n_pp

	# Handy things
	no_sp <- nrow(sim@params@species_params)
	no_w <- length(sim@params@w)
	idx <- 2:no_w
	# Hacky shortcut to access the correct element of a 2D array using 1D notation
	w_min_idx_array_ref <- (sim@params@species_params$w_min_idx-1) * no_sp + (1:no_sp)

	# sex ratio - DO SOMETHING LATER WITH THIS
	sex_ratio <- 0.5

	# Matrices for solver
	# Dynamics of background spectrum uses a semi-chemostat model (de Roos - ask Ken)
	A <- matrix(0,nrow=no_sp,ncol=no_w)
	B <- matrix(0,nrow=no_sp,ncol=no_w)
	S <- matrix(0,nrow=no_sp,ncol=no_w)

	# initialise n and nPP
	n <- sim@n[1,,]
	n_pp <- sim@n_pp[1,]

	for (i_time in 1:t_steps){

	    # Take this out. No n or npp so fmort could be precalculated
	    # Also f_mort is already called in z so not needed at all?
	    #f_mort <- getFMort(sim@params, effort[i_time,])
	    e_growth <- getEGrowth(sim@params, n=n, n_pp=n_pp)
	    z <- getZ(sim@params, n=n, n_pp=n_pp, effort=effort[i_time,])
	    rdd <- getRDD(sim@params, n=n, n_pp=n_pp, sex_ratio=sex_ratio)
	    m2_background <- getM2Background(sim@params, n=n, n_pp=n_pp)

	    # Iterate species one time step forward:
	    # See Ken's PDF
	    A[,idx] <- sweep(-e_growth[,idx-1]*dt, 2, sim@params@dw[idx], "/")
	    B[,idx] <- 1 + sweep(e_growth[,idx]*dt,2,sim@params@dw[idx],"/") + z[,idx]*dt
	    S[,idx] <- n[,idx]
	    # Boundary condition upstream end (recruitment)
	    B[w_min_idx_array_ref] <- 1+e_growth[w_min_idx_array_ref]*dt/sim@params@dw[sim@params@species_params$w_min_idx]+z[w_min_idx_array_ref]*dt
	    # Update first size group of n
	    n[w_min_idx_array_ref] <- (n[w_min_idx_array_ref] + rdd*dt/sim@params@dw[sim@params@species_params$w_min_idx]) / B[w_min_idx_array_ref]
	    # Invert matrix
	    for (i in 1:no_sp)
	        for (j in (sim@params@species_params$w_min_idx[i]+1):no_w)
	            n[i,j] <- (S[i,j] - A[i,j]*n[i,j-1]) / B[i,j]

	    # Dynamics of background spectrum uses a semi-chemostat model (de Roos - ask Ken)
	    tmp <- (sim@params@rr_pp * sim@params@cc_pp / (sim@params@rr_pp + m2_background))
	    n_pp <- tmp - (tmp - n_pp) * exp(-(sim@params@rr_pp+m2_background)*dt)
	    # save n and F if necessary
	    if (((i_time*dt) %% t_save)==0){
	        sim@n[(i_time*dt) / t_save +1,,] <- n
	        sim@n_pp[(i_time*dt) / t_save +1,] <- n_pp
	    }
	}
	# and end
	return(sim)
    }
)

#' Calculate initial population abundances for the community populations
#'
#' This function uses the model parameters and other parameters to calculate an initial population abundance for the
#' community populations. 
#' The returned population can be passed to the \code{project} method.
#'
#' @param params The model parameters. An object of type \code{MizerParams}
#' @param n0_mult Multiplier for the abundance at size 0. Default value is 1e10
#' @param slope0 The estimated initial slope of the size spectra. Default value is -(2/3) - 0.5
#' @export
#' @examples
#' data(species_params_gears)
#' params <- MizerParams(species_params_gears)
#' init_n <- getInitialN(params)
getInitialN<- function(params, n0_mult = 1e10, slope0 = -2/3 - 0.5){
    if (!is(params,"MizerParams"))
	stop("params argument must of type MizerParams")
#    if (slope0 > 0)
#	stop("slope0 must be less than or equal to 0")
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)

    n0 <- params@species_params$w_inf^(0.7-3+0.5+1)*n0_mult # From original NS model. Ask Ken!
    initial_n <- array(NA, dim=c(no_sp,no_w))
    dimnames(initial_n) <- dimnames(params@intake_max)
    initial_n[] <- unlist(tapply(params@w,1:no_w,function(wx,n0,w0,slope0)
	n0 * (wx/w0)^slope0, n0=n0, w0=min(params@w), slope0=slope0))
    # set densities at w > w_inf to 0
    initial_n[unlist(tapply(params@w,1:no_w,function(wx,w_inf) w_inf<wx, w_inf=params@species_params$w_inf))] <- 0
    # Also any densities at w < w_min set to 0
    initial_n[unlist(tapply(params@w,1:no_w,function(wx,w_min)w_min>wx, w_min=params@species_params$w_min))] <- 0    
    return(initial_n)
}
