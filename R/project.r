# Project method for the size based modelling package mizer

# Copyright 2012 Finlay Scott and Julia Blanchard. Distributed under the GPL 2 or later
# Maintainer: Finlay Scott, CEFAS

# project can dispatch with effort being different classes (missing, numeric, array)
# Meaty project dispatches with effort being an array 

# Soundtrack: M.B.A disc 2 - various artists
#' project method for the size based modelling
#'
#' Runs the size-based model simulation and projects the size based model through time.
#' \code{project} is called using an object of type \code{MizerParams} and an object that contains the effort of the fishing gears through time. 
#' The method returns an object of type \code{\link{MizerSim}} which can then be explored with a range of summary and plotting methods.
#'
#' @param object A \code{MizerParams} object
#' @param effort The effort of each fishing gear through time. See notes below. 
#' @param t_max The maximum time the projection runs for params. The default value is 100. However, this argument is not needed if an array is used for the \code{effort} argument, in which case this argument is ignored. See notes below.
#' @param dt Time step of the solver. The default value is 1.
#' @param t_save The frequency with which the output is stored. The default value is 1.
#' @param initial_n The initial populations of the species. See the notes below.
#' @param initial_n_pp The initial population of the background spectrum. It should be a numeric vector of the same length as the \code{w_full} slot of the \code{MizerParams} argument. By default the \code{cc_pp} slot of the \code{\link{MizerParams}} argument is used.
#'
#' @return An object of type of \code{MizerSim}
#' @note 
#' The \code{effort} argument specifies the level of fishing effort during the simulation. It can be specified in three different ways:
#' \itemize{
#' \item A single numeric value. This specifies the effort of all fishing gears which is constant through time (i.e. all the gears have the same constant effort).
#' \item A numerical vector which has the same length as the number of fishing gears. The values in the vector specify the constant fishing effort of each of the fishing gears, i.e. the effort is constant through time but each gear may have a different fishing effort.
#' \item A numerical array with dimensions time step x gear. This specifies the fishing effort of each gear at each time step. The order of gears in the array should be same as the order of gears in the \code{MizerParams} argument.
#'}
#'
#' If effort is specified as an array then the \code{t_max} argument is ignored and the maximum simulation time is calculated using the number of time steps in the effort array and the value of the \code{dt} argument.
#'
#' The \code{initial_n} argument is a matrix with dimensions species x size. The order of species must be the same as in the \code{MizerParams} argument. If the initial population is not specified, the argument is set by default by the \code{get_initial_n} function which is set up for a North Sea model.
#' @return An object of type \code{MizerSim}.
#' @export
#' @docType methods
#' @seealso \code{\link{MizerParams}}
#' @rdname project-methods
#' @examples
#' # Data set with different fishing gears
#' data(species_params_gears)
#' data(inter)
#' params <- MizerParams(species_params_gears, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # With constant fishing effort which is different for each gear
#' sim <- project(params, t_max = 20, effort = c(1,0.5,0,2))
#' # With fishing effort that varies through time for each gear
#' fishing_effort <- array(NA,dim=c(20,4))
#' fishing_effort[,1] <- seq(from=0,to=1,length=20)
#' fishing_effort[,2] <- seq(from=2,to=1,length=20)
#' fishing_effort[,3] <- seq(from=1,to=1,length=20)
#' fishing_effort[,4] <- seq(from=1,to=0.75,length=20)
#' sim <- project(params, effort = fishing_effort)
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

#' @rdname project-methods
#' @aliases project,MizerParams,array-method
setMethod('project', signature(object='MizerParams', effort='array'),
    function(object, effort, t_save=1, dt=1, initial_n=get_initial_n(object), initial_n_pp=object@cc_pp,  ...){
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
	# If no w_min_idx column in species_params, add one
	if (!("w_min_idx" %in% names(sim@params@species_params)))
	    sim@params@species_params$w_min_idx <- 1
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
    # We want the first time step only but cannot use drop as there may only be a single species
    n <- array(sim@n[1,,],dim=dim(sim@n)[2:3])
    dimnames(n) <- dimnames(sim@n)[2:3]
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
	    A[,idx] <- sweep(-e_growth[,idx-1,drop=FALSE]*dt, 2, sim@params@dw[idx], "/")
	    B[,idx] <- 1 + sweep(e_growth[,idx,drop=FALSE]*dt,2,sim@params@dw[idx],"/") + z[,idx,drop=FALSE]*dt
	    S[,idx] <- n[,idx,drop=FALSE]
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
#' The default parameters are based on the North Sea model.
#' The returned population can be passed to the \code{project} method.
#'
#' @param params The model parameters. An object of type \code{MizerParams}
#' @param n0_mult Multiplier for the abundance at size 0. Default value is 1e10
#' @param slope0 The estimated initial slope of the size spectra. Default value is -(2/3) - 0.5
#' @export
#' @return A matrix (species x size) of population abundances.
#' @examples
#' data(species_params_gears)
#' params <- MizerParams(species_params_gears)
#' init_n <- get_initial_n(params)
get_initial_n<- function(params, n0_mult = 1e10, slope0 = -2/3 - 0.5){
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
