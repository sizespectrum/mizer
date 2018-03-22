#' Methods used for projecting
#'
#' The methods defined in the file project_methods calculate the various
#' quantities needed to project the size-spectra forward in time, using the
#' model described in section 3 of the mizer vignette.
#'
#' @section List of Methods:
#' In this list we relate the methods in this file to the quantities named in
#' the mizer vignette.
#' \tabular{llll}{
#'   Method name \tab Expression \tab Description \tab Section in vignette\cr
#'   \code{\link{getPhiPrey}} \tab \eqn{E_{a.i}(w)} \tab Available energy \tab 3.2 \cr
#'   \code{\link{getFeedingLevel}} \tab \eqn{f_i(w)} \tab Feeding level \tab 3.3 \cr
#'   \code{\link{getPredRate}} \tab \eqn{\phi_i(w_p/w) (1-f_i(w)) \gamma_i w^q N_i(w) dw} \tab Predation \tab 3.7 \cr
#'   \code{\link{getM2}} \tab \eqn{\mu_{p.i}(w)} \tab Predation mortality \tab 3.7 \cr
#'   \code{\link{getM2Background}} \tab \eqn{\mu_{p}(w)} \tab Predation mortality on background \tab 3.8 \cr
#'   \code{\link{getFMortGear}} \tab \eqn{F_{g,i}(w)} \tab Fishing mortality by gear \tab 8.3 \cr
#'   \code{\link{getFMort}} \tab \eqn{\mu_{f.i}(w)} \tab Total fishing mortality \tab 8.3 \cr
#'   \code{\link{getZ}} \tab \eqn{\mu_{i}(w)} \tab Total mortality \tab 3.7 \cr
#'   \code{\link{getEReproAndGrowth}} \tab \eqn{E_{r.i}(w)} \tab Energy put into growth and reproduction \tab 3.4 \cr
#'   \code{\link{getESpawning}} \tab \eqn{\psi_i(w)E_{r.i}(w)} \tab Energy put reproduction\tab 3.5 \cr
#'   \code{\link{getEGrowth}} \tab \eqn{g_i(w)} \tab Energy put growth \tab 3.4 \cr
#'   \code{\link{getRDI}} \tab \eqn{R_{p.i}} \tab Egg production \tab 3.5 \cr
#'   \code{\link{getRDD}} \tab \eqn{R_i} \tab Recruitment \tab 3.6 \cr
#' }
#'
#' @name project_methods
NULL

# Copyright 2012 Finlay Scott and Julia Blanchard.
# Distributed under the GPL 3 or later 
# Maintainer: Gustav Delius, University of York, <gustav.delius@york.ac.uk>

#' getPhiPrey method for the size based model
#' 
#' Calculates the amount \eqn{E_{a,i}(w)} of food exposed to each predator by
#' predator size. This method is used by the \code{\link{project}} method for
#' performing simulations.
#' @param object An \linkS4class{MizerParams} object
#' @param n A matrix of species abundances (species x size)
#' @param n_pp A vector of the background abundance by size
#' @param ... Other arguments (currently unused)
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
#' getPhiPrey(params,n,n_pp)
#' }
setGeneric('getPhiPrey', function(object, n, n_pp,...)
    standardGeneric('getPhiPrey'))

#' @rdname getPhiPrey
setMethod('getPhiPrey', signature(object='MizerParams', n = 'matrix', n_pp='numeric'),
    function(object, n, n_pp, ...){
#        cat("In getPhiPrey\n")
	# Check n dims
	if(dim(n)[1] != dim(object@interaction)[1])
	    stop("n does not have the right number of species (first dimension)")
	if(dim(n)[2] != length(object@w))
	    stop("n does not have the right number of size groups (second dimension)")
	if(length(n_pp) != length(object@w_full))
	    stop("n_pp does not have the right number of size groups")

	# The object@w vector only gives weights from the egg size up to the max fish size. 
	# However object@w_full gives more smaller weights and idx_sp are the index 
	# values of object@w_full such that (object@w_full)[idx_sp]=object@w
	idx_sp <- (length(object@w_full) - length(object@w) + 1):length(object@w_full)

	fishEaten <- matrix(0, nrow = dim(n)[1], ncol=length(object@w_full))
	# Looking at Equation (3.4), for available energy in the mizer vignette, 
	# we have, for our predator species i, that fishEaten[k] equals 
	# the sum over all species j of fish, of theta_{i,j}*N_j(wFull[k])        
	fishEaten[, idx_sp] <- object@interaction %*% n
	# The vector f2 equals everything inside integral (3.4) except the feeding 
	# kernel phi_i(w_p/w). 
	# We work in log-space so an extra multiplier w_p is introduced.
	f2 <- sweep(sweep(fishEaten, 2, n_pp, "+"), 2, object@w_full^2, "*")
	# Eq (3.4) is then a convolution integral in terms of f2[w_p] and phi[w_p/w].
	# We approximate the integral by the trapezoidal method. Using the
	# convolution theorem we can evaluate the resulting sum via fast fourier
	# transform.
	# mvfft() does a Fourier transform of each column of its argument, but
	# we need the Fourier transforms of each row, so we need to apply mvfft()
	# to the transposed matrices and then transpose again at the end.
	fullEnergy <- Re(t(mvfft(t(object@ft_pred_kernel_e) * mvfft(t(f2)), inverse=TRUE)))/length(object@w_full)
	# Due to numerical errors we might get negative entries. They should be 0
	fullEnergy[fullEnergy<0] <- 0

	return(fullEnergy[, idx_sp, drop=FALSE])
})


# Feeding level
# The amount of food consumed by a predator, by each predator size

#' getFeedingLevel method for the size based model
#' 
#' Calculates the amount of food \eqn{f_i(w)} consumed by a predator by predator
#' size based on food availability, search volume and maximum intake. This
#' method is used by the \code{\link{project}} method for performing
#' simulations.
#' @param object A \code{MizerParams} or \code{MizerSim} object
#' @param n A matrix of species abundance (species x size). Only used if 
#'   \code{object} argument is of type \code{MizerParams}.
#' @param n_pp A vector of the background abundance by size. Only used if 
#'   \code{object} argument is of type \code{MizerParams}.
#' @param phi_prey The PhiPrey matrix (optional) of dimension no. species x no. 
#'   size bins. If not passed in, it is calculated internally using the 
#'   \code{\link{getPhiPrey}} method. Only used if \code{object} argument is of type 
#'   \code{MizerParams}.
#' @param time_range Subset the returned fishing mortalities by time. The time 
#'   range is either a vector of values, a vector of min and max time, or a 
#'   single value. Default is the whole time range. Only used if the 
#'   \code{object} argument is of type \code{MizerSim}.
#' @param drop should extra dimensions of length 1 in the output be dropped, 
#'   simplifying the output. Defaults to TRUE.
#' @param ... Other arguments (currently unused).
#'   
#' @note If a \code{MizerParams} object is passed in, the method returns a two 
#' dimensional array (predator species x predator size) based on the abundances 
#' also passed in.
#' 
#' If a \code{MizerSim} object is passed in, the method returns a three
#' dimensional array (time step x predator species x predator size) with the
#' feeding level calculated at every time step in the simulation.
#' @seealso \code{\link{getPhiPrey}}
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
setGeneric('getFeedingLevel', function(object, n, n_pp, phi_prey, ...)
    standardGeneric('getFeedingLevel'))

#' getFeedingLevel method for a \code{MizerParams} object with already calculated \code{phi_prey} matrix.
#' @rdname getFeedingLevel
setMethod('getFeedingLevel', signature(object='MizerParams', n = 'matrix', n_pp='numeric', phi_prey='matrix'),
    function(object, n, n_pp, phi_prey, ...){
    # Check dims of phi_prey
        if (!all(dim(phi_prey) == c(nrow(object@species_params),length(object@w)))){
            stop("phi_prey argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
        }
        # encountered food = available food * search volume
        encount <- object@search_vol * phi_prey
        # calculate feeding level
        f <- encount/(encount + object@intake_max)
        return(f)
    }
)

#' getFeedingLevel method for a \code{MizerParams} object without the \code{phi_prey} matrix argument.
#' @rdname getFeedingLevel
setMethod('getFeedingLevel', signature(object='MizerParams', n = 'matrix', n_pp='numeric', phi_prey='missing'),
    function(object, n, n_pp, ...){
        phi_prey <- getPhiPrey(object, n=n, n_pp=n_pp)
        # encountered food = available food * search volume
        #encount <- object@search_vol * phi_prey
	    # calculate feeding level
        #f <- encount/(encount + object@intake_max)
        f <- getFeedingLevel(object=object, n=n, n_pp=n_pp, phi_prey=phi_prey)
	    return(f)
    }
)

#' getFeedingLevel method for a \code{MizerSim} object.
#' @rdname getFeedingLevel
setMethod('getFeedingLevel', signature(object='MizerSim', n = 'missing', 
                                       n_pp='missing', phi_prey='missing'),
    function(object, time_range=dimnames(object@n)$time, drop=FALSE, ...){
        time_elements <- get_time_elements(object,time_range)
        feed_time <- aaply(which(time_elements), 1, function(x){
            # Necessary as we only want single time step but may only have 1 species which makes using drop impossible
            n <- array(object@n[x,,],dim=dim(object@n)[2:3])
            dimnames(n) <- dimnames(object@n)[2:3]
			feed <- getFeedingLevel(object@params, n=n, n_pp = object@n_pp[x,])
			return(feed)}, .drop=drop)
        return(feed_time)
    }
)

# Predation rate
# Soundtrack: Nick Drake - Pink Moon

#' \code{getPredRate} method for the size based model
#' 
#' Calculates the predation rate of each predator species at size on prey size. 
#' In formulas \deqn{\int\phi_i(w_p/w) (1-f_i(w)) \gamma_i w^q N_i(w) dw}
#' This method is used by the \code{\link{project}} method for performing
#' simulations. In the simulations, it is combined with the interaction matrix
#' (see \code{\link{MizerParams}}) to calculate the realised predation mortality
#' (see \code{\link{getM2}}).
#' @param object A \code{MizerParams} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the background abundance by size.
#' @param feeding_level The current feeding level (optional). A matrix of size
#'   no. species x no. size bins. If not supplied, is calculated internally
#'   using the \code{getFeedingLevel()} method.
#'   
#' @return A two dimensional array (predator species x prey size), 
#'   where the prey size runs over community plus background spectrum.
#' @export
#' @seealso \code{\link{project}}, \code{\link{getM2}}, \code{\link{getFeedingLevel}} and \code{\link{MizerParams}}
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
setGeneric('getPredRate', function(object, n, n_pp, feeding_level)
    standardGeneric('getPredRate'))

#' \code{getPredRate} method with \code{feeding_level} argument.
#' @rdname getPredRate
# Called from project ->
setMethod('getPredRate', signature(object='MizerParams', n = 'matrix', 
                                   n_pp='numeric', feeding_level = 'matrix'),
    function(object, n, n_pp, feeding_level){
        if (!all(dim(feeding_level) == c(nrow(object@species_params),length(object@w)))){
            stop("feeding_level argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
        }
        noSpecies <- dim(object@interaction)[1]
        # Get indices of w_full that give w
        idx_sp <- (length(object@w_full) - length(object@w) + 1):length(object@w_full)
        # get period used in spectral integration
        no_P <- length(object@ft_pred_kernel_p[1,])
        # We express the intermediate values as a a convolution integral involving
        # two objects: Q[i,] and ft_pred_kernel_p[i,]. 
        # Here Q[i,] is all the integrand of (3.12) except the
        # feeding kernel and theta, and we sample it from 0 to P, but it is only 
        # non-zero from fishEggSize to X, where P = X + beta + 3*sigma, and X is the max 
        # fish size in the log space
        
        Q <- matrix(0, nrow = noSpecies, ncol = no_P )
        # We fill the middle of each row of Q with the proper values
        Q[, idx_sp] <- sweep((1-feeding_level)*object@search_vol*n, 2, object@w, "*")
        
        # We do our spectral integration in parallel over the different species 
        mortLonger <- Re(t(mvfft(t(object@ft_pred_kernel_p) * mvfft(t(Q)), inverse=TRUE)))/no_P
        # Unfortunately due to numerical errors some entries might be negative
        # So we have to set them to zero. Is this the fastest way to do that?
        mortLonger[mortLonger<0] <- 0
        # We drop some of the final columns to get our output
        return(mortLonger[, 1:length(object@w_full), drop = FALSE])
    }
)

#' \code{getPredRate} method without \code{feeding_level} argument.
#' @rdname getPredRate
setMethod('getPredRate', signature(object='MizerParams', n = 'matrix', n_pp='numeric', feeding_level = 'missing'),
    function(object, n, n_pp){
        feeding_level <- getFeedingLevel(object, n=n, n_pp=n_pp)
        pred_rate <- getPredRate(object=object, n=n, n_pp=n_pp, feeding_level = feeding_level)
        return(pred_rate)
    }
)

#############################################################

# getM2
# This uses the predation rate which is also used in M2background
# Too much overlap? Inefficient? Same thing is calculated twice

#' getM2 method for the size based model
#'
#' Calculates the total predation mortality rate \eqn{\mu_{p,i}(w_p)} on each prey
#' species by prey size. This method is used by the \code{\link{project}} method
#' for performing simulations.
#' @param object A \code{MizerParams} or \code{MizerSim} object.
#' @param n A matrix of species abundance (species x size). Only used if
#'   \code{object} argument is of type \code{MizerParams}.
#' @param n_pp A vector of the background abundance by size. Only used if
#'   \code{object} argument is of type \code{MizerParams}.
#' @param pred_rate An array of predation rates of dimension no. sp x no.
#'   community size bins x no. of size bins in whole spectra (i.e. community +
#'   background, the w_full slot). The array is optional. If it is not provided
#'   it is calculated by the \code{getPredRate()} method.
#' @param time_range Subset the returned fishing mortalities by time. The time
#'   range is either a vector of values, a vector of min and max time, or a
#'   single value. Default is the whole time range. Only used if the
#'   \code{object} argument is of type \code{MizerSim}.
#' @param drop Only used when object is of type \code{MizerSim}. Should
#'   dimensions of length 1 in the output be dropped, simplifying the output.
#'   Defaults to TRUE
#' @param ... Other arguments (currently unused).
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
#' # Get M2 at one time step
#' n <- sim@@n[21,,]
#' n_pp <- sim@@n_pp[21,]
#' getM2(params,n,n_pp)
#' # Get M2 at all saved time steps
#' getM2(sim)
#' # Get M2 over the time 15 - 20
#' getM2(sim, time_range = c(15,20))
#' }
setGeneric('getM2', function(object, n, n_pp, pred_rate, ...)
    standardGeneric('getM2'))

#' \code{getM2} method for \code{MizerParams} object with \code{pred_rate} argument.
#' @rdname getM2
setMethod('getM2', signature(object='MizerParams', n = 'missing', 
                             n_pp='missing', pred_rate = 'array'),
    function(object, pred_rate){
        if ((!all(dim(pred_rate) == c(nrow(object@species_params),length(object@w_full)))) | (length(dim(pred_rate))!=2)){
            stop("pred_rate argument must have 2 dimensions: no. species (",nrow(object@species_params),") x no. size bins in community + background (",length(object@w_full),")")
        }
        # Get indexes such that w_full[idx_sp[k]]==w[[k]]
        idx_sp <- (length(object@w_full) - length(object@w) + 1):length(object@w_full)
        
        m2 <- (t(object@interaction) %*% pred_rate)[, idx_sp, drop=FALSE]
        return(m2)
    }
)

#' \code{getM2} method for \code{MizerParams} object without \code{pred_rate} argument.
#' @rdname getM2
setMethod('getM2', signature(object='MizerParams', n = 'matrix', 
                             n_pp='numeric', pred_rate = 'missing'),
    function(object, n, n_pp){
      feeding_level <- getFeedingLevel(object, n=n, n_pp=n_pp)
      
      pred_rate <- getPredRate(object= object, n = n, 
                                   n_pp=n_pp, feeding_level = feeding_level)
      
      # Get indexes such that w_full[idx_sp[k]]==w[[k]]
      idx_sp <- (length(object@w_full) - length(object@w) + 1):length(object@w_full)
      
      m2 <- (t(object@interaction) %*% pred_rate)[, idx_sp, drop=FALSE]
      return(m2)
    }
)

#' \code{getM2} method for \code{MizerSim} object.
#' @rdname getM2
setMethod('getM2', signature(object='MizerSim', n = 'missing', n_pp='missing', pred_rate = 'missing'),
	function(object, time_range=dimnames(object@n)$time, drop=TRUE, ...){
		time_elements <- get_time_elements(object,time_range)
		m2_time <- aaply(which(time_elements), 1, function(x){
            n <- array(object@n[x,,],dim=dim(object@n)[2:3])
            dimnames(n) <- dimnames(object@n)[2:3]
			m2 <- getM2(object@params, n=n, n_pp = object@n_pp[x,])
			return(m2)
		}, .drop=drop)
	return(m2_time)
	}
)


#' getM2Background method for the size based model
#'
#' Calculates the predation mortality rate \eqn{\mu_p(w)} on the background spectrum
#' by prey size. Used by the \code{project} method for running size based
#' simulations.
#' @param object A \code{MizerParams} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the background abundance by size.
#' @param pred_rate An array of predation rates of dimension no. sp x no.
#'   community size bins x no. of size bins in whole spectra (i.e. community +
#'   background, the w_full slot). The array is optional. If it is not provided
#'   it is calculated by the \code{getPredRate()} method.
#'
#' @return A vector of predation mortalities by background prey size.
#' @seealso \code{\link{project}} and \code{\link{getM2}}.
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get M2 of the background spectrum at one time step
#' n <- sim@@n[21,,]
#' n_pp <- sim@@n_pp[21,]
#' getM2Background(params,n,n_pp)
#' }
setGeneric('getM2Background', function(object, n, n_pp, pred_rate)
    standardGeneric('getM2Background'))

#' \code{getM2Background} method with \code{pred_array} argument.
#' @rdname getM2Background
setMethod('getM2Background', signature(object='MizerParams', n = 'matrix', 
                                       n_pp='numeric', pred_rate='array'),
    function(object, n, n_pp, pred_rate){
        if ((!all(dim(pred_rate) == c(nrow(object@species_params),length(object@w_full)))) | (length(dim(pred_rate))!=2)){
            stop("pred_rate argument must have 2 dimensions: no. species (",nrow(object@species_params),") x no. size bins in community + background (",length(object@w_full),")")
        }
        # Since pred_rate here is the result of calling getPredRate, we have that M2background equals 
        # the sum of rows of pred_rate
        M2background <- colSums(pred_rate)
        return(M2background)
    }
)

#' \code{getM2Background} method without \code{pred_array} argument.
#' @rdname getM2Background
setMethod('getM2Background', signature(object='MizerParams', n = 'matrix', 
                                       n_pp='numeric',  pred_rate='missing'),
    function(object, n, n_pp, pred_rate){
        pred_rate <- getPredRate(object,n=n,n_pp=n_pp)
        # Since pred_rate here is the result of calling getPredRate, we have that M2background equals 
        # the sum of rows of pred_rate
        M2background <- getM2Background(object, n=n, n_pp=n_pp, pred_rate=pred_rate)
        return(M2background)
    }
)

# getFMortGear
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
#' @param ... Other arguments (currently unused).
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
setGeneric('getFMortGear', function(object, effort, ...)
    standardGeneric('getFMortGear'))

#' \code{getFMortGear} method for \code{MizerParams} object with constant effort.
#' @rdname getFMortGear
# Effort is a single value or a numeric vector.
# Effort has no time time dimension
setMethod('getFMortGear', signature(object='MizerParams', effort = 'numeric'),
	function(object, effort, ...){
        no_gear <- dim(object@catchability)[1]
        # If a single value, just repeat it for all gears
        if(length(effort) == 1)
            effort <- rep(effort, no_gear)
        if (length(effort) != no_gear)
            stop("Effort must be a single value or a vector as long as the number of gears\n")
        # Streamlined for speed increase - note use of recycling
        out <- object@selectivity
        out[] <- effort * c(object@catchability) * c(object@selectivity)
        return(out)
    }
)

#' \code{getFMortGear} method for \code{MizerParams} object with time changing effort.
#' @rdname getFMortGear
# Always returns a 4D array: time x gear x species x size
setMethod('getFMortGear', signature(object='MizerParams', effort = 'matrix'),
    function(object, effort, ...){
	no_gear <- dim(object@catchability)[1]
	if (dim(effort)[2] != no_gear)
	    stop("Effort array must have a single value or a vector as long as the number of gears for each time step\n")
    # Make the output array - note that we put time as last dimension and then aperm before returning 
    # This is because of the order of the values when we call the other getFMortGear method
    # Fill it up with by calling the other method and passing in each line of the effort matrix
    out <- array(NA, dim=c(dim(object@selectivity), dim(effort)[1]), dimnames= c(dimnames(object@selectivity), list(time = dimnames(effort)[[1]])))
    out[] <- apply(effort, 1, function(x) getFMortGear(object, x))
    out <- aperm(out, c(4,1,2,3))
    return(out)
    }
)

# Returns the fishing mortality: time * gear * species * size
#' \code{getFMortGear} method for \code{MizerSim} object.
#' @rdname getFMortGear
setMethod('getFMortGear', signature(object='MizerSim', effort='missing'),
    function(object,effort, time_range=dimnames(object@effort)$time, ...){
        time_elements <- get_time_elements(object,time_range, slot_name="effort")
        f_mort_gear <- getFMortGear(object@params, object@effort, ...)
        return(f_mort_gear[time_elements,,,,drop=FALSE])
})

# Total fishing mortality from all gears
# species x size and maybe also by time if effort is time based

#' Get the total fishing mortality rate from all fishing gears by time, species and
#' size.
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
#' @param ... Other arguments passed to \code{getFMortGear} method.
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
#' If the object argument is of class \code{MizerSim} then the effort slot of the \code{MizerSim} object is used and the \code{effort} argument is not used.
#' @export
#' @seealso \code{getFMortGear}, \code{project}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # Get the total fishing mortality when effort is constant for all gears and time:
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
#' # Get the total fishing mortality using the effort already held in a MizerSim
#' object.
#' sim <- project(params, t_max = 20, effort = 0.5)
#' getFMort(sim)
#' getFMort(sim, time_range = c(10,20))
#' }
setGeneric('getFMort', function(object, effort, ...)
    standardGeneric('getFMort'))

#' \code{getFMort} method for \code{MizerParams} object with constant effort.
#' @rdname getFMort
# Called from project -> getZ -> 
setMethod('getFMort', signature(object='MizerParams', effort='numeric'),
    function(object, effort, ...){
	fMortGear <- getFMortGear(object, effort, ...)
	fMort <- colSums(fMortGear)
	return(fMort)
})

#' \code{getFMort} method for \code{MizerParams} object with time changing effort.
#' @rdname getFMort
setMethod('getFMort', signature(object='MizerParams', effort='matrix'),
    function(object, effort, ...){
	fMortGear <- getFMortGear(object, effort, ...)
	fMort <- apply(fMortGear, c(1,3,4), sum)
	return(fMort)
})

#' \code{getFMort} method for \code{MizerSim} object.
#' @rdname getFMort
setMethod('getFMort', signature(object='MizerSim', effort='missing'),
    function(object, effort, time_range=dimnames(object@effort)$time, drop=TRUE, ...){
    	time_elements <- get_time_elements(object,time_range, slot_name="effort")
    	fMort <- getFMort(object@params, object@effort, ...)
    	return(fMort[time_elements,,,drop=drop])
    }
)


# get total Z
#' getZ method for the size based model
#'
#' Calculates the total mortality rate \eqn{\mu_i(w)} on each species by size from
#' predation mortality (M2), background mortality (M) and fishing mortality for
#' a single time step.
#' @param object A \code{MizerParams} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the background abundance by size.
#' @param effort A numeric vector of the effort by gear or a single numeric
#'   effort value which is used for all gears.
#' @param m2 A two dimensional array of predation mortality (optional). Has
#'   dimensions no. sp x no. size bins in the community. If not supplied is
#'   calculated using the \code{getM2()} method.
#'
#' @return A two dimensional array (prey species x prey size). 
#'
#' @export
#' @seealso \code{\link{getM2}}, \code{\link{getFMort}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the total mortality at a particular time step
#' getZ(params,sim@@n[21,,],sim@@n_pp[21,],effort=0.5)
#' }
setGeneric('getZ', function(object, n, n_pp, effort, m2)
    standardGeneric('getZ'))

#' \code{getZ} method with \code{m2} argument.
#' @rdname getZ
# Called from project()
setMethod('getZ', signature(object='MizerParams', n = 'matrix', 
                            n_pp = 'numeric', effort='numeric', m2 = 'matrix'),
    function(object, n, n_pp, effort, m2){
        if (!all(dim(m2) == c(nrow(object@species_params),length(object@w)))){
            stop("m2 argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
        }
        return(m2 + object@mu_b + getFMort(object, effort = effort))
    }
)

#' \code{getZ} method without \code{m2} argument.
#' @rdname getZ
setMethod('getZ', signature(object='MizerParams', n = 'matrix', 
                            n_pp = 'numeric', effort='numeric', m2 = 'missing'),
    function(object, n, n_pp, effort){
        m2 <- getM2(object, n=n, n_pp=n_pp)
        z <- getZ(object, n=n, n_pp=n_pp, effort=effort, m2=m2)
        return(z)
    }
)


# Energy after metabolism and movement
#' getEReproAndGrowth method for the size based model
#'
#' Calculates the energy rate available by species and size for reproduction and
#' growth after metabolism and movement have been accounted for: \eqn{E_{r.i}(w)}.
#' Used by the \code{project} method for performing simulations.
#' @param object A \code{MizerParams} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the background abundance by size.
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
setGeneric('getEReproAndGrowth', function(object, n, n_pp, feeding_level)
    standardGeneric('getEReproAndGrowth'))

#' \code{getEReproAndGrowth} method with \code{feeding_level} argument.
#' @rdname getEReproAndGrowth
setMethod('getEReproAndGrowth', signature(object='MizerParams', n = 'matrix', 
                                          n_pp = 'numeric', feeding_level='matrix'),
    function(object, n, n_pp, feeding_level){
        if (!all(dim(feeding_level) == c(nrow(object@species_params),length(object@w)))){
            stop("feeding_level argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
        }
        # assimilated intake
        e <- sweep(feeding_level * object@intake_max,1,object@species_params$alpha,"*", check.margin=FALSE)
        # Subtract basal metabolism and activity 
        e <- e - object@std_metab - object@activity
        e[e<0] <- 0 # Do not allow negative growth
        return(e)
    }
)

#' \code{getEReproAndGrowth} method without \code{feeding_level} argument.
#' @rdname getEReproAndGrowth
setMethod('getEReproAndGrowth', signature(object='MizerParams', n = 'matrix', 
                                          n_pp = 'numeric', feeding_level='missing'),
    function(object, n, n_pp){
        feeding_level <- getFeedingLevel(object, n=n, n_pp=n_pp)
        e <- getEReproAndGrowth(object, n=n, n_pp=n_pp, feeding_level=feeding_level)
        return(e)
    }
)

# Energy left for reproduction
# assimilated food intake, less metabolism and activity, split between reproduction and growth

#' getESpawning method for the size based model
#'
#' Calculates the energy rate available by species and size for reproduction after
#' metabolism and movement have been accounted for: \eqn{\psi_i(w)E_{r.i}(w)}.
#' Used by the \code{project} method for performing simulations.
#' @param object A \code{MizerParams} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the background abundance by size.
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
#' getESpawning(params,sim@@n[21,,],sim@@n_pp[21,])
#' }
setGeneric('getESpawning', function(object, n, n_pp, e)
    standardGeneric('getESpawning'))

#' \code{getESpawning} method with \code{e} argument.
#' @rdname getESpawning
setMethod('getESpawning', signature(object='MizerParams', n = 'matrix', 
                                    n_pp = 'numeric', e = 'matrix'),
    function(object, n, n_pp, e){
        if (!all(dim(e) == c(nrow(object@species_params),length(object@w)))){
            stop("e argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
        }
        e_spawning <- object@psi * e 
        return(e_spawning)
    }
)

#' \code{getESpawning} method without \code{e} argument.
#' @rdname getESpawning
setMethod('getESpawning', signature(object='MizerParams', n = 'matrix', 
                                    n_pp = 'numeric', e = 'missing'),
    function(object, n, n_pp){
	    e <- getEReproAndGrowth(object,n=n,n_pp=n_pp)
        e_spawning <- getESpawning(object, n=n, n_pp=n_pp, e=e)
	    return(e_spawning)
    }
)

#' getEGrowth method for the size based model
#'
#' Calculates the energy rate \eqn{g_i(w)} available by species and size for growth
#' after metabolism, movement and reproduction have been accounted for. Used by
#' the \code{\link{project}} method for performing simulations.
#' @param object A \linkS4class{MizerParams} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the background abundance by size.
#' @param e The energy available for reproduction and growth (optional, although
#'   if specified, e_spawning must also be specified). A matrix of size no.
#'   species x no. size bins. If not supplied, is calculated internally using
#'   the \code{\link{getEReproAndGrowth}} method.
#' @param e_spawning The energy available for spawning (optional, although if
#'   specified, e must also be specified). A matrix of size no. species x no.
#'   size bins. If not supplied, is calculated internally using the
#'   \code{\link{getESpawning}} method.
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
setGeneric('getEGrowth', function(object, n, n_pp, e_spawning, e)
    standardGeneric('getEGrowth'))

#' \code{getEGrowth} method with \code{e_spawning} and \code{e} arguments.
#' @rdname getEGrowth 
setMethod('getEGrowth', signature(object='MizerParams', n = 'matrix', 
                                  n_pp = 'numeric', e_spawning='matrix', e='matrix'),
    function(object, n, n_pp, e_spawning, e){
        if (!all(dim(e_spawning) == c(nrow(object@species_params),length(object@w)))){
            stop("e_spawning argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
        }
        if (!all(dim(e) == c(nrow(object@species_params),length(object@w)))){
            stop("e argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
        }
        # Assimilated intake less activity and metabolism
        # energy for growth is intake - energy for growth
        e_growth <- e - e_spawning
        return(e_growth)
    }
)

#' \code{getEGrowth} method without \code{e_spawning} and \code{e} arguments.
#' @rdname getEGrowth
setMethod('getEGrowth', signature(object='MizerParams', n = 'matrix', 
                                  n_pp = 'numeric', e_spawning='missing', e='missing'),
    function(object, n, n_pp){
        # Assimilated intake less activity and metabolism
        e <- getEReproAndGrowth(object,n=n,n_pp=n_pp)
        e_spawning <- getESpawning(object,n=n,n_pp=n_pp)
        # energy for growth is intake - energy for growth
        e_growth <- getEGrowth(object, n=n, n_pp=n_pp, e_spawning=e_spawning, e=e)
        return(e_growth)
    }
)

#' getRDI method for the size based model
#'
#' Calculates the density independent recruitment (total egg production)
#' \eqn{R_{p.i}} before density dependence, by species. Used by the
#' \code{project} method for performing simulations.
#' @param object A \code{MizerParams} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the background abundance by size.
#' @param e_spawning The energy available for spawning (optional). A matrix of
#'   size no. species x no. size bins. If not supplied, is calculated internally
#'   using the \code{\link{getESpawning}} method.
#' @param sex_ratio Proportion of the population that is female. Default value
#'   is 0.5.
#' @param ... Other arguments (currently unused).
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
setGeneric('getRDI', function(object, n, n_pp, e_spawning, ...)
    standardGeneric('getRDI'))

#' \code{getRDI} method with \code{e_spawning} argument.
#' @rdname getRDI
setMethod('getRDI', signature(object='MizerParams', n = 'matrix', 
                              n_pp = 'numeric', e_spawning='matrix'),
    function(object, n, n_pp, e_spawning, sex_ratio = 0.5){
        if (!all(dim(e_spawning) == c(nrow(object@species_params),length(object@w)))){
            stop("e_spawning argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
        }
        e_spawning_pop <- (e_spawning*n) %*% object@dw
        rdi <- sex_ratio*(e_spawning_pop * object@species_params$erepro)/object@w[object@species_params$w_min_idx] 
        return(rdi)
    }
)

#' \code{getRDI} method without \code{e_spawning} argument.
#' @rdname getRDI
setMethod('getRDI', signature(object='MizerParams', n = 'matrix', 
                              n_pp = 'numeric', e_spawning='missing'),
    function(object, n, n_pp, sex_ratio = 0.5){
        e_spawning <- getESpawning(object, n=n, n_pp=n_pp)
        rdi <- getRDI(object, n=n, n_pp=n_pp, e_spawning=e_spawning, sex_ratio=sex_ratio)
        return(rdi)
    }
)

#' getRDD method for the size based model
#'
#' Calculates the density dependent recruitment (total egg production) \eqn{R_i}
#' for each species. This is the flux entering the smallest size class of each
#' species. The density dependent recruitment is the density independent
#' recruitment after it has been put through the density dependent
#' stock-recruitment relationship function. This method is used by the
#' \code{project} method for performing simulations.
#' @param object An \code{MizerParams} object
#' @param n A matrix of species abundance (species x size)
#' @param n_pp A vector of the background abundance by size
#' @param rdi A matrix of density independent recruitment (optional) with
#'   dimensions no. sp x 1. If not specified rdi is calculated internally using
#'   the \code{\link{getRDI}} method.
#' @param sex_ratio Proportion of the population that is female. Default value
#'   is 0.5
#' @param ... Other arguments (currently unused).
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
setGeneric('getRDD', function(object, n, n_pp, rdi, ...)
    standardGeneric('getRDD'))

#' \code{getRDD} method with \code{rdi} argument.
#' @rdname getRDD
setMethod('getRDD', signature(object='MizerParams', n = 'matrix', 
                              n_pp = 'numeric', rdi='matrix'),
    function(object, n, n_pp, rdi, sex_ratio = 0.5){
        if (!all(dim(rdi) == c(nrow(object@species_params),1))){
            stop("rdi argument must have dimensions: no. species (",nrow(object@species_params),") x 1")
        }
        rdd <- object@srr(rdi = rdi, species_params = object@species_params)
        return(rdd)
})

#' \code{getRDD} method without \code{rdi} argument.
#' @rdname getRDD
setMethod('getRDD', signature(object='MizerParams', n = 'matrix', n_pp = 'numeric', rdi='missing'),
    function(object, n, n_pp, sex_ratio = 0.5){
    	rdi <- getRDI(object, n=n, n_pp=n_pp, sex_ratio = sex_ratio)
    	rdd <- getRDD(object, n=n, n_pp=n_pp, rdi=rdi, sex_ratio=sex_ratio)
    	return(rdd)
    }
)


# get_time_elements
# internal function to get the array element references of the time dimension
# for the time based slots of a MizerSim object
# time_range can be character or numeric
# Necessary to include a slot_name argument because the effort and abundance
# slots have different time dimensions
get_time_elements <- function(sim,time_range,slot_name="n"){
    if (!is(sim,"MizerSim"))
        stop("First argument to get_time_elements function must be of class MizerSim")
    time_range <- range(as.numeric(time_range))
    # Check that time range is even in object
    sim_times <- as.numeric(dimnames(sim@effort)$time)
    sim_time_range <- range(sim_times)
    if ((time_range[1] < sim_time_range[1]) | (time_range[2] > sim_time_range[2]))
	    stop("Time range is outside the time range of the model")
    time_elements <- (sim_times >= time_range[1]) & (sim_times <= time_range[2])
    names(time_elements) <- dimnames(sim@effort)$time
    return(time_elements)
}

