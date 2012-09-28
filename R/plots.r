# Plotting methods for the MizerSim class

# Copyright 2012 Finlay Scott and Julia Blanchard. Distributed under the GPL 2 or later
# Maintainer: Finlay Scott, CEFAS

# Biomass through time
# Pretty easy - user could do it themselves by hand for fine tuning if necessary
 
#' Plot the biomass of each species through time
#'
#' After running a projection, the biomass of each species can plotted against time. The biomass is calculated within user defined limits (see \code{\link{getBiomass}}
#' This plot is pretty to do by hand in case you want to fiddle about with colours and linetypes etc. Just look at the source code for details.
#' 
#' @param object An object of class \code{MizerSim}
#' @param min_w minimum weight of species to be used in the calculation
#' @param max_w maximum weight of species to be used in the calculation
#' @param min_l minimum length of species to be used in the calculation
#' @param max_l maximum length of species to be used in the calculation
#'
#' @return An ggplot2 object
#' @export
#' @docType methods
#' @rdname plotBiomass-methods
setGeneric('plotBiomass', function(object, ...)
    standardGeneric('plotBiomass'))
#' @rdname plotBiomass-methods
#' @aliases plotBiomass,MizerSim-method
setMethod('plotBiomass', signature(object='MizerSim'),
    function(object, ...){
	b <- getBiomass(object, ...)
	bm <- melt(b)
	p <- ggplot(bm) + geom_line(aes(x=time,y=value, colour=sp, linetype=sp)) + scale_y_continuous(trans="log10", name="Biomass") + scale_x_continuous(name="Time") #+ scale_linetype_discrete(name="Species")
	print(p)
	return(p)
    })


