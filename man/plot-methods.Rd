\docType{methods}
\name{plot,MizerSim,missing-method}
\alias{plot,MizerSim,missing-method}
\title{Summary plot for \code{MizerSim} objects}
\usage{
\S4method{plot}{MizerSim,missing}(x, ...)
}
\arguments{
  \item{object}{An object of class \code{MizerSim}}

  \item{...}{For additional arguments see the documentation
  for \code{\link{plotBiomass}},
  \code{\link{plotFeedingLevel}},\code{\link{plotSpectra}},\code{\link{plotM2}}
  and \code{\link{plotFMort}}.}

  \item{x}{the coordinates of points in the plot.
  Alternatively, a single plotting structure, function or
  \emph{any \R object with a \code{plot} method} can be
  provided.}
}
\value{
A viewport object
}
\description{
After running a projection, produces 5 plots in the same
window: feeding level, abundance spectra, predation
mortality and fishing mortality of each species by size;
and biomass of each species through time. This method just
uses the other plotting methods and puts them all in one
window.
}
\examples{
\dontrun{
data(NS_species_params_gears)
data(inter)
params <- MizerParams(NS_species_params_gears, inter)
sim <- project(params, effort=1, t_max=20, t_save = 2)
plot(sim)
plot(sim, time_range = 10:20) # change time period for size-based plots
plot(sim, min_w = 10, max_w = 1000) # change size range for biomass plot
}
}

