\docType{methods}
\name{plot}
\alias{plot}
\alias{plot,MizerSim,missing-method}
\title{Summary plot for \class{MizerObject}s}
\arguments{
  \item{object}{An object of class \code{MizerSim}}

  \item{...}{For additional arguments see the documentation
  for \code{\link{plotBiomass}},
  \code{\link{plotFeedingLevel}},\code{\link{plotSpectra}},\code{\link{plotM2}}
  and \code{\link{plotFMort}}.}
}
\value{
  A viewport object
}
\description{
  After running a projection, produces 5 plots in the same
  window: feeding level, abundance spectra, predation
  mortality and fishing mortality of each species by size;
  and biomass of each species through time. This method
  just uses the other plotting methods and puts them all in
  one window.
}
\examples{
\dontrun{
data(species_params_gears)
data(inter)
params <- MizerParams(species_params_gears, inter)
sim <- project(params, effort=1, t_max=20, t_save = 2)
plot(sim)
plot(sim, time_range = 10:20) # change time period for size-based plots
plot(sim, min_l = 10, max_l = 25) # change size range for biomass plot
}
}

