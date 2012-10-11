\docType{methods}
\name{plotNSpectra}
\alias{plotNSpectra}
\alias{plotNSpectra,MizerSim-method}
\title{Plot the abundance spectra of each species and the background population}
\arguments{
  \item{object}{An object of class \code{MizerSim}}

  \item{min_w}{Minimum weight to be plotted (useful for
  truncating the background spectrum). Default value is a
  tenth of the minimum size value of the community.}

  \item{time_range}{The time range (either a vector of
  values, a vector of min and max time, or a single value)
  to average the abundances over. Default is the final time
  step.}
}
\value{
  A ggplot2 object
}
\description{
  After running a projection, the spectra of the abundance
  of each species and the background population can be
  plotted. The abundance is averaged over the specified
  time range (a single value for the time range can be
  used).
}
\examples{
data(species_params_gears)
data(inter)
params <- MizerParams(species_params_gears, inter)
sim <- project(params, effort=1, t_max=20, t_save = 2)
plotNSpectra(sim)
plotNSpectra(sim, min_w = 1e-6)
plotNSpectra(sim, time_range = 10:20)
}

