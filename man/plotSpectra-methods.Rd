\docType{methods}
\name{plotSpectra}
\alias{plotSpectra}
\alias{plotSpectra,MizerSim-method}
\title{Plot the abundance spectra of each species and the background population}
\arguments{
  \item{object}{An object of class \code{MizerSim}}

  \item{time_range}{The time range (either a vector of
  values, a vector of min and max time, or a single value)
  to average the abundances over. Default is the final time
  step.}

  \item{min_w}{Minimum weight to be plotted (useful for
  truncating the background spectrum). Default value is a
  tenth of the minimum size value of the community.}

  \item{biomass}{A boolean value. Should the biomass
  spectrum (TRUE) be plotted or the abundance in numbers
  (FALSE). Default is TRUE.}
}
\value{
  A ggplot2 object
}
\description{
  After running a projection, the spectra of the abundance
  of each species and the background population can be
  plotted. The abundance is averaged over the specified
  time range (a single value for the time range can be
  used). The abundance can be in terms of numbers or
  biomass, depending on the \code{biomass} argument.
}
\examples{
\dontrun{
data(species_params_gears)
data(inter)
params <- MizerParams(species_params_gears, inter)
sim <- project(params, effort=1, t_max=20, t_save = 2)
plotSpectra(sim)
plotSpectra(sim, min_w = 1e-6)
plotSpectra(sim, time_range = 10:20)
plotSpectra(sim, time_range = 10:20, biomass = FALSE)
}
}

