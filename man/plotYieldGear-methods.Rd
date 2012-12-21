\docType{methods}
\name{plotYieldGear}
\alias{plotYieldGear}
\alias{plotYieldGear,MizerSim-method}
\title{Plot the total yield of each species by gear through time}
\arguments{
  \item{object}{An object of class \code{MizerSim}}
}
\value{
  A ggplot2 object
}
\description{
  After running a projection, the total yield of each
  species by fishing gear can be plotted against time.
}
\examples{
\dontrun{
data(species_params_gears)
data(inter)
params <- MizerParams(species_params_gears, inter)
sim <- project(params, effort=1, t_max=20, t_save = 2)
plotYieldGear(sim)
}
}

