\docType{methods}
\name{plotYield}
\alias{plotYield}
\alias{plotYield,MizerSim-method}
\title{Plot the total yield of each species through time}
\arguments{
  \item{object}{An object of class \code{MizerSim}}
}
\value{
  A ggplot2 object
}
\description{
  After running a projection, the total yield of each
  species across all fishing gears can be plotted against
  time.
}
\examples{
\dontrun{
data(species_params_gears)
data(inter)
params <- MizerParams(species_params_gears, inter)
sim <- project(params, effort=1, t_max=20, t_save = 2)
plotYield(sim)
}
}

