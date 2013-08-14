\docType{methods}
\name{plotM2}
\alias{plotM2}
\alias{plotM2,MizerSim-method}
\alias{plotM2-method}
\title{Plot M2 of each species by size}
\arguments{
  \item{object}{An object of class \code{MizerSim}}

  \item{time_range}{The time range (either a vector of
  values, a vector of min and max time, or a single value)
  to average the abundances over. Default is the final time
  step.}

  \item{print_it}{Display the plot, or just return the
  ggplot2 object}
}
\value{
  A ggplot2 object
}
\description{
  After running a projection, plot M2 of each species by
  size. M2 is averaged over the specified time range (a
  single value for the time range can be used to plot a
  single time step).
}
\examples{
\dontrun{
data(species_params_gears)
data(inter)
params <- MizerParams(species_params_gears, inter)
sim <- project(params, effort=1, t_max=20, t_save = 2)
plotM2(sim)
plotM2(sim, time_range = 10:20)
}
}
\seealso{
  \code{\link{getM2}}

  \code{\link{getM2}}
}

