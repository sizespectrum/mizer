\docType{methods}
\name{plotFeedingLevel}
\alias{plotFeedingLevel}
\alias{plotFeedingLevel,MizerSim-method}
\title{Plot the feeding level of each species by size}
\arguments{
  \item{object}{An object of class \code{MizerSim}.}

  \item{time_range}{The time range (either a vector of
  values, a vector of min and max time, or a single value)
  to average the abundances over. Default is the final time
  step.}
}
\value{
  A ggplot2 object
}
\description{
  After running a projection, plot the feeding level of
  each species by size. The feeding level is averaged over
  the specified time range (a single value for the time
  range can be used).
}
\examples{
\dontrun{
data(species_params_gears)
data(inter)
params <- MizerParams(species_params_gears, inter)
sim <- project(params, effort=1, t_max=20, t_save = 2)
plotFeedingLevel(sim)
plotFeedingLevel(sim, time_range = 10:20)
}
}
\seealso{
  \code{\link{getFeedingLevel}}
}

