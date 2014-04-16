\docType{methods}
\name{getN}
\alias{getN}
\alias{getN,MizerSim-method}
\alias{getN-method}
\title{Calculate the total abundance in terms of numbers of species within a size range}
\usage{
getN(object, ...)

\S4method{getN}{MizerSim}(object, ...)
}
\arguments{
  \item{object}{An object of class \code{MizerSim}.}

  \item{min_w}{minimum weight of species to be used in the
  calculation.}

  \item{max_w}{maximum weight of species to be used in the
  calculation.}

  \item{min_l}{minimum length of species to be used in the
  calculation.}

  \item{max_l}{maximum length of species to be used in the
  calculation.}
}
\value{
An array containing the total numbers (time x species)
}
\description{
Calculates the total numbers through time of the species in
the \code{MizerSim} class within user defined size limits.
The default option is to use the whole size range You can
specify minimum and maximum weight or lengths for the
species. Lengths take precedence over weights (i.e. if both
min_l and min_w are supplied, only min_l will be used)
}
\examples{
\dontrun{
data(NS_species_params_gears)
data(inter)
params <- MizerParams(NS_species_params_gears, inter)
# With constant fishing effort for all gears for 20 time steps
sim <- project(params, t_max = 20, effort = 0.5)
getN(sim)
getN(sim, min_w = 10, max_w = 1000)
}
}

