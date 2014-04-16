\docType{methods}
\name{getMeanMaxWeight}
\alias{getMeanMaxWeight}
\alias{getMeanMaxWeight,MizerSim-method}
\alias{getMeanMaxWeight-method}
\title{Calculate the mean maximum weight of the community}
\usage{
getMeanMaxWeight(object, ...)

\S4method{getMeanMaxWeight}{MizerSim}(object,
  species = 1:nrow(object@params@species_params), measure = "both", ...)
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

  \item{species}{numeric or character vector of species to
  include in the calculation.}

  \item{measure}{The measure to return. Can be 'numbers',
  'biomass' or 'both'}
}
\value{
A matrix or vector containing the mean maximum weight of
the community through time
}
\description{
Calculates the mean maximum weight of the community through
time. This can be calculated by numbers or biomass. The
calculation is the sum of the w_inf * abundance of each
species, divided by the total abundance community, where
abundance is either in biomass or numbers. You can specify
minimum and maximum weight or length range for the species.
Lengths take precedence over weights (i.e. if both min_l
and min_w are supplied, only min_l will be used). You can
also specify the species to be used in the calculation.
}
\examples{
\dontrun{
data(NS_species_params_gears)
data(inter)
params <- MizerParams(NS_species_params_gears, inter)
sim <- project(params, effort=1, t_max=10)
getMeanMaxWeight(sim)
getMeanMaxWeight(sim, species=c("Herring","Sprat","N.pout"))
getMeanMaxWeight(sim, min_w = 10, max_w = 5000)
}
}

