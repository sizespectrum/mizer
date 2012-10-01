\docType{methods}
\name{plotBiomass}
\alias{plotBiomass}
\alias{plotBiomass,MizerSim-method}
\title{Plot the biomass of each species through time}
\arguments{
  \item{object}{An object of class \code{MizerSim}}

  \item{min_w}{minimum weight of species to be used in the
  calculation}

  \item{max_w}{maximum weight of species to be used in the
  calculation}

  \item{min_l}{minimum length of species to be used in the
  calculation}

  \item{max_l}{maximum length of species to be used in the
  calculation}
}
\value{
  An ggplot2 object
}
\description{
  After running a projection, the biomass of each species
  can plotted against time. The biomass is calculated
  within user defined limits (see \code{\link{getBiomass}}
  This plot is pretty to do by hand in case you want to
  fiddle about with colours and linetypes etc. Just look at
  the source code for details.
}

