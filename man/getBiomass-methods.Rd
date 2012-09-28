\docType{methods}
\name{getBiomass}
\alias{getBiomass}
\title{Calculate the biomass of species within a size range}
\arguments{
  \item{object}{An object of class \code{MizerSim}}

  \item{other}{params to set up size range (min_l, max_l,
  min_w, max_w)}

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
  An array containing the biomass (time x species)
}
\description{
  Calculates the total biomass through time of the species
  in the \code{MizerSim} class within user defined size
  limits. The default option is to use the whole size range
}

