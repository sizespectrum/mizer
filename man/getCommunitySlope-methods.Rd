\docType{methods}
\name{getCommunitySlope}
\alias{getCommunitySlope}
\alias{getCommunitySlope,MizerSim-method}
\alias{getCommunitySlope-method}
\title{Calculate the slope of the community abundance}
\arguments{
  \item{object}{An object of class \code{MizerSim}.}

  \item{species}{Numeric or character vector of species to
  include in the calculation.}

  \item{biomass}{Boolean. If TRUE (default), the abunance
  is based on biomass, if FALSE the abundance is based on
  numbers.}

  \item{min_w}{Minimum weight of species to be used in the
  calculation.}

  \item{max_w}{Maximum weight of species to be used in the
  calculation.}

  \item{min_l}{Minimum length of species to be used in the
  calculation.}

  \item{max_l}{Maximum length of species to be used in the
  calculation.}
}
\value{
  A data frame with time step, slope, intercept and R2
  values.
}
\description{
  Calculates the slope of the community abundance through
  time by performing a linear regression on the logged
  total numerical abundance at weight and logged weights
  (natural logs, not log to base 10, are used). You can
  specify minimum and maximum weight or length range for
  the species. Lengths take precedence over weights (i.e.
  if both min_l and min_w are supplied, only min_l will be
  used). You can also specify the species to be used in the
  calculation.
}
\examples{
data(species_params_gears)
data(inter)
params <- MizerParams(species_params_gears, inter)
sim <- project(params, effort=1, t_max=40, dt = 1, t_save = 1)
# Slope based on biomass, using all species and sizes
slope_biomass <- getCommunitySlope(sim)
# Slope based on numbers, using all species and sizes
slope_numbers <- getCommunitySlope(sim, biomass=FALSE)
# Slope based on biomass, using all species and sizes between 10cm and 100cm
slope_biomass <- getCommunitySlope(sim, min_l = 10, max_l = 100)
# Slope based on biomass, using only demersal species and sizes between 10cm and 100cm
dem_species <- c("Dab","Whiting","Sole","Gurnard","Plaice","Haddock", "Cod","Saithe")
slope_biomass <- getCommunitySlope(sim, species = dem_species, min_l = 10, max_l = 100)
}

