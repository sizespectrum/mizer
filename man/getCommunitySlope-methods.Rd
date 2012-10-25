\docType{methods}
\name{getCommunitySlope}
\alias{getCommunitySlope}
\alias{getCommunitySlope,MizerSim-method}
\title{Calculate the slope of the community abundance}
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

  \item{species}{numeric or character vector of species to
  include in the calculation.}
}
\value{
  A data frame with time step, slope, intercept and R2
  values.
}
\description{
  Calculates the slope of the community abundance through
  time by perfomring a linear regression on the logged
  total numerical abundance at weight and logged weights.
  You can specify minimum and maximum weight or length
  range for the species. Lengths take precedence over
  weights (i.e. if both min_l and min_w are supplied, only
  min_l will be used). You can also specify the species to
  be used in the calculation.
}

