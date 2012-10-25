\docType{methods}
\name{getMeanWeight}
\alias{getMeanWeight}
\alias{getMeanWeight,MizerSim-method}
\title{Calculate the mean weight of the community}
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
  include in the calculation}
}
\value{
  A vector containing the mean weight of the community
  through time
}
\description{
  Calculates the mean weight of the community through time.
  This is simply the total biomass of the community divided
  by the abundance in numbers. You can specify minimum and
  maximum weight or length range for the species. Lengths
  take precedence over weights (i.e. if both min_l and
  min_w are supplied, only min_l will be used). You can
  also specify the species to be used in the calculation.
}

