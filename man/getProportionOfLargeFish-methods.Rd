\docType{methods}
\name{getProportionOfLargeFish}
\alias{getProportionOfLargeFish}
\alias{getProportionOfLargeFish,MizerSim-method}
\title{Calculate the proportion of large fish}
\arguments{
  \item{object}{An object of class \code{MizerSim}}

  \item{species}{numeric or character vector of species to
  include in the calculation}

  \item{min_w}{minimum weight of species to be used in the
  calculation.}

  \item{max_w}{maximum weight of species to be used in the
  calculation.}

  \item{min_l}{minimum length of species to be used in the
  calculation.}

  \item{max_l}{maximum length of species to be used in the
  calculation.}

  \item{threshold_w}{the size used as the cutoff between
  large and small fish. Default value is 100.}

  \item{threshold_l}{the size used as the cutoff between
  large and small fish.}

  \item{biomass_proportion}{a boolean value. If TRUE the
  proportion calculated is based on biomass, if FALSE it is
  based on numbers of individuals. Default is TRUE.}
}
\value{
  An array containing the proportion of large fish through
  time
}
\description{
  Calculates the proportion of large fish through time in
  the \code{MizerSim} class within user defined size
  limits. The default option is to use the whole size
  range. You can specify minimum and maximum size ranges
  for the species and also the threshold size for large
  fish. Sizes can be expressed as weight or size. Lengths
  take precedence over weights (i.e. if both min_l and
  min_w are supplied, only min_l will be used). You can
  also specify the species to be used in the calculation.
  This method can be used to calculate the Large Fish
  Index. The proportion is based on either abundance or
  biomass.
}

