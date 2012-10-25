\name{getInitialN}
\alias{getInitialN}
\title{Calculate initial population abundances for the community populations}
\usage{
  getInitialN(params, n0_mult = 1e+10, slope0 = -2/3 - 0.5)
}
\arguments{
  \item{params}{The model parameters. An object of type
  \code{MizerParams}}

  \item{n0_mult}{Multiplier for the abundance at size 0.
  Default value is 1e10}

  \item{slope0}{The estimated initial slope of the size
  spectra. Default value is -(2/3) - 0.5}
}
\description{
  This function uses the model parameters and other
  parameters to calculate an initial population abundance
  for the community populations. The returned population
  can be passed to the \code{project} method.
}
\examples{
data(species_params_gears)
params <- MizerParams(species_params_gears)
init_n <- getInitialN(params)
}

