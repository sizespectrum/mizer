\name{get_initial_n}
\alias{get_initial_n}
\title{Calculate initial population abundances for the community populations}
\usage{
  get_initial_n(params, n0_mult = 1e+10, a = 0.35,
    slope0 = -2/3 - 0.5)
}
\arguments{
  \item{params}{The model parameters. An object of type
  \code{MizerParams}}

  \item{n0_mult}{Multiplier for the abundance at size 0.
  Default value is 1e10}

  \item{slope0}{The estimated initial slope of the size
  spectra. Default value is -(2/3) - 0.5}
}
\value{
  A matrix (species x size) of population abundances.
}
\description{
  This function uses the model parameters and other
  parameters to calculate an initial population abundance
  for the community populations. The default parameters are
  based on the North Sea model. The returned population can
  be passed to the \code{project} method.
}
\examples{
data(species_params_gears)
params <- MizerParams(species_params_gears)
init_n <- get_initial_n(params)
}

