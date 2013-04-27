\docType{methods}
\name{getM2Background}
\alias{getM2Background}
\alias{getM2Background,MizerParams,matrix,numeric-method}
\title{getM2Background method for the size based model}
\arguments{
  \item{object}{A \code{MizerParams} object.}

  \item{n}{A matrix of species abundance (species x size).}

  \item{n_pp}{A vector of the background abundance by
  size.}
}
\value{
  A vector of predation mortalities by background prey
  size.
}
\description{
  Calculates the predation mortality on the background
  spectrum by prey size. Used by the \code{project} method
  for running size based simulations.
}
\examples{
data(species_params_gears)
data(inter)
params <- MizerParams(species_params_gears, inter)
# With constant fishing effort for all gears for 20 time steps
sim <- project(params, t_max = 20, effort = 0.5)
# Get M2 of the background spectrum at one time step
n <- sim@n[21,,]
n_pp <- sim@n_pp[21,]
getM2Background(params,n,n_pp)
}
\seealso{
  \code{\link{project}} and \code{\link{getM2}}.
}

