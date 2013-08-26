\docType{methods}
\name{getPhiPrey}
\alias{getPhiPrey}
\alias{getPhiPrey,MizerParams,matrix,numeric-method}
\alias{getPhiPrey-method}
\title{getPhiPrey method for the size based model}
\arguments{
  \item{object}{An \code{MizerParams} object}

  \item{n}{A matrix of species abundances (species x size)}

  \item{n_pp}{A vector of the background abundance by size}
}
\value{
  A two dimensional array (predator species x predator
  size)
}
\description{
  Calculates the amount of food exposed to each predator by
  predator size. This method is used by the
  \code{\link{project}} method for performing simulations.
}
\examples{
data(NS_species_params_gears)
data(inter)
params <- MizerParams(NS_species_params_gears, inter)
# With constant fishing effort for all gears for 20 time steps
sim <- project(params, t_max = 20, effort = 0.5)
n <- sim@n[21,,]
n_pp <- sim@n_pp[21,]
getPhiPrey(params,n,n_pp)
}
\seealso{
  \code{\link{project}}
}

