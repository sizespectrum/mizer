\docType{methods}
\name{getRDD}
\alias{getRDD}
\alias{getRDD,MizerParams,matrix,numeric-method}
\title{getRDD method for the size based model}
\arguments{
  \item{object}{An \code{MizerParams} object}

  \item{n}{A matrix of species abundance (species x size)}

  \item{n_pp}{A vector of the background abundance by size}

  \item{sex_ratio}{Proportion of the population that is
  female. Default value is 0.5}
}
\value{
  A numeric vector the length of the number of species
}
\description{
  Calculates the density dependent recruitment (flux
  entering the smallest size class of each species) for
  each species. The density independent recruitment after
  it has been put through the density dependence function.
}

