\docType{methods}
\name{getPredRate}
\alias{getPredRate}
\alias{getPredRate,MizerParams,matrix,numeric-method}
\title{getPredRate method for the size based model}
\arguments{
  \item{object}{An \code{MizerParams} object}

  \item{n}{A matrix of species abundance (species x size)}

  \item{n_pp}{A vector of the background abundance by size}
}
\value{
  A three dimensional array (predator species x predator
  size x prey size)
}
\description{
  Calculates the predation rate on each prey species by
  predator and prey size
}

