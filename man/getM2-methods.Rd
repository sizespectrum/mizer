\docType{methods}
\name{getM2}
\alias{getM2}
\alias{getM2,MizerParams,matrix,numeric-method}
\title{getM2 method for the size based model}
\arguments{
  \item{object}{An \code{MizerParams} object}

  \item{n}{A matrix of species abundance (species x size)}

  \item{n_pp}{A vector of the background abundance by size}
}
\value{
  A two dimensional array (prey species x prey size)
}
\description{
  Calculates the predation mortality on each prey species
  by prey size
}

