\docType{methods}
\name{getEForReproAndGrowth}
\alias{getEForReproAndGrowth}
\alias{getEForReproAndGrowth,MizerParams,matrix,numeric-method}
\title{getEForReproAndGrowth method for the size based model}
\arguments{
  \item{object}{An \code{MizerParams} object}

  \item{n}{A matrix of species abundance (species x size)}

  \item{n_pp}{A vector of the background abundance by size}
}
\value{
  A two dimensional array (prey species x prey size)
}
\description{
  Calculates the energy available by species and size for
  reproduction growth after metabolism and movement have
  been accounted for
}

