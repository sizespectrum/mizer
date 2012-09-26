\docType{methods}
\name{getESpawning}
\alias{getESpawning}
\alias{getESpawning,MizerParams,matrix,numeric-method}
\title{getESpawning method for the size based model}
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
  reproduction after metabolism and movement have been
  accounted for
}

