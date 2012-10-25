\docType{methods}
\name{getPhiPrey}
\alias{getPhiPrey}
\alias{getPhiPrey,MizerParams,matrix,numeric-method}
\title{getPhiPrey method for the size based model}
\arguments{
  \item{object}{An \code{MizerParams} object}

  \item{n}{A matrix of species abundance (species x size)}

  \item{n_pp}{A vector of the background abundance by size}
}
\value{
  A two dimensional array (predator species x predator
  size)
}
\description{
  Calculates the amount of food exposed to each predator by
  predator size
}

