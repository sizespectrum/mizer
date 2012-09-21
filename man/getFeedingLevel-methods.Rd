\docType{methods}
\name{getFeedingLevel}
\alias{getFeedingLevel}
\alias{getFeedingLevel,MizerParams,matrix,numeric-method}
\title{getFeedingLevel method for the size based model}
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
  Calculates the amount of food consumed by a predator by
  predator size based on food availability, search volume
  and maximum intake
}

