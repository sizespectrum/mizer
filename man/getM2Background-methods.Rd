\docType{methods}
\name{getM2Background}
\alias{getM2Background}
\alias{getM2Background,MizerParams,matrix,numeric-method}
\title{getM2Background method for the size based model}
\arguments{
  \item{object}{An \code{MizerParams} object}

  \item{n}{A matrix of species abundance (species x size)}

  \item{n_pp}{A vector of the background abundance by size}
}
\value{
  A vector (prey size)
}
\description{
  Calculates the predation mortality on the background
  spectrum by prey size
}

