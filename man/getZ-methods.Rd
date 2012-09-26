\docType{methods}
\name{getZ}
\alias{getZ}
\alias{getZ,MizerParams,matrix,numeric,numeric-method}
\title{getZ method for the size based model}
\arguments{
  \item{object}{An \code{MizerParams} object}

  \item{n}{A matrix of species abundance (species x size)}

  \item{n_pp}{A vector of the background abundance by size}

  \item{effort}{A numeric vector of the effort by gear or a
  single numeric effort value which is used for all gears}
}
\value{
  A two dimensional array (prey species x prey size)
}
\description{
  Calculates the total mortality on each prey species by
  prey size
}

