\docType{class}
\name{MizerSim-class}
\alias{MizerSim-class}
\title{MizerSim}
\description{
A class that holds the results of projecting a
\code{\link{MizerParams}} object through time.
}
\details{
\code{MizerSim} objects are created by using the
\code{\link{project}} method on an object of type
\code{MizerParams}.

There are several plotting methods available to explore the
contents of a \code{MizerSim} object. See the package
vignette for more details.
}
\section{Slots}{
  \describe{ \item{\code{params}:}{An object of type
  \code{MizerParams}.} \item{\code{n}:}{Array that stores
  the projected community population abundances by time,
  species and size} \item{\code{effort}:}{Array that stores
  the fishing effort through time by time and gear.}
  \item{\code{n_pp}:}{Array that stores the projected
  background population by time and size} }
}
\seealso{
\code{\link{project}} \code{\link{MizerParams}}
}

