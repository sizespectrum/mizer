\docType{methods}
\name{project}
\alias{project}
\alias{project,MizerParams,array-method}
\alias{project,MizerParams,missing-method}
\alias{project,MizerParams,numeric-method}
\title{project method for the size based modelling}
\arguments{
  \item{object}{A \code{MizerParams} object}

  \item{effort}{The effort of each fishing gear through
  time. Can be an array (time x gear), a vector (each gear
  has a constant effort through time), or a single numeric
  value (each gear has the same effort, constant through
  time).}

  \item{t_max}{The maximum time the projection runs for
  params. Not needed if an array is used for the
  \code{effort} argument.}

  \item{dt}{Time step of the solver}

  \item{t_save}{Store the output at every \code{t_save}}
}
\value{
  An object of type of \code{MizerSim}
}
\description{
  Projects the size based model through time.
  \code{project()} is called using an object of type
  \code{MizerParams} and an object that stores the effort
  of the fishing gears through time.
}
\note{
  details and whatnot about setting up t from effort
}

