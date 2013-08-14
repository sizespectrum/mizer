\docType{methods}
\name{project}
\alias{project}
\alias{project,MizerParams,array-method}
\alias{project,MizerParams,missing-method}
\alias{project,MizerParams,numeric-method}
\alias{project-method}
\title{project method for the size based modelling}
\arguments{
  \item{object}{A \code{MizerParams} object}

  \item{effort}{The effort of each fishing gear through
  time. See notes below.}

  \item{t_max}{The maximum time the projection runs for
  params. The default value is 100. However, this argument
  is not needed if an array is used for the \code{effort}
  argument, in which case this argument is ignored. See
  notes below.}

  \item{dt}{Time step of the solver. The default value is
  0.1.}

  \item{t_save}{The frequency with which the output is
  stored. The default value is 1.}

  \item{initial_n}{The initial populations of the species.
  See the notes below.}

  \item{initial_n_pp}{The initial population of the
  background spectrum. It should be a numeric vector of the
  same length as the \code{w_full} slot of the
  \code{MizerParams} argument. By default the \code{cc_pp}
  slot of the \code{\link{MizerParams}} argument is used.}
}
\value{
  An object of type of \code{MizerSim}

  An object of type \code{MizerSim}.
}
\description{
  Runs the size-based model simulation and projects the
  size based model through time. \code{project} is called
  using an object of type \code{MizerParams} and an object
  that contains the effort of the fishing gears through
  time. The method returns an object of type
  \code{\link{MizerSim}} which can then be explored with a
  range of summary and plotting methods.
}
\note{
  The \code{effort} argument specifies the level of fishing
  effort during the simulation. It can be specified in
  three different ways: \itemize{ \item A single numeric
  value. This specifies the effort of all fishing gears
  which is constant through time (i.e. all the gears have
  the same constant effort). \item A numerical vector which
  has the same length as the number of fishing gears. The
  values in the vector specify the constant fishing effort
  of each of the fishing gears, i.e. the effort is constant
  through time but each gear may have a different fishing
  effort. \item A numerical array with dimensions time step
  x gear. This specifies the fishing effort of each gear at
  each time step. The order of gears in the array should be
  same as the order of gears in the \code{MizerParams}
  argument. }

  If effort is specified as an array then the \code{t_max}
  argument is ignored and the maximum simulation time is
  calculated using the number of time steps in the effort
  array and the value of the \code{dt} argument.

  The \code{initial_n} argument is a matrix with dimensions
  species x size. The order of species must be the same as
  in the \code{MizerParams} argument. If the initial
  population is not specified, the argument is set by
  default by the \code{get_initial_n} function which is set
  up for a North Sea model.
}
\examples{
# Data set with different fishing gears
data(species_params_gears)
data(inter)
params <- MizerParams(species_params_gears, inter)
# With constant fishing effort for all gears for 20 time steps
sim <- project(params, t_max = 20, effort = 0.5)
# With constant fishing effort which is different for each gear
sim <- project(params, t_max = 20, effort = c(1,0.5,0,2))
# With fishing effort that varies through time for each gear
fishing_effort <- array(NA,dim=c(20,4))
fishing_effort[,1] <- seq(from=0,to=1,length=20)
fishing_effort[,2] <- seq(from=2,to=1,length=20)
fishing_effort[,3] <- seq(from=1,to=1,length=20)
fishing_effort[,4] <- seq(from=1,to=0.75,length=20)
sim <- project(params, effort = fishing_effort)
}
\seealso{
  \code{\link{MizerParams}}
}

