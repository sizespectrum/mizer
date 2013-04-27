\docType{methods}
\name{getFMort}
\alias{getFMort}
\alias{getFMort,MizerParams,matrix-method}
\alias{getFMort,MizerParams,numeric-method}
\alias{getFMort,MizerSim,missing-method}
\title{Get the total fishing mortality from all fishing gears by time, species and size}
\arguments{
  \item{object}{A \code{MizerParams} object or a
  \code{MizerSim} object}

  \item{effort}{The effort of each fishing gear. Only
  needed if the object argument is of class
  \code{MizerParams}. See notes below.}

  \item{time_range}{Subset the returned fishing mortalities
  by time. The time range is either a vector of values, a
  vector of min and max time, or a single value. Default is
  the whole time range. Only used if the \code{object}
  argument is of type \code{MizerSim}.}
}
\value{
  An array. If the effort argument has a time dimension,
  the output array has three dimensions (time x species x
  size). If the effort argument does not have a time
  dimension, the output array has two dimensions (species x
  size).
}
\description{
  Calculates the fishing mortality from all gears by
  species and size at each time step in the \code{effort}
  argument. The total fishing mortality is just the sum of
  the fishing mortalities imposed by each gear.
}
\note{
  Here: fishing mortality = catchability x selectivity x
  effort.

  The \code{effort} argument is only used if a
  \code{MizerParams} object is passed in. The \code{effort}
  argument can be a two dimensional array (time x gear), a
  vector of length equal to the number of gears (each gear
  has a different effort that is constant in time), or a
  single numeric value (each gear has the same effort that
  is constant in time). The order of gears in the
  \code{effort} argument must be the same the same as in
  the \code{MizerParams} object.

  If the object argument is of class \code{MizerSim} then
  the effort slot of the \code{MizerSim} object is used and
  the \code{effort} argument is not used.
}
\examples{
data(species_params_gears)
data(inter)
params <- MizerParams(species_params_gears, inter)
# Get the total fishing mortality when effort is constant for all gears and time:
getFMort(params, effort = 1)
# Get the total fishing mortality when effort is different between the four gears but constant in time:
getFMort(params, effort = c(0.5,1,1.5,0.75))
# Get the total fishing mortality when effort is different between the four gears and changes with time:
effort <- array(NA, dim = c(20,4))
effort[,1] <- seq(from=0, to = 1, length=20)
effort[,2] <- seq(from=1, to = 0.5, length=20)
effort[,3] <- seq(from=1, to = 2, length=20)
effort[,4] <- seq(from=2, to = 1, length=20)
getFMort(params, effort=effort)
# Get the total fishing mortality using the effort already held in a MizerSim object.
sim <- project(params, t_max = 20, effort = 0.5)
getFMort(sim)
getFMort(sim, time_range = c(10,20))
}
\seealso{
  \code{getFMortGear}, \code{project}
}

