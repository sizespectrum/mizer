\docType{methods}
\name{MizerSim}
\alias{MizerSim}
\alias{MizerSim,MizerParams-method}
\alias{MizerSim-method}
\title{Constructor for the \code{MizerSim} class}
\usage{
MizerSim(object, ...)

\S4method{MizerSim}{MizerParams}(object, t_dimnames = NA, t_max = 100,
  t_save = 1, ...)
}
\arguments{
  \item{object}{a MizerParams object}

  \item{...}{other arguments including:

  \itemize{ \item{\code{t_dimnames} Numeric vector that is
  used for the time dimensions of the slots. Default = NA.}
  \item{\code{t_max} The maximum time step of the
  simulation. Only used if t_dimnames = NA. Default value =
  100.} \item{\code{t_save} How often should the results of
  the simulation be stored. Only used if t_dimnames = NA.
  Default value = 1.} }}
}
\value{
An object of type \code{MizerSim}
}
\description{
A constructor for the \code{MizerSim} class. This is used
by the \code{project} method to create \code{MizerSim}
objects of the right dimensions. It is not necessary for
users to use this constructor.
}
\examples{
\dontrun{
data(NS_species_params_gears)
data(inter)
params <- MizerParams(NS_species_params_gears, inter)
sim <- project(params)
}
}
\seealso{
\code{\link{project}} \code{\link{MizerParams}}
\code{\link{MizerSim-class}}
}

