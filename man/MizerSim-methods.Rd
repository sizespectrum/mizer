\docType{methods}
\name{MizerSim}
\alias{MizerSim}
\alias{MizerSim,MizerParams-method}
\title{Constructor for the \code{MizerSim} class}
\arguments{
  \item{object}{a MizerParams object}

  \item{...}{other arguments including:

  \itemize{ \item{\code{t_max} The maximum time step of the
  simulation.} \item{\code{t_save} How often should the
  results of the simulation be stored.} }}
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
data(species_params_gears)
data(inter)
params <- MizerParams(species_params_gears, inter)
sim <- project(params)
}
}
\seealso{
  \code{\link{project}} \code{\link{MizerParams}}
  \code{\link{MizerSim-class}}
}

