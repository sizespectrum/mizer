\docType{methods}
\name{summary}
\alias{summary}
\alias{summary,MizerParams-method}
\alias{summary,MizerSim-method}
\title{Summary method}
\description{
  Outputs a general summary of the structure and content of
  the object
}
\examples{
params <- MizerParams(object=3, species_names = c("cod", "haddock", "whiting"))
summary(params)
params <- MizerParams(object=3, species_names = c("cod", "haddock", "whiting"))
sim <- project(params, effort=1, t_max=5)
summary(sim)
}

