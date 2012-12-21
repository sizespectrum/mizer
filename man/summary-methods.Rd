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
data(species_params_gears)
data(inter)
params <- MizerParams(species_params_gears,inter)
summary(params)
data(species_params_gears)
data(inter)
params <- MizerParams(species_params_gears,inter)
sim <- project(params, effort=1, t_max=5)
summary(sim)
}

