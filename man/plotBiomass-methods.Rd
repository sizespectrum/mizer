\docType{methods}
\name{plotBiomass}
\alias{plotBiomass}
\alias{plotBiomass,MizerSim-method}
\alias{plotBiomass-method}
\title{Plot the biomass of each species through time}
\arguments{
  \item{object}{An object of class \code{MizerSim}.}

  \item{min_w}{Minimum weight of species to be used in the
  calculation.}

  \item{max_w}{Maximum weight of species to be used in the
  calculation.}

  \item{min_l}{Minimum length of species to be used in the
  calculation.}

  \item{max_l}{Maximum length of species to be used in the
  calculation.}

  \item{print_it}{Display the plot, or just return the
  ggplot2 object}
}
\value{
  A ggplot2 object
}
\description{
  After running a projection, the biomass of each species
  can be plotted against time. The biomass is calculated
  within user defined size limits (see
  \code{\link{getBiomass}}). This plot is pretty easy to do
  by hand. It just gets the biomass using the
  \code{\link{getBiomass}} method and plots using the
  ggplot2 package. You can then fiddle about with colours
  and linetypes etc. Just look at the source code for
  details.
}
\examples{
\dontrun{
data(NS_species_params_gears)
data(inter)
params <- MizerParams(NS_species_params_gears, inter)
sim <- project(params, effort=1, t_max=20, t_save = 2)
plotBiomass(sim)
plotBiomass(sim, min_l = 10, max_l = 25)
}
}
\seealso{
  \code{\link{getBiomass}}
}

