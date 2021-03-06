% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plotBiomass}
\alias{plotBiomass}
\alias{plotlyBiomass}
\title{Plot the biomass of species through time}
\usage{
plotBiomass(
  sim,
  species = NULL,
  start_time,
  end_time,
  y_ticks = 6,
  ylim = c(NA, NA),
  total = FALSE,
  background = TRUE,
  highlight = NULL,
  return_data = FALSE,
  ...
)

plotlyBiomass(
  sim,
  species = NULL,
  start_time,
  end_time,
  y_ticks = 6,
  ylim = c(NA, NA),
  total = FALSE,
  background = TRUE,
  highlight = NULL,
  ...
)
}
\arguments{
\item{sim}{An object of class \linkS4class{MizerSim}}

\item{species}{The species to be selected. Optional. By default all target
species are selected. A vector of species names, or a
numeric vector with the species indices, or a logical vector indicating for
each species whether it is to be selected (TRUE) or not.}

\item{start_time}{The first time to be plotted. Default is the beginning
of the time series.}

\item{end_time}{The last time to be plotted. Default is the end of the
time series.}

\item{y_ticks}{The approximate number of ticks desired on the y axis}

\item{ylim}{A numeric vector of length two providing lower and upper limits
for the y axis. Use NA to refer to the existing minimum or maximum. Any
values below 1e-20 are always cut off.}

\item{total}{A boolean value that determines whether the total biomass from
all species is plotted as well. Default is FALSE.}

\item{background}{A boolean value that determines whether background species
are included. Ignored if the model does not contain background species.
Default is TRUE.}

\item{highlight}{Name or vector of names of the species to be highlighted.}

\item{return_data}{A boolean value that determines whether the formatted data
used for the plot is returned instead of the plot itself. Default value is FALSE}

\item{...}{
  Arguments passed on to \code{\link[=get_size_range_array]{get_size_range_array}}
  \describe{
    \item{\code{min_w}}{Smallest weight in size range. Defaults to smallest weight in
the model.}
    \item{\code{max_w}}{Largest weight in size range. Defaults to largest weight in the
model.}
    \item{\code{min_l}}{Smallest length in size range. If supplied, this takes
precedence over \code{min_w}.}
    \item{\code{max_l}}{Largest length in size range. If supplied, this takes precedence
over \code{max_w}.}
  }}
}
\value{
A ggplot2 object
}
\description{
After running a projection, the biomass of each species can be plotted
against time. The biomass is calculated within user defined size limits
(min_w, max_w, min_l, max_l, see \code{\link[=getBiomass]{getBiomass()}}).
}
\examples{
\donttest{
# Set up example MizerParams and MizerSim objects
sim <- project(NS_params, effort = 1, t_max = 20, t_save = 0.2, progress_bar = FALSE)

plotBiomass(sim)
plotBiomass(sim, species = c("Cod", "Haddock"), total = TRUE)
plotBiomass(sim, min_w = 10, max_w = 1000)
plotBiomass(sim, start_time = 10, end_time = 15)
plotBiomass(sim, y_ticks = 3)

# Returning the data frame
fr <- plotBiomass(sim, return_data = TRUE)
str(fr)
}
}
\seealso{
\link{plotting_functions}, \code{\link[=getBiomass]{getBiomass()}}

Other plotting functions: 
\code{\link{animateSpectra}()},
\code{\link{plot,MizerSim,missing-method}},
\code{\link{plotDiet}()},
\code{\link{plotFMort}()},
\code{\link{plotFeedingLevel}()},
\code{\link{plotGrowthCurves}()},
\code{\link{plotPredMort}()},
\code{\link{plotSpectra}()},
\code{\link{plotYieldGear}()},
\code{\link{plotYield}()},
\code{\link{plotting_functions}}
}
\concept{plotting functions}
