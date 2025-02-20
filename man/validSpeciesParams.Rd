% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validSpeciesParams.R
\name{validSpeciesParams}
\alias{validSpeciesParams}
\alias{validGivenSpeciesParams}
\title{Validate species parameter data frame}
\usage{
validSpeciesParams(species_params)

validGivenSpeciesParams(species_params)
}
\arguments{
\item{species_params}{The user-supplied species parameter data frame}
}
\value{
For \code{validSpeciesParams()}: A valid species parameter data frame with
additional parameters with default values.

For \code{validGivenSpeciesParams()}: A valid species parameter data frame
without additional parameters.
}
\description{
These functions check the validity of a species parameter frame and, where
necessary, make corrections. \code{validGivenSpeciesParams()} only checks and
corrects the given species parameters but does not add default values for
species parameters that were not provided. \code{validSpeciesParams()} first calls
\code{validGivenSpeciesParams()} but then goes further by adding default values
for species parameters that were not provided.
}
\details{
\code{validGivenSpeciesParams()} checks the validity of the given species
parameter It throws an error if
\itemize{
\item the \code{species} column does not exist or contains duplicates
\item the maximum size is not specified for all species
}

If a weight-based parameter is missing but the corresponding length-based
parameter is given, as well as the \code{a} and \code{b} parameters for length-weight
conversion, then the weight-based parameters are added. If both length and
weight are given, then weight is used and a warning is issued if the two are
inconsistent.

If a \code{w_inf} column is given but no \code{w_max} then the value from \code{w_inf} is
used. This is for backwards compatibility. But note that the von Bertalanffy
parameter \code{w_inf} is not the maximum size of the largest individual, but the
asymptotic size of an average individual.

Some inconsistencies in the size parameters are resolved as follows:
\itemize{
\item Any \code{w_mat} that is not smaller than \code{w_max} is set to \code{w_max / 4}.
\item Any \code{w_mat25} that is not smaller than \code{w_mat} is set to NA.
\item Any \code{w_min} that is not smaller than \code{w_mat} is set to \code{0.001} or
\code{w_mat /10}, whichever is smaller.
\item Any \code{w_repro_max} that is not larger than \code{w_mat} is set to \code{4 * w_mat}.
}

The row names of the returned data frame will be the species names.
If \code{species_params} was provided as a tibble it is converted back to an
ordinary data frame.

The function tests for some typical misspellings of parameter names, like
wrong capitalisation or missing underscores and issues a warning if it
detects such a name.

\code{validSpeciesParams()} first calls \code{validGivenSpeciesParams()} but then
goes further by adding default values for species parameters that were not
provided. The function sets default values if any of the following species
parameters are missing or NA:
\itemize{
\item \code{w_repro_max} is set to \code{w_max}
\item \code{w_mat} is set to \code{w_max/4}
\item \code{w_min} is set to \code{0.001}
\item \code{alpha} is set to \code{0.6}
\item \code{interaction_resource} is set to \code{1}
\item \code{n} is set to \code{3/4}
\item \code{p} is set to \code{n}
}

Note that the species parameters returned by these functions are not
guaranteed to produce a viable model. More checks of the parameters are
performed by the individual rate-setting functions (see \code{\link[=setParams]{setParams()}} for the
list of these functions).
}
\seealso{
\code{\link[=species_params]{species_params()}}, \code{\link[=validGearParams]{validGearParams()}}, \code{\link[=validParams]{validParams()}}, \code{\link[=validSim]{validSim()}}
}
\concept{helper}
