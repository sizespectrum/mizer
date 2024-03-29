% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setFishing.R
\name{validGearParams}
\alias{validGearParams}
\title{Check validity of gear parameters and set defaults}
\usage{
validGearParams(gear_params, species_params)
}
\arguments{
\item{gear_params}{Gear parameter data frame}

\item{species_params}{Species parameter data frame}
}
\value{
A valid gear parameter data frame
}
\description{
The function returns a valid gear parameter data frame that can be used
by \code{setFishing()} or it gives an error message.
}
\details{
The gear_params data frame is allowed to have zero rows, but if it has
rows, then the following requirements apply:
\itemize{
\item There must be columns \code{species} and \code{gear} and any species - gear pair is
allowed to appear at most once. Any species that appears must also appear
in the \code{species_params} data frame.
\item There must be a \code{sel_func} column. If a selectivity function is not
supplied, it will be set to "knife_edge".
\item There must be a \code{catchability} column. If a catchability is not supplied,
it will be set to 1.
\item All the parameters required by the selectivity functions must be provided.
}

If gear_params is empty, then this function tries to find the necessary
information in the species_params data frame. This restricts each species
to be fished by only one gear. Defaults are used for information that can
not be found in the species_params dataframe, as follows:
\itemize{
\item If there is no \code{gear} column or it is NA then a new gear named after the
species is introduced.
\item If there is no \code{sel_func} column or it is NA then \code{knife_edge} is used.
\item If there is no \code{catchability} column or it is NA then this is set to 1.
\item If the selectivity function is \code{knife_edge} and no \code{knife_edge_size} is
provided, it is set to \code{w_mat}.
}

The row names of the returned data frame are of the form
"species, gear".

When \code{gear_params} is \code{NULL} and there is no gear information in
\code{species_params}, then a gear called \code{knife_edge_gear} is set up with a
\code{knife_edge} selectivity for each species and a \code{knive_edge_size} equal to
\code{w_mat}. Catchability is set to 0.3 for all species.
}
\seealso{
\code{\link[=gear_params]{gear_params()}}
}
\concept{helper}
