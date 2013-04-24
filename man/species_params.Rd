\docType{data}
\name{species_params}
\alias{species_params}
\alias{species_params_gears}
\alias{species_params_nogears}
\title{Example parameters data sets}
\description{
  There are two example data sets,
  \code{species_params_gears} and
  \code{species_params_nogears}. They are both data frames
  and are identical except that one has an additional
  column specifying the fishing gear that operates on each
  species. The data frames contain all the necessary
  information to be used by the \code{\link{MizerParams}}
  constructor. The data set without a fishing gear column
  will set up a model in which each species is fished by a
  separate gear.
}
\references{
  The North Sea paper (Blanchard et al)
}
\seealso{
  \code{MizerParams}
}
\keyword{datasets}
\keyword{species_params_nogears}

