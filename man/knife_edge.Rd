\name{knife_edge}
\alias{knife_edge}
\title{Size based knife-edge selectivity function}
\usage{
  knife_edge(w, w_limit, is_min = TRUE)
}
\arguments{
  \item{w}{The size of the individual.}

  \item{w_limit}{The size at which the knife-edge
  operates.}

  \item{is_min}{TRUE or FALSE. Sizes equal to or greater
  than (TRUE) or less than (FALSE) are selected.}
}
\description{
  A knife-edge selectivity function. The direction of the
  knife-edge is determined by the \code{is_min} argument.
  If \code{is_min} is TRUE, then all sizes equal to or
  greater than \code{w_limit} are selected. If
  \code{is_min} is FALSE, then all sizes equal to or less
  than \code{w_limit} are selected.
}

