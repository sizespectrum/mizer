\name{knife_edge}
\alias{knife_edge}
\title{Size based knife-edge selectivity function}
\usage{
  knife_edge(w, knife_edge_size, knife_is_min = TRUE)
}
\arguments{
  \item{w}{The size of the individual.}

  \item{knife_edge_size}{The size at which the knife-edge
  operates.}

  \item{knife_is_min}{TRUE or FALSE. Sizes equal to or
  greater than (TRUE) or less than (FALSE) are selected.}
}
\description{
  A knife-edge selectivity function. The direction of the
  knife-edge is determined by the \code{knife_is_min}
  argument. If \code{knife_is_min} is TRUE, then all sizes
  equal to or greater than \code{knife_edge_size} are
  selected. If \code{knife_is_min} is FALSE, then all sizes
  equal to or less than \code{knife_edge_size} are
  selected.
}

