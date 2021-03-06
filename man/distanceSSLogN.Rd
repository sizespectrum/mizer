% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/steady.R
\name{distanceSSLogN}
\alias{distanceSSLogN}
\title{Measure distance between current and previous state in terms of fish abundances}
\usage{
distanceSSLogN(params, current, previous)
}
\arguments{
\item{params}{MizerParams}

\item{current}{A named list with entries \code{n}, \code{n_pp} and \code{n_other}
describing the current state}

\item{previous}{A named list with entries \code{n}, \code{n_pp} and \code{n_other}
describing the previous state}
}
\value{
The sum of squares of the difference in the logs of the (nonzero)
fish abundances n:
\code{sum((log(current$n) - log(previous$n))^2)}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#experimental}{\figure{lifecycle-experimental.svg}{options: alt='[Experimental]'}}}{\strong{[Experimental]}}

Calculates the sum squared difference between log(N) in current and previous
state. This function can be used in \code{\link[=projectToSteady]{projectToSteady()}} to decide when
sufficient convergence to steady state has been achieved.
}
\seealso{
Other distance functions: 
\code{\link{distanceMaxRelRDI}()}
}
\concept{distance functions}
\concept{helper}
