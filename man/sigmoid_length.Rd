% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selectivity_funcs.R
\name{sigmoid_length}
\alias{sigmoid_length}
\title{Length based sigmoid selectivity function}
\usage{
sigmoid_length(w, l25, l50, species_params, ...)
}
\arguments{
\item{w}{Vector of sizes.}

\item{l25}{the length which gives a selectivity of 25\%.}

\item{l50}{the length which gives a selectivity of 50\%.}

\item{species_params}{A list with the species params for the current species.
Used to get at the length-weight parameters \code{a} and \code{b}.}

\item{...}{Unused}
}
\value{
Vector of selectivities at the given sizes.
}
\description{
A sigmoid shaped selectivity function. Based on two parameters \code{l25} and
\code{l50} which determine the length at which 25\% and 50\% of the stock is
selected respectively. As the size-based model is weight based, and this
selectivity function is length based, it uses the
length-weight parameters \code{a} and \code{b} to convert between length and weight.
}
\seealso{
Other selectivity functions: 
\code{\link{double_sigmoid_length}()},
\code{\link{knife_edge}()},
\code{\link{sigmoid_weight}()}
}
\concept{selectivity functions}
