% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selectivity_funcs.R
\name{sigmoid_weight}
\alias{sigmoid_weight}
\title{Weight based sigmoidal selectivity function}
\usage{
sigmoid_weight(w, sigmoidal_weight, sigmoidal_sigma, ...)
}
\arguments{
\item{w}{Vector of sizes.}

\item{sigmoidal_weight}{The weight at which the knife-edge operates.}

\item{sigmoidal_sigma}{The width of the selection function.}

\item{...}{Unused}
}
\value{
Vector of selectivities at the given sizes.
}
\description{
A sigmoidal selectivity function with 50\% selectivity at
weight \code{sigmoidal_weight} \eqn{=w_{\text{sigmoid}}} and width \code{sigmoidal_sigma} \eqn{=\sigma}.
\deqn{S(w) = \left(1 + \left(\frac{w}{w_{\text{sigmoid}}}\right)^{-\sigma}\right)^{-1}}{S(w) = (1 + (w/sigmoidal_weight)^{-sigmoidal_sigma})^{-1}}
}
\seealso{
\code{\link[=gear_params]{gear_params()}} for setting the selectivity parameters.

Other selectivity functions: 
\code{\link{double_sigmoid_length}()},
\code{\link{knife_edge}()},
\code{\link{sigmoid_length}()}
}
\concept{selectivity functions}
