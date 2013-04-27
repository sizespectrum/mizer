# Selectivity functions for size based model
# First argument has to be w

# Copyright 2012 Finlay Scott and Julia Blanchard. Distributed under the GPL 2 or later
# Maintainer: Finlay Scott, CEFAS

# Based on length
#sigmoidLength <- function(L,L25,L50)
#{
#  SR <- L50 - L25
#  S1 <- L50*log(3)/SR
#  S2 <- S1 / L50
#  return(1 / (1 + exp(S1 - S2*L)))
#}

#' Length based sigmoid selectivity function
#'
#' A sigmoid shaped selectivity function. Based on two parameters \code{l25} and \code{l50} which determine the length at which 25\% and 50\% of the stock is selected respectively. As the size-based model is weight based, and this selectivity function is length based, it is also necessary to supply the length-weight parameters \code{a} and \code{b}.
#'
#' @param w the size of the individual.
#' @param l25 the length which gives a selectivity of 25\%.
#' @param l50 the length which gives a selectivity of 50\%.
#' @param a the multiplier of the length-weight function.
#' @param b the exponent of the length-weight function.
#' @export
sigmoid_length <- function(w,l25,l50,a,b)
{
    l <- (w/a)^(1/b)
    sr <- l50 - l25
    s1 <- l50*log(3)/sr
    s2 <- s1 / l50
    return(1 / (1 + exp(s1 - s2*l)))
}

