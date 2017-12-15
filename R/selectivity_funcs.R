# Selectivity functions for size based model
# First argument to the function has to be w

# Copyright 2012 Finlay Scott, Julia Blanchard and Ken Andersen.
# Distributed under the GPL 3 or later 
# Maintainer: Gustav Delius, University of York, <gustav.delius@york.ac.uk>

#' Length based sigmoid selectivity function
#'
#' A sigmoid shaped selectivity function. Based on two parameters \code{l25} and
#' \code{l50} which determine the length at which 25\% and 50\% of the stock is
#' selected respectively. As the size-based model is weight based, and this
#' selectivity function is length based, it is also necessary to supply the
#' length-weight parameters \code{a} and \code{b}.
#'
#' @param w the size of the individual.
#' @param l25 the length which gives a selectivity of 25\%.
#' @param l50 the length which gives a selectivity of 50\%.
#' @param a the multiplier of the length-weight function.
#' @param b the exponent of the length-weight function.
#' @export
sigmoid_length <- function(w,l25,l50,a,b) {
    l <- (w/a)^(1/b)
    sr <- l50 - l25
    s1 <- l50*log(3)/sr
    s2 <- s1 / l50
    return(1 / (1 + exp(s1 - s2*l)))
}

#' Size based knife-edge selectivity function
#'
#' A knife-edge selectivity function where only sizes greater or equal to
#' \code{knife_edge_size} are selected.
#'
#' @param w The size of the individual.
#' @param knife_edge_size The size at which the knife-edge operates.
#' @export
knife_edge <- function(w, knife_edge_size) {
    sel <- rep(0, length(w))
    sel[w >= knife_edge_size] <- 1
    return(sel)
} 
