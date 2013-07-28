# Selectivity functions for size based model
# First argument to the function has to be w

# Copyright 2012 Finlay Scott, Julia Blanchard and Ken Andersen.
# Distributed under the GPL 2 or later
# Maintainer: Finlay Scott, CEFAS

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

#' Size based knife-edge selectivity function
#'
#' A knife-edge selectivity function. The direction of the knife-edge is determined by the \code{knife_is_min} argument. If \code{knife_is_min} is TRUE, then all sizes equal to or greater than \code{knife_edge_size} are selected. If \code{knife_is_min} is FALSE, then all sizes equal to or less than \code{knife_edge_size} are selected.
#'
#' @param w The size of the individual.
#' @param knife_edge_size The size at which the knife-edge operates.
#' @param knife_is_min TRUE or FALSE. Sizes equal to or greater than (TRUE) or less than (FALSE) are selected.
#' @export
knife_edge <- function(w, knife_edge_size, knife_is_min = TRUE){
    sel <- rep(0, length(w))
    if (knife_is_min == TRUE){
        sel[w >= knife_edge_size] <- 1
    }
    if (knife_is_min == FALSE){
        sel[w <= knife_edge_size] <- 1
    }
    return(sel)
} 
