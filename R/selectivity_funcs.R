# Selectivity functions for size based model
# First argument to the function has to be w

# Copyright 2012 Finlay Scott, Julia Blanchard and Ken Andersen.
# Distributed under the GPL 3 or later 
# Maintainer: Gustav Delius, University of York, <gustav.delius@york.ac.uk>

#' Length based sigmoid selectivity function
#'
#' A sigmoid shaped selectivity function. Based on two parameters \code{l25} and
#' \code{l50} which determine the length at which 25% and 50% of the stock is
#' selected respectively. As the size-based model is weight based, and this
#' selectivity function is length based, it uses the
#' length-weight parameters `a` and `b` to convert between length and weight.
#'
#' @param w the size of the individual.
#' @param l25 the length which gives a selectivity of 25%.
#' @param l50 the length which gives a selectivity of 50%.
#' @param species_params A list with the species params for the current species.
#'   Used to get at the length-weight parameters `a` and `b`
#' @param ... Unused
#' @export
sigmoid_length <- function(w, l25, l50, species_params, ...) {
    assert_that(is.numeric(w) && is.numeric(l25) && is.numeric(l50))
    assert_that(l25 > 0 && l25 < l50)
    a <- species_params[["a"]]
    b <- species_params[["b"]]
    if (is.null(a) || is.null(b)) {
        stop("The selectivity function needs the weight-length parameters ",
             "`a` and `b` to be provided in the species_params data frame.")
    }
    l <- (w / a)^(1 / b)
    sr <- l50 - l25
    s1 <- l50 * log(3) / sr
    s2 <- s1 / l50
    return(1 / (1 + exp(s1 - s2 * l)))
}

#' Length based double-sigmoid selectivity function
#'
#' A hump-shaped selectivity function with a sigmoidal rise and an independent
#' sigmoidal drop-off. This drop-off is what distinguishes this from the
#' function [sigmoid_length()] and it is intended to model the escape
#' of large individuals from the fishing gear. 
#' 
#' The selectivity is obtained as the product of two sigmoidal curves, one
#' rising and one dropping. The sigmoidal rise is based on the two parameters
#' \code{l25} and \code{l50} which determine the length at which 25% and 50%
#' of the stock is selected respectively. The sigmoidal drop-off is based on the
#' two parameters \code{l50_right} and \code{l25_right} which determine the
#' length at which the selectivity curve has dropped back to 50% and 25%
#' respectively.
#' 
#' As the size-based model is weight based, and this
#' selectivity function is length based, it uses the
#' length-weight parameters `a` and `b` to convert between length and weight.
#'
#' @param w the size of the individual.
#' @param l25 the length which gives a selectivity of 25%.
#' @param l50 the length which gives a selectivity of 50%.
#' @param l50_right the length which gives a selectivity of 50%.
#' @param l25_right the length which gives a selectivity of 25%.
#' @param species_params A list with the species params for the current species.
#'   Used to get at the length-weight parameters `a` and `b`
#' @param ... Unused
#' @export
double_sigmoid_length <- function(w, l25, l50, l50_right, l25_right,
                                  species_params, ...) {
    assert_that(is.numeric(w) && 
                    is.numeric(l25) && is.numeric(l50) &&
                    is.numeric(l50_right) && is.numeric(l25_right))
    assert_that(l25 > 0, l25 < l50, l50_right < l25_right)
    
    a <- species_params[["a"]]
    b <- species_params[["b"]]
    if (is.null(a) || is.null(b)) {
      stop("The selectivity function needs the weight-length parameters ",
           "`a` and `b` to be provided in the species_params data frame.")
    }
    
    l <- (w / a)^(1 / b)
    
    sr <- l50 - l25
    s1 <- l50 * log(3) / sr
    s2 <- s1 / l50
    
    sr_right <- l50_right - l25_right
    s1_right <- l50_right * log(3) / sr_right
    s2_right <- s1_right / l50_right
    
    return(1 / (1 + exp(s1 - s2 * l)) / (1 + exp(s1_right - s2_right * l)))
}

#' Weight based knife-edge selectivity function
#'
#' A knife-edge selectivity function where weights greater or equal to
#' `knife_edge_size` are selected.
#'
#' @param w The size of the individual.
#' @param knife_edge_size The weight at which the knife-edge operates.
#' @param ... Unused
#' @export
knife_edge <- function(w, knife_edge_size, ...) {
    sel <- rep(0, length(w))
    sel[w >= knife_edge_size] <- 1
    return(sel)
} 

#' Weight based sigmoidal selectivity function
#'
#' A sigmoidal selectivity function with 50% selectivity at
#' weight `sigmoidal_weight` and width `sigmoidal_sigma`.
#'
#' @param w The size of the individual.
#' @param sigmoidal_weight The weight at which the knife-edge operates.
#' @param sigmoidal_sigma The width of the selection function
#' @param ... Unused
#' @export
sigmoid_weight <- function(w, sigmoidal_weight, sigmoidal_sigma, ...) {
  return((1 + (w / sigmoidal_weight) ^ (-sigmoidal_sigma)) ^ (-1) )
} 


