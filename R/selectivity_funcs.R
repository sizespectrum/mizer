# Selectivity functions for size based model
# First argument to the function has to be w

# Copyright 2012 Finlay Scott, Julia Blanchard and Ken Andersen.
# Distributed under the GPL 3 or later 
# Maintainer: Gustav Delius, University of York, <gustav.delius@york.ac.uk>

#' Length based sigmoid selectivity function
#'
#' A sigmoid shaped selectivity function. Based on two parameters \code{l25} and
#' \code{l50} which determine the length at which 25% and 50% of the stock is
#' selected respectively.
#'
#' You would not usually call this function directly. Instead, set the `sel_func`
#' column in [gear_params()] to `"sigmoid_length"` and provide the `l25` and
#' `l50` values as additional columns. [setFishing()] will then call this
#' function automatically when calculating the selectivity array.
#'
#' The selectivity is given by the logistic function
#' \deqn{S(l) = \frac{1}{1 + \exp\left(\log(3)\frac{l50 -l}{l50 - l25}\right)}}{S(l) = 1/(1 + exp(log(3)*(l50 -l) / (l50 - l25)))}
#' As the mizer model is weight based, and this
#' selectivity function is length based, it uses the
#' length-weight parameters `a` and `b` to convert between length and weight
#' \deqn{l = \left(\frac{w}{a}\right)^{1/b}}{l = (w/a)^(1/b)}
#'
#' @param w Vector of sizes.
#' @param l25 the length which gives a selectivity of 25%.
#' @param l50 the length which gives a selectivity of 50%.
#' @param species_params A list with the species params for the current species.
#'   Used to get at the length-weight parameters `a` and `b`.
#' @param ... Unused
#' @return Vector of selectivities at the given sizes.
#' @export
#' @seealso [gear_params()] for setting the selectivity parameters.
#' @family selectivity functions
#' @examples
#' # Selectivity at weight given l25 = 10 cm, l50 = 15 cm
#' # using length-weight parameters a = 0.01, b = 3
#' sp <- list(a = 0.01, b = 3)
#' w <- c(1, 10, 100, 500, 1000)
#' sigmoid_length(w, l25 = 10, l50 = 15, species_params = sp)
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
#' function [sigmoid_length()] and it is intended to model the escape of large
#' individuals from the fishing gear.
#'
#' You would not usually call this function directly. Instead, set the `sel_func`
#' column in [gear_params()] to `"double_sigmoid_length"` and provide the
#' `l25`, `l50`, `l50_right` and `l25_right` values as additional columns.
#' [setFishing()] will then call this function automatically when calculating
#' the selectivity array.
#'
#' The selectivity is obtained as the product of two sigmoidal curves, one
#' rising and one dropping. The sigmoidal rise is based on the two parameters
#' \code{l25} and \code{l50} which determine the length at which 25% and 50% of
#' the stock is selected respectively. The sigmoidal drop-off is based on the
#' two parameters \code{l50_right} and \code{l25_right} which determine the
#' length at which the selectivity curve has dropped back to 50% and 25%
#' respectively. The selectivity is given by the function \deqn{S(l) =
#' \frac{1}{1 + \exp\left(\log(3)\frac{l50 -l}{l50 - l25}\right)}\frac{1}{1 +
#' \exp\left(\log(3)\frac{l50_{right} -l}{l50_{right} -
#' l25_{right}}\right)}}{S(l) = 1/(1 + exp(log(3)*(l50 -l) / (l50 - l25)))/(1
#' + exp(log(3)*(l50 -l) / (l50 - l25)))}
#'
#' As the size-based model is weight based, and this selectivity function is
#' length based, it uses the length-weight parameters `a` and `b` to convert
#' between length and weight. \deqn{l = \left(\frac{w}{a}\right)^{1/b}}{l =
#' (w/a)^(1/b)}
#' 
#' @param w Vector of sizes.
#' @param l25 the length which gives a selectivity of 25%.
#' @param l50 the length which gives a selectivity of 50%.
#' @param l50_right the length which gives a selectivity of 50%.
#' @param l25_right the length which gives a selectivity of 25%.
#' @param species_params A list with the species params for the current species.
#'   Used to get at the length-weight parameters `a` and `b`
#' @param ... Unused
#' @return Vector of selectivities at the given sizes.
#'   Requires `l25 < l50 < l50_right < l25_right`.
#' @export
#' @seealso [gear_params()] for setting the selectivity parameters.
#' @family selectivity functions
#' @examples
#' # Hump-shaped selectivity: rises from l25=10 to l50=15,
#' # then drops back to 50% at l50_right=40 and 25% at l25_right=50
#' sp <- list(a = 0.01, b = 3)
#' w <- c(1, 10, 100, 500, 1000)
#' double_sigmoid_length(w, l25 = 10, l50 = 15,
#'                       l50_right = 40, l25_right = 50,
#'                       species_params = sp)
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
#' `knife_edge_size` are fully selected and no fish smaller than this size
#' are selected.
#'
#' You would not usually call this function directly. Instead, set the `sel_func`
#' column in [gear_params()] to `"knife_edge"` and provide `knife_edge_size` as
#' an additional column. [setFishing()] will then call this function
#' automatically when calculating the selectivity array.
#'
#' @param w Vector of sizes.
#' @param knife_edge_size The weight at which the knife-edge operates.
#' @param ... Unused
#' @return Vector of selectivities at the given sizes.
#' @export
#' @seealso [gear_params()] for setting the `knife_edge_size` parameter.
#' @family selectivity functions
#' @examples
#' knife_edge(w = c(1, 10, 100, 1000), knife_edge_size = 100)
knife_edge <- function(w, knife_edge_size, ...) {
    sel <- rep(0, length(w))
    sel[w >= knife_edge_size] <- 1
    return(sel)
} 

#' Weight based sigmoidal selectivity function
#'
#' A sigmoidal selectivity function with 50% selectivity at
#' weight `sigmoidal_weight` \eqn{=w_{\text{sigmoid}}} and width `sigmoidal_sigma` \eqn{=\sigma}.
#' \deqn{S(w) = \left(1 + \left(\frac{w}{w_{\text{sigmoid}}}\right)^{-\sigma}\right)^{-1}}{S(w) = (1 + (w/sigmoidal_weight)^{-sigmoidal_sigma})^{-1}}
#'
#' You would not usually call this function directly. Instead, set the `sel_func`
#' column in [gear_params()] to `"sigmoid_weight"` and provide `sigmoidal_weight`
#' and `sigmoidal_sigma` as additional columns. [setFishing()] will then call
#' this function automatically when calculating the selectivity array.
#'
#' @param w Vector of sizes.
#' @param sigmoidal_weight The weight at which selectivity is 50%.
#' @param sigmoidal_sigma The width of the selection function.
#' @param ... Unused
#' @return Vector of selectivities at the given sizes.
#' @export
#' @seealso [gear_params()] for setting the selectivity parameters.
#' @family selectivity functions
#' @examples
#' sigmoid_weight(w = c(1, 10, 100, 1000),
#'                sigmoidal_weight = 100, sigmoidal_sigma = 3)
sigmoid_weight <- function(w, sigmoidal_weight, sigmoidal_sigma, ...) {
  return((1 + (w / sigmoidal_weight) ^ (-sigmoidal_sigma)) ^ (-1))
} 
