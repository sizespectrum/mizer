#' Check whether two objects are different
#' 
#' Check whether two objects are numerically different, ignoring all attributes
#' 
#' @param a First object
#' @param b Second object
#' 
#' @return TRUE or FALSE
different <- function(a, b) {
    !isTRUE(all.equal(a, b, check.attributes = FALSE, scale = 1, tolerance = 0))
}