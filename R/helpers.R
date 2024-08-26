#' Check whether two objects are different
#' 
#' Check whether two objects are numerically different, ignoring all attributes.
#' 
#' We use this helper function in particular to see if a new value for a slot
#' in MizerParams is different from the existing value in order to give the
#' appropriate messages.
#' 
#' @param a First object
#' @param b Second object
#' 
#' @return TRUE or FALSE
#' @concept helper
different <- function(a, b) {
    !isTRUE(all.equal(a, b, check.attributes = FALSE, scale = 1, 
                      tolerance = 10 * .Machine$double.eps))
}

#' Length-weight conversion
#' 
#' For each species, convert between length and weight using the relationship
#' \deqn{w_i = a_i l_i^{b_i}}{w_i = a_i l_i^b_i} or 
#' \deqn{l_i = (w_i / a_i)^{1/b_i}}{l_i = (w_i / a_i)^{1/b_i}
#' where `a` and `b` are taken from the species parameter data frame and 
#' \eqn{i}{i} is the species index.
#' 
#' This is useful for converting a length-based species parameter to a 
#' weight-based species parameter.
#' 
#' If any `a` or `b` parameters are missing the default values `a = 0.01` and
#' `b = 3` are used for the missing values.
#' 
#' @param l Lengths in cm. Either a single number used for all species or a
#'   vector with one number for each species.
#' @param w Weights in grams. Either a single number used for all species or a
#'   vector with one number for each species.
#' @param params A species parameter data frame or a MizerParams object.
#' @return A vector with one entry for each species. `l2w()` returns a vector
#' of weights in grams and `w2l()` returns a vector of lengths in cm.
#' @export
#' @concept helper
l2w <- function(l, params) {
    assert_that(is.numeric(l))
    sp <- params
    if (is(params, "MizerParams")) {
        sp <- params@species_params
    }
    if (!is.data.frame(sp)) {
        stop("The second argument must be either a MizerParams object or a
             species paramter data frame.")
    }
    if (length(l) != 1 && length(l) != nrow(sp)) {
        stop("The length of 'l' should be one or equal to the number of species.")
    }
    sp <- sp %>%
        set_species_param_default("a", 0.01,
                                  "Using default values for 'a' parameter.") %>%
        set_species_param_default("b", 3,
                                  "Using default values for 'a' parameter.")
    
    sp[["a"]] * l^sp[["b"]]
}

#' @rdname l2w
#' @export
w2l <- function(w, params) {
    assert_that(is.numeric(w))
    sp <- params
    if (is(params, "MizerParams")) {
        sp <- params@species_params
    }
    if (!is.data.frame(sp)) {
        stop("The second argument must be either a MizerParams object or a
             species paramter data frame.")
    }
    if (length(w) != 1 && length(w) != nrow(sp)) {
        stop("The length of 'w' should be one or equal to the number of species.")
    }
    sp <- sp %>%
        set_species_param_default("a", 0.01,
                                  "Using default values for 'a' parameter.") %>%
        set_species_param_default("b", 3,
                                  "Using default values for 'a' parameter.")
    
    (w / sp[["a"]])^(1 / sp[["b"]])
}
