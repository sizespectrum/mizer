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
#' \deqn{l_i = (w_i / a_i)^{1/b_i}}{l_i = (w_i / a_i)^{1/b_i}}
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
#' @param species_params A species parameter data frame or a MizerParams object.
#' @return A vector with one entry for each species. `l2w()` returns a vector
#' of weights in grams and `w2l()` returns a vector of lengths in cm.
#' @export
#' @concept helper
l2w <- function(l, species_params) {
    assert_that(is.numeric(l))
    sp <- species_params
    if (is(species_params, "MizerParams")) {
        sp <- species_params@species_params
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
w2l <- function(w, species_params) {
    assert_that(is.numeric(w))
    sp <- species_params
    if (is(species_params, "MizerParams")) {
        sp <- species_params@species_params
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

#' Helper function to calculate the steady state abundance using the upwind-difference scheme
#'
#' @param growth A numeric vector of growth rates.
#' @param mort A numeric vector of mortality rates.
#' @param dw A numeric vector of the size step.
#' @param idx A numeric vector of indices.
#' @param N0 The initial egg density.
#' @return A numeric vector representing the steady state abundances.
#' @keywords internal
get_steady_state_n <- function(params, g, mu, N0) {
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)
    n <- matrix(0, nrow = no_sp, ncol = no_w,
                dimnames = list(params@species_params$species, dimnames(params@initial_n)[[2]]))

    # Use get_transport_coefs to compute the coefficients with dt = 1
    # and no recruitment flux (since we handle the boundary manually)
    coefs <- get_transport_coefs(params, n, g, mu, dt = 1,
                                 recruitment_flux = rep(0, no_sp))
    
    a <- coefs$a
    # For steady state, the diagonal term \tilde{B} is B - 1
    b <- coefs$b - 1
    c <- coefs$c
    S <- coefs$S
    
    # Boundary conditions at the start of the size spectrum:
    # A_j = 0, B_j = 1, C_j = 0, S_j = N0
    j_start <- params@w_min_idx
    idxs_start <- cbind(1:no_sp, j_start)
    a[idxs_start] <- 0
    b[idxs_start] <- 1
    c[idxs_start] <- 0
    S[idxs_start] <- N0
    
    # Boundary conditions for sizes w > w_max
    # We want population to be 0 for these sizes
    w_idx_mat <- matrix(1:no_w, nrow = no_sp, ncol = no_w, byrow = TRUE)
    w_max_idx <- params@species_params$w_max_idx
    if (is.null(w_max_idx)) {
        w_max_idx <- sapply(1:no_sp, function(i) {
            sum(params@w <= params@species_params$w_max[i])
        })
    }
    
    mask_above <- w_idx_mat > w_max_idx
    if (any(mask_above)) {
        a[mask_above] <- 0
        b[mask_above] <- 1
        c[mask_above] <- 0
        S[mask_above] <- 0
        
        # We also set c_j = 0 at the w_max_idx boundary to prevent flux
        # out of the modelled spectrum from affecting the steady state below
        idxs_max <- cbind(1:no_sp, w_max_idx)
        c[idxs_max] <- 0
    }
    
    # project_n_loop expects a,b,c,S with same dimensions
    # Provide j_start to the c++ loop
    n <- project_n_loop(n, a, b, c, S, j_start)
    
    return(n)
}
