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
get_steady_state_n <- function(growth, mort, dw, diffusion = rep(0, length(dw)),
                               idx, N0) {
    if (any(diffusion[idx] > 0)) {
        # Steady state solution of the upwind-difference scheme used in project
        # We solve the system A_j * N_{j-1} + B_j * N_j + C_j * N_{j+1} = 0
        # The N_j are the densities at the size classes j in min_idx:max_idx
        min_idx <- min(idx)
        max_idx <- max(idx) + 1
        n <- max_idx - min_idx + 1
        
        # We need to subset the rate arrays to the range min_idx:max_idx
        # However, the coefficients at j depend on j-1, j, j+1.
        # So we need to be careful with indexing.
        # We will construct a, b, c vectors of length n.
        # The j-th element of these vectors corresponds to the size class min_idx + j - 1.
        
        # Ranges for the relevant size classes
        j_range <- min_idx:max_idx
        
        # Diffusion coefficient D_i(w)
        d <- diffusion
        # Growth rate g_i(w)
        g <- growth
        # Mortality rate mu_i(w)
        mu <- mort
        
        # Initialize vectors
        a <- numeric(n)
        b <- numeric(n)
        c <- numeric(n)
        rhs <- numeric(n)
        
        # Calculate coefficients for the inner points
        # The indices into the full arrays are j
        # The indices into the small arrays are k = j - min_idx + 1
        
        # We only need to form the equations for j from min_idx to max_idx.
        # But for j = min_idx we have the boundary condition N = N0.
        
        # Boundary condition at min_idx
        # N_{min_idx} = N0
        b[1] <- 1
        rhs[1] <- N0
        # a[1] and c[1] are 0
        
        # Now loop or vectorize for the rest
        # We iterate k from 2 to n.
        # This corresponds to j from min_idx + 1 to max_idx.
        if (n > 1) {
            k <- 2:n
            j <- j_range[k]
            
            # Using formulas from transport.R, divided by dt
            # a_j = - 1/dw_j * (g_{j-1} + D_{j-1} / (2 * dw_{j-1}))
            term_diff_minus_1 <- 0.5 * d[j - 1] / dw[j - 1]
            a[k] <- - (g[j - 1] + term_diff_minus_1) / dw[j]
            
            # c_j = - 1/dw_j * (D_{j+1} / (2 * dw_j))
            # Note: At the last bin j=max_idx, we assume N_{j+1} = 0 (or flux is handled)
            # If j < length(dw), we compute c normally.
            # If j == length(dw), term_diff_plus_1 involves d[length+1]??
            # In project_n/transport, c[no_w] is set to 0.
            # Here max_idx could be the last bin.
            c[k] <- 0 # Default to 0
            
            # Create a mask for valid j+1
            valid_c <- j < length(dw)
            if (any(valid_c)) {
                # Only compute for valid j
                # Indices in k that are valid
                k_valid <- k[valid_c]
                j_valid <- j[valid_c]
                term_diff_plus_1 <- 0.5 * d[j_valid + 1] / dw[j_valid]
                c[k_valid] <- - term_diff_plus_1 / dw[j_valid]
            }
            
            # b_j = mu_j + 1/dw_j * (g_j + D_j / (2 * dw_j) + D_j / (2 * dw_{j-1}))
            term_diff_j <- 0.5 * d[j] / dw[j]
            term_diff_j_minus_1 <- 0.5 * d[j] / dw[j - 1]
            b[k] <- mu[j] + (g[j] + term_diff_j + term_diff_j_minus_1) / dw[j]
        }
        
        # Solve
        n_exact <- thomas_solve(a, b, c, rhs)
        return(n_exact)
    }
    
    # Steady state solution of the upwind-difference scheme used in project
    n_exact <- c(1, cumprod(growth[idx] / ((growth + mort * dw)[idx + 1])))
    if (!missing(N0)) {
        n_exact <- N0 * n_exact
    }
    return(n_exact)
}
