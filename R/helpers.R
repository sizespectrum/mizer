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

#' Bin average of a power law over geometric bins
#'
#' Computes the exact average of the power law \eqn{w^d} over each bin
#' \eqn{[w_j, w_{j+1}]}, i.e.
#' \deqn{\overline{w^d}_j = \frac{1}{\Delta w_j}\int_{w_j}^{w_{j+1}} w^d\, dw.}
#'
#' The integral has a closed form, so the result is exact (not merely second
#' order):
#' \deqn{\overline{w^d}_j = \frac{w_{j+1}^{d+1} - w_j^{d+1}}{(d+1)\,\Delta w_j},
#'   \quad d \neq -1,}
#' \deqn{\overline{w^d}_j = \frac{\ln(w_{j+1}/w_j)}{\Delta w_j},
#'   \quad d = -1.}
#'
#' This is used by the bin-averaged (second-order) code paths that need the
#' average of a power-law rate over each bin, for example [setExtMort()] and
#' the resource capacity and rate in [setResource()]. The grid does not need to
#' be geometric; only the left bin edges `w` and the bin widths `dw` are used,
#' with \eqn{w_{j+1} = w_j + \Delta w_j}.
#'
#' An optional upper cutoff `w_max` handles a knife-edge truncation of the power
#' law (for example the resource carrying capacity, which is \eqn{\kappa
#' w^{-\lambda}} below `w_pp_cutoff` and zero above it). The bin straddling the
#' cutoff then receives the *partial* bin-average — the power-law average over
#' the part of the bin below `w_max`, divided by the full bin width — and bins
#' entirely above the cutoff get zero.
#'
#' @param w Numeric vector of left bin edges \eqn{w_j}.
#' @param dw Numeric vector of bin widths \eqn{\Delta w_j} (same length as `w`).
#' @param d Single numeric exponent of the power law.
#' @param w_max Optional upper cutoff. The power law is taken to be zero above
#'   `w_max`, with the straddling bin getting the partial average. Defaults to
#'   `Inf` (no cutoff), in which case the result is identical to the uncut
#'   formula.
#'
#' @return A numeric vector (same length as `w`) of bin averages of
#'   \eqn{w^d} (truncated at `w_max` when supplied).
#' @concept helper
#' @keywords internal
power_law_bin_average <- function(w, dw, d, w_max = Inf) {
    # Upper integration edge of each bin, clipped at the cutoff. The pmax keeps
    # bins lying entirely above the cutoff from contributing a negative (they
    # collapse to a zero-width interval and so average to zero).
    w_next <- pmax(pmin(w + dw, w_max), w)
    if (d == -1) {
        log(w_next / w) / dw
    } else {
        (w_next^(d + 1) - w^(d + 1)) / ((d + 1) * dw)
    }
}

#' Geometric bin centres of the size grid
#'
#' Internal helper for the second-order plotting code. A finite-volume cell
#' average \eqn{N_j = (1/\Delta w_j)\int_{w_j}^{w_{j+1}} N\,dw} does not live at
#' the left bin edge \eqn{w_j} but at the geometric bin centre
#' \deqn{w^*_j = \sqrt{w_j\,w_{j+1}} = w_j\sqrt{\beta},}
#' where \eqn{\beta = w_{j+1}/w_j} is the (constant) bin ratio of the
#' logarithmic grid. This is the log-symmetric, second-order-correct location at
#' which to plot a bin-averaged quantity (it is exact for the community spectrum
#' \eqn{N\propto w^{-2}}). It is a uniform half-bin shift to the right on the log
#' axis, the same for the consumer grid `w` and the full prey/resource grid
#' `w_full`.
#'
#' @param params A MizerParams object.
#' @param w_full If `TRUE`, return the centres of the full (resource) grid
#'   `params@w_full`; otherwise the consumer grid `params@w`.
#' @return A numeric vector of geometric bin centres, one per grid node.
#' @concept helper
#' @keywords internal
bin_midpoints <- function(params, w_full = FALSE) {
    # The grid is geometric, so the bin ratio beta is constant across the whole
    # (consumer and resource) grid.
    beta <- params@w_full[2] / params@w_full[1]
    w <- if (w_full) params@w_full else params@w
    w * sqrt(beta)
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
                                  "Using default values for 'b' parameter.")

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
                                  "Using default values for 'b' parameter.")

    (w / sp[["a"]])^(1 / sp[["b"]])
}

#' Calculate steady state abundance
#'
#' This function calculates the steady state abundance by solving the
#' transport equation with given growth and mortality rates. It sets up a
#' tri-diagonal system and solves it.
#'
#' @param params A MizerParams object
#' @param g A matrix of growth rates (species x size)
#' @param mu A matrix of mortality rates (species x size)
#' @param D A matrix of diffusion rates (species x size)
#' @param N0 A vector with the abundance at the smallest size for each species
#' @param max_iterations Maximum number of Picard iterations used when a flux
#'   limiter is active.
#' @param tol Relative convergence tolerance for the Picard iteration.
#' @param relax Under-relaxation factor in (0, 1] for the Picard iteration when a
#'   flux limiter is active.
#' @return A matrix with the steady state abundance
#'
#' @details
#' The spatial discretisation of the advective flux is read from the
#' `flux` entry of the `second_order_w` slot of `params`. With a second-order
#' flux scheme active the steady state must match the one that [project()] converges
#' to. Because the limiter depends on the solution, the steady state is then
#' found by an under-relaxed Picard iteration: the limiter is frozen at the
#' current iterate, the resulting tridiagonal system is solved, and the iterate
#' is updated towards that solution, repeating until it converges. (At `dt = 1`
#' the limited operator is not diagonally dominant, so the plain fixed-point map
#' only stalls; under-relaxation makes it converge.)
#' @concept helper
get_steady_state_n <- function(params, g, mu, D, N0,
                               max_iterations = 500, tol = 1e-10,
                               relax = 0.3) {
    # Advective-flux scheme read from the model's second_order_w slot.
    flux_limiter <- flux_limiter_scheme(params)
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)
    n <- matrix(0, nrow = no_sp, ncol = no_w,
                dimnames = list(params@species_params$species, dimnames(params@initial_n)[[2]]))

    j_start <- params@w_min_idx
    idxs_start <- cbind(1:no_sp, j_start)

    for (iteration in seq_len(max_iterations)) {
        # Coefficients with dt = 1 and no recruitment flux (the boundary is
        # handled manually below). The flux limiter is frozen at the current
        # iterate `n`; with flux_limiter = "none" this is the upwind operator and
        # the loop exits after a single, exact solve.
        coefs <- get_transport_coefs(params, n, g, mu, dt = 1,
                                     recruitment_flux = rep(0, no_sp),
                                     d = D, flux_limiter = flux_limiter)

        a <- coefs$a
        # For steady state, the diagonal term \tilde{B} is B - 1
        b <- coefs$b - 1
        c <- coefs$c
        # Steady-state right-hand side: zero except the fixed egg density at the
        # lower boundary. (We build this directly rather than reusing coefs$S,
        # which equals the frozen iterate `n`.)
        S <- matrix(0, nrow = no_sp, ncol = no_w, dimnames = dimnames(n))

        # Boundary conditions at the start of the size spectrum:
        # A_j = 0, B_j = 1, C_j = 0, S_j = N0
        a[idxs_start] <- 0
        b[idxs_start] <- 1
        c[idxs_start] <- 0
        S[idxs_start] <- N0

        n_solve <- project_n_loop(matrix(0, nrow = no_sp, ncol = no_w,
                                         dimnames = dimnames(n)),
                                  a, b, c, S, j_start)

        # Without a limiter the system is linear, so one solve is exact.
        if (flux_limiter == "none") {
            return(n_solve)
        }

        # Under-relax the update and keep the iterate non-negative (the limited
        # operator is not an M-matrix, so a bare solve can dip slightly below 0).
        n_new <- (1 - relax) * n + relax * n_solve
        n_new[n_new < 0] <- 0
        change <- max(abs(n_new - n)) / max(abs(n_new))
        n <- n_new
        if (is.finite(change) && change < tol) {
            break
        }
    }

    return(n)
}
