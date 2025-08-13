#' Set initial abundances to single-species steady state abundances
#'
#' `r lifecycle::badge("experimental")`
#' This first calculates growth and death rates that arise from the current
#' initial abundances. Then it uses these growth and death rates to
#' determine the steady-state abundances of the selected species.
#' 
#' If a non-zero emigration rate is set for the species and is chosen so high
#' that the species cannot sustain it, then an error is thrown.
#'
#' If the species parameters `d_over_g` is set, then the diffusion rate is
#' calculated from the growth rate at the smallest size and the `d_over_g`
#' parameter. If `d_over_g` is not set, then the diffusion is set to zero.
#'
#' The result of applying this function is of course not a multi-species steady
#' state, because after changing the abundances of the selected species the
#' growth and death rates will have changed.
#'
#' @param params A MizerParams object
#' @param species The species to be selected. Optional. By default all target
#'   species are selected. A vector of species names, or a numeric vector with
#'   the species indices, or a logical vector indicating for each species
#'   whether it is to be selected (TRUE) or not.
#' @param keep A string determining which quantity is to be kept constant. The
#'   choices are "egg" which keeps the egg density constant, "biomass" which
#'   keeps the total biomass of the species constant and "number" which keeps
#'   the total number of individuals constant.
#' @return A MizerParams object in which the initial abundances of the selected
#'   species are changed to their single-species steady state abundances.
#' @export
steadySingleSpecies <-
    function(params, species = NULL,
             keep = c("egg", "biomass", "number")) {
        species <- valid_species_arg(params, species)

        # Prepare to keep biomass or number in case it was requested
        keep <- match.arg(keep)
        biomass <- getBiomass(params, use_cutoff = TRUE)
        number <- getN(params)

        # Get the growth and mortality rates used for the steady state calculation
        w <- w(params)
        growth_all <- getEGrowth(params)
        mort_all <- getMort(params)

        # Loop through all species and calculate their steady state abundances
        for (sp in species) {
            growth <- growth_all[sp, ]
            mort <- mort_all[sp, ]
            emigration <- emigration(params)[sp, ]
            n <- params@species_params[sp, "n"]

            w_min_idx <- params@w_min_idx[sp]
            w_max_idx <- sum(w <= params@species_params[sp, "w_max"])
            idx <- w_min_idx:w_max_idx

            # Check that species can grow to maturity at least
            w_mat_idx <- sum(w <= params@species_params[sp, "w_mat"])
            if (any(growth[w_min_idx:w_mat_idx] == 0)) {
                stop(sp, " cannot grow to maturity")
            }

            # Keep egg density constant
            N0 <- params@initial_n[sp, w_min_idx]

            d_over_g <- params@species_params[sp, "d_over_g"]
            if (is.null(d_over_g)) {
                # If d_over_g is not given, we use the old method
                idx <- w_min_idx:(w_max_idx - 1)
                params@initial_n[sp, ] <- 0
                params@initial_n[sp, w_min_idx:w_max_idx] <-
                    get_steady_state_n(growth, mort, params@dw, idx, N0)
            } else {
                message("Determining steady state with diffusion")
                sol <- solve_ode_steady_state(growth[idx], mort[idx],
                                              emigration[idx],
                                              d_over_g, N0, w[idx], n)
                if (any(sol <= 0)) {
                    stop(sp, " can not sustain the level of emigration.")
                }
                params@initial_n[sp, idx] <- sol
            }

        }

        if (any(is.infinite(params@initial_n))) {
            stop("Candidate steady state holds infinities")
        }
        if (any(is.na(params@initial_n) | is.nan(params@initial_n))) {
            stop("Candidate steady state holds non-numeric values")
        }

        if (keep == "biomass") {
            factor <- biomass / getBiomass(params, use_cutoff = TRUE)
            params@initial_n <- params@initial_n * factor
        }
        if (keep == "number") {
            factor <- number / getN(params)
            params@initial_n <- params@initial_n * factor
        }

        params@time_modified <- lubridate::now()
        params
    }

# Helper function to solve steady state ODE
solve_ode_steady_state <- function(growth, mort, emigration,
                                   d_over_g, N0, w, n) {
    N <- length(w) - 2  # Number of internal points
    if (length(mort) != N + 2 || length(growth) != N + 2) {
        stop("Growth and mortality vectors must have the same length.")
    }
    h <- log(w[2] / w[1])  # Size step

    # Determine diffusion rate so that at offspring size we have
    # d(w)=d_over_g * g(w) * w
    g_0 <- growth[1] / w[1]^n
    d_0 <- d_over_g * g_0
    diffusion <- d_0 * w^(n + 1)

    # Rescalings to convert from w to x = log(w/w_0)
    n0 <- N0 * w[1]  # Initial condition for the ODE
    dtilde <- diffusion / w^2
    gtilde <- growth / w - 0.5 * dtilde
    etilde <- emigration * w

    # Coefficients for the tridiagonal matrix
    U <- (dtilde / 2)[3:(N+1)]    # Upper diagonal
    L <- (dtilde / 2 + h * gtilde)[2:N]     # Lower diagonal
    D <- (-dtilde - h * gtilde - h^2 * mort)[2:(N+1)]  # Main diagonal

    # Solve using sparse LU decomposition
    # Create sparse tridiagonal matrix
    A <- Matrix::bandSparse(N, k = c(-1, 0, 1), diagonals = list(L, D, U))
    # Right-hand side
    b <- h^2 * etilde[2:(N+1)]
    b[1] <- b[1] - L[1] * n0

    solution <- numeric(N + 2)
    solution[1] <- n0  # Set the first value to n0
    solution[2:(N + 1)] <- Matrix::solve(A, b)
    solution[N + 2] <- 0  # Boundary condition at x_max

    # Convert back to original size space
    solution <- solution / w

    return(solution)
}
