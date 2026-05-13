#' Calculate initial population abundances
#'
#' This function uses the model parameters and other parameters to calculate
#' initial values for the species number densities. These initial
#' abundances are calculated from the steady state of a simplified version of
#' the model.
#'
#' @param params The model parameters. An object of type
#'   \linkS4class{MizerParams}.
#' @param a `r lifecycle::badge("deprecated")` This argument was used only by
#'   the retired defaults edition 1 calculation and is now ignored.
#' @param n0_mult `r lifecycle::badge("deprecated")` This argument was used
#'   only by the retired defaults edition 1 calculation and is now ignored.
#' @export
#' @concept helper
#' @return An `ArraySpeciesBySize` object (species x size) of population abundances.
#' @examples
#' init_n <- get_initial_n(NS_params)
get_initial_n <- function(params, n0_mult = NULL, a = 0.35) {
    if (!is(params,"MizerParams"))
        stop("params argument must of type MizerParams")
    if (!missing(a)) {
        lifecycle::deprecate_warn(
            "3.0.0",
            "get_initial_n(a)",
            details = "`a` was used only by the retired defaults edition 1 calculation and is now ignored."
        )
    }
    if (!is.null(n0_mult)) {
        lifecycle::deprecate_warn(
            "3.0.0",
            "get_initial_n(n0_mult)",
            details = "`n0_mult` was used only by the retired defaults edition 1 calculation and is now ignored."
        )
    }
    no_sp <- nrow(params@species_params)

    p <- params
    p@initial_n[] <- 0
    p@initial_n_pp <- p@resource_params$kappa *
        p@w_full ^ (-p@resource_params$lambda)
    p@interaction[] <- 0
    income <- getEReproAndGrowth(p) + p@metab
    mort <- getFMort(p)
    growth <- getEGrowth(p)
    diffusion <- getDiffusion(p)
    N0_vec <- numeric(no_sp)

    for (i in seq_len(no_sp)) {
        # At small sizes the income should be A w^n. Determine A
        # Use w_min_idx + 1 in case user has implemented reduced growth
        # for the smallest size class (see e.g. #241)
        iw <- p@w_min_idx[i] + 1
        A <- income[i, iw] / (p@w[iw] ^ p@species_params[[i, "n"]])

        mort[i, ] <- mort[i, ] + 0.4 * A * p@w ^ (p@species_params[[i, "n"]] - 1)

        # We start with an arbitrary population at the smallest size class
        N0_vec[i] <- 1
    }

    n_exact_matrix <- get_steady_state_n(p, growth, mort, diffusion, N0_vec)

    for (i in seq_len(no_sp)) {
        idxs <- p@w_min_idx[i]:(min(which(c(growth[i, ], 0) <= 0)) - 1)
        # Steady state solution of the upwind-difference scheme used in project
        p@initial_n[i, idxs] <- n_exact_matrix[i, idxs]
    }
    p <- matchBiomasses(p)
    return(ArraySpeciesBySize(p@initial_n, value_name = "Number density",
                             params = p))
}
