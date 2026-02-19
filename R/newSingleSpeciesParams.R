#' Set up parameters for a single species in a power-law background
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' This functions creates a \code{MizerParams} object with a single
#' species. This species is embedded in a fixed power-law community spectrum
#' \deqn{N_c(w) = \kappa w^{-\lambda}}
#' This community provides the food income for the species. Cannibalism is
#' switched off. The predation mortality arises only from the predators in the
#' power-law community and it is assumed that the predators in the community
#' have the same feeding parameters as the foreground species. The function has
#' many arguments, all of which have default values.
#'
#' @details
#' In addition to setting up the parameters, this function also sets up an
#' initial condition that is close to steady state, under the assumption of
#' no fishing.
#'
#' @param species_name A string with a name for the species. Will be used in
#'   plot legends.
#' @param w_max Maximum size of species
#' @param w_min Egg size of species
#' @param eta Ratio between maturity size \code{w_mat} and maximum size
#'   \code{w_max}. Default is 10^(-0.6), approximately 1/4. Ignored if
#'   \code{w_mat} is supplied explicitly.
#' @param w_mat Maturity size of species. Default value is
#'   \code{eta * w_max}.
#' @param k_vb `r lifecycle::badge("deprecated")` The von Bertalanffy growth
#'   parameter.
#' @param w_inf `r lifecycle::badge("deprecated")` The argument has been
#'   renamed to `w_max`.
#' @inheritParams newTraitParams
#' @export
#' @importFrom lifecycle deprecated
#' @return An object of type \code{MizerParams}
#' @family functions for setting up models
#' @examples
#' params <- newSingleSpeciesParams()
#' sim <- project(params, t_max = 5, effort = 0)
#' plotSpectra(sim)
newSingleSpeciesParams <-
    function(species_name = "Target species",
             w_max = 100,
             w_min = 0.001,
             eta = 10^(-0.6),
             w_mat = w_max * eta,
             no_w = log10(w_max / w_min) * 20 + 1,
             n = 3 / 4,
             p = n,
             lambda = 2.05,
             kappa = 0.005,
             alpha = 0.4,
             h = 30,
             beta = 100,
             sigma = 1.3,
             f0 = 0.6,
             fc = 0.25,
             ks = NA,
             gamma = NA,
             ext_mort_prop = 0,
             reproduction_level = 0,
             R_factor = deprecated(),
             w_inf = deprecated(),
             k_vb = deprecated()) {
        if (lifecycle::is_present(w_inf)) {
            lifecycle::deprecate_warn(
                when = "2.4.0.0",
                what = "newSingleSpeciesParams(w_inf)",
                with = "newSingleSpeciesParams(w_max)"
            )
            w_max <- w_inf
            w_mat <- w_max * eta
            no_w <- log10(w_max / w_min) * 20 + 1
        }
        if (lifecycle::is_present(k_vb)) {
            lifecycle::deprecate_warn(
                when = "2.4.0.0",
                what = "newSingleSpeciesParams(k_vb)",
                details = "Specify 'h' instead."
            )
        }
    assert_that(is.string(species_name), length(species_name) == 1)

    ## Much of the following code is copied from newTraitParams

    ## Check validity of parameters ----
    if (ext_mort_prop >= 1 || ext_mort_prop < 0) {
        stop("ext_mort_prop can not take the value ", ext_mort_prop,
             " because it should be the proportion of the total mortality",
             " coming from sources other than predation.")
    }
    if (w_min <= 0) {
        stop("The egg size w_min must be greater than zero.")
    }
    no_w <- round(no_w)
    if (no_w < 1) {
        stop("The number of size bins no_w must be a positive integer")
    }
    if (no_w < log10(w_max/w_min)*5) {
        no_w <- round(log10(w_max / w_min) * 5 + 1)
        message(paste("Increased no_w to", no_w, "so that there are 5 bins ",
                      "for an interval from w and 10w."))
    }
    if (no_w > 10000) {
        message("Running a simulation with ", no_w,
                " size bins is going to be very slow.")
    }
    if (w_min >= w_mat) {
        stop("The egg size w_min must be smaller than ",
             "the maturity size w_mat")
    }
    if (w_mat >= w_max) {
        stop("The maturity size w_mat must be ",
             "smaller the maximum size w_max")
    }
    if (!all(c(n, lambda, kappa, alpha, h, beta, sigma, f0) > 0)) {
        stop("The parameters n, lambda, kappa, alpha, h, beta, sigma ",
             "and f0, if supplied, need to be positive.")
    }
    if (!is.na(fc) && (fc < 0 || fc > f0)) {
        stop("The critical feeding level must lie between 0 and f0")
    }
    if (!is.na(gamma)) {  # If gamma is supplied, f0 is ignored
        f0 <- NA
    }
    if (lifecycle::is_present(R_factor)) {
        lifecycle::deprecate_soft("2.3.0", "newTraitParams(R_factor)",
                                  "newTraitParams(reproduction_level)",
                                  "Set `reproduction_level = 1 / R_factor`.")
        reproduction_level <- 1 / R_factor
    }

    ## Build Params Object ----
    erepro <- 0.1  # Will be changed later to achieve coexistence
    species_params <- data.frame(
        species = species_name,
        w_min = w_min,
        w_max = w_max,
        w_mat = w_mat,
        w_min_idx = 1,
        h = h,
        gamma = gamma,
        ks = ks,
        f0 = f0,
        fc = fc,
        beta = beta,
        sigma = sigma,
        z0 = 0,
        alpha = alpha,
        erepro = erepro,
        stringsAsFactors = FALSE
    )
    params <-
        suppressMessages(newMultispeciesParams(
            species_params,
            min_w = w_min,
            no_w = no_w,
            max_w = w_max,
            w_pp_cutoff = Inf,
            lambda = lambda,
            kappa = kappa,
            n = n,
            p = p,
            resource_dynamics = "resource_constant"
        ))
    # No cannibalism
    params@interaction[] <- 0

    w <- params@w
    dw <- params@dw
    h <- params@species_params[["h"]]
    ks <- params@species_params$ks
    f0 <- get_f0_default(params)

    ## Construct steady state solution ----

    # Get constants for steady-state solution
    # Predation mortality rate coefficient
    mu0 <- get_power_law_mort(params)
    # Add backgound mortality rate
    mu0 <- mu0 / (1 - ext_mort_prop)
    hbar <- alpha * h * f0 - ks
    if (hbar < 0) {
        stop("The feeding level is not sufficient to maintain the fish.")
    }
    pow <- mu0 / hbar / (1 - n)
    if (pow < 1) {
        message("The ratio of death rate to growth rate is too small, leading to
                an accumulation of fish at their largest size.")
    }

    initial_n <- params@psi  # get array with correct dimensions and names
    initial_n[, ] <- 0
    mumu <- mu0 * w^(n - 1)  # Death rate
    params@mu_b[] <- mumu
    comment(params@mu_b) <- "power-law"
    i_inf <- sum(params@w <= w_max)  # index of maximum size
    idx <- 1:(i_inf - 1)
    idxs <- 1:i_inf
    gg <- hbar * w^n * (1 - params@psi[1, ])  # Growth rate
    # Steady state solution of the upwind-difference scheme used in project
    initial_n[1, idxs] <- get_steady_state_n(gg, mumu, dw, params@diffusion[1, ], idx)

    # The resource was already set up by newMultispeciesParams()
    initial_n_pp <- params@initial_n_pp

    # Normalise abundance so that the maximum ratio between the species
    # abundance and community abundance is 1/2
    fish <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)
    # index in fish spectrum
    imax <- which.max(initial_n[1, ] / initial_n_pp[fish])
    # corresponding resource index
    pmax <- imax + length(params@w_full) - length(params@w)
    initial_n <- initial_n / initial_n[1, imax] * initial_n_pp[pmax] / 2
    params@initial_n <- initial_n

    ## Set reproduction to meet boundary condition ----
    params@species_params$erepro <- params@species_params$erepro *
        getRequiredRDD(params) / getRDI(params)
    params@given_species_params$erepro <- params@species_params$erepro

    params <- setBevertonHolt(params, reproduction_level = reproduction_level)

    return(params)
}
