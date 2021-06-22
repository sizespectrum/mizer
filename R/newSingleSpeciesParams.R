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
#' Although the steady state is often stable without imposing a stock
#' recruitment relationship, the function can set a Beverton-Holt type stock
#' recruitment relationship that imposes a maximal reproduction rate that is a
#' multiple of the recruitment rate at steady state. That multiple is set by the
#' argument \code{R_factor}.
#'
#' @param species_name A string with a name for the species. Will be used in
#'   plot legends.
#' @param w_inf Asymptotic size of species
#' @param w_min Egg size of species
#' @param eta Ratio between maturity size \code{w_mat} and asymptotic size
#'   \code{w_inf}. Default is 10^(-0.6), approximately 1/4. Ignored if
#'   \code{w_mat} is supplied explicitly.
#' @param w_mat Maturity size of species. Default value is
#'   \code{eta * w_inf}.
#' @param k_vb The von Bertalanffy growth parameter.
#' @inheritParams newTraitParams
#' @param version A string specifying the version of mizer. If you want to make
#'   sure that your code will still set up the model with the exact same default
#'   values even if it is run with a future versions of mizer that would choose
#'   defaults differently, then set this argument to "v2.2.0".
#' @export
#' @return An object of type \code{MizerParams}
#' @family functions for setting up models
#' @examples
#' params <- newSingleSpeciesParams()
#' sim <- project(params, t_max = 5, effort = 0)
#' plotSpectra(sim)
newSingleSpeciesParams <- 
    function(species_name = "Target species",
             w_inf = 100,
             w_min = 0.001,
             eta = 10^(-0.6),
             w_mat = w_inf * eta,
             no_w = log10(w_inf / w_min) * 20 + 1,
             n = 3/4,
             p = n,
             lambda = 2.05,
             kappa = 0.005,
             alpha = 0.4,
             k_vb = 1,
             beta = 100,
             sigma = 1.3,
             f0 = 0.6,
             fc = 0.25,
             ks = NA,
             gamma = NA,
             ext_mort_prop = 0,
             R_factor = 4,
             version) {
    assert_that(is.string(species_name), length(species_name) == 1)
    no_sp <- 1
    ## Much of the following code is copied from newTraitParams

    ## Check validity of parameters ----
    if (ext_mort_prop >= 1 || ext_mort_prop < 0) {
        stop("ext_mort_prop can not take the value ", ext_mort_prop,
             " because it should be the proportion of the total mortality",
             " coming from sources other than predation.")
    }
    if (R_factor <= 1) {
        message("R_factor needs to be larger than 1. Setting R_factor=1.01")
        R_factor <- 1.01
    }
    if (w_min <= 0) {
        stop("The egg size w_min must be greater than zero.")
    }
    no_w <- round(no_w)
    if (no_w < 1) {
        stop("The number of size bins no_w must be a positive integer")
    }
    if (no_w < log10(w_inf/w_min)*5) {
        no_w <- round(log10(w_inf / w_min) * 5 + 1)
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
    if (w_mat >= w_inf) {
        stop("The maturity size w_mat must be ",
             "smaller the maximum size w_inf")
    }
    if (!all(c(n, lambda, kappa, alpha, k_vb, beta, sigma, f0) > 0)) {
        stop("The parameters n, lambda, kappa, alpha, k_vb, beta, sigma ",
             "and f0, if supplied, need to be positive.")
    }
    if (!is.na(fc) && (fc < 0 || fc > f0)) {
        stop("The critical feeding level must lie between 0 and f0")
    }
    if (!is.na(gamma)) {  # If gamma is supplied, f0 is ignored
        f0 <- NA
    }

    ## Build Params Object ----
    erepro <- 0.1  # Will be changed later to achieve coexistence
    species_params <- data.frame(
        species = species_name,
        w_min = w_min,
        w_inf = w_inf,
        w_mat = w_mat,
        w_min_idx = 1,
        k_vb =  k_vb,
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
            max_w = w_inf,
            w_pp_cutoff = w_inf,
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
    h <- params@species_params$h
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
    i_inf <- sum(params@w <= w_inf)  # index of asymptotic size
    idx <- 1:(i_inf - 1)
    idxs <- 1:i_inf
    gg <- hbar * w^n * (1 - params@psi[1, ])  # Growth rate
    # Steady state solution of the upwind-difference scheme used in project
    initial_n[1, idxs] <- c(1, cumprod(gg[idx] / ((gg + mumu * dw)[idx + 1])))

    # The resource was already set up by newMultispeciesParams()
    initial_n_pp <- params@initial_n_pp

    # Normalise abundance so that the maximum ratio between the species
    # abundance and community abundance is 1/2
    fish <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)
    imax <- which.max(initial_n[1, ] / initial_n_pp[fish]) # index in fish spectrum
    pmax <- imax + length(params@w_full) - length(params@w) # corresponding resource index
    initial_n <- initial_n / initial_n[1, imax] * initial_n_pp[pmax] / 2
    params@initial_n <- initial_n

    ## Set reproduction to meet boundary condition ----
    params@species_params$erepro <- params@species_params$erepro *
        get_required_reproduction(params) / getRDI(params) 
    
    params <- setBevertonHolt(params, R_factor = R_factor)

    return(params)
}
