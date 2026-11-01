#' Steady-state abundance per recruit
#'
#' `r lifecycle::badge("experimental")`
#' Calculates the steady-state abundance density that a recruitment of one
#' individual per year produces, with growth, mortality and diffusion held
#' fixed.
#'
#' With the rates held fixed, the species decouple and the steady-state
#' abundance \eqn{N_i(w)} of species \eqn{i} is proportional to its rate of
#' recruitment \eqn{R_i}. The abundance per recruit is the constant of
#' proportionality, \eqn{N_i(w) / R_i}. Integrated over a size range, it gives
#' the expected time that a recruit spends in that range. Weighted with the
#' rate of egg production, it gives the lifetime reproductive output
#' [getR0()].
#'
#' The steady state is found by solving the same discretised transport
#' equation that [project()] solves, with the scheme selected by
#' [second_order_w()] and including size diffusion, so it is exact for the
#' model's numerics. The abundance is held at zero above each species' `w_max`,
#' as in [project()].
#'
#' The rates default to those in the initial state of `params`. A species that
#' is absent from that state still has a growth and mortality rate, the ones a
#' few individuals of it would experience, so it has an abundance per recruit
#' too. To fill in a size range that a model does not resolve with its
#' steady-state profile, scale the abundance per recruit by the recruitment.
#'
#' @param params A MizerParams object
#' @param e_growth An array (species x size) with the energy available for
#'   growth, as calculated by [getEGrowth()]. Defaults to the growth rate in
#'   the initial state of `params`.
#' @param mort An array (species x size) with the mortality rate, as
#'   calculated by [getMort()]. Defaults to the mortality rate in the initial
#'   state of `params`.
#' @param diffusion An array (species x size) with the diffusion rate, as
#'   calculated by [getDiffusion()]. Defaults to the diffusion rate in the
#'   initial state of `params`.
#' @param ... Unused.
#' @return An [ArraySpeciesBySize] object (species x size) with the abundance
#'   density per recruit, in units of year/g.
#' @seealso [getR0()], [steadySingleSpecies()]
#' @export
#' @examples
#' n1 <- getNPerRecruit(NS_params)
#' plot(n1)
#'
#' # The steady state for the current rates and recruitment
#' n <- n1 * getRDD(NS_params)
getNPerRecruit <- function(params, ...) {
    UseMethod("getNPerRecruit")
}

#' @rdname getNPerRecruit
#' @export
getNPerRecruit.MizerParams <- function(params,
                                       e_growth = getEGrowth(params),
                                       mort = getMort(params),
                                       diffusion = getDiffusion(params),
                                       ...) {
    params <- validParams(params)
    e_growth <- as_species_by_size(e_growth, params, "e_growth")
    mort <- as_species_by_size(mort, params, "mort")
    diffusion <- as_species_by_size(diffusion, params, "diffusion")
    no_sp <- nrow(params@species_params)

    # With the rates fixed, the steady state of each species is proportional
    # to its egg density, so solve with unit egg density ...
    n <- get_steady_state_n(params, e_growth, mort, diffusion,
                            N0 = rep(1, no_sp))
    # ... and divide by the recruitment that holds that solution steady. The
    # coefficients are recalculated at the solution because the flux limiter
    # depends on it.
    coefs <- get_transport_coefs(params, n, e_growth, mort, dt = 1,
                                 recruitment_flux = numeric(no_sp),
                                 d = diffusion,
                                 flux_limiter = flux_limiter_scheme(params))
    n <- n / required_recruitment(params, n, coefs)
    dimnames(n) <- dimnames(params@initial_n)

    ArraySpeciesBySize(n, value_name = "Abundance per recruit",
                       units = "year/g", type = "density", params = params,
                       representation = "average")
}

#' Lifetime reproductive output
#'
#' `r lifecycle::badge("experimental")`
#' Calculates the lifetime reproductive output \eqn{R_0} of each species: the
#' number of eggs that one recruit produces over its life, with growth,
#' mortality, diffusion and the investment into reproduction held fixed.
#'
#' \eqn{R_0} decides whether a species can increase when it is rare: at low
#' abundance a species grows if \eqn{R_0 > 1} and declines if \eqn{R_0 < 1}.
#' It is the density-independent rate of egg production (see [getRDI()]) of
#' the abundance per recruit \eqn{n_1(w)} (see [getNPerRecruit()]):
#' \deqn{R_0 = \frac{\epsilon}{2 w_{min}} \int n_1(w) E_r(w) \psi(w) \, dw.}{R_0 = (\epsilon/(2 w_min)) \int n_1(w) E_r(w) \psi(w) dw.}
#' It ignores the stock-recruitment relationship, which is exact at low
#' abundance, where density dependence does not act.
#'
#' Both the abundance per recruit and the integral are calculated with the
#' discretisation that [project()] uses, selected by [second_order_w()], and
#' include size diffusion. A hand-written calculation that follows the
#' first-order upwind recursion for the abundance per recruit is only right for
#' the default scheme without diffusion.
#'
#' @section Relation to the steady state:
#' In a steady state the recruitment is the one returned by
#' [getRequiredRDD()], so there
#' \deqn{R_0 = \frac{R_{di}}{R_{dd}},}{R_0 = R_di / R_dd,}
#' the ratio of `getRDI(params)` to `getRequiredRDD(params)`. With
#' density-independent reproduction (`R_max = Inf`) every steady state
#' therefore has \eqn{R_0 = 1}. With Beverton-Holt reproduction the
#' recruitment is reduced by the factor \eqn{1 - r}, where \eqn{r} is the
#' reproduction level (see [reproduction_level()]), so a steady state has
#' \eqn{R_0 = 1 / (1 - r)}, not 1.
#'
#' The ratio `getRDI(params) / getRequiredRDD(params)` is \eqn{R_0} only in a
#' steady state, and it is 0/0 for a species that is absent. `getR0()` needs
#' neither, because it holds the rates fixed and solves for the abundance per
#' recruit. Because \eqn{R_0} is proportional to the reproductive efficiency
#' `erepro`, dividing `erepro` by \eqn{R_0} makes the steady state for the
#' current rates a steady state of a model with density-independent
#' reproduction.
#'
#' @section Rates from another state:
#' The rate arguments default to the rates in the initial state of `params`,
#' but each can be taken from a different state. For example
#' `getR0(params, mort = getMort(other_params))` combines the growth and
#' reproduction that the food in `params` supports with the mortality in
#' `other_params`. This separates the effect of food from that of predation on
#' whether a species can invade.
#'
#' The rate of egg production is calculated with the model's RDI function, as
#' in [getRDI()]. This assumes that it is proportional to the abundance, as
#' [mizerRDI()] is.
#'
#' @inheritParams getNPerRecruit
#' @inheritParams valid_species_arg
#' @param e_repro An array (species x size) with the energy invested into
#'   reproduction, as calculated by [getERepro()]. Defaults to the investment
#'   in the initial state of `params`.
#' @return A named vector with the lifetime reproductive output of each
#'   selected species.
#' @seealso [getNPerRecruit()], [getRDI()], [getRequiredRDD()],
#'   [reproduction_level()]
#' @export
#' @examples
#' # NS_params is close to a Beverton-Holt steady state, so R0 is close to
#' # 1 / (1 - r)
#' getR0(NS_params)
#' 1 / (1 - reproduction_level(NS_params))
#'
#' # R0 with the same food but with twice the mortality
#' getR0(NS_params, mort = 2 * getMort(NS_params))
#'
#' # R0 of Cod when it is absent, where RDI / required recruitment is 0/0
#' params <- NS_params
#' initialN(params)["Cod", ] <- 0
#' getR0(params, species = "Cod")
getR0 <- function(params, ...) {
    UseMethod("getR0")
}

#' @rdname getR0
#' @export
getR0.MizerParams <- function(params, species = NULL,
                              e_growth = getEGrowth(params),
                              mort = getMort(params),
                              diffusion = getDiffusion(params),
                              e_repro = getERepro(params), ...) {
    params <- validParams(params)
    species <- valid_species_arg(params, species)
    e_growth <- as_species_by_size(e_growth, params, "e_growth")
    mort <- as_species_by_size(mort, params, "mort")
    diffusion <- as_species_by_size(diffusion, params, "diffusion")
    e_repro <- as_species_by_size(e_repro, params, "e_repro")

    n <- unclass_rate(getNPerRecruit(params, e_growth = e_growth, mort = mort,
                                     diffusion = diffusion))
    # R0 is the egg production of the abundance per recruit, calculated by the
    # model's RDI function as in getRDI().
    if (params@rates_funcs$RDI == "mizerRDI") {
        r0 <- projectRDI(params, n = n, n_pp = params@initial_n_pp,
                         n_other = params@initial_n_other, t = 0,
                         e_repro = e_repro, e_growth = e_growth, mort = mort,
                         diffusion = diffusion)
    } else {
        f <- get(params@rates_funcs$RDI)
        r0 <- f(params, n = n, n_pp = params@initial_n_pp,
                n_other = params@initial_n_other, t = 0,
                e_repro = e_repro, e_growth = e_growth, mort = mort,
                diffusion = diffusion)
    }
    names(r0) <- params@species_params$species
    r0[species]
}

#' Check a rate supplied as an argument and return it as a plain matrix
#'
#' @param x The supplied rate, which must be a numeric array (species x size)
#'   of the same shape as `initialN(params)`.
#' @param params A MizerParams object.
#' @param arg The name of the argument, for the error message.
#' @return `x` as a plain matrix with the dimnames of `initialN(params)`.
#' @noRd
as_species_by_size <- function(x, params, arg) {
    if (!is.numeric(x) || !identical(dim(x), dim(params@initial_n))) {
        stop("`", arg, "` must be a numeric array (species x size) with ",
             nrow(params@initial_n), " rows and ", ncol(params@initial_n),
             " columns, the same shape as `initialN(params)`.")
    }
    x <- unclass_rate(x)
    dimnames(x) <- dimnames(params@initial_n)
    x
}
