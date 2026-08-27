# Summary and indicator functions for mizer package

# Copyright 2012 Finlay Scott and Julia Blanchard.
# Copyright 2018 Gustav Delius and Richard Southwell.
# Distributed under the GPL 3 or later
# Maintainer: Gustav Delius, University of York, <gustav.delius@york.ac.uk>

# Soundtrack: The Definitive Lead Belly

# Summary functions ####

#' Description of summary functions
#'
#' Mizer provides a range of functions to summarise the results of a simulation.
#'
#' A list of available summary functions is given in the table below.
#' \tabular{lll}{
#'   Function \tab Returns \tab Description \cr
#'   [getDiet()] \tab Three dimensional array (predator x size x prey) \tab Diet of predator at size, resolved by prey species \cr
#'   [getTrophicLevel()] \tab `ArraySpeciesBySize` (species x size) \tab Trophic level of individuals at size, accounting for ontogenetic diet shifts \cr
#'   [getTrophicLevelBySpecies()] \tab Named vector (species) \tab Consumption-rate-weighted mean trophic level of each species \cr
#'   [getSSB()] \tab Two dimensional array (time x species) \tab Total Spawning Stock Biomass (SSB) of each species through time where SSB is calculated as the sum of weight of all mature individuals. \cr
#'   [getBiomass()] \tab Two dimensional array (time x species) \tab Total biomass of each species through time. \cr
#'   [getN()] \tab Two dimensional array (time x species) \tab Total abundance of each species through time. \cr
#'   [getFeedingLevel()] \tab Three dimensional array (time x species x size) \tab Feeding level of each species by size through time. \cr
#'   \code{\link{getM2}} \tab Three dimensional array (time x species x size) \tab The predation mortality imposed on each species by size through time. \cr
#'   [getFMort()] \tab Three dimensional array (time x species x size) \tab Total fishing mortality on each species by size through time. \cr
#'   [getFMortGear()] \tab Four dimensional array (time x gear x species x size) \tab Fishing mortality on each species by each gear at size through time. \cr
#'   [getYieldGear()] \tab Three dimensional array (time x gear x species) \tab Total yield by gear and species through time. \cr
#'   [getYield()] \tab Two dimensional array (time x species) \tab Total yield of each species across all gears through time. \cr
#'   [sizeIntegral()] \tab Named vector (species) or two dimensional array (time x species) \tab Any integral over the size spectrum, from which all of the above are built. Use it to write your own summary function. \cr
#' }
#'
#' @section Writing your own summary function:
#' The entry point is [sizeIntegral()]. It selects the size range, applies the
#' bin-averaging appropriate to the model's [second_order_w()] setting and wraps
#' the result in the right mizer array class, so a summary function built on it
#' is automatically consistent with the quadrature the model is actually using.
#' Pass the whole weighting factor \eqn{K(w)} evaluated on the size grid, but
#' neither the bin widths `params@dw` nor any bin-averaging of your own:
#' `sizeIntegral()` supplies both.
#'
#' If your quantity involves the predation kernel, take the kernel from
#' [encounter_kernel()] rather than from [pred_kernel()], and pair it with the
#' plain point prey weight `params@w_full * params@dw_full`. That weight is a
#' normalisation which the kernel construction is built to cancel, not a
#' quadrature weight, so it is the one place where you must *not* bin-average.
#' Pairing the point-sampled [pred_kernel()] with a bin-averaged prey weight
#' applies the prey-bin integral twice.
#'
#' @seealso [indicator_functions], [plotting_functions]
#' @name summary_functions
NULL

#' Get diet of predator at size, resolved by prey species
#'
#' Calculates the rate at which a predator of a particular species and size
#' consumes biomass of each prey species, resource, and other components of the
#' ecosystem. Returns either the rates in grams/year or the proportion of the
#' total consumption rate.
#'
#' The rates \eqn{D_{ij}(w)} at which a predator of species \eqn{i}
#' and size \eqn{w} consumes biomass from prey species \eqn{j} are
#' calculated from the predation kernel \eqn{\phi_i(w, w_p)},
#' the search volume \eqn{\gamma_i(w)}, the feeding level \eqn{f_i(w)}, the
#' species interaction matrix \eqn{\theta_{ij}} and the prey abundance density
#' \eqn{N_j(w_p)}:
#' \deqn{
#' D_{ij}(w, w_p) = (1-f_i(w)) \gamma_i(w) \theta_{ij}
#' \int N_j(w_p) \phi_i(w, w_p) w_p dw_p.
#' }
#' The prey index \eqn{j} runs over all species and the resource.
#'
#' Extra columns are added for the external encounter rate and for any extra
#' ecosystem components in your model for which you have defined an encounter
#' rate function. These encounter rates are multiplied by \eqn{1-f_i(w)} to give
#' the rate of consumption of biomass from these extra components.
#'
#' This function performs the same integration as [getEncounter()] but does not
#' aggregate over prey species, and multiplies by \eqn{1-f_i(w)} to get the
#' consumed biomass rather than the available biomass. Outside the range of
#' sizes for a predator species the returned rate is zero. Summing the result
#' of `getDiet(proportion = FALSE)` over prey therefore reproduces
#' `getEncounter(params) * (1 - getFeedingLevel(params))`, whichever quadrature
#' scheme the model uses (see [second_order_w()]).
#'
#' @param object A \linkS4class{MizerParams} or \linkS4class{MizerSim} object.
#' @param proportion If TRUE (default) the function returns the diet as a
#'   proportion of the total consumption rate. If FALSE it returns the
#'   consumption rate in grams per year.
#' @param ... Additional arguments that depend on the class of `object`.
#'
#'   **For a \linkS4class{MizerParams} object:**
#'   \describe{
#'     \item{`n`}{A matrix of species abundances (species x size). Defaults to
#'       the initial abundances stored in `object`.}
#'     \item{`n_pp`}{A vector of the resource abundance by size. Defaults to the
#'       initial resource abundance stored in `object`.}
#'     \item{`n_other`}{A named list of the abundances of other dynamical
#'       components. Defaults to the initial values stored in `object`.}
#'   }
#'
#'   **For a \linkS4class{MizerSim} object:**
#'   \describe{
#'     \item{`time_range`}{The time range over which to return the diet. Either
#'       a vector of values, a vector of min and max time, or a single value.
#'       Defaults to the whole time range of the simulation.}
#'     \item{`drop`}{If `TRUE` then any dimension of length 1 is removed from
#'       the returned array.}
#'   }
#'
#' @return
#' * `MizerParams`: An array (predator species x predator size x
#'   (prey species + resource + other components)). Dimnames are "predator",
#'   "w", and "prey".
#' * `MizerSim`: A four-dimensional array (time x predator species x predator
#'   size x prey) with the diet at each selected saved time step. If
#'   `drop = TRUE` then dimensions of length 1 are removed.
#'
#' @export
#' @family summary functions
#' @concept summary_function
#' @seealso [plotDiet()]
#' @examples
#' diet <- getDiet(NS_params)
#' str(diet)
#' \donttest{
#' # For a MizerSim the diet is returned at each saved time step
#' sim <- project(NS_params, t_max = 20, effort = 0.5)
#' # Diet at the saved time steps over years 15 - 20
#' diet <- getDiet(sim, time_range = c(15, 20))
#' str(diet)
#' }
getDiet <- function(object, proportion = TRUE, ...) {
    UseMethod("getDiet")
}
#' @export
getDiet.MizerParams <- function(object,
                    proportion = TRUE,
                    n = initialN(object),
                    n_pp = initialNResource(object),
                    n_other = initialNOther(object), ...) {
    # The code is based on that for getEncounter()
    params <- validParams(object)
    species <- params@species_params$species
    no_sp <- length(species)
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)
    no_other <- length(params@other_encounter)
    other_names <- names(params@other_encounter)
    assert_that(identical(dim(n), c(no_sp, no_w)),
                length(n_pp) == no_w_full)
    diet <- array(0, dim = c(no_sp, no_w, no_sp + 2 + no_other),
                  dimnames = list("predator" = species,
                                  "w" = dimnames(params@initial_n)$w,
                                  "prey" = c(species,
                                             "Resource", "External",
                                             other_names)))
    # idx_sp are the index values of object@w_full such that
    # object@w_full[idx_sp] = object@w
    idx_sp <- (no_w_full - no_w + 1):no_w_full

    # The prey vectors below carry the plain point weight `w * dw`, exactly as
    # in mizerEncounter(). This is deliberate also when second-order
    # bin-averaging is switched on: there the prey-bin integral has already been
    # folded into the kernel by setPredKernel(), which builds
    # `ft_pred_kernel_e` as the kernel integrated over the prey bin and divides
    # out the `w * dw` that the prey vector supplies. Weighting the prey vector
    # by the bin-averaged `w` as well would apply that quadrature twice and
    # inflate the diet by (1 + beta) / 2 (issue #474). Summed over prey, the
    # diet must reproduce `getEncounter() * (1 - getFeedingLevel())` under both
    # schemes.

    # If the user has set a custom kernel we can not use fft.
    if (!is.null(comment(params@pred_kernel))) {
        # pred_kernel is predator species x predator size x prey size
        # We want to multiply this by the prey abundance, which is
        # prey species by prey size, sum over prey size. We use matrix
        # multiplication for this. Then we multiply 1st and 3rd
        ae <- matrix(params@pred_kernel[, , idx_sp, drop = FALSE],
                     ncol = no_w) %*%
            t(sweep(n, 2, params@w * params@dw, "*"))
        diet[, , 1:no_sp] <- ae
        # Eating the resource
        diet[, , no_sp + 1] <- rowSums(sweep(
            params@pred_kernel, 3, params@dw_full * params@w_full * n_pp, "*"),
            dims = 2)
    } else {
        prey <- matrix(0, nrow = no_sp + 1, ncol = no_w_full)
        prey[1:no_sp, idx_sp] <- sweep(n, 2, params@w * params@dw, "*")
        prey[no_sp + 1, ] <- n_pp * params@w_full * params@dw_full
        ft <- array(rep(params@ft_pred_kernel_e, times = no_sp + 1) *
                        rep(mvfft(t(prey)), each = no_sp),
                    dim = c(no_sp, no_w_full, no_sp + 1))
        # We now have an array predator x wave number x prey
        # To Fourier transform back we need a matrix of wave number x everything
        ft <- matrix(aperm(ft, c(2, 1, 3)), nrow = no_w_full)
        ae <- array(Re(mvfft(ft, inverse = TRUE) / no_w_full),
                    dim = c(no_w_full, no_sp, no_sp + 1))
        ae <- ae[idx_sp, , , drop = FALSE]
        ae <- aperm(ae, c(2, 1, 3))
        # Due to numerical errors we might get negative or very small entries that
        # should be 0
        ae[ae < 1e-18] <- 0
        diet[, , 1:(no_sp + 1)] <- ae
    }
    # Multiply by interaction matrix, including resource, and then by
    # search volume
    inter <- cbind(params@interaction, params@species_params$interaction_resource)
    diet[, , 1:(no_sp + 1)] <- sweep(sweep(diet[, , 1:(no_sp + 1), drop = FALSE],
                                           c(1, 3), inter, "*"),
                                     c(1, 2), params@search_vol, "*")
    # Add diet from external sources
    diet[, , no_sp + 2] <- params@ext_encounter
    # Add diet from other components
    for (i in seq_along(params@other_encounter)) {
        diet[, , no_sp + 2 + i] <-
            do.call(params@other_encounter[[i]],
                    list(params = params,
                         n = n, n_pp = n_pp, n_other = n_other,
                         component = names(params@other_encounter)[[i]]))
    }

    # Correct for satiation and keep only entries corresponding to fish sizes
    f <- getFeedingLevel(params, n, n_pp)
    fish_mask <- n > 0
    diet <- sweep(diet, c(1, 2), (1 - f) * fish_mask, "*")
    if (proportion) {
        total <- rowSums(diet, dims = 2)
        diet <- sweep(diet, c(1, 2), total, "/")
        diet[is.nan(diet)] <- 0
    }
    return(diet)
}

#' @export
getDiet.MizerSim <- function(object, proportion = TRUE,
                             time_range, drop = FALSE, ...) {
    sim <- object
    time_elements <- get_sim_rate_time_elements(sim, time_range)

    # Compute the diet at each selected saved time step from that step's
    # abundances and stack the results into a time x predator x size x prey
    # array. This mirrors the other `MizerSim` rate getters, which evaluate the
    # quantity per saved step (see `get_species_size_rate_from_sim()`); the
    # plotting functions then average the computed quantity over a time range.
    diet_time <- plyr::aaply(which(time_elements), 1, function(time_idx) {
        slice <- get_sim_rate_slice(sim, time_idx)
        getDiet(sim@params, n = slice$n, n_pp = slice$n_pp,
                n_other = slice$n_other, proportion = proportion, ...)
    }, .drop = FALSE)
    names(dimnames(diet_time))[[1]] <- "time"

    diet_time[, , , , drop = drop]
}


#' Get trophic level of individuals at size
#'
#' `r lifecycle::badge("experimental")`
#' Calculates the trophic level of individuals of each species at each size,
#' assuming the system is in a steady state. The trophic level of an individual
#' is defined as 1 more than the consumption-rate-weighted average trophic level
#' of all the prey it has consumed during its lifetime up to the current size.
#' The resource is given a size-dependent trophic level (see below).
#'
#'
#' @details
#' In the traditional non-size-resolved approach, all individuals of a species
#' have the same diet composition \eqn{D_{ij}}, defined as the proportion of
#' total biomass intake of species \eqn{i} that comes from species \eqn{j}.
#' The trophic levels then satisfy
#' \deqn{T_i = 1 + \sum_j D_{ij}\,T_j,}
#' which is solved as a linear system \eqn{(I - D)\,\mathbf{T} = \mathbf{1}}.
#'
#' In mizer, diet composition changes as an individual grows, so we must
#' integrate over the individual's lifetime. Assuming a steady state so that
#' the growth rate \eqn{g_i(w)} and prey densities depend only on size and not
#' on time, we can replace the integral over time since birth by an integral
#' over weight using \eqn{dt = dw / g_i(w)}. The trophic level
#' \eqn{T_i(w)} of an individual of species \eqn{i} at weight \eqn{w} is
#' then
#' \deqn{
#' T_i(w) = 1 + \frac{
#'   \int_{w_0}^{w} \frac{1}{g_i(w')} \sum_j \int r_{ij}(w', w_p)\, T_j(w_p)\, dw_p\, dw'
#' }{
#'   \int_{w_0}^{w} \frac{1}{g_i(w')} \sum_j \int r_{ij}(w', w_p)\, dw_p\, dw'
#' },
#' }
#' where \eqn{w_0} is the egg size and \eqn{r_{ij}(w, w_p)} is the rate at
#' which a predator of species \eqn{i} at weight \eqn{w} consumes biomass from
#' prey species \eqn{j} at weight \eqn{w_p}:
#' \deqn{
#' r_{ij}(w, w_p) = \theta_{ij}\,\gamma_i(w)\,(1 - f_i(w))\,\phi_i(w/w_p)\,
#'   N_j(w_p)\,w_p.
#' }
#' The sum over \eqn{j} runs over all species and the resource. The resource is
#' assigned a size-dependent trophic level
#' \deqn{
#'   T_R(w) = \max\left(1,\; 1 + \frac{\log(w / w_R)}{\log(\beta_R)}\right),
#' }
#' where \eqn{w_R} is an average size of primary producers (which therefore have
#' trophic level 1) and \eqn{\beta_R} is an average predator/prey mass ratio for
#' the resource (for example zooplankton). This adds one trophic level for each
#' factor of \eqn{\beta_R} increase in resource size, with a floor at 1 so that
#' the resource trophic level never drops below the primary-producer level.
#' Both the numerator and the denominator (which equals the total biomass
#' consumed over the predator's lifetime from egg size to current weight
#' \eqn{w}) therefore include the resource.
#'
#' This equation can be viewed as a linear system
#' \eqn{(I - D)\,\mathbf{T} = \mathbf{1}} in which the entries of
#' \eqn{\mathbf{T}} are indexed by \eqn{(i, w)} and the matrix \eqn{D} encodes
#' the lifetime-integrated diet composition. The system is solved iteratively
#' from small to large sizes, exploiting the fact that prey are typically much
#' smaller than the predator (large predator-to-prey mass ratio), so that the
#' trophic levels of all relevant prey sizes are already known when computing
#' \eqn{T_i(w)}.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param n A matrix of species abundances (species x size). Defaults to
#'   the initial abundances stored in `params`.
#' @param n_pp A vector of the resource abundance by size. Defaults to the
#'   initial resource abundance stored in `params`.
#' @param n_other A named list of the abundances of other dynamical
#'   components. Defaults to the initial values stored in `params`.
#' @param w_R An average size (in grams) of primary producers in the resource
#'   spectrum, used to set the size-dependent resource trophic level. Defaults
#'   to `1e-10`.
#' @param beta_R An average predator/prey mass ratio for the resource spectrum,
#'   used to set the size-dependent resource trophic level. Must be greater than
#'   `1`. Defaults to `1000`.
#' @param ... Unused
#'
#' @return An `ArraySpeciesBySize` object (species x size) with the trophic
#'   level of individuals at each size. Entries below the egg size of each
#'   species are \code{NA}.
#'
#' @export
#' @family summary functions
#' @concept summary_function
#' @seealso [getTrophicLevelBySpecies()]
#' @examples
#' tl <- getTrophicLevel(NS_params)
#' plot(tl)
getTrophicLevel <- function(params,
                            n = initialN(params),
                            n_pp = initialNResource(params),
                            n_other = initialNOther(params),
                            w_R = 1e-10,
                            beta_R = 1000,
                            ...) {
    UseMethod("getTrophicLevel")
}

#' @export
getTrophicLevel.MizerParams <- function(params,
                                        n = initialN(params),
                                        n_pp = initialNResource(params),
                                        n_other = initialNOther(params),
                                        w_R = 1e-10,
                                        beta_R = 1000,
                                        ...) {
    params <- validParams(params)
    assert_that(is.number(w_R), w_R > 0,
                is.number(beta_R), beta_R > 1)
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)
    idx_sp <- (no_w_full - no_w + 1):no_w_full

    # Pre-compute rates
    encounter <- getEncounter(params, n, n_pp, n_other)  # no_sp x no_w
    feeding_level <- getFeedingLevel(params, n, n_pp, n_other)  # no_sp x no_w
    growth <- getEGrowth(params, n, n_pp, n_other)  # no_sp x no_w
    # Total consumption = (1 - f) * E, used for denominator accumulator
    consumption <- (1 - feeding_level) * encounter  # no_sp x no_w

    # Full predation kernel array (no_sp x no_w x no_w_full). The numerator
    # below is a trophic-level-weighted copy of the encounter integral, and its
    # ratio to `consumption` is only a trophic level if the two use the same
    # quadrature. encounter_kernel() therefore returns the kernel that
    # mizerEncounter() effectively uses: the point-sampled kernel on the default
    # path, and the prey-bin-integrated kernel when second-order bin-averaging
    # is on (issue #474).
    pred_kernel <- encounter_kernel(params)

    # prey_mass_tl[j, p] = N_j(w_p) * T_j(w_p) * w_p * dw_p
    # The prey-biomass weight is the plain point weight `w * dw`, matching
    # mizerEncounter(). Under second-order the prey-bin integral lives in
    # `pred_kernel`, so bin-averaging this weight as well would apply it twice.
    # Initialised with T_j = 1 everywhere; updated as trophic levels are computed
    prey_mass_tl <- sweep(n, 2, params@w * params@dw, "*")  # no_sp x no_w

    # Size-dependent trophic level of the resource:
    #   T_R(w) = max(1, 1 + log(w / w_R) / log(beta_R))
    # where w_R is an average primary-producer size and beta_R an average
    # predator/prey mass ratio for the resource. Unlike the species trophic
    # levels, this is fixed by the formula, so the resource's whole contribution
    # to the numerator can be precomputed before the size loop. This mirrors the
    # `phi_prey_background` term in mizerEncounter().
    tl_R <- pmax(1, 1 + log(params@w_full / w_R) / log(beta_R))  # no_w_full
    # TL-weighted resource biomass per bin (length no_w_full), weighted exactly
    # as the resource is weighted in mizerEncounter().
    resource_mass_tl <- tl_R * n_pp * params@w_full * params@dw_full
    # Resource contribution to the TL-weighted encounter for every predator/size:
    #   ae_R[i, k] = sum_p kernel[i, k, p] * resource_mass_tl[p]
    ae_R <- rowSums(sweep(pred_kernel, 3, resource_mass_tl, "*"), dims = 2)
    E_tl_resource <- params@search_vol *
        (params@species_params$interaction_resource * ae_R)  # no_sp x no_w

    # Cumulative numerator A_i and denominator B_i (integrals weighted by 1/g dw)
    cumA <- numeric(no_sp)
    cumB <- numeric(no_sp)

    # Output matrix: NA below egg size of each species
    tl <- matrix(NA_real_, nrow = no_sp, ncol = no_w,
                 dimnames = dimnames(params@initial_n))

    # Iterate from smallest to largest size, building up trophic levels
    for (k in seq_len(no_w)) {
        # Trophic-level-weighted encounter for all predator species at size w[k]:
        # E_tl[i] = gamma_i(w_k) * sum_j theta_ij * sum_p kernel[i,k,p] * N_tl[j,p] * w_p * dw_p
        pred_kernel_k <- matrix(pred_kernel[, k, idx_sp], nrow = no_sp)
        ae_k <- pred_kernel_k %*% t(prey_mass_tl)  # no_sp x no_sp
        E_tl <- params@search_vol[, k] * rowSums(params@interaction * ae_k)
        # Add the TL-weighted resource contribution (resource trophic level is
        # given by the formula, so this term is precomputed).
        E_tl <- E_tl + E_tl_resource[, k]

        # Update trophic level for each species active at this size
        for (i in seq_len(no_sp)) {
            if (k < params@w_min_idx[i]) next
            g_ik <- growth[i, k]
            if (!is.finite(g_ik) || g_ik <= 0) {
                # No growth: carry forward previous trophic level
                tl[i, k] <- if (k > params@w_min_idx[i]) tl[i, k - 1L] else 1
                next
            }
            weight <- params@dw[k] / g_ik
            cumA[i] <- cumA[i] + (1 - feeding_level[i, k]) * E_tl[i] * weight
            cumB[i] <- cumB[i] + consumption[i, k] * weight
            tl[i, k] <- if (cumB[i] > 0) 1 + cumA[i] / cumB[i] else 1
        }

        # Update prey_mass_tl for size k with newly computed trophic levels
        active_k <- k >= params@w_min_idx
        tl_k <- tl[, k]
        tl_k[is.na(tl_k)] <- 1
        prey_mass_tl[active_k, k] <-
            n[active_k, k] * tl_k[active_k] * params@w[k] * params@dw[k]
    }

    return(ArraySpeciesBySize(tl, value_name = "Trophic level", params = params))
}


#' Get mean trophic level of each species
#'
#' `r lifecycle::badge("experimental")`
#' Calculates the consumption-rate-weighted mean trophic level of each species,
#' defined as
#' \deqn{
#'   T_i = \frac{\int r_i(w)\,N_i(w)\,T_i(w)\,dw}
#'              {\int r_i(w)\,N_i(w)\,dw},
#' }
#' where \eqn{r_i(w) = (1 - f_i(w))\,E_i(w)} is the consumption rate of an
#' individual of species \eqn{i} at weight \eqn{w}, \eqn{N_i(w)} is the
#' abundance density, and \eqn{T_i(w)} is the size-resolved trophic level
#' from [getTrophicLevel()]. As in [getTrophicLevel()], the resource is given a
#' size-dependent trophic level controlled by the `w_R` and `beta_R` arguments.
#'
#' @inheritParams getTrophicLevel
#'
#' @return A named vector with the mean trophic level for each species.
#'
#' @export
#' @family summary functions
#' @concept summary_function
#' @seealso [getTrophicLevel()]
#' @examples
#' getTrophicLevelBySpecies(NS_params)
getTrophicLevelBySpecies <- function(params,
                                     n = initialN(params),
                                     n_pp = initialNResource(params),
                                     n_other = initialNOther(params),
                                     w_R = 1e-10,
                                     beta_R = 1000,
                                     ...) {
    UseMethod("getTrophicLevelBySpecies")
}

#' @export
getTrophicLevelBySpecies.MizerParams <- function(params,
                                                  n = initialN(params),
                                                  n_pp = initialNResource(params),
                                                  n_other = initialNOther(params),
                                                  w_R = 1e-10,
                                                  beta_R = 1000,
                                                  ...) {
    params <- validParams(params)
    tl <- getTrophicLevel(params, n = n, n_pp = n_pp, n_other = n_other,
                          w_R = w_R, beta_R = beta_R)
    tl[is.na(tl)] <- 0
    encounter <- getEncounter(params, n, n_pp, n_other)
    feeding_level <- getFeedingLevel(params, n, n_pp, n_other)
    # Consumption rate per individual times abundance density
    consumption_n <- (1 - feeding_level) * encounter * n  # no_sp x no_w
    # Weighted mean trophic level: integral of consumption_n * T * dw / integral of consumption_n * dw
    numerator <- (consumption_n * tl) %*% params@dw
    denominator <- consumption_n %*% params@dw
    tl_by_sp <- numerator[, 1] / denominator[, 1]
    tl_by_sp[denominator[, 1] == 0] <- NA_real_
    return(tl_by_sp)
}

#' Calculate the SSB of species
#'
#' Calculates the spawning stock biomass (SSB) for each species. For a
#' `MizerSim` object this is returned for every saved time; for a
#' `MizerParams` object it is calculated from the initial state. SSB is the
#' total mass of all mature individuals.
#'
#' @param object An object of class `MizerParams` or `MizerSim`.
#'
#' @return If called with a MizerParams object, a named vector with the SSB in
#'   grams for each species in the model. If called with a MizerSim object, a
#'   `ArrayTimeBySpecies` object (time x species) containing the SSB in grams
#'   at each time step for all species.
#' @export
#' @family summary functions
#' @concept summary_function
#' @examples
#' ssb <- getSSB(NS_sim)
#' ssb[c("1972", "2010"), c("Herring", "Cod")]
getSSB <- function(object) {
    UseMethod("getSSB")
}
#' @export
getSSB.MizerSim <- function(object) {
    params <- object@params
    # The weighting is the product maturity * w, which sizeIntegral() bin-averages
    # as a single weighting factor.
    sizeIntegral(object, weighting = sweep(params@maturity, 2, params@w, "*"),
                 value_name = "Spawning stock biomass", units = "g")
}
#' @export
getSSB.MizerParams <- function(object) {
    sizeIntegral(object,
                 weighting = sweep(object@maturity, 2, object@w, "*"))
}


#' Calculate the total biomass of each species within a size range at each time
#' step.
#'
#' Calculates the total biomass through time within user defined size limits.
#' The default option is to use the size range starting at the size specified
#' by the `biomass_cutoff` species parameter, if it is set, or else the full
#' size range of each species. You can specify minimum
#' and maximum weight or length range for the species. Lengths take precedence
#' over weights (i.e. if both min_l and min_w are supplied, only min_l will be
#' used).
#'
#' @details
#' When no size range arguments are provided, the function checks if the
#' `biomass_cutoff` column exists in the species parameters. If it does,
#' those values are used as the minimum weight for each species. For species
#' with NA values in `biomass_cutoff`, the default minimum weight (smallest
#' weight in the model) is used.
#'
#' @param object An object of class `MizerParams` or `MizerSim`.
#' @param use_cutoff If TRUE, the `biomass_cutoff` column in the
#'   species parameters is used as the minimum weight for each species (ignoring any
#'   size range arguments in `...`). If FALSE (default), the specified size range
#'   arguments are used, if provided, or the full size range of the species is used.
#' @inheritDotParams get_size_range_array -params
#'
#' @return If called with a MizerParams object, a named vector with the biomass
#'   in grams for each species in the model. If called with a MizerSim object,
#'   an `ArrayTimeBySpecies` object (time x species) containing the biomass in
#'   grams at each time step for all species.
#' @export
#' @family summary functions
#' @concept summary_function
#' @examples
#' biomass <- getBiomass(NS_sim)
#' biomass["1972", "Herring"]
#' biomass <- getBiomass(NS_sim, min_w = 10, max_w = 1000)
#' biomass["1972", "Herring"]
#'
#' # If species_params contains a `biomass_cutoff`` column, it can be used
#' # as the minimum weight when use_cutoff = TRUE
#' species_params(NS_sim@params)$biomass_cutoff <- 10
#' biomass <- getBiomass(NS_sim, use_cutoff = TRUE)  # Uses biomass_cutoff as min_w
#' biomass["1972", "Herring"]
getBiomass <- function(object, use_cutoff = FALSE, ...) {
    UseMethod("getBiomass")
}
#' @export
getBiomass.MizerSim <- function(object, use_cutoff = FALSE, ...) {
    params <- object@params
    if (use_cutoff && "biomass_cutoff" %in% names(params@species_params)) {
        sizeIntegral(object, weighting = params@w,
                     min_w = cutoff_min_w(params),
                     value_name = "Biomass", units = "g")
    } else {
        sizeIntegral(object, weighting = params@w, ...,
                     value_name = "Biomass", units = "g")
    }
}
#' @export
getBiomass.MizerParams <- function(object, use_cutoff = FALSE, ...) {
    params <- object
    if (use_cutoff && "biomass_cutoff" %in% names(params@species_params)) {
        sizeIntegral(params, weighting = params@w,
                     min_w = cutoff_min_w(params))
    } else {
        sizeIntegral(params, weighting = params@w, ...)
    }
}


#' Calculate the number of individuals within a size range
#'
#' Calculates the number of individuals within user-defined size limits. The
#' default option is to use the whole size range. You can specify minimum and
#' maximum weight or lengths for the species. Lengths take precedence over
#' weights (i.e. if both min_l and min_w are supplied, only min_l will be used)
#'
#' @param object An object of class `MizerParams` or `MizerSim`.
#' @inheritDotParams get_size_range_array -params
#'
#' @return If called with a MizerParams object, a named vector with the numbers
#'   for each species in the model. If called with a MizerSim object, a
#'   `ArrayTimeBySpecies` object (time x species) containing the numbers at
#'   each time step for all species.
#' @export
#' @family summary functions
#' @concept summary_function
#' @examples
#' numbers <- getN(NS_sim)
#' numbers["1972", "Herring"]
#' # The above gave a huge number, because that included all the larvae.
#' # The number of Herrings between 10g and 1kg is much smaller.
#' numbers <- getN(NS_sim, min_w = 10, max_w = 1000)
#' numbers["1972", "Herring"]
getN <- function(object, ...) {
    UseMethod("getN")
}
#' @export
getN.MizerSim <- function(object, ...) {
    sizeIntegral(object, ..., value_name = "Abundance")
}
#' @export
getN.MizerParams <- function(object, ...) {
    sizeIntegral(object, ...)
}


#' Calculate the rate at which biomass of each species is fished by each gear
#'
#' This yield rate is given in grams per year. It is calculated at each time
#' step saved in the MizerSim object.
#'
#' For details of how the yield rate is defined see the help page of
#' [getYield()].
#'
#' @param object An object of class `MizerParams` or `MizerSim`.
#'
#' @return If called with a MizerParams object, an array (gear x species) with
#'   the yield rate in grams per year from each gear for each species in the
#'   model. If called with a MizerSim object, an array (time x gear x species)
#'   containing the yield rate at each time step.
#' @export
#' @family summary functions
#' @concept summary_function
#' @seealso [getYield()]
#' @examples
#' yield <- getYieldGear(NS_sim)
#' dim(yield)
#' yield["1972", , "Herring"]
getYieldGear <- function(object) {
    UseMethod("getYieldGear")
}
#' @export
getYieldGear.MizerSim <- function(object) {
    # The weighting is the product F * w (time x gear x sp x w), which
    # sizeIntegral() bin-averages as a single weighting factor.
    f_gear <- getFMortGear(object)
    sizeIntegral(object,
                 weighting = sweep(f_gear, length(dim(f_gear)),
                                   object@params@w, "*"))
}
#' @export
getYieldGear.MizerParams <- function(object) {
    f_gear <- getFMortGear(object)  # gear x sp x w
    sizeIntegral(object,
                 weighting = sweep(f_gear, length(dim(f_gear)), object@w, "*"))
}

#' Calculate the rate at which biomass of each species is fished
#'
#' This yield rate is given in grams per year. It is calculated at each time
#' step saved in the MizerSim object.
#'
#' @details
#' The yield rate \eqn{y_i(t)} for species \eqn{i} at time \eqn{t} is defined as
#' \deqn{y_i(t)=\int\mu_{f.i}(w, t)N_i(w, t)w dw}
#' where \eqn{\mu_{f.i}(w, t)} is the fishing mortality of an individual of
#' species \eqn{i} and weight \eqn{w} at time \eqn{t} and \eqn{N_i(w, t)} is the
#' abundance density of such individuals.  The factor of \eqn{w} converts the
#' abundance density into a biomass density and the integral aggregates the
#' contribution from all sizes.
#'
#' The total catch in a time period from \eqn{t_1} to  \eqn{t_2} is the integral
#' of the yield rate over that period:
#' \deqn{C = \int_{t_1}^{t2}y_i(t)dt}
#' In practice, as the yield rate is only available
#' at the saved times, one can only approximate this integral by averaging over
#' the available yield rates during the time period and multiplying by the time
#' period. The less the yield changes between the saved values, the more
#' accurate this approximation is. So the approximation can be improved by
#' saving simulation results at smaller intervals, using the `t_save` argument
#' to [project()]. But this is only a concern if abundances change quickly
#' during the time period of interest.
#'
#' @param object An object of class `MizerParams` or `MizerSim`.
#'
#' @return If called with a `MizerParams` object, a named numeric vector with
#'   the yield rate in grams per year for each species in the model. If called
#'   with a `MizerSim` object, an `ArrayTimeBySpecies` object (time x species)
#'   containing the yield rate in grams per year at each saved time step.
#' @export
#' @family summary functions
#' @concept summary_function
#' @seealso [getYieldGear()]
#' @examples
#' yield <- getYield(NS_sim)
#' yield[c("1972", "2010"), c("Herring", "Cod")]
#'
#' # Running simulation for another year, saving intermediate time steps
#' params <- finalParams(NS_sim)
#' sim <- project(params, t_save = 0.1, t_max = 1,
#'                t_start = 2010, progress_bar = FALSE)
#' # The yield rate for Herring decreases during the year
#' getYield(sim)[, "Herring"]
#' # We approximate the total catch in the year by averaging over the year
#' sum(getYield(sim)[1:10, "Herring"] / 10)
getYield <- function(object) {
    UseMethod("getYield")
}
#' @export
getYield.MizerSim <- function(object) {
    # The weighting is the product F * w (time x sp x w), which sizeIntegral()
    # bin-averages as a single weighting factor.
    f <- getFMort(object, drop = FALSE)
    sizeIntegral(object,
                 weighting = sweep(f, length(dim(f)), object@params@w, "*"),
                 value_name = "Yield rate", units = "g/year")
}
#' @export
getYield.MizerParams <- function(object) {
    f <- getFMort(object, drop = FALSE)  # sp x w
    sizeIntegral(object,
                 weighting = sweep(f, length(dim(f)), object@w, "*"))
}


#' Get growth curves giving weight as a function of age
#'
#' @param object MizerSim or MizerParams object. If given a
#'   \linkS4class{MizerSim} object, uses the growth rates at the final time of a
#'   simulation to calculate the size at age. If given a
#'   \linkS4class{MizerParams} object, uses the initial growth rates instead.
#' @inheritParams valid_species_arg
#' @param max_age The age up to which to run the growth curve. Default is 20.
#' @param percentage Boolean value. If TRUE, the size is given as a percentage
#'   of the maximal size.
#'
#' @return An array (species x age) containing the weight in grams.
#' @export
#' @family summary functions
#' @examples
#' growth_curves <- getGrowthCurves(NS_params, species = c("Cod", "Haddock"))
#' str(growth_curves)
#'
#' library(ggplot2)
#' ggplot(melt(growth_curves)) +
#'   geom_line(aes(Age, value)) +
#'   facet_wrap(~ Species, scales = "free") +
#'   ylab("Size[g]") + xlab("Age[years]")
getGrowthCurves <- function(object,
                            species = NULL,
                            max_age = 20,
                            percentage = FALSE) {
    UseMethod("getGrowthCurves")
}
#' @export
getGrowthCurves.MizerSim <- function(object,
                            species = NULL,
                            max_age = 20,
                            percentage = FALSE) {
    params <- finalParams(object)
    getGrowthCurves(params, species, max_age, percentage)
}
#' @export
getGrowthCurves.MizerParams <- function(object,
                            species = NULL,
                            max_age = 20,
                            percentage = FALSE) {
    params <- object
    species <- valid_species_arg(params, species)
    # reorder list of species to coincide with order in params
    idx <- which(params@species_params$species %in% species)
    species <- params@species_params$species[idx]
    age <- seq(0, max_age, length.out = 50)
    ws <- array(dim = c(length(species), length(age)),
                dimnames = list(Species = species, Age = age))
    g <- getEGrowth(params)
    for (j in seq_along(species)) {
        i <- idx[j]
        w_max <- params@species_params$w_max[[i]]
        # Build the interpolation points, appending a point at `w_max` where
        # growth is zero. Drop any grid points at or above `w_max` first to
        # avoid duplicate x values, which would trigger a warning from
        # approxfun().
        keep <- params@w < w_max
        g_fn <- stats::approxfun(c(params@w[keep], w_max),
                                 c(g[i, keep], 0))
        myodefun <- function(t, state, parameters) {
            return(list(g_fn(state)))
        }
        ws[j, ] <- deSolve::ode(y = params@w[params@w_min_idx[i]],
                                times = age, func = myodefun)[, 2]
        if (percentage) {
            ws[j, ] <- ws[j, ] / params@species_params$w_max[i] * 100
        }
    }
    return(ws)
}

#' Get size range array
#'
#' Helper function that returns an array (species x size) of logical values
#' indicating whether that size bin is within the size limits specified by the
#' arguments. Either the size limits can be the same for all species or they
#' can be specified as vectors with one value for each species in the model.
#'
#' @param params MizerParams object
#' @param min_w Smallest weight in size range. Defaults to smallest weight in
#'   the model.
#' @param max_w Largest weight in size range. Defaults to largest weight in the
#'   model.
#' @param min_l Smallest length in size range. If supplied, this takes
#'   precedence over `min_w`.
#' @param max_l Largest length in size range. If supplied, this takes precedence
#'   over `max_w`.
#' @param ... Unused
#'
#' @return A logical array (species x size), with dimnames `sp` and `w`.
#'
#' @section Length to weight conversion:
#' If `min_l` is specified there is no need to specify `min_w` and so on.
#' However, if a length is specified (minimum or maximum) then it is necessary
#' for the species parameter data.frame to include the parameters `a` and `b`
#' that determine the relation between length \eqn{l} and weight \eqn{w} by
#' \deqn{w = a l^b.}
#'
#' It is possible to mix length and weight constraints, e.g. by supplying a
#' minimum weight and a maximum length, but this must be done the same for
#' all species. The default values are the minimum and
#' maximum weights of the spectrum, i.e., the full range of the size spectrum is
#' used.
#' @export
#' @concept helper
get_size_range_array <- function(params, min_w = min(params@w),
                                 max_w = max(params@w),
                                 min_l = NULL, max_l = NULL, ...) {
    no_sp <- nrow(params@species_params)
    if (!is.null(min_l) || !is.null(max_l)) {
        if (any(!c("a", "b") %in% names(params@species_params))) {
            stop("species_params slot must have columns 'a' and 'b' for ",
                 "length-weight conversion")
        }
        if (anyNA(params@species_params[["a"]]) ||
            anyNA(params@species_params[["b"]])) {
            stop("There must be no NAs in the species_params columns 'a' and 'b'.")
        }
    }
    if (!is.null(min_l)) {
        if (length(min_l) != 1 && length(min_l) != no_sp) {
            stop("min_l must be a single number or a vector with one value for each species.")
        }
        min_w <- params@species_params[["a"]] *
            min_l ^ params@species_params[["b"]]
    }
    if (!is.null(max_l)) {
        if (length(max_l) != 1 && length(max_l) != no_sp) {
            stop("max_l must be a single number or a vector with one value for each species.")
        }
        max_w <- params@species_params[["a"]] *
            max_l ^ params@species_params[["b"]]
    }
    if (length(min_w) == 1) {
        min_w <- rep(min_w, no_sp)
    }
    if (length(max_w) == 1) {
        max_w <- rep(max_w, no_sp)
    }
    if (length(min_w) != no_sp || length(max_w) != no_sp) {
        stop("min_w and max_w must be a single number of a vector with one
             value for each species.")
    }
    min_n <- t(sapply(min_w, function(x) params@w >= x))
    max_n <- t(sapply(max_w, function(x) params@w <= x))
    size_n <- min_n & max_n
    dimnames(size_n) <- list(sp = params@species_params$species,
                             w = signif(params@w, 3))
    return(size_n)
}


#### summary for MizerParams ####
#' @export
summary.MizerParams <- function(object, ...) {
    params <- validParams(object)
    cat("An object of class \"", as.character(class(params)), "\" \n", sep = "")
    # Display metadata if available
    md <- getMetadata(params)
    if (!is.null(md$title) && nchar(md$title) > 0)
        cat("Title: ", md$title, "\n", sep = "")
    if (!is.null(md$description) && nchar(md$description) > 0)
        cat("Description: ", md$description, "\n", sep = "")
    if (!is.null(md$authors)) {
        authors <- md$authors
        if (is.list(authors)) {
            author_names <- sapply(authors, function(a) {
                if (is.list(a)) a$name else as.character(a)
            })
            cat("Authors: ", paste(author_names, collapse = ", "), "\n", sep = "")
        } else {
            cat("Authors: ", paste(authors, collapse = ", "), "\n", sep = "")
        }
    }
    if (!is.null(md$doi) && length(md$doi) > 0)
        cat("DOI: ", paste(md$doi, collapse = ", "), "\n", sep = "")
    if (!is.null(md$url) && length(md$url) > 0)
        cat("URL: ", paste(md$url, collapse = ", "), "\n", sep = "")
    cat("mizer version: ", as.character(md$mizer_version), "\n", sep = "")
    if (!is.null(md$time_created))
        cat("Created: ", format(md$time_created), "\n", sep = "")
    if (!is.null(md$time_modified))
        cat("Modified: ", format(md$time_modified), "\n", sep = "")
    cat("Consumer size spectrum:\n")
    cat("\tminimum size:\t", signif(min(params@w)), "\n", sep = "")
    cat("\tmaximum size:\t", signif(max(params@w)), "\n", sep = "")
    cat("\tno. size bins:\t", length(params@w), "\n", sep = "")
    # Length of background?
    cat("Resource size spectrum:\n")
    cat("\tminimum size:\t", signif(min(params@w_full)), "\n", sep = "")
    cat("\tmaximum size:\t", signif(max(params@w_full[params@initial_n_pp > 0])),
        "\n", sep = "")
    cat("\tno. size bins:\t", length(params@w_full[params@initial_n_pp > 0]),
        "\t(", length(params@w_full), " size bins in total)\n", sep = "")
    # Whether the model is at its steady state is not visible from any of the
    # parameters above, and getting it wrong is the most common way a
    # calibration goes quietly wrong, so it is reported here.
    report <- tryCatch(steady_drift_report(params), error = function(e) NULL)
    residual <- if (is.null(report)) NA_real_ else report$drift
    header_shown <- FALSE
    if (is.finite(residual)) {
        cat("Steady state:\n")
        header_shown <- TRUE
        cat("\tbiomass drift:\t", signif(residual, 2), " /year\t",
            if (residual <= steady_residual_tol()) {
                "(at steady state)"
            } else {
                paste0("(not at steady state, largest in ",
                       steady_variable_txt(report, quote = FALSE),
                       " - run tuneSteadyState())")
            }, "\n", sep = "")
    }
    # The components are not in the drift above, so a model with a drifting
    # component reads as steady on that line alone. Listing them here is what
    # stops this block being misleading for exactly the models that need it.
    if (!is.null(report) && length(report$other) > 0) {
        if (!header_shown) cat("Steady state:\n")
        for (component in names(report$other)) {
            drift <- report$other[[component]]
            cat("\t", component, " drift:\t",
                if (is.finite(drift)) {
                    paste0(signif(drift, 2), " /year")
                } else {
                    "not measurable"
                },
                "\t(not included in the biomass drift)\n", sep = "")
        }
    }
    cat("Species details:\n")
    sel_params <- intersect(c("species", "w_inf", "w_mat", "w_min", "f0", "fc",
                              "age_mat", "beta", "sigma"),
                            names(params@species_params))
    sp <- params@species_params[, sel_params]
    rownames(sp) <- NULL
    print(sp)
    cat("\nFishing gear details:\n")
    cat(sprintf("%-13s %s %s", "Gear", "Effort", " Target species"), "\n",
        "----------------------------------\n")
    gears <- dimnames(params@catchability)$gear
    for (i in seq_len(dim(params@catchability)[1])) {
        splist <- dimnames(params@catchability)$sp[params@catchability[i, ] > 0]
        cat(sprintf("%-14s %1.2f   %s",
                gears[i], params@initial_effort[[gears[i]]],
                toString(splist)), "\n")
    }
    invisible(params)
}


#### summary for MizerSim ####
#' @export
summary.MizerSim <- function(object, ...) {
    cat("An object of class \"", as.character(class(object)), "\" \n", sep = "")
    cat("Parameters:\n")
    # Report the effort that was actually used during the simulation (stored in
    # the `effort` slot) rather than the model's `initial_effort`. We delegate
    # to the params summary on a copy whose `initial_effort` holds the effort
    # used, averaged over the saved timesteps for each gear.
    params <- object@params
    used <- object@effort
    gears <- dimnames(used)$gear
    params@initial_effort[gears] <- colMeans(used)
    summary(params)
    # Flag any gears whose effort was not constant over the saved timesteps, so
    # that the averaged value shown above is not mistaken for a fixed effort.
    varied <- apply(used, 2, function(x) !isTRUE(all.equal(min(x), max(x))))
    if (any(varied)) {
        ranges <- vapply(gears[varied], function(g) {
            sprintf("%s (%1.2f to %1.2f)", g, min(used[, g]), max(used[, g]))
        }, character(1))
        cat("\tNote: effort varied over time for ", toString(ranges),
            "; mean shown above.\n", sep = "")
    }
    cat("Simulation parameters:\n")
    # Need to store t_max and dt in a description slot? Or just in simulation
    # time parameters? Like a list?
    cat("\tTime period: ",
        min(as.numeric(dimnames(object@n)$time)), " to ",
        max(as.numeric(dimnames(object@n)$time)),
        "\n", sep = "")
    cat("\tOutput stored every ",
        as.numeric(dimnames(object@n)$time)[2] -
            as.numeric(dimnames(object@n)$time)[1], " years\n", sep = "")
    cat("\tTime step ", object@sim_params$dt, sep = "")
    cat("\tMethod: ", object@sim_params$method, sep = "")
    invisible(object)
}


# The default str() on an S4 object recurses through every slot. Many
# MizerParams slots hold arrays that carry a `params` attribute pointing back to
# the whole object, so the default output would recursively dump the entire
# model, flooding the console. We instead print a one-line summary per slot.
#' @export
str.MizerParams <- function(object, max.level = NA, ...) {
    cat("Formal class 'MizerParams' [package \"mizer\"] with ",
        length(slotNames(object)), " slots\n", sep = "")
    if (!is.na(max.level) && max.level <= 0) {
        return(invisible(NULL))
    }
    slot_max_level <- if (is.na(max.level)) NA else max.level - 1
    for (name in slotNames(object)) {
        val <- slot(object, name)
        val_str <- utils::capture.output(utils::str(val, max.level = slot_max_level, ...))
        if (length(val_str) > 0) {
            cat(" @ ", name, paste(rep(" ", max(0, 20 - nchar(name))), collapse = ""), ": ", val_str[1], "\n", sep = "")
            if (length(val_str) > 1) {
                cat(paste0("   ", val_str[-1], collapse = "\n"), "\n")
            }
        }
    }
    invisible(NULL)
}


# Like str.MizerParams: avoid the default S4 recursion, which would expand the
# embedded MizerParams (and its array slots' `params` back-references) in full.
#' @export
str.MizerSim <- function(object, max.level = NA, ...) {
    cat("Formal class 'MizerSim' [package \"mizer\"] with ",
        length(slotNames(object)), " slots\n", sep = "")
    if (!is.na(max.level) && max.level <= 0) {
        return(invisible(NULL))
    }
    slot_max_level <- if (is.na(max.level)) NA else max.level - 1
    for (name in slotNames(object)) {
        val <- slot(object, name)
        val_str <- utils::capture.output(utils::str(val, max.level = slot_max_level, ...))
        if (length(val_str) > 0) {
            cat(" @ ", name, paste(rep(" ", max(0, 20 - nchar(name))), collapse = ""), ": ", val_str[1], "\n", sep = "")
            if (length(val_str) > 1) {
                cat(paste0("   ", val_str[-1], collapse = "\n"), "\n")
            }
        }
    }
    invisible(NULL)
}
