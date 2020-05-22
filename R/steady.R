#' Set initial values to a steady state for the model
#' 
#' The steady state is found by running the dynamics while keeping reproduction
#' and other components constant until the size spectra no longer change (or
#' until time `t_max` is reached if earlier) Then the reproductive efficiencies
#' are set to the values that give the level of reproduction observed in that
#' steady state.
#' 
#' @param params A \linkS4class{MizerParams} object
#' @param t_max The maximum number of years to run the simulation. Default is 100.
#' @param t_per The simulation is broken up into shorter runs of `t_per` years,
#'   after each of which we check for convergence. Default value is 7.5. This
#'   should be chosen as an odd multiple of the timestep `dt` in order to be
#'   able to detect period 2 cycles.
#' @param tol The simulation stops when the relative change in the egg
#'   production RDI over `t_per` years is less than `tol` for every species.
#'   Default value is 1/100.
#' @param dt The time step to use in `project()`.
#' @param return_sim If TRUE, the function returns the MizerSim object holding
#'   the result of the simulation run. If FALSE (default) the function returns
#'   a MizerParams object with the "initial" slots set to the steady state.
#' @param progress_bar A shiny progress object to implement a progress bar in a
#'   shiny app. Default FALSE.
#' @export
#' @examples
#' \dontrun{
#' params <- newTraitParams()
#' params@species_params$gamma[5] <- 3000
#' params <- setSearchVolume(params)
#' params <- steady(params)
#' plotSpectra(params)
#' }
steady <- function(params, t_max = 100, t_per = 7.5, tol = 10^(-2),
                   dt = 0.1, return_sim = FALSE, progress_bar = TRUE) {
    assert_that(is(params, "MizerParams"),
                noNA(getRDD(params)))
    p <- params
    
    if (is(progress_bar, "Progress")) {
        # We have been passed a shiny progress object
        progress_bar$set(message = "Finding steady state", value = 0)
        proginc <- 1/ceiling(t_max/t_per)
    }
    
    # Force the reproduction to stay at the current level
    p@species_params$constant_reproduction <- getRDD(p)
    p@rates_funcs$RDD <- "constantRDD"
    old_rdi <- getRDI(p)
    rdi_limit <- old_rdi / 1e7
    # Force other componens to stay at current level
    old_other_dynamics <- p@other_dynamics
    for (res in names(p@other_dynamics)) {
        p@other_dynamics[[res]] <- "constant_other"
    }
    
    n <- p@initial_n
    n_pp <- p@initial_n_pp
    n_other <- p@initial_n_other
    sim <- p
    for (ti in (1:ceiling(t_max/t_per))) {
        # advance shiny progress bar
        if (is(progress_bar, "Progress")) {
            progress_bar$inc(amount = proginc)
        }
        if (return_sim) {
            sim <- project(sim, dt = dt, t_max = t_per, t_save = t_per,
                           initial_n = n, initial_n_pp = n_pp, 
                           initial_n_other = n_other)
        } else {
            sim <- project(p, dt = dt, t_max = t_per, t_save = t_per,
                           initial_n = n, initial_n_pp = n_pp, 
                           initial_n_other = n_other)
        }
        no_t <- dim(sim@n)[1]
        n[] <- sim@n[no_t, , ]
        n_pp[] <- sim@n_pp[no_t, ]
        n_other <- sim@n_other[no_t, ]
        new_rdi <- getRDI(p, n, n_pp, n_other)
        deviation <- max(abs((new_rdi - old_rdi)/old_rdi))
        if (any(new_rdi < rdi_limit)) {
            if (return_sim) {
                message("One of the species is going extinct.")
                break
            }
            extinct <- p@species_params$species[new_rdi < rdi_limit]
            stop(paste(extinct, collapse = ", "),
                 " are going extinct.")
        }
        if (deviation < tol) {
            break
        }
        old_rdi <- new_rdi
    }
    if (deviation >= tol) {
        warning("Simulation run in steady() did not converge after ", 
                ti * t_per,
                " years. Residual relative rate of change = ", deviation)
    } else {
        message("Steady state was reached before ", ti * t_per, " years.")
    }
    
    # Restore original RDD and other dynamics
    p@rates_funcs$RDD <- params@rates_funcs$RDD
    p@other_dynamics <- old_other_dynamics
    
    no_sp <- length(p@species_params$species)
    p@initial_n[] <- n
    p@initial_n_pp[] <- n_pp
    p@initial_n_other[] <- n_other
    
    # Retune the values of erepro so that we get the correct level of
    # reproduction
    p <- retune_erepro(p)
    
    if (return_sim) {
        sim@params <- p
        return(sim)
    } else {
        return(p)
    }
}

#' Retune reproduction efficiency to maintain initial egg abundances
#'
#' Sets the reproductive efficiency for all species so that the rate of egg
#' production exactly compensates for the loss from the first size class due
#' to growth and mortality. Turns off the external density dependence in the
#' reproduction rate by setting the `RDD` function to
#' [noRDD()]
#'
#' @inheritParams steady
#' @param species A vector of the names of the species to be affected or a
#'   boolean vector indicating for each species whether it is to be affected
#'   (TRUE) or not. By default all species are affected
#' @return A MizerParams object
#' @export
retune_erepro <- function(params, species = species_params(params)$species) {
    assert_that(is(params, "MizerParams"))
    
    no_sp <- nrow(params@species_params)
    if (is.logical(species)) {
        if (length(species) != no_sp) {
            stop("The boolean species argument has the wrong length")
        }
    } else {
        species <- dimnames(params@initial_n)$sp %in% species
        if (length(species) == 0) {
            warning("The species argument matches none of the species in the params object")
            return(params)
        }
    }
    mumu <- getMort(params)
    gg <- getEGrowth(params)
    rdi <- getRDI(params)
    eff <- params@species_params$erepro
    for (i in (1:no_sp)[species]) {
        gg0 <- gg[i, params@w_min_idx[i]]
        mumu0 <- mumu[i, params@w_min_idx[i]]
        DW <- params@dw[params@w_min_idx[i]]
        if (!rdi[i] == 0) {
            eff[i] <- params@species_params$erepro[i] *
                (params@initial_n[i, params@w_min_idx[i]] *
                     (gg0 + DW * mumu0)) / rdi[i]
        }
        else {
            eff[i] <- 0.1
        }
    }
    params@species_params$erepro <- eff
    return(setReproduction(params, RDD = "noRDD"))
}



# Helper function to keep other components constant
constant_other <- function(params, n_other, component, ...) {
    n_other[[component]]
}
