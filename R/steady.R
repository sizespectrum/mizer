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
#'   after each of which we check for convergence. Default value is 1.5. This
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
steady <- function(params, t_max = 100, t_per = 1.5, tol = 10^(-2),
                   dt = 0.1, return_sim = FALSE, progress_bar = TRUE) {
    params <- validParams(params)
    assert_that(noNA(getRDD(params)))
    if ((t_per < dt) || !isTRUE(all.equal((t_per - round(t_per / dt) * dt), 0))) {
        stop("t_per must be a positive multiple of dt")
    }
    t_dimnames <-  seq(0, t_max, by = t_per)
    
    if (is(progress_bar, "Progress")) {
        # We have been passed a shiny progress object
        progress_bar$set(message = "Finding steady state", value = 0)
        proginc <- 1/ceiling(t_max/t_per)
    }
    
    if (return_sim) {
        # create MizerSim object
        sim <- MizerSim(params, t_dimnames =  t_dimnames)
    }
    
    # Force the reproduction to stay at the current level
    params@species_params$constant_reproduction <- getRDD(params)
    old_rdd_fun <- params@rates_funcs$RDD
    params@rates_funcs$RDD <- "constantRDD"
    old_rdi <- getRDI(params)
    if (any(is.na(old_rdi)) || any(old_rdi <= 0)) {
        stop("The project function expects positive RDI for all species.")
    }
    rdi_limit <- old_rdi / 1e7
    # Force other componens to stay at current level
    old_other_dynamics <- params@other_dynamics
    for (res in names(params@other_dynamics)) {
        params@other_dynamics[[res]] <- "constant_other"
    }
    
    # get functions
    resource_dynamics_fn <- get(params@resource_dynamics)
    other_dynamics_fns <- lapply(params@other_dynamics, get)
    rates_fns <- lapply(params@rates_funcs, get)
    
    n_list <- list(n = params@initial_n,
                   n_pp = params@initial_n_pp,
                   n_other = params@initial_n_other)
    
    for (i in seq_along(t_dimnames)) {
        # advance shiny progress bar
        if (is(progress_bar, "Progress")) {
            progress_bar$inc(amount = proginc)
        }
        n_list <- project_simple(params, n = n_list$n, n_pp = n_list$n_pp,
                                 n_other = n_list$n_other, t = 0,
                                 dt = dt, steps = round(t_per / dt),
                                 effort = params@initial_effort,
                                 resource_dynamics_fn = resource_dynamics_fn,
                                 other_dynamics_fns = other_dynamics_fns,
                                 rates_fns = rates_fns)
        if (return_sim) {
            # Store result
            sim@n[i, , ] <- n_list$n
            sim@n_pp[i, ] <- n_list$n_pp
            sim@n_other[i, ] <- n_list$n_other
            sim@effort[i, ] <- params@initial_effort
        }
        
        new_rdi <- n_list$rates$rdi
        deviation <- max(abs((new_rdi - old_rdi)/old_rdi))
        if (any(new_rdi < rdi_limit)) {
            if (return_sim) {
                message("One of the species is going extinct.")
                break
            }
            extinct <- params@species_params$species[new_rdi < rdi_limit]
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
                i * t_per,
                " years. Residual relative rate of change = ", deviation)
    } else {
        message("Steady state was reached before ", i * t_per, " years.")
    }
    
    # Restore original RDD and other dynamics
    params@rates_funcs$RDD <- old_rdd_fun
    params@other_dynamics <- old_other_dynamics
    
    params@initial_n[] <- n_list$n
    params@initial_n_pp[] <- n_list$n_pp
    params@initial_n_other[] <- n_list$n_other
    
    # Retune the values of erepro so that we get the correct level of
    # reproduction
    params <- retune_erepro(params)
    
    if (return_sim) {
        sim@params <- params
        return(sim)
    } else {
        return(params)
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
#' @inheritParams valid_species_arg
#' @return A MizerParams object with updated values for the `erepro` column
#'   in the `species_params` data frame.
#' @export
retune_erepro <- function(params, species = species_params(params)$species) {
    assert_that(is(params, "MizerParams"))
    species <- valid_species_arg(params, species)
    if (length(species) == 0) return(params)

    mumu <- getMort(params)
    gg <- getEGrowth(params)
    rdi <- getRDI(params)
    eff <- params@species_params$erepro
    for (i in seq_len(nrow(params@species_params))[species]) {
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


#' Helper function to keep other components constant
#' 
#' @param params MizerParams object
#' @param n_other Abundances of other components
#' @param component Name of the component that is being updated
#' @param ... Unused
#' @export
constant_other <- function(params, n_other, component, ...) {
    n_other[[component]]
}

#' Helper function to assure validity of species argument
#' 
#' @param params A MizerParams object
#' @param species A vector of the names of the species to be affected or a
#'   boolean vector indicating for each species whether it is to be affected
#'   (TRUE) or not. By default all species are affected.
#' @export
valid_species_arg <- function(params, species) {
    no_sp <- nrow(params@species_params)
    if (is.logical(species)) {
        if (length(species) != no_sp) {
            stop("The boolean `species` argument has the wrong length")
        }
    } else {
        invalid <- setdiff(species, dimnames(params@initial_n)$sp)
        if (length(invalid) > 0) {
            warning("The following species do not exist: ", 
                    toString(invalid))
        }
        species <- dimnames(params@initial_n)$sp %in% species
        if (length(species) == 0) {
            warning("The species argument matches none of the species in the params object")
        }
    }
    species
}