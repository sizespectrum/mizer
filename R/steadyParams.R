#' Set initial values to a steady state for the model
#'
#' The steady state is found by running the dynamics while keeping reproduction,
#' resource and other components constant until the size spectra no longer
#' change much (or until time `t_max` is reached, if earlier).
#' 
#' If the model use Beverton-Holt reproduction then the reproduction parameters
#' are set to values that give the level of reproduction observed in that
#' steady state. The `preserve` argument can be used to specify which of the 
#' reproduction parameters should be preserved.
#' 
#' A precursor of this function, called [steady()] differs from this by not
#' keeping the resource constant while running the dynamics to steady state.
#' 
#' @param params A \linkS4class{MizerParams} object
#' @param t_max The maximum number of years to run the simulation. Default is 100.
#' @param t_per The simulation is broken up into shorter runs of `t_per` years,
#'   after each of which we check for convergence. Default value is 1.5. This
#'   should be chosen as an odd multiple of the timestep `dt` in order to be
#'   able to detect period 2 cycles.
#' @param dt The time step to use in `project()`.
#' @param tol The simulation stops when the relative change in the egg
#'   production RDI over `t_per` years is less than `tol` for every species.
#' @param return_sim If TRUE, the function returns the MizerSim object holding
#'   the result of the simulation run. If FALSE (default) the function returns
#'   a MizerParams object with the "initial" slots set to the steady state.
#' @param preserve `r lifecycle::badge("experimental")`
#'   Specifies whether the `reproduction_level` should be preserved (default)
#'   or the maximum reproduction rate `R_max` or the reproductive
#'   efficiency `erepro`. See [setBevertonHolt()] for an explanation
#'   of the `reproduction_level`.
#' @param progress_bar A shiny progress object to implement a progress bar in a
#'   shiny app. Default FALSE.
#' @export
#' @examples
#' \dontrun{
#' params <- newTraitParams()
#' species_params(params)$gamma[5] <- 3000
#' params <- steadyParams(params)
#' plotSpectra(params)
#' }
steadyParams <- 
    function(params, t_max = 100, t_per = 1.5, dt = 0.1,
             tol = 0.1 * dt, return_sim = FALSE, 
             preserve = c("reproduction_level", "erepro", "R_max"),
             progress_bar = TRUE) {
    params <- validParams(params)
    
    if (params@rates_funcs$RDD == "BevertonHoltRDD") {
        preserve <- match.arg(preserve)
        old_reproduction_level <- getReproductionLevel(params)
        old_R_max <- params@species_params$R_max
        old_erepro <- params@species_params$erepro
    }
    
    # Force the reproduction to stay at the current level
    params@species_params$constant_reproduction <- getRDD(params)
    old_rdd_fun <- params@rates_funcs$RDD
    params@rates_funcs$RDD <- "constantRDD"
    
    # Force resource to stay at current level
    old_resource_dynamics <- resource_dynamics(params)
    resource_dynamics(params) <- "resource_constant"
    
    # Force other components to stay at current level
    old_other_dynamics <- params@other_dynamics
    for (res in names(params@other_dynamics)) {
        params@other_dynamics[[res]] <- "constant_other"
    }
    
    object <- projectToSteady(params,
                              distance_func = distanceMaxRelRDI,
                              t_per = t_per,
                              t_max = t_max,
                              dt = dt,
                              tol = tol,
                              return_sim = return_sim,
                              progress_bar = progress_bar)
    if (return_sim) {
        params <- object@params
    } else {
        params <- object
    }
    # Restore original RDD and dynamics
    params@rates_funcs$RDD <- old_rdd_fun
    params@other_dynamics <- old_other_dynamics
    params@species_params$constant_reproduction <- NULL
    resource_dynamics(params) <- old_resource_dynamics
    
    if (params@rates_funcs$RDD == "BevertonHoltRDD") {
        if (preserve == "reproduction_level") {
            params <- setBevertonHolt(params, 
                                      reproduction_level = old_reproduction_level)
        } else if (preserve == "R_max") {
            params <- setBevertonHolt(params, 
                                      R_max = old_R_max)
        } else {
            params <- setBevertonHolt(params, erepro = old_erepro)
        }
    }
    
    # Set resource carrying capacity
    rr <- resource_rate(params)
    cc <- (getResourceMort(params) + rr) / rr * initialNResource(params)
    cc[rr == 0] <- 0
    resource_capacity(params) <- cc
    
    if (return_sim) {
        object@params <- params
        return(object)
    } else {
        params@time_modified <- lubridate::now()
        return(params)
    }
}
