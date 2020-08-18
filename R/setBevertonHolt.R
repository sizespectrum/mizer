#' Set Beverton-Holt density dependence
#' 
#' Takes a MizerParams object (with arbitrary density dependence) and sets a
#' Beverton-Holt density-dependence with a maximum reproduction rate that is a
#' chosen factor `R_factor` higher than the initial-state reproduction rate. At
#' the same time it adjusts the reproductive efficiency `erepro` so that the
#' steady-state abundances do not change. Setting `R_factor = Inf` switches off
#' all density dependence.
#' 
#' @param params A MizerParams object
#' @param R_factor The factor by which the maximum reproduction rate should be
#'   higher than the initial-state reproduction rate
#' 
#' @return A MizerParams object
#' @export
setBevertonHolt <- function(params, R_factor) {
    assert_that(is(params, "MizerParams"),
                is.numeric(R_factor),
                length(R_factor) %in% c(1, nrow(params@species_params)),
                all(R_factor > 1))
    
    rdi <- getRDI(params)
    rdd <- getRDD(params)
    params@species_params$R_max <- R_factor * getRDD(params)
    
    # erepro needs to be changed to
    # compensate for using a Beverton Holt relationship
    # because RDD = (1-1/R_factor) RDI
    params@species_params$erepro <- 
        params@species_params$erepro / (1 - 1 / R_factor) *
        rdd / rdi
    
    return(setReproduction(params, RDD = "BevertonHoltRDD"))
}

#' Alias for setBevertonHolt
#' 
#' An alias provided for backward compatibility with mizer version <= 2.0.4
#' @inherit setBevertonHolt
#' @export
setRmax <- setBevertonHolt
