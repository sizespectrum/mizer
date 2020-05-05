#' Set up parameters for a general multispecies model
#'
#' Sets up a multi-species size spectrum model by filling all slots in the
#' \linkS4class{MizerParams} object based on user-provided or default
#' parameters. It does this by creating an empty MizerParams object with
#' [emptyParams()] and then filling the slots by passing its arguments
#' to [setParams()]. There is a long list of arguments, but almost
#' all of them have sensible default values. All arguments are described in more
#' details in the sections below the list.
#' 
#' @inheritParams emptyParams
#' @inheritParams setInteraction
#' @inheritParams setPredKernel
#' @inheritParams setSearchVolume
#' @inheritParams setMaxIntakeRate
#' @inheritParams setMetabolicRate
#' @inheritParams setExtMort
#' @inheritParams setReproduction
#' @inheritParams setFishing
#' @inheritParams setResource
#' @param min_w_pp The smallest size of the resource spectrum. By default this
#'   is set to the smallest value at which any of the consumers can feed.
#' @param n The allometric growth exponent. This can be overruled for individual
#'   species by including a `n` column in the `species_params`. 
#'
#' @return An object of type \linkS4class{MizerParams}
#' 
#' @section Species parameters:
#' The only essential argument is a data frame that contains the species
#' parameters. The data frame is arranged species by parameter, so each column
#' of the parameter data frame is a parameter and each row has the values of the
#' parameters for one of the species in the model.
#'
#' There are two essential columns that must be included in the species
#' parameter data.frame and that do not have default values: the 
#' `species` column that should hold strings with the names of the
#' species and the `w_inf` column with the asymptotic sizes of the species. 
#' 
#' The species_params dataframe also needs to contain the parameters needed
#' by any predation kernel function or size selectivity function. This will
#' be mentioned in the appropriate sections below.
#' 
#' For all other species parameters, mizer will calculate default values if they
#' are not included in the species parameter data frame. They will be
#' automatically added when the `MizerParams` object is created. For these
#' parameters you can also specify values for only some species and leave the
#' other entries as NA and the missing values will be set to the defaults.
#' 
#' All the parameters will be mentioned in the following sections.
#' @inheritSection emptyParams Changes to species params
#' @inheritSection emptyParams Size grid
#' @inheritSection setParams Units in mizer
#' @inheritSection setInteraction Setting interactions
#' @inheritSection setPredKernel Setting predation kernel
#' @inheritSection setSearchVolume Setting search volume
#' @inheritSection setMaxIntakeRate Setting maximum intake rate
#' @inheritSection setMetabolicRate Setting metabolic rate
#' @inheritSection setExtMort Setting external mortality rate
#' @inheritSection setReproduction Setting reproduction
#' @inheritSection setFishing Setting fishing
#' @inheritSection setResource Setting resource dynamics
#' 
#' @export
#' @family functions for setting up models
#' @examples
#' \dontrun{
#' params <- newMultispeciesParams(NS_species_params_gears, inter)
#' }
newMultispeciesParams <- function(
    species_params,
    interaction = NULL,
    no_w = 100,
    min_w = 0.001,
    max_w = NA,
    min_w_pp = NA,
    # setPredKernel()
    pred_kernel = NULL,
    # setSearchVolume()
    search_vol = NULL,
    # setMaxIntakeRate()
    intake_max = NULL,
    # setMetabolicRate()
    metab = NULL,
    p = 0.7,
    # setExtMort
    z0 = NULL,
    z0pre = 0.6,
    z0exp = n - 1,
    # setReproduction
    maturity = NULL,
    repro_prop = NULL,
    RDD = "BevertonHoltRDD",
    # setResource
    resource_rate = NULL,
    resource_capacity = NULL,
    n = 2 / 3,
    r_pp = 10,
    kappa = 1e11,
    lambda = 2.05,
    w_pp_cutoff = 10,
    resource_dynamics = "resource_semichemostat",
    # setFishing
    gear_params = data.frame(),
    selectivity = NULL,
    catchability = NULL,
    initial_effort = NULL) {
    no_sp <- nrow(species_params)
    
    ## For backwards compatibility, allow r_max instead of R_max
    if (!("R_max" %in% names(species_params)) &&
        "r_max" %in% names(species_params)) {
        names(species_params)[names(species_params) == "r_max"] <- "R_max"
    }
    
    ## Create MizerParams object ----
    params <- emptyParams(species_params,
                          gear_params,
                          no_w = no_w, 
                          min_w = min_w,  
                          max_w = max_w, 
                          min_w_pp = min_w_pp)
    
    ## Fill the slots ----
    params <- params %>% 
        set_species_param_default("n", n) %>% 
        set_species_param_default("p", p)
    params <- set_species_param_default(params, "q", 
                                        lambda - 2 + params@species_params$n)
    if (is.null(interaction)) {
        interaction <- matrix(1, nrow = no_sp, ncol = no_sp)
    }
    params <-
        setParams(params,
                  # setInteraction
                  interaction = interaction,
                  # setPredKernel()
                  pred_kernel = pred_kernel,
                  # setSearchVolume()
                  search_vol = search_vol,
                  # setMaxIntakeRate()
                  intake_max = intake_max,
                  # setMetabolicRate()
                  metab = metab,
                  # setExtMort
                  z0 = z0,
                  z0pre = z0pre,
                  z0exp = z0exp,
                  # setReproduction
                  maturity = maturity,
                  repro_prop = repro_prop,
                  RDD = RDD,
                  # setResource
                  resource_rate = resource_rate,
                  resource_capacity = resource_capacity,
                  r_pp = r_pp,
                  kappa = kappa,
                  lambda = lambda,
                  n = n,
                  w_pp_cutoff = w_pp_cutoff,
                  resource_dynamics = resource_dynamics,
                  # setFishing
                  gear_params = gear_params,
                  selectivity = selectivity,
                  catchability = catchability,
                  initial_effort = initial_effort)
    
    params@initial_n <- get_initial_n(params)
    params@initial_n_pp <- params@cc_pp
    params@A <- rep(1, nrow(species_params))
    
    return(params)
}

#' Set or change any model parameters
#' 
#' This is a convenient wrapper function calling each of the following
#' functions
#' \itemize{
#' \item [setPredKernel()]
#' \item [setSearchVolume()]
#' \item [setInteraction()]
#' \item [setMaxIntakeRate()]
#' \item [setMetabolicRate()]
#' \item [setExtMort()]
#' \item [setReproduction()]
#' \item [setFishing()]
#' \item [setResource()]
#' }
#' See the Details section below for a discussion of how to use this function.
#' 
#' @param params A \linkS4class{MizerParams} object
#' @inheritParams setInteraction
#' @inheritDotParams setPredKernel
#' @inheritDotParams setSearchVolume
#' @inheritDotParams setMaxIntakeRate
#' @inheritDotParams setMetabolicRate
#' @inheritDotParams setExtMort
#' @inheritDotParams setReproduction
#' @inheritDotParams setFishing
#' @inheritDotParams setResource
#' 
#' @return A \linkS4class{MizerParams} object
#' 
#' @details 
#' Usually, if you are happy with the way mizer calculates its model functions
#' from the species parameters and only want to change the values of some
#' species parameters, you would make those changes in the `species_params` data
#' frame contained in the `params` object and then call the `setParams()`
#' function to effect the change. Note that just changing the species parameters
#' by themselves is not changing the model until you call `setParams()` or the
#' appropriate one of its sub-functions. Here is an example which assumes that
#' you have have a MizerParams object `params` in which you just want to change
#' one parameter of the third species:
#' ```
#' params@species_params$gamma[[3]] <- 1000
#' params <- setParams(params)
#' ```
#' Because of the way the R language works, `setParams` does not make the
#' changes to the `params` object that you pass to it but instead returns a new
#' params object. So to affect the change you call the function in the form
#' `params <- setParams(params, ...)`.
#' 
#' If you are not happy with the assumptions that mizer makes by default about
#' the shape of the model functions, for example if you want to change one of
#' the allometric scaling assumptions, you can do this by providing your
#' choice as an array in the appropriate argument to `setParams()`. The
#' sections below discuss all the model functions that you can change this way.
#' 
#' This function will use the species parameters in the `params` object to reset
#' the values of all the model functions that you do not specify explicitly when
#' calling this function, unless you have protected the corresponding slots with
#' a comment. If you have changed any of the model functions in the
#' `params` object previously and now want to make changes to a different slot,
#' you will want to call the appropriate change function individually. So in the
#' above example you would have used `params <- setSearchVolume(params)`
#' instead of `params <- setParams(params)`. 
#' 
#' If you have added a comment to a slot of the params object, then setParams()
#' and its subfunctions will not recalculate the value for that slot from the
#' species parameters. For example
#' ```
#' comment(params@search_vol) <- "This should not change"
#' species_params(params)$gamma <- 10
#' ```
#' will just issue a warning "The search volume has been commented and therefore
#' will not be recalculated from the species parameters". You can remove the
#' comment, and therefore allow recalculation of the slot, with
#' `comment(params@search_vol) <- NULL`.
#' 
#' @section Units in mizer:
#' Mizer uses grams to measure weight, centimetres to measure lengths, and
#' years to measure time.
#' 
#' Mizer is agnostic about whether abundances are given as 
#' 1. numbers per area, 
#' 2. numbers per volume or
#' 3. total numbers for the entire study area. 
#' 
#' You should make the choice most convenient for your application and then
#' stick with it. If you make choice 1 or 2 you will also have to choose a unit
#' for area or volume. Your choice will then determine the units for some of
#' the parameters. This will be mentioned when the parameters are discussed in
#' the sections below.
#' 
#' You choice will also affect the units of the quantities you may want to
#' calculate with the model. For example, the yield will be in grams/year/m^2 in
#' case 1 if you choose m^2 as your measure of area, in grams/year/m^3 in case 2
#' if you choose m^3 as your unit of volume, or simply grams/year in case 3. The
#' same comment applies for other measures, like total biomass, which will be
#' grams/area in case 1, grams/volume in case 2 or simply grams in case 3. When
#' mizer puts units on axes, for example in `plotBiomass`, it will simply
#' put grams, as appropriate for case 3.
#' 
#' You can convert between these choices. For example, if you use case 1, you
#' need to multiply with the area of the ecosystem to get the total quantity. 
#' If you work with case 2, you need to multiply by both area and the thickness 
#' of the productive layer. In that respect, case 2 is a bit cumbersome.
#' 
#' @inheritSection setInteraction Setting interactions
#' @inheritSection setPredKernel Setting predation kernel
#' @inheritSection setSearchVolume Setting search volume
#' @inheritSection setMaxIntakeRate Setting maximum intake rate
#' @inheritSection setMetabolicRate Setting metabolic rate
#' @inheritSection setExtMort Setting external mortality rate
#' @inheritSection setReproduction Setting reproduction
#' @inheritSection setFishing Setting fishing
#' @inheritSection setResource Setting resource dynamics
#' @export
#' @family functions for setting parameters
#' @examples
#' \dontrun{
#' params <- newTraitParams()
#' params@species_params$gamma[3] <- 1000
#' params <- setParams(params)
#' }
setParams <- function(params, interaction = NULL, ...) {
    validObject(params)
    params <- setResource(params, ...)
    params <- setInteraction(params, interaction)
    params <- setPredKernel(params, ...)
    params <- setMaxIntakeRate(params, ...)
    params <- setMetabolicRate(params, ...)
    params <- setExtMort(params, ...)
    # setSearchVolume() should be called only after 
    # setMaxIntakeRate() and setPredKernel()
    params <- setSearchVolume(params, ...)
    params <- setReproduction(params, ...)
    params <- setFishing(params, ...)
    return(params)
}
