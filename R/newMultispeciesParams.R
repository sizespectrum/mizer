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
#' @param info_level Controls the amount of information messages that are shown
#'   when the function sets default values for parameters. 
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
#' If you are not happy with any of the species parameter values used you can
#' always change them later with [species_params<-()].
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
#' params <- newMultispeciesParams(NS_species_params_gears, NS_interaction)
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
    ext_mort = NULL,
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
    gear_params = NULL,
    selectivity = NULL,
    catchability = NULL,
    initial_effort = NULL,
    info_level = 3,
    z0 = deprecated()) {
    if (lifecycle::is_present(z0)) {
        lifecycle::deprecate_warn("2.2.3", "newMultispeciesParams(z0)", "newMultispeciesParams(ext_mort)")
        ext_mort <- z0
    }
    
    # Define a signal handler that collects the information signals
    # into the `infos` list.
    infos <- list()
    collect_info <- function(cnd) {
        if (cnd$level >= info_level) {
            infos[[cnd$var]] <<- cnd$message
        }
    }
    # Register this signal handler
    withCallingHandlers(
        info_about_default = collect_info, {
    no_sp <- nrow(species_params)
    species_params <- validSpeciesParams(species_params)
    gear_params <- validGearParams(gear_params, species_params)
    
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
                  ext_mort = ext_mort,
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
    })
    if (length(infos) > 0) {
        message(paste(infos, collapse = "\n"))
    }
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
#' \item [setMetadata()]
#' }
#' See the Details section below for a discussion of how to use this function.
#' 
#' @param params A \linkS4class{MizerParams} object
#' @inheritParams setInteraction
#' @inheritDotParams setPredKernel -reset
#' @inheritDotParams setSearchVolume -reset
#' @inheritDotParams setMaxIntakeRate -reset
#' @inheritDotParams setMetabolicRate -reset
#' @inheritDotParams setExtMort -reset
#' @inheritDotParams setReproduction -reset
#' @inheritDotParams setFishing -reset
#' @inheritDotParams setResource -reset
#' @inheritDotParams setMetadata
#' 
#' @return A \linkS4class{MizerParams} object
#' 
#' @details 
#' If you are not happy with the assumptions that mizer makes by default about
#' the shape of the model functions, for example if you want to change one of
#' the allometric scaling assumptions, you can do this by providing your
#' choice as an array in the appropriate argument to `setParams()`. The
#' sections below discuss all the model functions that you can change this way.
#' 
#' Because of the way the R language works, `setParams` does not make the
#' changes to the `params` object that you pass to it but instead returns a new
#' params object. So to affect the change you call the function in the form
#' `params <- setParams(params, ...)`.
#' 
#' Usually, if you are happy with the way mizer calculates its model functions
#' from the species parameters and only want to change the values of some
#' species parameters, you would make those changes in the `species_params` data
#' frame contained in the `params` object using [species_params<-()]. 
#' Here is an example which assumes that
#' you have have a MizerParams object `params` in which you just want to change
#' the `gamma` parameter of the third species:
#' ```
#' species_params(params)$gamma[[3]] <- 1000
#' ```
#' Internally that will actually call `setParams()` to recalculate any of the
#' other parameters that are affected by the change in the species parameter.
#' 
#' `setParams()` will use the species parameters in the `params` object to
#' recalculate the values of all the model functions except those for which you have set custom
#' values.
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
#' @inheritSection setMetadata Setting model metadata
#' @export
#' @family functions for setting parameters
# The reason we list `interaction` explicitly rather than including it in
# the `...` is for backwards compatibitlity. It used to be the second argument.
setParams <- function(params, interaction = NULL, ...) {
    params <- suppressWarnings(validParams(params))
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
    params <- setMetadata(params, ...)
    
    colours <- params@species_params$linecolour
    if (!is.null(colours)) {
        names(colours) <- params@species_params$species
        params <- setColours(params, colours)
    }
    linetypes <- params@species_params$linetype
    if (!is.null(linetypes)) {
        names(linetypes) <- params@species_params$species
        params <- setLinetypes(params, linetypes)
    }
    
    validObject(params)
    params
}
