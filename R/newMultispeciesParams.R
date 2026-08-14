#' Set up parameters for a general multispecies model
#'
#' Sets up a multi-species size spectrum model by filling all slots in the
#' \linkS4class{MizerParams} object based on user-provided or default
#' parameters. There is a long list of arguments, but almost
#' all of them have sensible default values. The only required argument is
#' the `species_params` data frame. All arguments are described in more
#' details in the sections below the list.
#'
#' @inheritParams emptyParams
#' @inheritParams setInteraction
#' @inheritParams setPredKernel
#' @inheritParams setSearchVolume
#' @inheritParams setMaxIntakeRate
#' @inheritParams setMetabolicRate
#' @inheritParams setExtMort
#' @inheritParams setExtEncounter
#' @inheritParams setReproduction
#' @inheritParams setFishing
#' @inheritParams setResource
#' @param r_pp `r lifecycle::badge("deprecated")`. Use `resource_rate` argument
#'   instead.
#' @param kappa The coefficient \eqn{\kappa} of the resource carrying capacity
#'   power law \eqn{c_R(w) = \kappa\, w^{-\lambda}}, which also sets the initial
#'   resource abundance. See [resource_params()].
#' @param min_w_pp The smallest size of the resource spectrum. By default this
#'   is set to the smallest value at which any of the consumers can feed.
#' @param n The allometric growth exponent. This can be overruled for individual
#'   species by including a `n` column in the `species_params`.
#' @param p The allometric metabolic exponent. This can be overruled for
#'   individual species by including a `p` column in the `species_params`.
#' @param z0pre If `z0` is missing from the species parameter data frame,
#'   calculate it as `z0pre * w_inf ^ z0exp`. Defaults to 0.6.
#' @param z0exp The exponent used with `z0pre` to calculate a missing `z0`.
#'   Defaults to `n - 1`.
#' @param second_order_w `r lifecycle::badge("experimental")` Selects the
#'   second-order numerical scheme for the new model. Accepts the same values as
#'   the [second_order_w()] setter: a single logical (`TRUE` switches on both
#'   second-order flux and bin-averaging), a single flux scheme name
#'   (`"upwind"`, `"van_leer"` or `"centred"`), or a named vector with entries
#'   `flux` and/or `bin_average`. The `bin_average` choice is applied *before*
#'   the resource and abundance power laws are constructed, so they are built
#'   bin-averaged from the start (unlike setting [second_order_w()] on an
#'   existing object). The `flux` scheme governs time projection only, so the
#'   robust first-order upwind scheme is used for the construction-time
#'   steady-state solve and the chosen scheme is then activated for the returned
#'   model. Defaults to `FALSE` (the first-order behaviour of previous mizer).
#' @param info_level Controls the amount of information messages that are shown
#'   when the function sets default values for parameters. Higher levels lead
#'   to more messages, `info_level = 0` gives silence. The default is taken
#'   from the `mizer_info_level` option, see [default_info_level()].
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
#' species and the `w_inf` column with the von Bertalanffy asymptotic sizes of
#' the species in grams. (You could alternatively specify the corresponding
#' length in cm in an `l_inf` column.) The computational upper size boundary
#' `w_max` is not essential; if it is missing it defaults to `1.5 * w_inf`. For
#' backwards compatibility, if `w_inf` is missing it is taken from the
#' `w_repro_max` or `w_max` column instead.
#'
#' The `species_params dataframe` also needs to contain the parameters needed
#' by any predation kernel function (size selectivity function). This will
#' be mentioned in the appropriate sections below.
#'
#' For all other species parameters, mizer will calculate default values if they
#' are not included in the species parameter data frame. They will be
#' automatically added when the `MizerParams` object is created. For these
#' parameters you can also specify values for only some species and leave the
#' other entries as NA and the missing values will be set to the defaults.
#' So the `species_params` data frame saved in the returned MizerParams object
#' will differ from the one you supply because it will have the missing
#' species parameters filled in with default values.
#'
#' If you are not happy with any of the species parameter values used you can
#' always change them later with [species_params<-()].
#'
#' All the parameters will be mentioned in the following sections.
#' @inheritSection emptyParams Size grid
#' @inheritSection setParams Units in mizer
#' @inheritSection setInteraction Setting interaction matrix
#' @inheritSection setPredKernel Setting predation kernel
#' @inheritSection setSearchVolume Setting search volume
#' @inheritSection setMaxIntakeRate Setting maximum intake rate
#' @inheritSection setMetabolicRate Setting metabolic rate
#' @inheritSection setExtMort Setting external mortality rate
#' @inheritSection setExtEncounter Setting external encounter rate
#' @inheritSection setExtDiffusion Setting external diffusion rate
#' @inheritSection setReproduction Setting reproduction
#' @inheritSection setFishing Setting fishing
#' @inheritSection setResource Setting resource dynamics
#'
#' @section Setting initial values:
#' The initial values for the species number densities are set using the
#' function `get_initial_n()`. These are quite arbitrary and not very close to
#' the steady state abundances. We intend to improve this in the future.
#'
#' The initial resource number density \eqn{N_R(w)} is set to a power law with
#' coefficient `kappa` (\eqn{\kappa}) and exponent `-lambda` (\eqn{-\lambda}):
#' \deqn{N_R(w) = \kappa\, w^{-\lambda}}{c_R(w) = \kappa w^{-\lambda}}
#' for all \eqn{w} less than `w_pp_cutoff` and zero for sizes at or above
#' `w_pp_cutoff`.
#'
#' @export
#' @family functions for setting up models
#' @examples
#' params <- newMultispeciesParams(NS_species_params)
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
    # setExtEncounter
    ext_encounter = NULL,
    # setReproduction
    maturity = NULL,
    repro_prop = NULL,
    RDD = "BevertonHoltRDD",
    # setResource
    kappa = 1e11,
    n = 2 / 3,
    resource_rate = 10,
    resource_capacity = kappa,
    lambda = 2.05,
    w_pp_cutoff = 10,
    resource_dynamics = "resource_semichemostat",
    # setFishing
    gear_params = NULL,
    selectivity = NULL,
    catchability = NULL,
    initial_effort = NULL,
    second_order_w = FALSE,
    info_level = default_info_level(),
    z0 = deprecated(),
    r_pp = deprecated()) {

    if (lifecycle::is_present(r_pp)) {
        lifecycle::deprecate_warn("1.0.0", "newMultispeciesParams(r_pp)",
                                  "newMultispeciesParams(resource_rate)")
        resource_rate <- r_pp
    }
    if (lifecycle::is_present(z0)) {
        lifecycle::deprecate_warn("2.2.3", "newMultispeciesParams(z0)",
                                  "newMultispeciesParams(ext_mort)")
        ext_mort <- z0
    }

    # Collect the information signals raised while the model is built and
    # report them together at the end.
    with_info_level(info_level = info_level, {
    no_sp <- nrow(species_params)
    z0_missing <- if ("z0" %in% names(species_params)) {
        is.na(species_params$z0)
    } else {
        rep(TRUE, no_sp)
    }
    assert_that(is.number(z0pre), z0pre >= 0,
                is.number(z0exp))

    species_params <- set_species_param_default(species_params, "n", n)
    species_params <- set_species_param_default(species_params, "p", p)
    given_species_params <- validGivenSpeciesParams(species_params)

    species_params <- validSpeciesParams(species_params)
    gear_params <- validGearParams(gear_params, species_params)

    ## Create MizerParams object ----
    params <- emptyParams(given_species_params,
                          gear_params,
                          no_w = no_w,
                          min_w = min_w,
                          max_w = max_w,
                          min_w_pp = min_w_pp)

    # Resolve the requested second-order scheme. The `bin_average` choice is
    # applied now, before the rest of the parameters are computed, so that the
    # bin-averaged resource constructions (initial abundance, capacity, rate)
    # and the encounter/mortality quadratures all honour it. We set the slot
    # directly rather than via the `second_order_w<-` setter to avoid its extra
    # setParams() call here, since setParams() is run below anyway. The `flux`
    # scheme is deferred until the end (see below): it does not affect any
    # precomputed array, only time projection, and the construction-time
    # steady-state solve for the initial abundances is only robust with the
    # first-order upwind scheme.
    target_second_order_w <-
        validate_second_order_w(params@second_order_w, second_order_w)
    params@second_order_w[["bin_average"]] <-
        target_second_order_w[["bin_average"]]

    # Fill the slots ----
    if (is.null(interaction)) {
        interaction <- matrix(1, nrow = no_sp, ncol = no_sp)
    }

    params@initial_n_pp[] <- resource_power_law(params, kappa, lambda,
                                                w_max = w_pp_cutoff)
    params@resource_params$kappa <- kappa
    params@resource_params$lambda <- lambda
    # Needed by rate defaults, including the construction default for `z0`,
    # before `setResource()` installs the complete resource parameters below.
    params@resource_params$n <- n
    params@resource_params$w_pp_cutoff <- w_pp_cutoff

    params <- setParams(
        params,
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
        # setExtEncounter
        ext_encounter = ext_encounter,
        # setReproduction
        maturity = maturity,
        repro_prop = repro_prop,
        RDD = RDD,
        # setFishing. The `gear_params` are not passed here: they are already
        # stored in the params object by emptyParams() above and setFishing()
        # reads them from there.
        selectivity = selectivity,
        catchability = catchability,
        initial_effort = initial_effort)

    # `z0pre` and `z0exp` are construction arguments. `setParams()` has
    # installed the ordinary `z0` default at the same point in the species
    # table as before; now replace only the values that were missing from the
    # user's table with the requested construction default and update external
    # mortality. Do not add these calculated values to `given_species_params`.
    if (is.null(ext_mort) && any(z0_missing)) {
        z0_default <- z0pre * params@species_params$w_inf^z0exp
        params@species_params$z0[z0_missing] <- z0_default[z0_missing]
        params <- setExtMort(params)
    }

    params <- setResource(
        params,
        resource_rate = resource_rate,
        resource_capacity = resource_capacity,
        resource_dynamics = resource_dynamics,
        lambda = lambda,
        n = n,
        w_pp_cutoff = w_pp_cutoff,
        balance = FALSE)

    params@initial_n[] <- get_initial_n(params)
    # TODO: The next line can be removed after release of mizer 3.0
    params@A <- rep(1, nrow(species_params))
    })
    # Now that the initial abundances have been computed with the robust upwind
    # solver, switch on the requested advective-flux scheme for projection.
    params@second_order_w[["flux"]] <- target_second_order_w[["flux"]]
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
#' \item [setExtEncounter()]
#' \item [setExtDiffusion()]
#' \item [setReproduction()]
#' \item [setFishing()]
#' }
#' Note that [setResource()] is **not** among them: the resource rate, capacity
#' and dynamics are not changed by `setParams()` and have to be set with
#' [setResource()]. Passing a resource argument to `setParams()` gives an error
#' rather than being silently ignored. See the Details section below for a
#' discussion of how to use this function.
#'
#' @param object A \linkS4class{MizerParams} object
#' @inheritParams setInteraction
#' @inheritDotParams setPredKernel -reset
#' @inheritDotParams setSearchVolume -reset
#' @inheritDotParams setMaxIntakeRate -reset
#' @inheritDotParams setMetabolicRate -reset -p
#' @inheritDotParams setExtMort -reset
#' @inheritDotParams setExtEncounter -reset
#' @inheritDotParams setExtDiffusion -reset
#' @inheritDotParams setReproduction -reset
#' @inheritDotParams setFishing -reset
#' @param reset If set to TRUE then all the rate arrays that `setParams()` sets
#'   are recalculated from the species parameters, even if they had previously
#'   been overwritten with custom values. The default is FALSE, in which case
#'   arrays that have been set manually are left alone. This is passed on to
#'   each of the setter functions listed above, so it thaws all of them at once.
#'   To thaw only one of them, call that setter with `reset = TRUE` instead.
#' @param info_level Controls the amount of information messages that are shown.
#'   Higher levels lead to more messages, `info_level = 0` gives silence. The
#'   default is taken from the `mizer_info_level` option, see
#'   [default_info_level()]. Note that the report that a change cannot take
#'   effect because the rate array it feeds has been set manually is a warning
#'   rather than a message, but it too is silenced by `info_level = 0`.
#'
#' @return A \linkS4class{MizerParams} object
#'
#' @details
#' If you are not happy with the assumptions that mizer makes by default about
#' the size-dependence of parameters, for example if you want to change one of
#' the allometric scaling assumptions, you can do this by providing your
#' choice as an array in the appropriate argument to `setParams()`. The
#' sections below discuss all the model functions that you can change this way.
#'
#' Because of the way the R language works, `setParams` does not make the
#' changes to the `params` object that you pass to it but instead returns a new
#' params object. So to affect the change you call the function in the form
#' `params <- setParams(params, ...)`.
#'
#' Usually, if you are happy with the way mizer calculates the size-dependent parameters
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
#' recalculate the values of all the parameter arrays except those for which you
#' have set custom values.
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
#' Your choice will also affect the units of the quantities you may want to
#' calculate with the model. For example, the yield will be in grams/year/m^2 in
#' case 1 if you choose m^2 as your measure of area, in grams/year/m^3 in case 2
#' if you choose m^3 as your unit of volume, or simply grams/year in case 3. The
#' same comment applies for other measures, like total biomass, which will be
#' grams/area in case 1, grams/volume in case 2 or simply grams in case 3. When
#' mizer puts units on axes in plots, it will choose the units appropriate for
#' case 3. So for example in [plotBiomass()] it gives the unit as grams.
#'
#' You can convert between these choices. For example, if you use case 1, you
#' need to multiply with the area of the ecosystem to get the total quantity.
#' If you work with case 2, you need to multiply by both area and the thickness
#' of the productive layer. In that respect, case 2 is a bit cumbersome. The
#' function [scaleModel()] is useful to change the units you are using.
#'
#' @inheritSection setInteraction Setting interaction matrix
#' @inheritSection setPredKernel Setting predation kernel
#' @inheritSection setSearchVolume Setting search volume
#' @inheritSection setMaxIntakeRate Setting maximum intake rate
#' @inheritSection setMetabolicRate Setting metabolic rate
#' @inheritSection setExtMort Setting external mortality rate
#' @inheritSection setExtEncounter Setting external encounter rate
#' @inheritSection setExtDiffusion Setting external diffusion rate
#' @inheritSection setReproduction Setting reproduction
#' @inheritSection setFishing Setting fishing
#' @export
#' @family functions for setting parameters
# The reason we list `interaction` explicitly rather than including it in
# the `...` is for backwards compatibility. It used to be the second argument.
setParams <- function(object, interaction = NULL,
                      info_level = default_info_level(), ..., reset = FALSE) {
    UseMethod("setParams")
}
#' @export
setParams.MizerParams <- function(object, interaction = NULL,
                                  info_level = default_info_level(), ...,
                                  reset = FALSE) {
    # The setters below all declare `... Unused`, so without this check any
    # misspelled or misplaced argument would be silently ignored.
    check_set_params_dots(match.call(expand.dots = FALSE)[["..."]])
    assert_that(is.flag(reset))

    # Collect the information signals raised by the individual setters and
    # report them together at the end.
    with_info_level(info_level = info_level, {
    params <- validParams(object)

    params <- setInteraction(params, interaction)
    params <- setPredKernel(params, ..., reset = reset)
    params <- setMaxIntakeRate(params, ..., reset = reset)
    params <- setMetabolicRate(params, ..., reset = reset)
    params <- setExtMort(params, ..., reset = reset)
    params <- setExtEncounter(params, ..., reset = reset)
    params <- setExtDiffusion(params, ..., reset = reset)
    # setSearchVolume() should be called only after
    # setMaxIntakeRate() and setPredKernel()
    params <- setSearchVolume(params, ..., reset = reset)
    params <- setReproduction(params, ..., reset = reset)
    params <- setFishing(params, ..., reset = reset)

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
    })
    params
}

# The setter functions that `setParams()` calls. `setInteraction()` is not in
# the list because its `interaction` argument is a named argument of
# `setParams()` rather than part of the `...`.
set_params_setters <- c("setPredKernel", "setMaxIntakeRate", "setMetabolicRate",
                        "setExtMort", "setExtEncounter", "setExtDiffusion",
                        "setSearchVolume", "setReproduction", "setFishing")

# Names of the arguments that `setParams()` passes on to its setters. Collected
# from the formals of the setters themselves so that this stays correct when a
# setter gains or loses an argument.
set_params_arg_names <- function() {
    args <- unlist(lapply(set_params_setters,
                          function(f) names(formals(get(f)))))
    # `params` is supplied by `setParams()` and `reset` is one of its own
    # arguments.
    setdiff(unique(args), c("params", "object", "reset", "..."))
}

# Where a user is likely to have meant an argument that `setParams()` does not
# have. Anything not listed here is reported as an unknown argument.
set_params_elsewhere <- c(
    resource_rate     = "`setResource()`",
    resource_capacity = "`setResource()`",
    resource_level    = "`setResource()`",
    resource_dynamics = "`setResource()`",
    r_pp              = "`setResource()`",
    kappa             = "`setResource()`",
    lambda            = "`setResource()`",
    w_pp_cutoff       = "`setResource()`",
    balance           = "`setResource()`",
    n                 = "`species_params<-()` or `setResource()`",
    gear_params       = "`gear_params<-()`",
    species_params    = "`species_params<-()`",
    interaction_resource = "`species_params<-()`",
    second_order_w    = "`second_order_w<-()`",
    use_predation_diffusion = "`use_predation_diffusion<-()`"
)

#' Check the arguments passed to `setParams()` in its `...`
#'
#' Each of the setters that [setParams()] calls declares `... Unused`, so an
#' argument that none of them recognises would otherwise be accepted in silence.
#'
#' @param dots The `...` of the call to [setParams()], as returned by
#'   `match.call(expand.dots = FALSE)[["..."]]`. Using the unevaluated call
#'   rather than `list(...)` keeps the arguments from being evaluated twice.
#' @return NULL, invisibly. Called for the error it throws.
#' @noRd
check_set_params_dots <- function(dots) {
    if (length(dots) == 0) {
        return(invisible(NULL))
    }
    names <- names(dots)
    if (is.null(names) || any(names == "")) {
        stop("All arguments to `setParams()` after the first must be named.")
    }
    unknown <- setdiff(names, set_params_arg_names())
    if (length(unknown) == 0) {
        return(invisible(NULL))
    }
    elsewhere <- intersect(unknown, names(set_params_elsewhere))
    unrecognised <- setdiff(unknown, elsewhere)
    msg <- paste0("`setParams()` does not have ",
                  if (length(unknown) == 1) "an argument " else "arguments ",
                  paste0("`", unknown, "`", collapse = ", "), ".")
    if (length(elsewhere) > 0) {
        msg <- c(msg,
                 paste0("Use ", set_params_elsewhere[elsewhere],
                        " to set `", elsewhere, "`."))
    }
    if (length(unrecognised) > 0) {
        msg <- c(msg,
                 paste0("Check for a typo in ",
                        paste0("`", unrecognised, "`", collapse = ", "),
                        ", or set it as a species parameter with ",
                        "`species_params<-()`."))
    }
    stop(paste(msg, collapse = "\n"))
}
