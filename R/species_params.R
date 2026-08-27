#' Species parameters
#'
#' These functions allow you to get or set the species-specific parameters
#' stored in a MizerParams object.
#'
#'
#' There are a lot of species parameters and we will list them all below, but
#' most of them have sensible default values. The only required columns are
#' `species` for the species name and `w_inf` for its von Bertalanffy
#' asymptotic size. However if you have information about the values of other
#' parameters then you should provide them.
#'
#' Three species parameters describe maximum sizes and play distinct roles:
#'
#' * `w_inf` is the von Bertalanffy asymptotic size of an average individual.
#'   It is the required maximum-size parameter and is used to set default values
#'   for `w_max`, `w_repro_max` and `w_mat`.
#' * `w_repro_max` is the size at which a typical mature individual invests all
#'   of its available energy into reproduction, see [setReproduction()]. It is
#'   not a hard ceiling on size and defaults to `w_inf`.
#' * `w_max` is purely a computational boundary: it sets the upper end of the
#'   size grid and the range of plots. It defaults to `1.5 * w_inf`. For
#'   backwards compatibility, if `w_inf` is not supplied it is taken from
#'   `w_repro_max` or `w_max` instead.
#'
#' Mizer distinguishes between the species parameters that you have given
#' explicitly and the species parameters that have been calculated by mizer or
#' set to default values. You can retrieve the given species parameters with
#' `given_species_params()` and the calculated ones with
#' `calculated_species_params()`. You get all species_params with
#' `species_params()`.
#'
#' When you change species parameters with `species_params<-()`, mizer
#' automatically detects which parameters you have changed. It records these
#' changed parameters in `given_species_params` so that they are protected
#' against being overwritten by future recalculations. It then re-calculates
#' the quantities that depend on the changed parameters. Changes to observation,
#' direct-runtime or other custom columns that base mizer does not use to build
#' a cached quantity do not trigger that recalculation. Unknown columns on an
#' extension object retain the conservative recalculation path because an
#' extension setter may use them.
#'
#' There are some species parameters that are used to set up the
#' size-dependent parameters that are used in the mizer model:
#'
#' * `gamma` and `q` are used to set the search volume, see [setSearchVolume()].
#' * `h` and `n` are used to set the maximum intake rate, see [setMaxIntakeRate()].
#' * `k`, `ks` and `p` are used to set activity and basic metabolic rate,
#'   see [setMetabolicRate()].
#' * `z0`, `z_ext` and `d` are used to set the external mortality rate, see
#'   [setExtMort()].
#' * `E_ext` and `n` are used to set the external encounter rate, see
#'   [setExtEncounter()].
#' * `D_ext` and `n` are used to set the external diffusion rate, see
#'   [setExtDiffusion()].
#' * `w_mat`, `w_mat25`, `w_repro_max` and `m` are used to set the allocation to
#'   reproduction, see [setReproduction()].
#' * `pred_kernel_type` specifies the shape of the predation kernel. The default
#'   is a "lognormal", for other options see the "Setting predation kernel"
#'   section in the help for [setPredKernel()].
#' * `beta` and `sigma` are parameters of the lognormal predation kernel, see
#'   [lognormal_pred_kernel()]. The Gaussian-mixture kernel instead uses the
#'   list-columns `kernel_p`, `kernel_mean`, and `kernel_sd`, see
#'   [gaussian_mixture_pred_kernel()]. There will be other parameters if you are
#'   using other predation kernel functions.
#'
#' When you change one of the above species parameters using
#' `species_params<-()` or `given_species_params<-()`, the new value will be
#' used to update the corresponding size-dependent rates automatically, unless
#' you have set those size-dependent rates manually, in which case the
#' corresponding species parameters will be ignored. Mizer warns you when that
#' happens, because the value is then in the species parameter table without
#' having any effect on the model. The warning names the call that puts the
#' rate back under the control of the species parameters.
#'
#' There are some species parameters that are used directly in the model
#' rather than being used for setting up size-dependent parameters:
#'
#' * `alpha` is the assimilation efficiency, the proportion of the consumed
#'   biomass that can be used for growth, metabolism and reproduction, see
#'   the help for [getEReproAndGrowth()].
#' * `w_min` is the egg size.
#' * `interaction_resource` sets the interaction strength with the resource,
#'   see "Predation encounter" section in the help for [getEncounter()].
#' * `erepro` is the reproductive efficiency, the proportion of the energy
#'   invested into reproduction that is converted to egg biomass, see
#'   [getRDI()].
#' * `R_max` is the parameter in the Beverton-Holt density dependence added to
#'   the reproduction, see [setBevertonHolt()]. There will be other such
#'   parameters if you use other density dependence functions, see the
#'   "Density dependence" section in the help for [setReproduction()].
#'
#' Two parameters are used only by functions that need to convert between
#' weight and length:
#'
#' * `a` and `b` are the parameters in the allometric weight-length
#'   relationship \eqn{w = a l ^ b}.
#'
#' If you have supplied the `a` and `b` parameters, then you can replace weight
#' parameters like `w_inf`, `w_max`, `w_mat`, `w_mat25`, `w_repro_max` and
#' `w_min` by their corresponding length parameters `l_inf`, `l_max`, `l_mat`,
#' `l_mat25`, `l_repro_max` and `l_min`.
#'
#' You can also keep both, and change either of them later. Mizer keeps the two
#' consistent by the rule that the one you gave last wins, and if you gave both
#' at the same time the weight wins. So on a model set up with lengths you can
#' still set `w_mat` with `species_params<-()` and mizer will update `l_mat` to
#' match, and if you set `l_mat` it will update `w_mat` as always. When you
#' supply a length and a weight together that do not agree, mizer uses the
#' weight and warns you that it has changed the length to match.
#'
#' The rule is applied when a species parameter data frame is put into a model.
#' A data frame that you have taken out of a model and are editing on its own
#' is left exactly as you write it: no conversions, no checks and no warnings
#' until you assign it back with `species_params<-()` or
#' `given_species_params<-()`, which is when mizer can tell which values you
#' changed. A data frame that was never in a model, for example one you pass to
#' [validSpeciesParams()], carries no such history, so a length and a weight
#' that disagree there count as given at the same time and the weight wins.
#'
#' The parameters that are only used to calculate default values for other
#' parameters are:
#'
#' * `f0` is the feeding level and is used to get a default value for the
#'   coefficient of the search volume `gamma`, see [get_gamma_default()].
#' * `fc` is the critical feeding level below which the species can not
#'   maintain itself. This is used to get a default value for the coefficient
#'   `ks` of the metabolic rate, see [get_ks_default()].
#' * `age_mat` is the age at maturity and is used to get a default value for
#'   the coefficient `h` of the maximum intake rate, see [get_h_default()].
#' * If `age_mat` is not supplied, mizer used the von Bertalanffy parameters
#'   `k_vb`, `w_inf` and `t0` as well as the weight-length exponent `b` to
#'   determine it. This is unreliable and is therefore not recommended.
#'
#' Changing these parameters with `species_params<-()` will trigger a
#' recalculation of the downstream parameters, provided they are not protected
#' by being explicitly given.
#'
#' There are other species parameters that are used in tuning the model to
#' observations:
#'
#' * `biomass_observed` and `biomass_cutoff` allow you to specify for each
#'   species the total observed biomass above some cutoff size. This is
#'   used by [calibrateBiomass()] and [matchBiomasses()].
#'
#' The total annual fisheries yield is not a species parameter but a gear
#' parameter, because it is observed for each gear separately, see
#' [gear_params()]. For backwards compatibility mizer still accepts a
#' `yield_observed` column in the species parameter data frame, see
#' [get_yield_observed()].
#'
#' Finally there are two species parameters that control the way the species are
#' represented in plots:
#'
#' * `linecolour` specifies the colour and can be any valid R colour value.
#' * `linetype` specifies the line type ("solid", "dashed", "dotted", "dotdash",
#'    "longdash", "twodash" or "blank")
#'
#' Other species-specific information that is related to how the species is
#' fished is specified in a gear parameter data frame, see [gear_params()].
#' However in the case where each species is caught by only a single gear,
#' this information can also optionally be provided as species parameters and
#' [newMultispeciesParams()] will transfer them to the `gear_params` data frame.
#' However changing these parameters later in the species parameter data frames
#' will have no effect.
#'
#' You are allowed to include additional columns in the species parameter
#' data frames. They will simply be ignored by mizer but will be stored in the
#' MizerParams object, in case your own code makes use of them.
#'
#' @param object A MizerParams object, a MizerSim object or a data frame
#' @param params A MizerParams object.
#' @param value A data frame with the new species parameters.
#' @param strict Whether to raise an error, rather than correct silently, for
#'   the inconsistencies that can be corrected. Used internally.
#' @param check_misspellings Whether to report column names that look like
#'   misspellings of standard species parameter names. `TRUE` by default.
#'   Mizer passes `FALSE` when it re-validates a table whose columns it has
#'   already checked, so that the report is made once, when the column is
#'   introduced, rather than again on every rebuild.
#' @param recalculate Whether `species_params<-()` should be allowed to
#'   re-derive calculated species parameters and rates that depend on a changed
#'   parameter. Defaults to `TRUE`; mizer still skips the rebuild when all
#'   changes are to columns with no cached dependants. See the section "Setting
#'   species parameters without recalculation" below before setting it to
#'   `FALSE`.
#' @param x An object to test with `is.species_params()` or
#'   `is.given_species_params()`.
#' @param ... Other arguments passed to methods.
#' @section Extracting a column with `$`:
#' `species_params(params)$w_mat` returns the column as a vector named after the
#' species. Unlike `$` on an ordinary data frame, it does **not** partially
#' match the column name. Partial matching is dangerous here because so many
#' species parameter names are prefixes of others: in a model without
#' length-weight parameters `species_params(params)$a` used to return the
#' `alpha` column and `$b` the `beta` column, complete with species names, so
#' code converting weights to lengths silently got the assimilation efficiency
#' and the preferred predator/prey mass ratio instead. Writing was never
#' partially matched, so reads and writes disagreed about which column `$b`
#' meant.
#'
#' A name that is not a column now gives `NULL`, so
#' `is.null(species_params(params)$foo)` is a reliable way of testing whether a
#' parameter is present. If the name would have partially matched a column
#' under the old behaviour you also get a warning naming that column, because
#' that is exactly the case where existing code changes its meaning. The same
#' holds for [gear_params()].
#'
#' @section Removing a species parameter column:
#' A column that is missing from the table you assign is one you no longer
#' supply, and mizer takes it out of the given species parameters. What happens
#' next depends on whether mizer knows how to produce the parameter itself: one
#' that it calculates comes straight back as a calculated value, while one that
#' it knows nothing about leaves the model altogether. Both setters work this
#' way, so either of these removes a custom column:
#'
#' ```
#' species_params(params)$my_col <- NULL
#' given_species_params(params)$my_col <- NULL
#' ```
#'
#' and either of these hands a parameter you had given back to mizer, exactly as
#' setting it to `NA` in the given species parameters does:
#'
#' ```
#' species_params(params)$gamma <- NULL
#' given_species_params(params)$gamma <- NULL
#' ```
#'
#' This is how an extension package withdraws a species parameter it added when
#' the user switches the extension off, without reaching into the
#' `species_params` slot. Mizer reports the removal at `info_level` 3, see
#' [setParams()].
#'
#' Because the whole table is compared, assigning a table with only some of the
#' model's columns withdraws all the others. There is not much room for
#' surprise: `species_params<-()` validates what you give it, so a table without
#' `species` and one of `w_inf`, `w_max` or `w_repro_max` is an error rather
#' than a partial update. Still, edit the table you got from the accessor rather
#' than building a new one from a handful of columns.
#'
#' @section Setting species parameters without recalculation:
#' `species_params(params, recalculate = FALSE) <- value` records the values you
#' changed among the given species parameters, so that they are not calculated
#' away later, and stores `value` as the species parameters. It then stops
#' there: the calculated species parameters are not re-derived from the given
#' ones, no missing parameters are filled in with their default values, and none
#' of the size-dependent rates are recalculated. Your species parameters are
#' stored as you supplied them, after the same checks and length-to-weight
#' conversions that writing into the `species_params` slot would trigger.
#'
#' This is for code that has worked out a species parameter *together with* the
#' rate array that the parameter determines, for example an optimiser that fits
#' `ks` and the matching `metab`, or `z_ext` and the matching `mu_b`. There the
#' recalculation is not just wasted work but would overwrite the rates the
#' caller has just set.
#'
#' The object you get back is only as consistent as you make it. Mizer will not
#' check that the species parameters you supplied agree with the rate arrays in
#' the object, nor that they agree with the other species parameters that are
#' normally derived from them. Unless you are setting the affected rates
#' yourself, use the default `recalculate = TRUE`.
#' @return `species_params()`: Data frame containing all species parameters
#'   currently stored in the model.
#'
#'   `species_params<-()`: Updates the `given_species_params` with any
#'   parameters you have changed, and recalculates the full species parameter
#'   table and model parameters when a changed column has cached dependants. A
#'   column that `value` does not have is one you no longer supply and is
#'   removed from the given species parameters, see the section "Removing a
#'   species parameter column" below. With `recalculate = FALSE` it only does
#'   the recording and stores the parameters you supplied, see the section
#'   "Setting species parameters without recalculation" below.
#'
#'   `given_species_params()`: Data frame containing the species parameter
#'   values that were supplied explicitly by the user.
#'
#'   `given_species_params<-()`: Replaces the authoritative table of parameters
#'   that are to count as explicit user input. Every non-`NA` entry in `value`
#'   is recorded as given, even when it is numerically equal to the value
#'   currently in `species_params()`. This lets you protect a calculated value
#'   against future recalculation. An `NA` entry, or removal of a column, hands
#'   a previously given parameter back to mizer's calculation -- and where there
#'   is no mizer calculation to hand it back to, removing the column removes the
#'   parameter from the model, see the section "Removing a species parameter
#'   column" below. Dependent quantities are recalculated only when the
#'   replacement can change them; merely marking the current value as given does
#'   not rebuild the model.
#'
#'   This setter also warns when a change you asked for cannot take effect,
#'   namely when the parameter is
#'   overridden by another one you have already given (`f0` by `gamma`, `fc` by
#'   `ks`, `age_mat` by `h`), when the rate array it feeds has been set by hand
#'   and so is no longer calculated, or when it is a gear parameter that mizer
#'   reads from [gear_params()] instead. `species_params<-()` stays quiet about
#'   all three. `given_species_params<-()` has no `recalculate` argument;
#'   where you need to record values without recalculating, use
#'   `species_params<-()` or [record_given_species_params()].
#'
#'   `calculated_species_params()`: Data frame containing only those species
#'   parameter entries that are not explicit user input. Columns that would
#'   consist entirely of `NA` values are dropped.
#' @export
#' @seealso [validSpeciesParams()], [setParams()]
#' @family functions for setting parameters
species_params <- function(object, ...) {
    UseMethod("species_params")
}

#' @rdname species_params
#' @usage NULL
#' @export
species_params.MizerParams <- function(object, ...) {
    object@species_params
}

#' @rdname species_params
#' @usage NULL
#' @export
species_params.MizerSim <- function(object, ...) {
    object@params@species_params
}

#' @rdname species_params
#' @usage NULL
#' @export
species_params.data.frame <- function(object, strict = FALSE,
                                      check_misspellings = TRUE, ...) {
    sp <- given_species_params(object, strict = strict,
                               check_misspellings = check_misspellings)
    if ("w_inf" %in% names(sp)) {
        sp <- set_species_param_default(sp, "w_max", 1.5 * sp$w_inf)
        sp <- set_species_param_default(sp, "w_repro_max", sp$w_inf)
        sp <- set_species_param_default(sp, "w_mat", sp$w_inf / 4)
    }
    # Only parameters that no single rate setter owns are defaulted here. A
    # parameter that exactly one `setX()` function reads is defaulted by that
    # function instead, so that each default has a single home. See the
    # "Where defaults live" section of the `default_parameters` vignette.
    sp <- set_species_param_default(sp, "w_min", 0.001)
    sp <- set_species_param_default(sp, "alpha", 0.6)
    sp <- set_species_param_default(sp, "n", 3/4)
    sp <- set_species_param_default(sp, "is_background", FALSE)
    sp <- set_species_param_default(sp, "a", 0.01)
    sp <- set_species_param_default(sp, "b", 3)
    class(sp) <- c("species_params", setdiff(class(sp), c("given_species_params", "species_params")))
    check_and_convert_species_params(sp)
}

#' @rdname species_params
#' @usage NULL
#' @export
species_params.species_params <- function(object, strict = FALSE, ...) {
    species_params.data.frame(object, strict = strict, ...)
}

#' @rdname species_params
#' @export
`species_params<-` <- function(object, recalculate = TRUE, value) {
    UseMethod("species_params<-")
}

#' @rdname species_params
#' @usage NULL
#' @export
`species_params<-.MizerParams` <- function(object, recalculate = TRUE, value) {
    # The columns the user actually handed us. `validSpeciesParams()` below
    # fills in defaults, and a default that the model does not already carry
    # would otherwise look like a column the user has just added and be
    # recorded as a given species parameter, freezing it against every later
    # recalculation.
    supplied_cols <- names(value)
    # Report a likely misspelling only for a column the model does not have
    # yet, so that it is reported when it is introduced and not again on every
    # later edit of an unrelated parameter. Both tables are consulted, the same
    # baseline `given_species_params<-()` uses, so that the two setters cannot
    # disagree about which columns are new.
    check_species_params_misspellings(
        setdiff(supplied_cols, union(names(object@species_params),
                                     names(object@given_species_params))))
    # Columns the model has and the assigned table does not. The user is no
    # longer supplying them, so they are dropped from the given species
    # parameters below: mizer then calculates afresh the ones it knows how to
    # calculate, and the ones it does not know simply go away. `species` cannot
    # be withdrawn; `validSpeciesParams()` rejects a table without it.
    withdrawn <- setdiff(names(object@species_params),
                         c(supplied_cols, "species"))
    # A length parameter whose weight has just been set has to follow the new
    # weight before anything converts the weight away again.
    value <- reconcile_length_weight(value, object@species_params)
    if (recalculate) {
        value <- validSpeciesParams(value, check_misspellings = FALSE)
    } else {
        # Only the checks and conversions that writing into the
        # `@species_params` slot triggers anyway. In particular no default
        # values are filled in: a default that is not already in the slot would
        # look like a change and be recorded as a given species parameter.
        if (!is.species_params(value)) {
            class(value) <- c("species_params",
                              setdiff(class(value),
                                      c("species_params",
                                        "given_species_params")))
        }
        value <- check_and_convert_species_params(value)
    }
    if (!all(value$species == object@species_params$species)) {
        stop("The species names in the new species parameter data frame do not match the species names in the model.")
    }

    # Find what changed compared to old species_params and record it among the
    # given species parameters. Only what the user actually supplied is user
    # input: a default that mizer filled in is not, and neither is the default
    # it may have just put back into a column the user has withdrawn.
    object@given_species_params <-
        record_given_species_params(object@given_species_params,
                                    value[intersect(names(value),
                                                    supplied_cols)],
                                    object@species_params)
    # The withdrawn columns leave the given species parameters, so that the
    # rebuild below no longer regenerates them from there.
    for (col in withdrawn) {
        object@given_species_params <-
            set_column(object@given_species_params, col, NULL)
    }
    if (!recalculate) {
        # Store the supplied parameters as they are. The calculated species
        # parameters are not re-derived and the rates are not recalculated, so
        # it is up to the caller to keep them consistent with the new values.
        object@species_params <- value
    } else {
        changed <- species_param_changes(value, object@species_params)
        if (!needs_species_recalculation(object, value, changed)) {
            object@species_params <- value
            object@time_modified <- lubridate::now()
        } else {
            # Deliberately quiet about the changes that will have no impact,
            # either because another given parameter takes precedence or
            # because the rate array they feed has been set by hand. Only
            # `given_species_params<-()` reports those, so that this setter stays
            # usable in scripts without producing diagnostics they did not ask
            # for. See the note on the two setters in the documentation of
            # `given_species_params<-()`.
            object <- rebuild_from_given(object, value)
        }
    }

    # Report only after the replacement succeeded, and only for columns that
    # validation did not restore to the given species parameters.
    signal_removed_species_params(
        setdiff(withdrawn, names(object@given_species_params)))
    object
}

# Rebuild the species parameters from the given species parameters and
# recalculate all rates. This is the tail that `species_params<-()` and
# `given_species_params<-()` share, so that the two setters cannot drift apart.
#
# `keep` is the species parameter table whose columns are carried over where
# they are not regenerated from the given species parameters: for
# `species_params<-()` the table the user supplied, for
# `given_species_params<-()` the model's current species parameters, which the
# user has not edited. Such columns are for example parameters set directly on
# the `@species_params` slot. The species parameters that a rate setter owns
# are excluded: `setParams()` below derives those afresh, and carrying the
# previously calculated value over would leave it looking like a given value
# and so stop it responding to the change. Where the user did supply such a
# parameter it has been recorded among the given species parameters and so is
# already part of the rebuilt table.
rebuild_from_given <- function(object, keep) {
    # No misspelling check: these columns came from the user through one of the
    # setters, which checked them as they arrived. Checking again here would
    # repeat the report on every rebuild, which is what #581 was.
    new_sp <- validSpeciesParams(object@given_species_params,
                                 check_misspellings = FALSE)
    extra_cols <- setdiff(names(keep),
                          c(names(new_sp), setter_owned_species_params()))
    for (col in extra_cols) {
        new_sp[[col]] <- keep[[col]]
    }
    object@species_params <- new_sp
    return(suppressMessages(setParams(object)))
}

# Which entries changed between two species parameter tables. A removed column
# is represented as changed for every species. This deliberately builds on
# `changed_entries()`, the package's single definition of entry-wise equality.
species_param_changes <- function(value, old) {
    no_sp <- nrow(old)
    cols <- union(names(value), names(old))
    changed <- lapply(cols, function(col) {
        if (!col %in% names(value)) {
            return(rep(TRUE, no_sp))
        }
        changed_entries(value[[col]], old[[col]], no_sp)
    })
    names(changed) <- cols
    changed[vapply(changed, any, logical(1))]
}

# Species parameters for which a numerical change can alter something cached
# in a MizerParams object. Standard columns not on the leaf list default to the
# safe choice of recalculating. Unknown columns on a base MizerParams object are
# leaves, but on an extension object they are conservatively treated as
# dependencies because an extension setter may read them.
recalculation_species_params <- function(params, value) {
    leaf <- c(
        # Read directly by summaries, calibration or plotting code
        "biomass_observed", "biomass_cutoff", "number_observed",
        "number_cutoff", "yield_observed", "catch_observed", "legend_name",
        "is_background",
        # Read directly by reproduction-density-dependence functions
        "erepro", "R_max", "constant_recruitment", "constant_reproduction",
        "ricker_b", "sheperd_b", "sheperd_c"
    )
    required <- setdiff(known_species_params_columns(), leaf)

    # A custom predation kernel can introduce arbitrary species parameters.
    # Its formal arguments are dependencies even when their names are unknown
    # to base mizer.
    types <- unique(c(params@species_params$pred_kernel_type,
                      value$pred_kernel_type))
    types <- types[!is.na(types)]
    for (type in types) {
        fun <- get0(paste0(type, "_pred_kernel"))
        if (is.function(fun)) {
            required <- union(required,
                              setdiff(names(formals(fun)), c("ppmr", "...")))
        }
    }

    if (usesExtensionDispatch(params)) {
        required <- union(required, names(value))
    }
    required
}

# Decide whether changed entries require the expensive rebuild through
# setParams(). For the given-parameter setter, `old_given` distinguishes a
# numerical change from a pure promotion of the current calculated value to an
# explicitly given value.
needs_species_recalculation <- function(params, value, changed,
                                        old_given = NULL) {
    no_sp <- nrow(value)
    # A column the new given table no longer has is one the user has withdrawn.
    # It never shows up among the changes, which are indexed by the columns the
    # new table does have, so it is looked for first: a withdrawal on its own
    # is a change to the model even when nothing else moved.
    removed <- if (is.null(old_given)) {
        character()
    } else {
        setdiff(names(old_given), names(value))
    }
    for (col in removed) {
        # Either the withdrawn value was feeding the model, or the column is
        # still standing in the species parameters and has to be taken out of
        # them too. Only the rebuild does either.
        if (any(explicit_species_param_entries(old_given[[col]], no_sp)) ||
                col %in% names(params@species_params)) {
            return(TRUE)
        }
    }
    if (length(changed) == 0) {
        # A no-op `species_params<-()` call has historically been a way to
        # rebuild an object after package code changed its slots directly.
        # Preserve that repair behaviour. A no-op replacement of the given
        # table, by contrast, changes no model value and needs no rebuild.
        return(is.null(old_given))
    }
    required <- recalculation_species_params(params, value)
    if (is.null(old_given)) {
        # A column the assigned table no longer has is one the user has
        # withdrawn. `species_params<-()` has already taken it out of the given
        # species parameters, and only the rebuild can work out what the model
        # looks like without it.
        if (any(!names(changed) %in% names(value))) {
            return(TRUE)
        }
        return(any(names(changed) %in% required))
    }

    for (col in names(changed)) {
        sel <- changed[[col]]
        new_explicit <- explicit_species_param_entries(value[[col]], no_sp)
        old_explicit <- explicit_species_param_entries(old_given[[col]], no_sp)

        # Handing any explicit value back to calculation needs a rebuild.
        if (any(sel & old_explicit & !new_explicit)) {
            return(TRUE)
        }
        if (!col %in% required) {
            next
        }

        # A newly explicit value that equals the current model value changes
        # provenance only. Every other change to a dependent parameter rebuilds.
        promotion <- sel & !old_explicit & new_explicit
        same_as_current <- !changed_entries(value[[col]],
                                            params@species_params[[col]], no_sp)
        if (any(sel & !(promotion & same_as_current))) {
            return(TRUE)
        }
    }
    FALSE
}

explicit_species_param_entries <- function(x, n) {
    if (is.null(x)) {
        return(rep(FALSE, n))
    }
    entry_has_value(x, n)
}

# Apply explicit entries from a given-parameter table to the full species
# parameter table without deriving anything else. Called only after the change
# classifier has established that no entry was handed back to calculation.
merge_given_species_params <- function(sp, given, changed) {
    no_sp <- nrow(sp)
    for (col in names(changed)) {
        if (col == "species") {
            next
        }
        sel <- changed[[col]] &
            explicit_species_param_entries(given[[col]], no_sp)
        if (!any(sel)) {
            next
        }
        new <- given[[col]]
        old <- sp[[col]]
        if (is.null(old)) {
            sp <- set_column(sp, col, new)
        } else if (!is.null(dim(old)) && length(dim(old)) >= 2) {
            old[sel, ] <- new[sel, , drop = FALSE]
            sp <- set_column(sp, col, old)
        } else {
            old[sel] <- new[sel]
            sp <- set_column(sp, col, old)
        }
    }
    sp
}

#' Record the species parameters that have changed
#'
#' `r lifecycle::badge("experimental")`
#' Compares the new species parameters in `value` against the old ones in
#' `old_sp` and records the entries that have changed in the given species
#' parameter data frame `given`. This is the change detection used by
#' `species_params<-()`, exported so that code which updates the species
#' parameters by other means can record its changes the same way.
#'
#' @details
#' Mizer distinguishes between the species parameters that were given
#' explicitly and those that it calculated itself, see [species_params()]. Only
#' the given ones are protected: whenever the species parameters are
#' recalculated -- which every use of `species_params<-()` triggers -- the
#' calculated ones are derived afresh from the given ones. So code that
#' computes a species parameter and writes it into the `species_params` slot
#' directly has its work silently undone by the next parameter change, unless
#' the value is also recorded among the given species parameters.
#'
#' The usual way to record a parameter is to set it with `species_params<-()`,
#' which also rebuilds the species parameters and recalculates all the rates
#' that depend on them. This function is the recording step on its own, for the
#' case where the caller has already updated the affected rates itself, for
#' example an optimiser that fits a species parameter and the rate array it
#' determines together. Rebuilding and recalculating would then be wasted work,
#' and can even undo the caller's own adjustment.
#'
#' Only the values that have actually changed are recorded. This matters:
#' recording an unchanged value would turn a calculated species parameter into
#' a given one and thereby stop it from responding to changes in the parameters
#' it is derived from. The comparison is made entry by entry, so a parameter is
#' protected only for the species whose value changed. `NA` is compared as a
#' value rather than as an unknown, so `NA` staying `NA` does not count as a
#' change. A column that is not present in `old_sp` at all is taken to be new
#' and is recorded in full.
#'
#' @param given The given species parameter data frame to record into, usually
#'   `given_species_params(params)`.
#' @param value A data frame holding the new species parameters.
#' @param old_sp A data frame holding the species parameters as they were
#'   before the change. Must have one row per species, in the same order as
#'   `value` and `given`.
#'
#' @return The updated `given` data frame.
#' @export
#' @seealso [species_params()], [given_species_params()]
#' @concept helper
#' @examples
#' params <- NS_params
#' sp_before <- species_params(params)
#' given_before <- given_species_params(params)
#'
#' # Set a species parameter and the rate it determines, without going through
#' # `species_params<-()` and its recalculation of every other rate.
#' params@species_params$ks[1] <- species_params(params)$ks[1] * 2
#' params@metab[1, ] <- params@metab[1, ] * 2
#'
#' # Record the change so that it is not recalculated away later
#' params@given_species_params <-
#'     record_given_species_params(given_species_params(params),
#'                                 species_params(params), sp_before)
#'
#' # Only the entry that changed has been recorded
#' given_species_params(params)$ks == given_before$ks
record_given_species_params <- function(given, value, old_sp) {
    assert_that(is.data.frame(given),
                is.data.frame(value),
                is.data.frame(old_sp))
    no_sp <- nrow(old_sp)
    if (nrow(value) != no_sp || nrow(given) != no_sp) {
        stop("`given`, `value` and `old_sp` must all have one row per species.")
    }

    common_cols <- intersect(names(value), names(old_sp))
    for (col in common_cols) {
        new_vals <- value[[col]]
        changed <- changed_entries(new_vals, old_sp[[col]], no_sp)

        if (any(changed)) {
            if (!col %in% names(given)) {
                given[[col]] <- if (is.list(new_vals)) {
                    vector("list", length(new_vals))
                } else {
                    NA
                }
            }
            given[[col]][changed] <- new_vals[changed]
        }
    }
    # A column that `old_sp` does not have at all is new and is recorded in
    # full. Assigned one at a time rather than with `cbind()`, which would
    # return a plain data frame and so strip the `given_species_params` class,
    # stopping `$` from naming the entries after the species.
    new_cols <- setdiff(names(value), names(old_sp))
    for (col in new_cols) {
        given[[col]] <- value[[col]]
    }

    given
}

#' @rdname species_params
#' @return `is.species_params()` returns `TRUE` if `x` is a `species_params`
#'   object, `FALSE` otherwise.
#' @export
is.species_params <- function(x) {
    inherits(x, "species_params")
}

# Recognised species_params column names, used by check_for_misspellings() to
# flag likely typos. This is not an exhaustive list of every possible column
# (users may add custom columns), but covers the standard parameters so that a
# near miss can be detected. Grouped roughly by purpose.
known_species_params_columns <- function() {
    c(# identity and sizes
      "species", "w_max", "w_mat", "w_mat25", "w_min", "w_inf",
      "w_repro_max", "w_min_idx",
      # length-based equivalents and length-weight parameters
      "l_max", "l_mat", "l_mat25", "l_min", "l_inf", "l_repro_max", "a", "b",
      # von Bertalanffy growth
      "k_vb", "t0", "age_mat",
      # physiology
      "h", "k", "ks", "gamma", "alpha", "beta", "sigma",
      "n", "p", "q", "m", "z0", "fc", "f0", "erepro",
      "d", "z_ext", "D_ext", "E_ext",
      # reproduction
      "R_max", "r_max", "constant_recruitment", "constant_reproduction",
      "ricker_b", "sheperd_b", "sheperd_c",
      # predation kernel
      "pred_kernel_type", "kernel_exp", "kernel_l_l", "kernel_u_l",
      "kernel_l_r", "kernel_u_r", "kernel_p", "kernel_mean", "kernel_sd",
      "ppmr_min", "ppmr_max",
      # fishing
      "gear", "sel_func", "catchability", "knife_edge_size",
      "yield_observed", "catch_observed",
      # interactions
      "interaction_resource", "interaction_p",
      # observations
      "biomass_observed", "biomass_cutoff", "number_observed", "number_cutoff",
      # flags and plotting
      "is_background", "linecolour", "linetype", "legend_name")
}

# The species parameters that a rate-setting function derives or defaults
# itself, grouped by the setter that owns each one. See the "Where defaults
# live" section of the `default_parameters` vignette.
#
# `setParams()` fills all of these in again from the given species parameters,
# so code that rebuilds the species parameters must drop the values calculated
# last time rather than carry them over. Carrying a calculated value over would
# make it look like an explicit one and stop it from responding to a change in
# the parameters it is derived from. A value that the user did supply is
# recorded among the given species parameters and so is unaffected.
setter_owned_species_params <- function() {
    c(# setMaxIntakeRate()
      "h",
      # setSearchVolume()
      "gamma", "q",
      # setMetabolicRate()
      "ks", "k", "p",
      # setExtMort()
      "z0", "z_ext", "d",
      # setExtEncounter()
      "E_ext",
      # setExtDiffusion()
      "D_ext",
      # setInteraction()
      "interaction_resource",
      # setPredKernel()
      "pred_kernel_type", "beta", "sigma",
      # setReproduction()
      "w_mat25", "m", "erepro", "R_max")
}

# Familiar abbreviations / capitalisation mistakes that should always be flagged
# even when further than the fuzzy-match threshold from a recognised name.
curated_species_params_misspellings <- function() {
    c("wmin", "wmax", "wmat", "wmat25", "w_mat_25", "Rmax",
      "Species", "Gamma", "Beta", "Sigma", "Alpha",
      "W_min", "W_max", "W_mat", "e_repro", "Age_mat", "w_max_mat")
}

# The misspelling check for species parameter columns. The three constants it
# needs are always the same, and the `var` is what handlers key on, so it is
# named in exactly one place rather than repeated at every call site.
check_species_params_misspellings <- function(cols) {
    check_for_misspellings(cols, known_species_params_columns(),
                           "species parameter", "species_params",
                           curated_species_params_misspellings())
}

# The size parameters that can be given either as a weight or as a length. The
# first entry of each pair is the weight, the second the length it converts
# from.
length_weight_mappings <- function() {
    list(
        c("w_mat", "l_mat"),
        c("w_mat25", "l_mat25"),
        c("w_repro_max", "l_repro_max"),
        c("w_inf", "l_inf"),
        c("w_max", "l_max"),
        c("w_min", "l_min")
    )
}

# Which entries of a species parameter column have changed, comparing entry by
# entry. This is mizer's single definition of "this species parameter changed":
# `record_given_species_params()` uses it to decide what to record and
# `given_species_params<-()` uses it to decide what to report on, so that the
# two cannot drift apart.
#
# `NA` is compared as a value rather than as an unknown: `NA` staying `NA` is
# no change, but a value replaced by `NA` is one, because that is how the user
# hands a parameter back to mizer's calculation. A column that did not exist
# before is compared against `NA`, so adding one that holds only `NA` changes
# nothing.
#
# `new` is the column as it is now and must have one entry per species;
# anything else is taken to say nothing about what changed. `prev` is the
# column as it was, `NULL` where there was none.
changed_entries <- function(new, prev, n) {
    if (is.null(new) || column_length(new) != n) {
        return(rep(FALSE, n))
    }
    if (is.null(prev)) {
        prev <- NA
    } else if (column_length(prev) != n) {
        # Nothing to compare against, so everything counts as changed.
        return(rep(TRUE, n))
    }
    # `==` is only reliable for plain atomic vectors. It is undefined (and
    # errors) for list columns or columns holding S4/other objects, and for a
    # matrix column it returns a matrix rather than one logical per species. In
    # those cases fall back to a per-species `identical()` comparison.
    simple <- is.atomic(new) && is.atomic(prev) &&
        is.null(dim(new)) && is.null(dim(prev))
    if (simple) {
        ch <- !((new == prev) | (is.na(new) & is.na(prev)))
        ch[is.na(ch)] <- TRUE
        return(rep_len(ch, n))
    }
    get_row <- function(x, i) {
        d <- dim(x)
        if (!is.null(d) && length(d) >= 2) {
            x[i, ]
        } else if (length(x) == 1) {
            # The stand-in for a column that did not exist before
            x[[1]]
        } else {
            x[[i]]
        }
    }
    !vapply(seq_len(n),
            function(i) identical(get_row(new, i), get_row(prev, i)),
            logical(1))
}

# The number of species a species parameter column holds values for. A matrix
# column has one row per species, any other column one entry per species.
column_length <- function(x) {
    d <- dim(x)
    if (!is.null(d) && length(d) >= 2) nrow(x) else length(x)
}

# Which entries of a species parameter column hold a value. Only a plain atomic
# column can say that it does not; a list or matrix column is taken to hold a
# value for every species.
entry_has_value <- function(x, n) {
    if (is.atomic(x) && is.null(dim(x)) && length(x) == n) {
        return(!is.na(x))
    }
    rep(TRUE, n)
}

# Set a column without triggering the reactive validation
set_column <- function(x, col, values) {
    saved_class <- class(x)
    class(x) <- "data.frame"
    x[[col]] <- values
    class(x) <- saved_class
    x
}

# Fill in the length-weight parameters that a given species parameter table
# leaves out
#
# `a` and `b` are defaulted rather than given, so a table of given species
# parameters need not contain them even though the model always does. Without
# them none of the length/weight rules can be applied, so the model's values
# are put in where the table does not supply them. The result is only used
# while the length and weight parameters are brought into line and is undone by
# `restore_length_weight_params()` afterwards, so that the values mizer
# defaulted do not end up recorded as given.
#
# Returns `given` unchanged when the model has nothing to contribute.
fill_length_weight_params <- function(given, sp) {
    if (!is.data.frame(sp) || !is.data.frame(given) ||
            nrow(sp) != nrow(given) ||
            !all(c("a", "b") %in% names(sp))) {
        return(given)
    }
    for (col in c("a", "b")) {
        model_value <- sp[[col]]
        given_value <- given[[col]]
        given <- set_column(given, col,
                            if (is.null(given_value) ||
                                    length(given_value) != nrow(given)) {
                                model_value
                            } else {
                                ifelse(is.na(given_value), model_value,
                                       given_value)
                            })
    }
    given
}

# Put the length-weight parameters back the way the caller supplied them,
# undoing `fill_length_weight_params()`. `supplied` is the `a` and `b` columns
# of the table before they were filled in, or no columns at all where the
# caller supplied neither.
restore_length_weight_params <- function(given, supplied) {
    for (col in c("a", "b")) {
        if (col %in% names(supplied)) {
            given <- set_column(given, col, supplied[[col]])
        } else {
            given <- set_column(given, col, NULL)
        }
    }
    given
}

# Bring a length and weight parameter pair into line after a change
#
# A size can be given either as a weight or as the length it converts to. When
# both are present they have to agree, and which of the two determines the
# other follows a simple rule: the one that was given last wins, and if both
# were given at the same time the weight wins.
#
# So where the weight has just changed the length is derived from it, whether
# or not the length changed as well, and where only the length changed the
# weight is derived from the length as it always was.
#
# This is decided when a species parameter data frame is put into a model, by
# `species_params<-()` and `given_species_params<-()`: `old` is then the table
# of the same kind that the model currently holds, against which the incoming
# one is compared. Editing a species parameter data frame on its own changes
# nothing until it is assigned into a model, and a table that was edited by
# itself carries no record of the order in which its columns were changed, so a
# length and a weight that both differ from the model's count as given at the
# same time and the weight wins.
#
# When `old` is `NULL`, or does not line up row by row with `x`, nothing can be
# said about what changed and nothing is done here; the default of letting the
# weight win is then applied by `check_and_convert_species_params()`.
reconcile_length_weight <- function(x, old) {
    if (is.null(old) || !is.data.frame(old) || nrow(old) != nrow(x) ||
            !all(c("a", "b") %in% names(x))) {
        return(x)
    }
    changed <- function(new, prev) changed_entries(new, prev, nrow(x))
    # The conversions are done inline rather than through `l2w()` and `w2l()`,
    # which would repeat their argument checks for every pair. The guard above
    # has already established what those checks would look for.
    a <- x[["a"]]
    b <- x[["b"]]
    for (m in length_weight_mappings()) {
        pw <- m[1]
        pl <- m[2]
        if (!all(c(pw, pl) %in% names(x)) || !all(c(pw, pl) %in% names(old))) {
            next
        }
        w_changed <- changed(x[[pw]], old[[pw]])
        l_changed <- changed(x[[pl]], old[[pl]])
        # The weight was just set, so the length follows it. This also covers
        # the case where both were set at the same time, where the weight wins.
        to_length <- w_changed
        # Only the length was set, so it determines the weight.
        to_weight <- l_changed & !w_changed
        if (any(to_length)) {
            vl <- (x[[pw]] / a)^(1 / b)
            sel <- to_length & !is.na(vl)
            if (any(sel)) {
                new_l <- x[[pl]]
                new_l[sel] <- vl[sel]
                x <- set_column(x, pl, new_l)
            }
        }
        if (any(to_weight)) {
            vw <- a * x[[pl]]^b
            sel <- to_weight & !is.na(vw)
            if (any(sel)) {
                new_w <- x[[pw]]
                new_w[sel] <- vw[sel]
                x <- set_column(x, pw, new_w)
            }
        }
    }
    x
}

check_f0 <- function(f0, species = NULL) {
    if (is.logical(f0) && all(is.na(f0))) {
        f0 <- as.numeric(f0)
    }
    if (!is.numeric(f0)) {
        stop("The species parameter `f0` must be numeric.")
    }
    invalid_f0 <- is.nan(f0) |
        (!is.na(f0) & (!is.finite(f0) | f0 < 0 | f0 >= 1))
    if (any(invalid_f0)) {
        msg <- paste0("The species parameter `f0` must be finite and in the ",
                      "interval [0, 1).")
        if (!is.null(species)) {
            msg <- paste0(msg, " Invalid values for: ",
                          paste(species[invalid_f0], collapse = ", "), ".")
        }
        stop(msg)
    }
    invisible(NULL)
}

check_and_convert_species_params <- function(x) {
    if ("f0" %in% names(x)) {
        check_f0(x[["f0"]], x[["species"]])
    }

    # Auto convert length to weight if allometric parameters exist
    if (all(c("a", "b") %in% names(x))) {
        # Converted inline rather than through `l2w()` and `w2l()`, whose
        # argument checks would otherwise be repeated for each of the six
        # pairs. The guard above covers what they would check for.
        a <- x[["a"]]
        b <- x[["b"]]
        for (m in length_weight_mappings()) {
            pw <- m[1]
            pl <- m[2]
            if (pl %in% names(x)) {
                # The weight that the length converts to
                vw <- a * x[[pl]]^b
                if (!(pw %in% names(x))) {
                    x <- set_column(x, pw, vw)
                    next
                }
                # A weight that is not known is taken from its length. This is
                # done per species, so that one species with no length no
                # longer loses the weight it does have.
                fill <- !is.na(vw) & is.na(x[[pw]])
                if (any(fill)) {
                    new_w <- x[[pw]]
                    new_w[fill] <- vw[fill]
                    x <- set_column(x, pw, new_w)
                }
                # Anything still disagreeing was not resolved above, so
                # neither value was just set, or there was nothing to compare
                # against. Both were then given at the same time, as far as we
                # can tell, and the weight wins: the length is brought into
                # line with it. Say so, because the length the caller gave is
                # being changed.
                disagree <- !is.na(vw) & !is.na(x[[pw]]) &
                    abs(x[[pw]] - vw) > 1e-10
                if (any(disagree)) {
                    vl <- (x[[pw]] / a)^(1 / b)
                    sel <- disagree & !is.na(vl)
                    if (any(sel)) {
                        signal_info(pl, paste0(
                            "For the following species the value of `", pl,
                            "` is not consistent with the value of `", pw,
                            "`, so I am using `", pw, "` and setting `", pl,
                            "` to match it: ",
                            paste(x$species[sel], collapse = ", ")),
                            level = 1, severity = "warning",
                            unhandled = "show")
                        new_l <- x[[pl]]
                        new_l[sel] <- vl[sel]
                        x <- set_column(x, pl, new_l)
                    }
                }
            }
        }
    }

    # Check w_mat < w_inf consistency
    if (all(c("w_mat", "w_inf") %in% names(x))) {
        wrong <- !is.na(x$w_mat) & !is.na(x$w_inf) & x$w_mat >= x$w_inf
        if (any(wrong)) {
            signal_info("w_mat", paste0(
                "For the species ",
                paste(x$species[wrong], collapse = ", "),
                " the value for `w_mat` is not smaller than that of `w_inf`."),
                level = 1, severity = "warning", unhandled = "show")
        }
    }

    x
}

#' @export
`$.species_params` <- function(x, name) {
    out <- exact_column(x, name, "species_params")
    if (!is.null(out) && !is.data.frame(out) && name != "species") {
        names(out) <- rownames(x)
    }
    out
}


#' @export
`[.species_params` <- function(x, i, j, ..., drop = FALSE) {
    out <- NextMethod("[")
    if (is.data.frame(out)) {
        class(out) <- class(x)
    }
    out
}

#' @export
print.species_params <- function(x, ...) {
    cat("An object of class \"", class(x)[1], "\" containing parameters for ", nrow(x), " species:\n", sep = "")
    core_cols <- c("species", "w_inf", "w_mat", "h", "ks", "z0", "z_ext")
    cols_to_show <- intersect(core_cols, names(x))
    extra_cols <- setdiff(names(x), core_cols)

    if (length(cols_to_show) < length(core_cols) && length(extra_cols) > 0) {
        num_to_add <- min(length(core_cols) - length(cols_to_show), length(extra_cols))
        cols_to_show <- c(cols_to_show, extra_cols[1:num_to_add])
        extra_cols <- extra_cols[-(1:num_to_add)]
    }

    print(as.data.frame(x)[, cols_to_show, drop = FALSE], row.names = FALSE, ...)

    if (length(extra_cols) > 0) {
        cat("With", length(extra_cols), "other parameters:", paste(extra_cols, collapse = ", "), "\n")
    }
    invisible(x)
}

#' @export
summary.species_params <- function(object, ...) {
    cat("Summary of species_params:\n")
    cat("Number of species:", nrow(object), "\n")
    num_cols <- names(object)[vapply(object, is.numeric, logical(1))]
    cat("Parameter ranges:\n")
    for (col in num_cols) {
        vals <- object[[col]]
        vals <- vals[is.finite(vals)]
        if (length(vals) > 0) {
            cat("  ", col, ": min = ", min(vals), ", max = ", max(vals), "\n", sep = "")
        }
    }
    invisible(object)
}


#' @rdname species_params
#' @export
given_species_params <- function(object, ...) {
    UseMethod("given_species_params")
}

#' @rdname species_params
#' @usage NULL
#' @export
given_species_params.MizerParams <- function(object, ...) {
    object@given_species_params
}

#' @rdname species_params
#' @usage NULL
#' @export
given_species_params.MizerSim <- function(object, ...) {
    object@params@given_species_params
}

#' @rdname species_params
#' @usage NULL
#' @export
given_species_params.data.frame <- function(object, strict = FALSE,
                                            check_misspellings = TRUE, ...) {
    assert_that(is.data.frame(object))
    # Convert a tibble back to an ordinary data frame
    sp <- as.data.frame(object, stringsAsFactors = FALSE)

    if (check_misspellings) {
        check_species_params_misspellings(names(sp))
    }

    # check species
    if (!("species" %in% colnames(sp))) {
        stop("The species params dataframe needs a column 'species' with the species names")
    }
    sp$species <- as.character(sp$species)
    species_names <- as.character(sp$species)
    no_sp <- nrow(sp)
    if (length(unique(species_names)) != no_sp) {
        stop("The species parameter data frame has multiple rows for the same species")
    }
    sp$species <- species_names
    row.names(sp) <- species_names

    # Allow r_max instead of R_max
    if (!("R_max" %in% names(sp)) && "r_max" %in% names(sp)) {
        names(sp)[names(sp) == "r_max"] <- "R_max"
    }

    # Convert lengths to weights
    if (all(c("a", "b") %in% names(sp))) {
        sp <- sp %>%
            set_species_param_from_length("w_mat", "l_mat") %>%
            set_species_param_from_length("w_mat25", "l_mat25") %>%
            set_species_param_from_length("w_repro_max", "l_repro_max") %>%
            set_species_param_from_length("w_inf", "l_inf") %>%
            set_species_param_from_length("w_max", "l_max") %>%
            set_species_param_from_length("w_min", "l_min")
    }

    # check w_inf
    if (!("w_inf" %in% names(sp))) {
        if ("w_repro_max" %in% names(sp)) {
            sp$w_inf <- sp$w_repro_max
            signal_info("w_inf", "The species parameter data frame is missing a `w_inf` column. I am using the values from the `w_repro_max` column instead.",
                        level = 1)
        } else if ("w_max" %in% names(sp)) {
            sp$w_inf <- sp$w_max
            signal_info("w_inf", "The species parameter data frame is missing a `w_inf` column. I am using the values from the `w_max` column instead. ",
                        level = 1)
        } else if (strict) {
            stop("You need to specify the asymptotic size `w_inf` for all species.")
        }
    }
    if ("w_inf" %in% names(sp)) {
        missing <- is.na(sp$w_inf)
        if (any(missing) && strict) {
            stop("You need to specify the asymptotic size `w_inf` for all species.")
        }
        if (!is.numeric(sp$w_inf) && strict) {
            stop("`w_inf` contains non-numeric values.")
        }
    }

    # check w_mat
    if ("w_mat" %in% names(sp) && "w_inf" %in% names(sp)) {
        wrong <- !is.na(sp$w_mat) & !is.na(sp$w_inf) & sp$w_mat >= sp$w_inf
        if (any(wrong)) {
            signal_info("w_mat", paste0(
                "For the species ",
                paste(sp$species[wrong], collapse = ", "),
                " the value for `w_mat` is not smaller than that of `w_inf`.",
                " I have corrected that by setting it to 25% of `w_inf`."),
                level = 1, severity = "warning", unhandled = "show")
            sp$w_mat[wrong] <- sp$w_inf[wrong] / 4
        }

        # check w_mat25
        if ("w_mat25" %in% names(sp)) {
            wrong <- !is.na(sp$w_mat) & !is.na(sp$w_mat25) &
                sp$w_mat25 >= sp$w_mat
            if (any(wrong)) {
                msg <- paste0("For the species ",
                              paste(sp$species[wrong], collapse = ", "),
                              " the value for `w_mat25` is not smaller than ",
                              "that of `w_mat`. I have corrected that by ",
                              "marking it as missing so that its default will ",
                              "be used.")
                signal_info("w_mat25", msg, level = 1, severity = "warning",
                            unhandled = "show")
                sp$w_mat25[wrong] <- NA
                # The length has to go with it. Otherwise
                # `check_and_convert_species_params()` below fills the missing
                # weight back in from `l_mat25` and the rejected value returns,
                # leaving the message above untrue. What comes back is `w_mat`
                # to within rounding error, small enough to pass the
                # `w_mat25 < w_mat` assertion in `setReproduction()` and large
                # enough to collapse the maturity ogive to a step function.
                if ("l_mat25" %in% names(sp)) {
                    sp$l_mat25[wrong] <- NA
                }
            }
        }

        # check w_min
        if ("w_min" %in% names(sp)) {
            wrong <- !is.na(sp$w_min) & !is.na(sp$w_mat) & sp$w_min >= sp$w_mat
            if (any(wrong)) {
                sp$w_min[wrong] <- pmin(0.001, sp$w_mat[wrong] / 10)
                signal_info("w_min", paste0(
                    "For the species ",
                    paste(sp$species[wrong], collapse = ", "),
                    " the value for `w_min` is not smaller than that of `w_mat`.",
                    " I have reduced the values."),
                    level = 1, severity = "warning", unhandled = "show")
            }
        }
    }

    # check w_repro_max
    if ("w_repro_max" %in% names(sp) && "w_mat" %in% names(sp)) {
        wrong <- !is.na(sp$w_repro_max) & !is.na(sp$w_mat) & sp$w_repro_max <= sp$w_mat
        if (any(wrong)) {
            signal_info("w_repro_max", paste0(
                "For the species ",
                paste(sp$species[wrong], collapse = ", "),
                " the value for `w_repro_max` is smaller than that of `w_mat`.",
                " I have corrected that by setting it to 4 times `w_mat`."),
                level = 1, severity = "warning", unhandled = "show")
            sp$w_repro_max[wrong] <- 4 * sp$w_mat[wrong]
        }
    }

    class(sp) <- c("given_species_params", "species_params", setdiff(class(sp), c("given_species_params", "species_params")))
    check_and_convert_species_params(sp)
}

#' @rdname species_params
#' @usage NULL
#' @export
given_species_params.given_species_params <- function(object, strict = FALSE, ...) {
    given_species_params.data.frame(object, strict = strict, ...)
}

#' @rdname species_params
#' @return `is.given_species_params()` returns `TRUE` if `x` is a
#'   `given_species_params` object, `FALSE` otherwise.
#' @export
is.given_species_params <- function(x) {
    inherits(x, "given_species_params")
}

#' @rdname species_params
#' @export
`given_species_params<-` <- function(object, value) {
    UseMethod("given_species_params<-")
}

#' @rdname species_params
#' @usage NULL
#' @export
`given_species_params<-.MizerParams` <- function(object, value) {
    params <- object
    # Worked out here, before the length/weight rules below rewrite `value`.
    # Only the columns the model does not have yet, so that a likely
    # misspelling is reported when it is introduced and not again on every
    # later edit of an unrelated parameter. Reported further down, together
    # with the other diagnostics this setter gives.
    new_cols <- setdiff(names(value),
                        union(names(params@species_params),
                              names(params@given_species_params)))
    # The length/weight rules need `a` and `b`, which the given species
    # parameters need not contain, so the model's values stand in while they
    # are applied.
    supplied_ab <- if (is.data.frame(value)) {
        value[intersect(c("a", "b"), names(value))]
    } else {
        # `validGivenSpeciesParams()` below rejects this with its own message
        list()
    }
    value <- fill_length_weight_params(value, params@species_params)
    # A length parameter whose weight has just been set has to follow the new
    # weight before anything converts the weight away again. This is the same
    # rule that `species_params<-()` applies, here comparing the incoming given
    # species parameters against the model's.
    value <- reconcile_length_weight(value, params@given_species_params)
    value <- validGivenSpeciesParams(value, check_misspellings = FALSE)
    value <- restore_length_weight_params(value, supplied_ab)
    if (!all(value$species == params@species_params$species)) {
        stop("The species names in the new species parameter data frame do not match the species names in the model.")
    }

    # Which species changed in which column, for the reports below. All three
    # reports work from this one list, so they cannot disagree about what
    # counts as a change. A column that the given species parameters do not
    # have yet is compared against `NA`, so that adding an all-`NA` column is
    # no change while clearing a given value to `NA` is one: the user is then
    # handing that parameter back to mizer's calculation.
    old_given <- params@given_species_params
    # Columns the user has withdrawn. They never appear among the changes
    # below, which are indexed by the columns the new table has.
    withdrawn <- setdiff(names(old_given), names(value))
    no_sp <- nrow(value)
    changed <- lapply(names(value), function(col) {
        changed_entries(value[[col]], old_given[[col]], no_sp)
    })
    names(changed) <- names(value)
    changed <- changed[vapply(changed, any, logical(1))]

    # Of those changes, the ones that gave a species a value. Only a value that
    # is there can be overruled by another parameter, so this is the narrower
    # question `signal_ignored_changes()` asks: clearing a value to `NA` is a
    # change, but not one it has anything to say about.
    specified <- lapply(names(changed), function(col) {
        changed[[col]] & entry_has_value(value[[col]], no_sp)
    })
    names(specified) <- names(changed)

    # Warn about the changes that will have no impact, either because another
    # given parameter takes precedence or because the rate array they feed has
    # been set by hand. A withdrawn column is a change too, although it is not
    # present in the table that `changed` is built from.
    changed_names <- union(names(changed), withdrawn)
    with_info_level({
        check_species_params_misspellings(new_cols)
        signal_ignored_changes(old_given, specified)
        signal_gear_params_changes(changed_names)
        signal_frozen_changes(params, changed_names)
    })

    params@given_species_params <- value
    if (!needs_species_recalculation(params, value, changed,
                                     old_given = old_given)) {
        sp <- params@species_params
        params@species_params <- merge_given_species_params(sp, value, changed)
        params@time_modified <- lubridate::now()
    } else {
        # The user has edited the given species parameters and not the full
        # table, so it is the model's own species parameters whose columns are
        # carried over where they are not rebuilt from the given ones. Not the
        # withdrawn ones though: `validSpeciesParams()` and the rate setters put
        # back the ones mizer knows how to calculate, and one it cannot
        # calculate would otherwise survive here and be reported by
        # `calculated_species_params()` as a value mizer had produced.
        keep <- params@species_params
        for (col in withdrawn) {
            keep <- set_column(keep, col, NULL)
        }
        params <- rebuild_from_given(params, keep)
    }
    # Report only after the replacement succeeded, and only for columns that
    # validation did not restore to the given species parameters.
    signal_removed_species_params(
        setdiff(withdrawn, names(params@given_species_params)))
    params
}

#' @rdname species_params
#' @export
calculated_species_params <- function(params) {
    # Identifying common columns
    common_cols <- intersect(names(params@species_params),
                             names(params@given_species_params))
    # Copy df1 to new_df
    calculated <- params@species_params
    # remove the entries that are also in given_species_params
    for (col in common_cols) {
        calculated[[col]] <- replace_with_na(calculated[[col]],
                                             params@given_species_params[[col]])
    }
    # Removing columns that only contain NAs
    calculated <- calculated %>%
        select(where(~ !all(is.na(.))))

    calculated$species <- params@species_params$species
    calculated <- calculated[, c("species", setdiff(names(calculated), "species")), drop = FALSE]

    return(calculated)
}

# Function to replace overlapping entries with NA
replace_with_na <- function(x, y) {
    ifelse(is.na(y), x, NA)
}

#' Set a species parameter to a default value
#'
#' If the species parameter does not yet exist in the species parameter data
#' frame, then create it and fill it with the default. Otherwise use the default
#' only to fill in any NAs. Optionally gives a message if the parameter
#' did not already exist. The signal has class `info_about_default`.
#' @param object Either a MizerParams object or a species parameter data frame
#' @param parname A string with the name of the species parameter to set
#' @param default A single default value or a vector with one default value for
#'   each species
#' @param message A string with a message to be issued when the parameter did
#'   not already exist
#' @return The `object` with an updated column in the species params data frame.
#' @export
#' @concept helper
set_species_param_default <- function(object, parname, default,
                                      message = NULL) {
    if (is(object, "MizerParams")) {
        species_params <- object@species_params
    } else {
        species_params <- object
    }
    assert_that(is.data.frame(species_params))
    assert_that(is.string(parname))
    no_sp <- nrow(species_params)
    if (length(default) == 1) {
        default <- rep(default, no_sp)
    }
    assert_that(length(default) == no_sp)
    if (!(parname %in% colnames(species_params))) {
        if (!missing(message)) {
            signal_info(parname, message)
        }
        species_params[parname] <- default
    } else {
        # We do not like factors
        if (is.factor(species_params[[parname]])) {
            species_params[[parname]] <- as.character(species_params[[parname]])
        }
        missing <- is.na(species_params[[parname]])
        if (any(missing)) {
            species_params[missing, parname] <- default[missing]
        }
    }
    if  (is(object, "MizerParams")) {
        object@species_params <- species_params
        return(object)
    } else {
        return(species_params)
    }
}



#' Get default value for h
#'
#' Sets `h` so that the species reaches maturity size `w_mat` at the maturity
#' age `age_mat` if it feeds at feeding level `f0`.
#'
#' If `age_mat` is missing in the species parameter data frame, then it is
#' calculated from the von Bertalanffy growth curve parameters `k_vb` and
#' (optionally `t0`) taken from the species parameter data frame. This is not
#' reliable and a warning is issued.
#'
#' If no growth information is given at all for a species, the default is set
#' to `h = 30`.
#'
#' See the [Maximum Intake Rate Coefficient](
#' https://sizespectrum.org/mizer/articles/default_parameters.html#h-default)
#' section of the "Calculation of Default Parameter Values" vignette for the
#' mathematical derivation.
#'
#' @param params A MizerParams object or a species parameter data frame
#' @return A vector with the values of h for all species
#' @export
#' @concept helper
#' @family functions calculating defaults
get_h_default <- function(params) {
    if (is(params, "MizerParams")) {
        species_params <- params@species_params
    } else {
        species_params <- validSpeciesParams(params)
    }
    assert_that("n" %in% names(species_params))
    species_params <- set_species_param_default(species_params, "f0", 0.6)
    if (!("h" %in% colnames(species_params))) {
        species_params[["h"]] <- rep(NA, nrow(species_params))
    }
    missing <- is.na(species_params[["h"]])
    if (any(missing)) {
        # The following should be assured by `validSpeciesParams()`
        assert_that(is.numeric(species_params$f0),
                    noNA(species_params$alpha),
                    "alpha" %in% names(species_params))
        signal_info("h", "No h provided for some species, so using age at maturity to calculate it.")
        if (!isTRUE(all.equal(species_params$n[missing], species_params$p[missing],
                              check.attributes = FALSE))) {
            signal_info("h", "Because you have n != p, the default value for `h` is not very good.",
                        level = 1)
        }
        species_params <- species_params %>%
            set_species_param_default("fc", 0.2) %>%
            set_species_param_default(
                "age_mat", age_mat_vB(species_params),
                strwrap("Because the age at maturity is not known, I need to
                        fall back to using von Bertalanffy parameters, where
                        available.")
            )
        w_mat <- species_params$w_mat
        w_min <- species_params$w_min
        age_mat <- species_params$age_mat
        n <- species_params[["n"]]
        h <- (w_mat^(1 - n) - w_min^(1 - n)) / age_mat / (1 - n) /
            species_params$alpha / (species_params$f0 - species_params$fc)

        species_params[missing, "h"] <- h[missing]

        # If no acceptable default could be calculated, set h=30
        missing <- is.na(species_params[["h"]]) | species_params[["h"]] <= 0
        if (any(missing)) {
            signal_info("h", "For species where no growth information is available the parameter h has been set to h = 30.")
            species_params[missing, "h"] <- 30
        }
    }
    return(species_params[["h"]])
}


#' Measure the available energy in a power-law prey spectrum
#'
#' Puts the model into the reference state used for deriving the default
#' `gamma` and `f0`: no fish, and a resource spectrum given by the power law
#' \eqn{\kappa w^{-\lambda}}. It then measures the resulting encounter rate at
#' the largest size and divides out its power-law dependence on size, so that
#' the returned number for each species is the coefficient \eqn{A_i} in
#' \eqn{E_i(w) = A_i\, w^{2 + q_i - \lambda}}.
#'
#' The measurement uses [mizerEncounter()] rather than [getEncounter()]. What is
#' wanted here is a property of the species parameters, not of the model's
#' dynamics, so it must not go through the extension dispatch chain or through a
#' rate function that the user has registered with [setRateFunction()]. An
#' extension that scales the search volume would otherwise fold its own factor
#' into mizer's `gamma` default, and because that default then determines the
#' search volume, the factor would be re-applied on every rebuild (issue #577).
#'
#' The caller is responsible for having set the `q` and `gamma` species
#' parameters and the `search_vol` slot that the measurement is to use.
#'
#' @param params A MizerParams object
#' @return A vector with one number for each species
#' @noRd
measure_avail_energy <- function(params) {
    params@initial_n[] <- 0
    if (defaults_edition() < 2) {
        # See issue #238
        params@species_params$interaction_resource <- 1
    }
    params@initial_n_pp[] <- resource_power_law(
        params, params@resource_params$kappa,
        params@resource_params$lambda)
    encounter <- mizerEncounter(params, n = params@initial_n,
                                n_pp = params@initial_n_pp,
                                n_other = params@initial_n_other, t = 0)
    encounter[, length(params@w)] /
        params@w[length(params@w)] ^
        (2 + params@species_params[["q"]] - params@resource_params$lambda)
}

#' Get default value for gamma
#'
#' Fills in any missing values for gamma so that fish feeding on a resource
#' spectrum described by the power law \eqn{\kappa w^{-\lambda}} achieve a
#' feeding level \eqn{f_0}. Only for internal use.
#'
#' See the [Search Volume Coefficient](
#' https://sizespectrum.org/mizer/articles/default_parameters.html#gamma-default)
#' section of the "Calculation of Default Parameter Values" vignette for the
#' mathematical derivation.
#'
#' The available energy is measured with mizer's own encounter rate,
#' [mizerEncounter()]. It is a property of the species parameters and of the
#' resource power law rather than of the model's dynamics, so the measurement
#' deliberately does not go through the extension dispatch chain, nor through an
#' encounter function registered with [setRateFunction()].
#'
#' @param params A MizerParams object
#' @return A vector with the values of gamma for all species
#' @export
#' @concept helper
#' @family functions calculating defaults
get_gamma_default <- function(params) {
    assert_that(is(params, "MizerParams"))
    species_params <- params@species_params %>%
        set_species_param_default("f0", 0.6)
    if (!("gamma" %in% colnames(species_params))) {
        species_params$gamma <- rep(NA, nrow(species_params))
    }
    check_f0(species_params$f0, species_params$species)
    missing <- is.na(species_params$gamma)
    if (any(missing)) {
        assert_that(is.number(params@resource_params$lambda),
                    is.number(params@resource_params$kappa),
                    is.numeric(species_params$f0))
        signal_info("gamma", "Using f0, h, lambda, kappa and the predation kernel to calculate gamma.")
        if (!"h" %in% names(params@species_params) ||
            any(is.na(species_params[["h"]]))) {
            species_params[["h"]] <- get_h_default(params)
        }
        # Calculate available energy by setting search_volume
        # coefficient to 1. We write the search volume into the slot directly
        # instead of going through `setSearchVolume()`, because that setter
        # refuses to recalculate a search volume that the user has set by hand
        # and would leave us measuring the available energy with the user's
        # array instead of the unit-gamma one (issue #488).
        params@species_params$gamma <- 1
        params <- set_q_default(params)
        params@search_vol[] <- compute_search_vol(params)
        # and setting a power-law prey spectrum
        avail_energy <- measure_avail_energy(params)
        # Now set gamma so that this available energy leads to f0
        if (any(is.infinite(species_params[["h"]][missing]))) {
            stop("Cannot calculate default `gamma` for species with `h = Inf`. Please supply `gamma` explicitly.")
        }
        gamma_default <- (species_params[["h"]] / avail_energy) *
            (species_params$f0 / (1 - species_params$f0))
        # Only overwrite missing gammas with calculated values
        if (any(!is.finite(gamma_default[missing]))) {
            stop("Could not calculate gamma.")
        }
        species_params$gamma[missing] <- gamma_default[missing]
    }
    return(species_params$gamma)
}

#' Get default value for f0
#'
#' Fills in any missing values for f0 so that if the prey abundance was
#' described by the power law \eqn{\kappa w^{-\lambda}} then the encounter rate
#' coming from the given `gamma` parameter would lead to the feeding level
#' \eqn{f_0}. This is thus doing the inverse of [get_gamma_default()].
#' Only for internal use.
#'
#' For species for which no value for `gamma` is specified in the species
#' parameter data frame, the `f0` values is kept as provided in the species
#' parameter data frame or it is set to 0.6 if it is not provided.
#'
#' See the [Target Feeding Level](
#' https://sizespectrum.org/mizer/articles/default_parameters.html#f0-default)
#' section of the "Calculation of Default Parameter Values" vignette for the
#' mathematical derivation.
#'
#' The available energy is measured with mizer's own encounter rate,
#' [mizerEncounter()]. It is a property of the species parameters and of the
#' resource power law rather than of the model's dynamics, so the measurement
#' deliberately does not go through the extension dispatch chain, nor through an
#' encounter function registered with [setRateFunction()].
#'
#' @param params A MizerParams object
#' @return A vector with the values of f0 for all species
#' @export
#' @concept helper
#' @family functions calculating defaults
get_f0_default <- function(params) {
    assert_that(is(params, "MizerParams"))
    species_params <- params@species_params %>%
        set_species_param_default("f0", 0.6)
    if (!("gamma" %in% colnames(species_params))) {
        species_params$gamma <- rep(NA, nrow(species_params))
    }
    given <- !is.na(species_params$gamma)
    if (any(given)) {
        assert_that(is.number(params@resource_params$lambda),
                    is.number(params@resource_params$kappa),
                    is.numeric(species_params$gamma))
        if (!"h" %in% names(params@species_params) ||
            any(is.na(species_params[["h"]]))) {
            species_params[["h"]] <- get_h_default(params)
        }
        # This is the inverse of `get_gamma_default()` and so has to measure
        # the available energy with the search volume implied by the given
        # `gamma`, not with whatever is in the slot, which the user may have
        # set by hand (issue #488). Only the species with a given `gamma` get
        # their search volume replaced, the others keep theirs so that no NAs
        # enter the calculation.
        params <- set_q_default(params)
        params@search_vol[given, ] <- compute_search_vol(params)[given, ]
        # Calculate available energy by setting a power-law prey spectrum
        avail_energy <- measure_avail_energy(params)
        # Now set f0 so that this available energy leads to f0
        f0_default <- 1 / (species_params[["h"]] / avail_energy + 1)
        if (any(is.na(f0_default[given]))) {
            stop("Could not calculate f0.")
        }
        # Only overwrite f0 for species where gamma was given
        species_params$f0[given] <- f0_default[given]
    }
    return(species_params$f0)
}

#' Get default value for `ks`
#'
#' Fills in any missing values for `ks` so that the critical feeding level needed
#' to sustain the species is as specified in the `fc` column in the species
#' parameter data frame. If that column is not provided the default critical
#' feeding level \eqn{f_c = 0.2} is used.
#'
#' See the [Standard Metabolic Rate Coefficient](
#' https://sizespectrum.org/mizer/articles/default_parameters.html#ks-default)
#' section of the "Calculation of Default Parameter Values" vignette for the
#' mathematical derivation.
#'
#' @param params A MizerParams object
#' @return A vector with the values of ks for all species
#' @export
#' @concept helper
#' @family functions calculating defaults
get_ks_default <- function(params) {
    assert_that(is(params, "MizerParams"),
                "n" %in% names(params@species_params),
                "p" %in% names(params@species_params))
    if (!"h" %in% names(params@species_params) ||
        any(is.na(params@species_params[["h"]]))) {
        params@species_params[["h"]] <- get_h_default(params)
    }
    params <- set_species_param_default(params, "fc", 0.2)
    sp <- params@species_params
    if (!"ks" %in% names(sp)) {
        sp$ks <- rep(NA_real_, nrow(sp))
    }
    missing_ks <- is.na(sp$ks)
    if (any(missing_ks) && any(is.infinite(sp[["h"]][missing_ks]))) {
        stop("Cannot calculate default `ks` for species with `h = Inf`. Please supply `ks` explicitly.")
    }
    ks_default <- sp$fc * sp$alpha * sp[["h"]] * sp$w_mat^(sp[["n"]] - sp[["p"]])

    message <- ("No ks column so calculating from critical feeding level.")
    sp <- set_species_param_default(sp, "ks", ks_default, message)
    if (any(is.na(sp$ks) |  is.infinite(sp$ks))) {
        stop("Could not calculate default values for the missing species ",
             "parameter ks. Got: ", sp$ks)
    }
    return(sp$ks)
}
