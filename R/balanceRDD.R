#' Balance density-dependent reproduction at the initial state
#'
#' Changes the reproduction parameters so that, at the initial state, the
#' density-dependent reproduction rate equals the rate that is required to keep
#' the egg abundance constant. A model whose spectra are at steady state is
#' then also at steady state in its reproduction. The density-dependence
#' function itself is not changed. Mizer calls this function whenever it has
#' changed a model in a way that needs the reproduction to be rebalanced, for
#' example in [matchBiomasses()], [matchGrowth()], [tuneSteadyState()] and
#' [addSpecies()].
#'
#' # Balancing functions
#'
#' `balanceRDD()` looks for a function called `balance_<RDD>()`, where
#' `<RDD>` is the name of the density-dependence function of the model, the
#' name you give to `setReproduction(params, RDD = ...)` or to
#' `setRateFunction(params, "RDD", ...)`. It is looked up in the same way as
#' the density-dependence function itself, so a balancing function that is
#' defined next to its density-dependence function is found.
#'
#' Mizer provides balancing functions for its own density-dependence
#' functions. By default they keep `erepro` and change the parameter of the
#' density dependence:
#'
#' * `balance_BevertonHoltRDD()` changes `R_max`. It also accepts one of the
#'   arguments `erepro`, `R_max` or `reproduction_level`, which choose the
#'   Beverton-Holt curve as described in [setBevertonHolt()].
#' * `balance_RickerRDD()` changes `ricker_b`. It also accepts `erepro`.
#' * `balance_SheperdRDD()` changes `sheperd_b` and keeps `sheperd_c`. It also
#'   accepts `erepro`.
#' * `balance_constantRDD()` changes `constant_reproduction`.
#' * `balance_noRDD()` changes `erepro`, because [noRDD()] has no parameter of
#'   its own.
#'
#' When `erepro` cannot be kept because the density-independent reproduction
#' rate is smaller than the required rate, `erepro` is increased to the
#' smallest value that works, with a warning.
#'
#' # Density-dependence functions without a balancing function
#'
#' If there is no balancing function, mizer keeps the parameters of the
#' density dependence and solves numerically for the values of `erepro` that
#' give the required reproduction rate, and says that it has done so. This
#' works for any density dependence in which the reproduction rate of a species
#' responds to the egg production, including density dependence that couples
#' several species. Only `erepro` can be adjusted in this way, so
#' arguments in `...` are then not allowed.
#'
#' If no value of `erepro` gives the required rate, for example because the
#' density dependence ignores the egg production, the parameters are left
#' unchanged with a warning. The same happens when the density dependence
#' depends on time, which mizer checks by evaluating it at a few different
#' times: such a model has no steady state to balance it to.
#'
#' # Writing a balancing function
#'
#' To control how your own density-dependence function `myRDD()` is balanced,
#' write a function `balance_myRDD()` with the arguments `params`, `rdi` and
#' `rdd_target`, followed by any arguments that choose which parameters are
#' changed. It should return `params` with the parameters changed so that
#' `myRDD()` gives `rdd_target` when the density-independent reproduction rate
#' is `rdi`, scaled by any change you make to `erepro`. It must leave alone the
#' species for which `rdd_target` is `NA`. It may assume that the
#' density-independent reproduction rate is proportional to `erepro`, as it is
#' with [mizerRDI()].
#'
#' # Extension packages
#'
#' An extension package that changes the reproduction rate through a
#' [projectRDD()] method should also provide a `balanceRDD()` method. That
#' method undoes the change to the target and passes it on to the next method,
#' which in the end balances the density-dependence function named in the
#' model. For an extension that adds immigration to the recruitment:
#'
#' ```r
#' projectRDD.myExtension <- function(params, rdi, ...) {
#'     NextMethod() + other_params(params)$immigration
#' }
#' balanceRDD.myExtension <- function(params, rdi = getRDI(params),
#'                                    rdd_target = getRequiredRDD(params),
#'                                    reproduction_level = NULL, ...) {
#'     immigration <- other_params(params)$immigration
#'     if (!is.null(reproduction_level)) {
#'         reproduction_level <- 1 - (1 - reproduction_level) *
#'             rdd_target / (rdd_target - immigration)
#'     }
#'     NextMethod(rdd_target = rdd_target - immigration,
#'                reproduction_level = reproduction_level)
#' }
#' ```
#'
#' A `reproduction_level` is the reproduction level of the whole reproduction
#' rate, as given by [reproduction_level()], so a method that passes one on
#' has to convert it to the level of the part it passes it to, as above.
#'
#' Without a `balanceRDD()` method, the balancing function of the named
#' density-dependence function would not know about the change, so mizer
#' instead solves numerically for `erepro`, as described above, through the
#' whole of [projectRDD()].
#'
#' @param params A MizerParams object
#' @param rdi The density-independent reproduction rate of each species.
#'   Defaults to the rate at the initial state.
#' @param rdd_target The density-dependent reproduction rate of each species
#'   that the balancing should achieve. Defaults to the rate that keeps the
#'   initial egg abundance constant, see [getRequiredRDD()]. Species for which
#'   it is `NA` are left alone.
#' @param ... Arguments passed on to the balancing function. They choose which
#'   parameters are changed, see the section *Balancing functions*.
#' @param info_level Controls the amount of information messages and warnings
#'   that are shown. Higher levels lead to more messages, `info_level = 0`
#'   gives silence. The default is taken from the `mizer_info_level` option,
#'   see [default_info_level()].
#' @param erepro Reproductive efficiency for each species. See
#'   [setBevertonHolt()].
#' @param R_max Maximum reproduction rate. See [setBevertonHolt()].
#' @param reproduction_level The reproduction level, see
#'   [reproduction_level()].
#' @param R_factor Legacy alternative for specifying
#'   `reproduction_level = 1 / R_factor`.
#'
#' @return The MizerParams object with the reproduction parameters changed.
#' @seealso [setBevertonHolt()], [reproduction_level()], [getRequiredRDD()]
#' @export
#' @family functions for setting parameters
#' @examples
#' params <- setReproduction(NS_params, RDD = "RickerRDD")
#' params <- balanceRDD(params)
#' species_params(params)$ricker_b
#' all.equal(getRDD(params), getRequiredRDD(params))
balanceRDD <- function(params, rdi = getRDI(params),
                       rdd_target = getRequiredRDD(params), ...,
                       info_level = default_info_level()) {
    UseMethod("balanceRDD")
}

#' @export
balanceRDD.MizerParams <- function(params, rdi = getRDI(params),
                                   rdd_target = getRequiredRDD(params), ...,
                                   info_level = default_info_level()) {
    with_info_level(info_level = info_level, {
    no_sp <- nrow(params@species_params)
    if (length(rdi) != no_sp || length(rdd_target) != no_sp) {
        stop("`rdi` and `rdd_target` need to have one value for each species.",
             call. = FALSE)
    }
    storage.mode(rdd_target) <- "double"
    # A `NULL` argument has not been given. Callers use this to pass on an
    # argument that they may or may not have.
    args <- list(...)
    args <- args[!vapply(args, is.null, logical(1))]
    if (all(is.na(rdd_target))) {
        return(params)
    }

    if (rdd_depends_on_time(params, rdi)) {
        signal_info("RDD", paste0(
            "The density-dependent reproduction rate of this model depends ",
            "on time, so there is no steady state to balance it to. The ",
            "reproduction parameters have been left unchanged."),
            level = 1, severity = "warning", unhandled = "show")
        return(params)
    }

    rdd_name <- params@rates_funcs$RDD
    wrappers <- rdd_wrappers(params)
    if (length(wrappers$unbalanced) > 0) {
        classes <- paste0("`", wrappers$unbalanced, "`", collapse = ", ")
        if (length(wrappers$balanced) > 0) {
            signal_info("RDD", paste0(
                "The `projectRDD()` method for ", classes, " changes the ",
                "reproduction rate but there is no `balanceRDD()` method for ",
                "it, while other extensions do have one, so the reproduction ",
                "can not be balanced. The reproduction parameters have been ",
                "left unchanged."),
                level = 1, severity = "warning", unhandled = "show")
            return(params)
        }
        reason <- paste0("the `projectRDD()` method for ", classes,
                         " changes the reproduction rate and there is no ",
                         "`balanceRDD()` method for it")
        rdd_fun <- function(x) {
            projectRDD(params, rdi = x, species_params = params@species_params,
                       t = 0)
        }
    } else {
        balance_fn <- get0(paste0("balance_", rdd_name), mode = "function")
        if (!is.null(balance_fn)) {
            unsupported <- setdiff(names(args), names(formals(balance_fn)))
            if (length(unsupported) > 0) {
                stop("`balance_", rdd_name, "()` does not accept the ",
                     "argument", if (length(unsupported) > 1) "s", " ",
                     paste0("`", unsupported, "`", collapse = ", "),
                     ". To switch to Beverton-Holt reproduction instead, ",
                     "use `setBevertonHolt()`.", call. = FALSE)
            }
            # Call by name and with symbols, so that an error message shows
            # a readable call rather than the whole params object.
            call_args <- c(list(quote(params), rdi = quote(rdi),
                                rdd_target = quote(rdd_target)),
                           sapply(names(args),
                                  function(n) call("[[", quote(args), n),
                                  simplify = FALSE))
            return(do.call(paste0("balance_", rdd_name), call_args,
                           envir = environment()))
        }
        reason <- paste0("there is no `balance_", rdd_name, "()` function")
        # Only the named function: any extension that changes the rate has
        # already adjusted `rdi` and `rdd_target` before passing them on.
        rdd_fun <- function(x) {
            projectRDD.MizerParams(params, rdi = x,
                                   species_params = params@species_params,
                                   t = 0)
        }
    }
    if (length(args) > 0) {
        stop("Because ", reason, ", mizer can only balance the reproduction ",
             "by adjusting `erepro`, so the argument",
             if (length(args) > 1) "s", " ",
             paste0("`", names(args), "`", collapse = ", "),
             " can not be used. To switch to Beverton-Holt reproduction ",
             "instead, use `setBevertonHolt()`.", call. = FALSE)
    }
    balance_rdd_numerically(params, rdi, rdd_target, rdd_fun, reason)
    })
}

#' Balance reproduction by solving numerically for `erepro`
#'
#' The fallback of [balanceRDD()] for density dependence without a balancing
#' function. Solves `rdd_fun(rdi * erepro_new / erepro) == rdd_target` for the
#' new `erepro` of all species at once, so that density dependence that
#' couples species is handled too. The unknowns are the logarithms of the
#' factors by which `erepro` changes, which keeps `erepro` positive.
#'
#' @param params A MizerParams object.
#' @param rdi,rdd_target As in [balanceRDD()].
#' @param rdd_fun A function giving the density-dependent reproduction rate of
#'   all species from their density-independent rate.
#' @param reason Why the numerical solution is used, completing a sentence
#'   starting with "Because".
#' @return The MizerParams object with `erepro` changed, or unchanged with a
#'   warning if no solution was found.
#' @noRd
balance_rdd_numerically <- function(params, rdi, rdd_target, rdd_fun,
                                    reason) {
    species <- params@species_params$species
    sel <- !is.na(rdd_target)
    stop_if_no_eggs(species[sel & rdi == 0 & rdd_target > 0],
                    "value of `erepro`")
    # A species without egg production, such as one that is absent, has zero
    # reproduction whatever its parameters, so there is nothing to balance.
    sel <- sel & rdi > 0
    if (!any(sel)) {
        return(params)
    }
    target <- rdd_target[sel]
    scale <- ifelse(target > 0, target, 1)
    residual <- function(u) {
        x <- rdi
        x[sel] <- rdi[sel] * exp(u)
        (rdd_fun(x)[sel] - target) / scale
    }
    sol <- tryCatch(
        nleqslv::nleqslv(rep(0, sum(sel)), residual,
                         control = list(ftol = 1e-12, xtol = 1e-14,
                                        maxit = 200)),
        error = function(e) NULL)
    unsolved <- if (is.null(sol) || !all(is.finite(sol$x))) {
        rep(TRUE, sum(sel))
    } else {
        !is.finite(sol$fvec) | abs(sol$fvec) > 1e-8
    }
    if (any(unsolved)) {
        signal_info("erepro", paste0(
            "Because ", reason, ", mizer tried to balance the reproduction by ",
            "adjusting `erepro`, but no value of `erepro` gives the required ",
            "reproduction rate for the following species: ",
            paste(species[sel][unsolved], collapse = ", "),
            ". The reproduction parameters have been left unchanged, so the ",
            "model is not at steady state."),
            level = 1, severity = "warning", unhandled = "show")
        return(params)
    }

    erepro <- params@species_params$erepro
    erepro[sel] <- erepro[sel] * exp(sol$x)
    params <- set_erepro(params, which(sel), erepro[sel])
    signal_info("erepro", paste0(
        "Because ", reason, ", mizer has balanced the reproduction by ",
        "adjusting `erepro` numerically."),
        level = 2, unhandled = "show")
    params@time_modified <- lubridate::now()
    params
}

#' Set `erepro` and record it as given
#'
#' The balancing functions treat the `erepro` they calculate as given, because
#' the user chose to have it calculated this way and won't want it changed by
#' the next call to `setParams()` or `given_species_params<-()`. Also warns
#' about values greater than 1.
#'
#' @param params A MizerParams object.
#' @param sp_idx Indices of the species to change.
#' @param values The new values of `erepro` for those species.
#' @return The MizerParams object.
#' @noRd
set_erepro <- function(params, sp_idx, values) {
    if (!"erepro" %in% names(params@given_species_params)) {
        params@given_species_params$erepro <- params@species_params$erepro
    }
    params@species_params$erepro[sp_idx] <- values
    params@given_species_params$erepro[sp_idx] <- values
    warn_erepro_above_one(params@species_params$species[sp_idx][values > 1])
    params
}

#' Set a parameter of the density dependence and record it as given
#'
#' @param params A MizerParams object.
#' @param par The name of the species parameter.
#' @param sp_idx Indices of the species to change.
#' @param values The new values for those species.
#' @return The MizerParams object.
#' @noRd
set_rdd_param <- function(params, par, sp_idx, values) {
    if (!par %in% names(params@species_params)) {
        params@species_params[[par]] <- NA_real_
    }
    if (!par %in% names(params@given_species_params)) {
        params@given_species_params[[par]] <- NA_real_
    }
    params@species_params[[par]][sp_idx] <- values
    params@given_species_params[[par]][sp_idx] <- values
    params
}

#' Warn about species that need an `erepro` greater than 1
#' @param species Names of the species.
#' @noRd
warn_erepro_above_one <- function(species) {
    if (length(species) == 0) return(invisible())
    signal_info("erepro", paste0(
        "The following species require an unrealistic value greater ",
        "than 1 for `erepro`: ",
        paste(species, collapse = ", ")),
        level = 1, severity = "warning", unhandled = "show")
}

#' Stop for species that have no egg production but need reproduction
#' @param species Names of the species.
#' @param what What could not sustain them, for the message.
#' @noRd
stop_if_no_eggs <- function(species, what) {
    if (length(species) == 0) return(invisible())
    stop("The following species have no reproduction, so no ", what,
         " can sustain them: ", paste(species, collapse = ", "),
         ". Check that they have mature individuals in `initialN(params)`.")
}

#' Expand values given for some or all species
#'
#' A single unnamed value is used for all species, an unnamed vector must have
#' one value per species and a named vector gives values for the named species.
#' `NA` values are dropped.
#'
#' @param params A MizerParams object.
#' @param values The values.
#' @return A named vector holding the non-`NA` values, named by species.
#' @noRd
expand_species_values <- function(params, values) {
    no_sp <- nrow(params@species_params)
    if (length(values) == 1 && is.null(names(values))) {
        values <- rep(values, no_sp)
    }
    if (length(values) != no_sp && is.null(names(values))) {
        stop("You need to supply a vector of length ", no_sp,
             " or a single number or a named vector.")
    }
    if (is.null(names(values))) {
        names(values) <- params@species_params$species
    }
    values <- values[!is.na(values)]
    if (length(values) > 0 && !all(is.numeric(values))) {
        stop("You provided invalid non-numeric values.")
    }
    values
}

#' Find the extensions that change the reproduction rate
#'
#' An extension class changes the reproduction rate when it has a
#' [projectRDD()] method. It takes part in balancing when it also has a
#' [balanceRDD()] method.
#'
#' @param params A MizerParams object.
#' @return A list with the character vectors `all` (classes with a
#'   `projectRDD()` method), `unbalanced` (those without a `balanceRDD()`
#'   method) and `balanced` (classes with a `balanceRDD()` method).
#' @noRd
rdd_wrappers <- function(params) {
    classes <- setdiff(class(params), "MizerParams")
    has_method <- function(generic) {
        vapply(classes, function(cl) {
            !is.null(utils::getS3method(generic, cl, optional = TRUE,
                                        envir = environment(projectRDD)))
        }, logical(1))
    }
    wraps <- has_method("projectRDD")
    balances <- has_method("balanceRDD")
    list(all = classes[wraps],
         unbalanced = classes[wraps & !balances],
         balanced = classes[balances])
}

#' Whether the reproduction is balanced by `balance_BevertonHoltRDD()`
#'
#' True when the model's density-dependence function is [BevertonHoltRDD()] and
#' every extension that changes the reproduction rate takes part in balancing,
#' so that the Beverton-Holt parameters `erepro`, `R_max` and
#' `reproduction_level` can be preserved.
#'
#' @param params A MizerParams object.
#' @return `TRUE` or `FALSE`.
#' @noRd
uses_beverton_holt_balancer <- function(params) {
    identical(params@rates_funcs$RDD, "BevertonHoltRDD") &&
        length(rdd_wrappers(params)$unbalanced) == 0
}

#' Whether the reproduction rate depends on time
#'
#' Evaluates the reproduction rate at a few times for the same
#' density-independent rate. A rate that changes, or that can not be evaluated
#' at one of those times, depends on time. The times are not round numbers so
#' that seasonal patterns and short pulses are unlikely to be missed.
#'
#' @param params A MizerParams object.
#' @param rdi The density-independent reproduction rate.
#' @return `TRUE` or `FALSE`.
#' @noRd
rdd_depends_on_time <- function(params, rdi) {
    rdd_at <- function(t) {
        tryCatch(unname(projectRDD(params, rdi = rdi,
                                   species_params = params@species_params,
                                   t = t)),
                 error = function(e) NULL)
    }
    rdd0 <- rdd_at(0)
    # If the rate can not be evaluated at all, leave the error to the
    # balancing.
    if (is.null(rdd0)) return(FALSE)
    for (t in c(0.0371, 0.2917, 0.5, 0.7453)) {
        if (!identical(rdd_at(t), rdd0)) return(TRUE)
    }
    FALSE
}

#' The reproduction level of any density dependence
#'
#' One minus the elasticity of the reproduction rate with respect to the
#' density-independent reproduction rate, calculated by a central difference
#' through [projectRDD()]. The density-independent rate of all species is
#' changed by the same factor, so for density dependence that couples species
#' this is the response to a change in the egg production of all of them.
#'
#' @param params A MizerParams object.
#' @return A named vector with the reproduction level of each species.
#' @noRd
rdd_elasticity_level <- function(params) {
    rdi <- getRDI(params)
    f <- function(x) {
        unname(projectRDD(params, rdi = x,
                          species_params = params@species_params, t = 0))
    }
    h <- 1e-6
    rdd <- f(rdi)
    level <- 1 - (f(rdi * (1 + h)) - f(rdi * (1 - h))) / (2 * h * rdd)
    # Without reproduction the level is the limit for small egg production:
    # zero if the reproduction rate vanishes with it.
    none <- rdd == 0
    level[none] <- ifelse(rdi[none] == 0, 0, NA_real_)
    names(level) <- params@species_params$species
    level
}

#' @rdname balanceRDD
#' @export
balance_noRDD <- function(params, rdi, rdd_target) {
    species <- params@species_params$species
    sel <- !is.na(rdd_target)
    stop_if_no_eggs(species[sel & rdi == 0 & rdd_target > 0],
                    "value of `erepro`")
    sel <- which(sel & rdi > 0)
    if (length(sel) == 0) return(params)
    params <- set_erepro(params, sel, params@species_params$erepro[sel] *
                             rdd_target[sel] / rdi[sel])
    params@time_modified <- lubridate::now()
    params
}

#' @rdname balanceRDD
#' @export
balance_constantRDD <- function(params, rdi, rdd_target) {
    sel <- which(!is.na(rdd_target))
    if (length(sel) == 0) return(params)
    params <- set_rdd_param(params, "constant_reproduction", sel,
                            rdd_target[sel])
    params@time_modified <- lubridate::now()
    params
}

#' @rdname balanceRDD
#' @export
balance_RickerRDD <- function(params, rdi, rdd_target, erepro = NULL) {
    balance_overcompensating(params, rdi, rdd_target, erepro,
                             par = "ricker_b",
                             b_fun = function(rdi, target, sel) {
                                 log(rdi / target) / rdi
                             })
}

#' @rdname balanceRDD
#' @export
balance_SheperdRDD <- function(params, rdi, rdd_target, erepro = NULL) {
    shape <- params@species_params$sheperd_c
    if (is.null(shape)) {
        stop("The species_params dataframe must contain a column sheperd_c.")
    }
    balance_overcompensating(params, rdi, rdd_target, erepro,
                             par = "sheperd_b",
                             b_fun = function(rdi, target, sel) {
                                 ((rdi - target) / target)^(1 / shape[sel]) /
                                     rdi
                             })
}

#' Balance a density dependence with one parameter by keeping `erepro`
#'
#' Shared by the Ricker and Shepherd balancing functions. Both can only reduce
#' the reproduction rate below the density-independent rate, so where that
#' rate is too small `erepro` is increased until it equals the required rate
#' and the density dependence is switched off.
#'
#' @param params A MizerParams object.
#' @param rdi,rdd_target,erepro As in [balanceRDD()].
#' @param par The name of the parameter of the density dependence.
#' @param b_fun A function of the density-independent rate, the target rate
#'   and the selected species indices that gives the parameter value.
#' @return The MizerParams object.
#' @noRd
balance_overcompensating <- function(params, rdi, rdd_target, erepro,
                                     par, b_fun) {
    sp <- params@species_params
    species <- sp$species
    erepro_new <- sp$erepro
    affected <- !is.na(rdd_target)
    if (!is.null(erepro)) {
        values <- expand_species_values(params, erepro)
        idx <- match(valid_species_arg(params, names(values)), species)
        erepro_new[idx] <- values
        affected <- affected & seq_along(species) %in% idx
    }
    stop_if_no_eggs(species[affected & rdi == 0 & rdd_target > 0],
                    "value of `erepro`")
    sel <- which(affected & rdi > 0)
    if (length(sel) == 0) return(params)

    rdi_new <- rdi[sel] * erepro_new[sel] / sp$erepro[sel]
    target <- rdd_target[sel]
    wrong <- rdi_new < target
    if (any(wrong)) {
        erepro_new[sel][wrong] <- (erepro_new[sel] * target / rdi_new)[wrong]
        rdi_new[wrong] <- target[wrong]
        signal_info("erepro", paste0(
            "For the following species `erepro` has been increased to the ",
            "smallest possible value: ",
            paste0("erepro[", species[sel][wrong], "] = ",
                   signif(erepro_new[sel][wrong], 3), collapse = "; ")),
            level = 1, severity = "warning", unhandled = "show")
    }
    params <- set_erepro(params, sel, erepro_new[sel])
    params <- set_rdd_param(params, par, sel, b_fun(rdi_new, target, sel))
    params@time_modified <- lubridate::now()
    params
}
