# Setter factories for scanModel()
#
# Copyright 2026 Gustav Delius.
# Distributed under the GPL 3 or later.

#' Setters for scanning a model
#'
#' `r lifecycle::badge("experimental")`
#' These functions build the `set_func` that [scanModel()] uses to apply each
#' scan value to the model. Each returns a function of `(params, value)` that
#' returns a modified `MizerParams` object, carrying attributes that let
#' [scanModel()] label the axis and mark reference lines without being told.
#'
#' You are not restricted to these. Any function of `(params, value)` returning
#' a `MizerParams` will do, as long as it is **idempotent**: with
#' `continuation = TRUE` it is applied to the object it returned at the previous
#' scan value, so applying it twice must give the same thing as applying it
#' once. Setting a value is idempotent; appending something is not, which is why
#' [scanFishingMortality()] checks whether its gear is already there.
#'
#' \describe{
#'   \item{`scanEffort()`}{Scans the fishing effort. With `gear = NULL` the same
#'     effort is applied to every gear, which is what a bifurcation diagram over
#'     fishing effort needs.}
#'   \item{`scanFishingMortality()`}{Scans the fishing mortality on one species
#'     while leaving the fishing on every other species alone. It does this by
#'     adding a temporary gear that catches only the target species with
#'     catchability 1, so that its effort *is* the fishing mortality, and
#'     switching off the catchability of the gears it replaces. If several gears
#'     catch the species you can name the one whose mortality is to be varied,
#'     and the others go on fishing unchanged.}
#'   \item{`scanSpeciesParam()`}{Scans any species parameter. It writes to
#'     [given_species_params()], so the change propagates through to the rates
#'     that depend on it rather than being frozen. Note that scanning a
#'     parameter that mizer calculates for itself from others — `gamma` when
#'     `f0` is given, say — will not work, because the calculated value
#'     overwrites what you set.}
#' }
#'
#' @param gear For `scanEffort()`, the name of the gear whose effort is scanned,
#'   or NULL (default) to scan the effort of every gear together. For
#'   `scanFishingMortality()`, the name of the gear whose fishing mortality on
#'   the target species is scanned. Only needed when several gears catch the
#'   species; if NULL, the fishing mortality from all of them is replaced.
#' @param species The name of the target species.
#' @param parameter The name of the species parameter to scan.
#'
#' @return A function of `(params, value)` returning a `MizerParams` object.
#'
#' @family scan functions
#' @concept scan
#' @seealso [scanModel()]
#' @export
#' @examples
#' \donttest{
#' # The fishing mortality on Cod alone, leaving the other species alone
#' plot(scanModel(NS_params, scan_values = seq(0, 1.2, 0.3),
#'                set_func = scanFishingMortality("Cod"),
#'                value_func = getYield, species = "Cod"))
#' }
scanEffort <- function(gear = NULL) {
    force(gear)
    f <- function(params, value) {
        if (is.null(gear)) {
            initial_effort(params) <- value
        } else {
            effort <- initial_effort(params)
            missing <- setdiff(gear, names(effort))
            if (length(missing) > 0) {
                stop("There is no gear called ",
                     paste(missing, collapse = ", "), ".")
            }
            effort[gear] <- value
            initial_effort(params) <- effort
        }
        params
    }
    attr(f, "scan_name") <- if (is.null(gear)) {
        "Fishing effort"
    } else {
        paste0("Effort of ", paste(gear, collapse = ", "))
    }
    attr(f, "current_scan_value") <- function(params) {
        effort <- initial_effort(params)
        unname(if (is.null(gear)) effort[[1]] else mean(effort[gear]))
    }
    f
}

#' @rdname scanEffort
#' @export
scanFishingMortality <- function(species, gear = NULL) {
    force(species)
    force(gear)
    f <- function(params, value) {
        # Install the temporary gear only if it is not already there. This is
        # what makes the function idempotent, which it has to be because with
        # continuation it is handed the object it returned last time.
        # Reinstalling would append a second "tmp" row and rebuild the
        # selectivity array at every scan value.
        if (!("tmp" %in% dimnames(params@catchability)[[1]])) {
            params <- install_tmp_gear(params, species, gear)
        }
        effort <- initial_effort(params)
        effort[["tmp"]] <- value
        initial_effort(params) <- effort
        params
    }
    attr(f, "scan_name") <- paste0("Fishing mortality on ", species)
    attr(f, "scan_units") <- "1/year"
    attr(f, "current_scan_value") <- function(params) {
        gp <- gear_params(params)
        gp$gear <- as.character(gp$gear)
        sel <- select_gear_rows(gp, species, gear)
        sum(initial_effort(params)[gp$gear[sel]] * gp$catchability[sel])
    }
    attr(f, "reference_lines") <- function(params) {
        sp <- params@species_params
        if (!("F_MSY" %in% names(sp))) return(NULL)
        idx <- match(species, sp$species)
        if (is.na(idx)) return(NULL)
        value <- sp$F_MSY[[idx]]
        if (is.na(value)) NULL else c(F_MSY = value)
    }
    f
}

#' @rdname scanEffort
#' @export
scanSpeciesParam <- function(species, parameter) {
    force(species)
    force(parameter)
    f <- function(params, value) {
        idx <- match(species, params@species_params$species)
        if (is.na(idx)) {
            stop("There is no species called ", species, " in this model.")
        }
        given_species_params(params)[[parameter]][[idx]] <- value
        params
    }
    attr(f, "scan_name") <- paste0(parameter, " of ", species)
    attr(f, "current_scan_value") <- function(params) {
        idx <- match(species, params@species_params$species)
        if (is.na(idx) || !(parameter %in% names(params@species_params))) {
            stop("The model has no ", parameter, " for ", species, ".")
        }
        params@species_params[[parameter]][[idx]]
    }
    f
}

#' The gear params rows whose fishing mortality is to be varied
#'
#' @param gp The gear params data frame, with a character `gear` column.
#' @param species The target species.
#' @param gear The selected gear, or NULL for all gears catching the species.
#' @return An integer vector of row indices.
#' @keywords internal
select_gear_rows <- function(gp, species, gear = NULL) {
    sel <- which(gp$species == species)
    if (length(sel) == 0) {
        stop(species, " is not selected by any gear.")
    }
    if (!is.null(gear)) {
        gear <- as.character(gear)
        if (length(gear) != 1) {
            stop("You can only select a single gear.")
        }
        chosen <- sel[gp$gear[sel] == gear]
        if (length(chosen) == 0) {
            stop("The gear ", gear, " does not catch ", species, ".")
        }
        return(chosen)
    }
    sel
}

#' Add a gear that exerts the fishing mortality being scanned
#'
#' Copies the selectivity of the first gear catching the species onto a new gear
#' called `"tmp"` with catchability 1, so that the effort of that gear is the
#' fishing mortality it exerts, and switches off the catchability of the gears it
#' replaces. The fishing on every other species is untouched.
#'
#' @param params A MizerParams object.
#' @param species The target species.
#' @param gear The gear whose mortality is replaced, or NULL for all of the
#'   gears catching the species.
#' @return The MizerParams object with the extra gear.
#' @keywords internal
install_tmp_gear <- function(params, species, gear = NULL) {
    species <- valid_species_arg(params, species, error_on_empty = TRUE)
    if (length(species) != 1) {
        stop("You can only scan the fishing mortality on one species at a ",
             "time.")
    }
    gp <- gear_params(params)
    gp$gear <- as.character(gp$gear)
    sel <- select_gear_rows(gp, species, gear)
    if (is.null(gear) && length(sel) > 1) {
        signal_info("scan", paste0(
            "Several gears catch ", species, ". The fishing mortality from ",
            "all of them will be replaced. Use the `gear` argument if you ",
            "want to vary the fishing mortality from only one of them."),
            unhandled = "show")
    }
    gp_extra <- gp[sel[[1]], ]
    gp_extra$gear <- "tmp"
    gp_extra$catchability <- 1
    gp$catchability[sel] <- 0
    gear_params(params) <- rbind(gp, gp_extra)
    initial_effort(params)["tmp"] <- 1
    params
}
