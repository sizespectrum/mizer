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
#'     and the others go on fishing unchanged. The added gear is given a name
#'     the model is not already using, so a model that happens to have a gear
#'     called `"scan"` is not disturbed.}
#'   \item{`scanSpeciesParam()`}{Scans any species parameter. It assigns to
#'     [species_params()], so the value is recorded as a given one and the
#'     change propagates through to the rates that depend on it. The parameter
#'     has to be one the model already has; add the column first if it is not.}
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
    # The name of the gear this setter installs, chosen on first use so that it
    # cannot collide with a gear the model already has.
    scan_gear <- NULL
    f <- function(params, value) {
        # Install the gear unless this model already carries the installation
        # this setter makes. Whether it does has to be judged from the state of
        # the gears, not from the name: a gear that merely looks like ours
        # leaves the gears we are supposed to have switched off still fishing,
        # so the scanned mortality would add to them instead of replacing them.
        # Skipping the install when it has already been done is what makes the
        # setter idempotent under continuation, where it is handed the object it
        # returned last time and re-installing would append a second row and
        # rebuild the selectivity array at every scan value.
        if (!scan_gear_installed(params, scan_gear, species, gear)) {
            if (is.null(scan_gear) || scan_gear %in% gear_names(params)) {
                # Either this is the first use, or a gear of that name is
                # present but is not our installation -- a setter can be handed
                # more than one model, and the name belongs to that model.
                # Leave it alone and take a name that is free here.
                scan_gear <<- free_gear_name(params)
            }
            params <- install_tmp_gear(params, species, gear, scan_gear)
        }
        effort <- initial_effort(params)
        effort[[scan_gear]] <- value
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
        sp <- species_params(params)
        idx <- match(species, sp$species)
        if (is.na(idx)) {
            stop("There is no species called ", species, " in this model.")
        }
        # The full species params table, not `given_species_params()`, which
        # holds only what the user supplied. Writing a parameter that is absent
        # there would build a malformed column rather than set a value, and
        # assigning the full table back is also what runs the validation and
        # propagates the change to the rates that depend on it.
        if (!(parameter %in% names(sp))) {
            stop("This model has no species parameter called `", parameter,
                 "`. Add the column to the species parameters before scanning ",
                 "it.")
        }
        sp[[parameter]][[idx]] <- value
        species_params(params) <- sp
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

#' A gear name that the model is not already using
#'
#' @param params A MizerParams object.
#' @return A string naming a gear that does not exist in `params`.
#' @keywords internal
free_gear_name <- function(params) {
    taken <- gear_names(params)
    # The leading dot marks the gear as one mizer added rather than one the
    # user set up, which makes an accidental clash with a real gear unlikely.
    # One more candidate than there are gears, so at least one is always free.
    candidates <- paste0(".scan_", seq_len(length(taken) + 1L))
    setdiff(candidates, taken)[[1]]
}

#' Every gear name a model uses
#'
#' @param params A MizerParams object.
#' @return A character vector of gear names.
#' @keywords internal
gear_names <- function(params) {
    unique(c(dimnames(params@catchability)[[1]],
             as.character(gear_params(params)$gear)))
}

#' Has a fishing-mortality scan already been installed in this model?
#'
#' A gear of the right name proves nothing: it might be one the model already
#' had, and setting its effort would then leave the fishing the scan is supposed
#' to replace still switched on, so the scanned mortality would be added to the
#' existing mortality rather than replacing it. Nor is it enough for the gear to
#' *look* like the one [scanFishingMortality()] adds, for the same reason.
#'
#' What is checked is therefore the whole installation, exactly as
#' [install_tmp_gear()] leaves it: the name carried by exactly one row, which
#' catches the target species with catchability 1, **and** at least one original
#' gear still present with every gear it was supposed to replace switched off.
#'
#' @param params A MizerParams object.
#' @param gear_name The name of the gear to check.
#' @param species The target species.
#' @param gear The gear whose mortality the scan replaces, or NULL for all of
#'   the gears catching the species.
#' @return TRUE if this model already carries the installation.
#' @keywords internal
scan_gear_installed <- function(params, gear_name, species, gear = NULL) {
    if (is.null(gear_name)) return(FALSE)
    gp <- gear_params(params)
    gp$gear <- as.character(gp$gear)
    # Exactly one row may carry the name, and that row must be ours. A gear of
    # that name with further rows catching other species is not our
    # installation, and scanning its effort would move those species too.
    scan <- gp$gear == gear_name
    if (sum(scan) != 1) return(FALSE)
    if (!identical(as.character(gp$species[scan]), as.character(species)) ||
            !isTRUE(gp$catchability[scan] == 1)) {
        return(FALSE)
    }
    replaced <- !scan & gp$species == species
    if (!is.null(gear)) {
        replaced <- replaced & gp$gear == gear
    }
    any(replaced) && all(gp$catchability[replaced] == 0)
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
#' called `gear_name` with catchability 1, so that the effort of that gear is
#' the fishing mortality it exerts, and switches off the catchability of the
#' gears it replaces. The fishing on every other species is untouched. The
#' effort of the new gear is left at whatever [gear_params<-()] gives it; the
#' caller sets it to the value being scanned.
#'
#' @param params A MizerParams object.
#' @param species The target species.
#' @param gear The gear whose mortality is replaced, or NULL for all of the
#'   gears catching the species.
#' @param gear_name The name to give the new gear.
#' @return The MizerParams object with the extra gear.
#' @keywords internal
install_tmp_gear <- function(params, species, gear = NULL, gear_name) {
    species <- valid_species_arg(params, species, error_on_empty = TRUE)
    if (length(species) != 1) {
        stop("You can only scan the fishing mortality on one species at a ",
             "time.")
    }
    gp <- gear_params(params)
    gp$gear <- as.character(gp$gear)
    sel <- select_gear_rows(gp, species, gear)
    if (is.null(gear) && length(sel) > 1) {
        signal_info(
            "scan",
            paste0("Several gears catch ", species, ". The fishing mortality ",
                   "from all of them will be replaced. Use the `gear` ",
                   "argument if you want to vary the fishing mortality from ",
                   "only one of them."),
            unhandled = "show"
        )
    }
    gp_extra <- gp[sel[[1]], ]
    gp_extra$gear <- gear_name
    gp_extra$catchability <- 1
    gp$catchability[sel] <- 0
    gear_params(params) <- rbind(gp, gp_extra)
    params
}
