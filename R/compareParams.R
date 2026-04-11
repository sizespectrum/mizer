#' Compare two MizerParams objects and print out differences
#'
#' @param params1 First MizerParams object
#' @param params2 Second MizerParams object
#' @param ... Additional arguments passed to the method.
#'
#' @return Invisibly returns a character vector of difference messages, one
#'   element per difference. As a side effect, prints the differences in a
#'   human-readable format.
#' @export
#' @examples
#' params1 <- NS_params
#' params2 <- params1
#' species_params(params2)$w_mat[1] <- 10
#' compareParams(params1, params2)
compareParams <- function(params1, params2, ...)
    UseMethod("compareParams")

#' @export
compareParams.MizerParams <- function(params1, params2, ...) {
    params1 <- validParams(params1)
    params2 <- validParams(params2)

    result <- character()

    # species parameters ----
    result <- c(result,
                compareSpeciesParams(params1@species_params,
                                     params2@species_params),
                compareSpeciesParams(params1@given_species_params,
                                     params2@given_species_params,
                                     text = "given species parameters"))

    # resource parameters ----
    res1 <- params1@resource_params
    res1 <- res1[order(names(res1))]
    res2 <- params2@resource_params
    res2 <- res2[order(names(res2))]
    res_eq <- all.equal(res1, res2, scale = 1)
    if (!isTRUE(res_eq)) {
        msg <- paste("The following resource parameters differ:",
                     toString(res_eq))
        result <- c(result, msg)
    }

    # size grid ----
    if (length(params1@w) != length(params2@w)) {
        msg <- "The number of community size bins is different."
        result <- c(result, msg)
    } else {
        if (!isTRUE(all.equal(params1@w, params2@w, scale = 1))) {
            msg <- "The community size bins differ."
            result <- c(result, msg)
        }
    }
    if (length(params1@w_full) != length(params2@w_full)) {
        msg <- "The number of resource size bins is different."
        result <- c(result, msg)
    } else {
        if (!isTRUE(all.equal(params1@w_full, params2@w_full, scale = 1))) {
            msg <- "The resource size bins differ."
            result <- c(result, msg)
        }
    }

    # other slots ----
    for (sl in slotNames(params1)) {
        if (sl %in% c("w", "w_full", "species_params",
                      "given_species_params", "resource_params")) {
            next
        }
        s1 <- slot(params1, sl)
        s2 <- slot(params2, sl)
        eq <- all.equal(s1, s2, scale = 1)
        if (!isTRUE(eq)) {
            attr(s1, "comment") <- NULL
            attr(s2, "comment") <- NULL
            eq_no_comment <- all.equal(s1, s2, scale = 1)
            if (isTRUE(eq_no_comment)) {
                c1 <- comment(slot(params1, sl))
                c2 <- comment(slot(params2, sl))
                msg <- paste0('The "', sl, '" slots differ only in their ',
                              "comment attributes.",
                              "\n  params1: ", if (is.null(c1)) "<none>" else paste(c1, collapse = " "),
                              "\n  params2: ", if (is.null(c2)) "<none>" else paste(c2, collapse = " "))
            } else {
                msg <- paste("The", sl, "slots do not agree:",
                             toString(eq))
                if (is.array(s1)) {
                    detail <- summariseArrayDiff(s1, s2)
                    if (length(detail) > 0) {
                        msg <- paste0(msg, "\n", detail)
                    }
                }
            }
            result <- c(result, msg)
        }
    }
    if (length(result) == 0) {
        cat("No differences\n")
        return(invisible("No differences"))
    }
    cat(result, sep = "\n\n")
    cat("\n")
    invisible(result)
}

summariseArrayDiff <- function(a1, a2) {
    abs_diff <- abs(a2 - a1)
    dn <- dimnames(a1)
    if (!is.null(dn) && !is.null(dn[[1]])) {
        row_max <- apply(abs_diff, 1, max)
        affected <- row_max[row_max > 0]
        if (length(affected) > 0) {
            detail <- paste(names(affected), signif(affected, 3),
                            sep = ": ", collapse = ", ")
            return(paste0("  Max |diff|: ", detail))
        }
    } else {
        return(paste0("  Max |diff|: ", signif(max(abs_diff), 3),
                      " (", sum(abs_diff > 0), " of ", length(abs_diff),
                      " elements differ)"))
    }
    character()
}

compareSpeciesParams <- function(species_params1,
                                 species_params2,
                                 text = "species parameters") {

    result <- character()

    # species parameters ----
    sp_names1 <- names(species_params1)
    sp_names2 <- names(species_params2)
    sp_diff1 <- setdiff(sp_names1, sp_names2)
    if (length(sp_diff1) > 0) {
        msg <- paste0("params1 has the following additional ",
                      text, ": ", toString(sp_diff1))
        result <- c(result, msg)
    }
    sp_diff2 <- setdiff(sp_names2, sp_names1)
    if (length(sp_diff2) > 0) {
        msg <- paste0("params2 has the following additional ",
                      text, ": ", toString(sp_diff2))
        result <- c(result, msg)
    }
    sp_names <- intersect(sp_names1, sp_names2)

    sp_eq <- all.equal(species_params1[, sp_names],
                       species_params2[, sp_names], scale = 1)
    if (!isTRUE(sp_eq)) {
        msg <- paste("The following", text, "differ:",
               toString(sp_eq))
        result <- c(result, msg)
        # ctable <- compareDF::compare_df(params1@species_params[, sp_names],
        #                       params2@species_params[, sp_names],
        #                       keep_unchanged_cols = FALSE)
        # print(ctable$comparison_df)
    }

    result
}
