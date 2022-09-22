#' Compare two MizerParams objects and print out differences
#'
#'`r lifecycle::badge("experimental")`
#'
#' @param params1 First MizerParams object
#' @param params2 Second MizerParams object
#' @export
#' @examples
#' \dontrun{
#' sp1 <- NS_species_params
#' params1 <- newMultispeciesParams(sp1)
#' sp2 <- sp1
#' sp2$w_mat[1] <- 10
#' params2 <- newMultispeciesParams(sp2)
#' compareParams(params1, params2)
#' }
compareParams <- function(params1, params2) {
    validObject(params1)
    validObject(params1)
    assert_that(is(params1, "MizerParams"))
    assert_that(is(params2, "MizerParams"))

    result <- character()

    # species parameters ----
    sp_names1 <- names(params1@species_params)
    sp_names2 <- names(params2@species_params)
    sp_diff1 <- setdiff(sp_names1, sp_names2)
    if (length(sp_diff1) > 0) {
        msg <- paste("params1 has the following additional species parameters:",
                     toString(sp_diff1))
        result <- c(result, msg)
    }
    sp_diff2 <- setdiff(sp_names2, sp_names1)
    if (length(sp_diff2) > 0) {
        msg <- paste("params2 has the following additional species parameters:",
                     toString(sp_diff2))
        result <- c(result, msg)
    }
    sp_names <- intersect(sp_names1, sp_names2)

    sp_eq <- all.equal(params1@species_params[, sp_names],
                       params2@species_params[, sp_names], scale = 1)
    if (!isTRUE(sp_eq)) {
        msg <- paste("The following species parameters differ:",
               toString(sp_eq))
        result <- c(result, msg)
        # ctable <- compareDF::compare_df(params1@species_params[, sp_names],
        #                       params2@species_params[, sp_names],
        #                       keep_unchanged_cols = FALSE)
        # print(ctable$comparison_df)
    }

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
        if (sl %in% c("w", "w_full", "species_params", "resource_params")) {
            next
        }
        eq <- all.equal(slot(params1, sl), slot(params2, sl), scale = 1)
        if (!isTRUE(eq)) {
            msg <- paste("The", sl, "slots do not agree:",
                         toString(eq))
            result <- c(result, msg)
        }
    }
    result
}
