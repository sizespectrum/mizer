# Cache rate arrays once to avoid recomputing them in every test block.
enc_small <- getEncounter(NS_params_small)
pred_mort_small <- getPredMort(NS_params_small)
pred_rate_small <- getPredRate(NS_params_small)

test_that("ArraySpeciesBySize constructor works", {
    mat <- matrix(1:30, nrow = 3, ncol = 10)
    rownames(mat) <- NS_params_small@species_params$species
    colnames(mat) <- signif(NS_params_small@w[1:10], 3)

    rate <- ArraySpeciesBySize(mat, value_name = "Test rate", units = "g/year")

    expect_s3_class(rate, "ArraySpeciesBySize")
    expect_true(is.ArraySpeciesBySize(rate))
    expect_false(is.ArraySpeciesBySize(mat))
    expect_true(is.matrix(rate))
    expect_identical(dim(rate), dim(mat))
    expect_identical(attr(rate, "value_name"), "Test rate")
    expect_identical(attr(rate, "units"), "g/year")
})

test_that("ArraySpeciesBySize constructor validates input", {
    expect_error(ArraySpeciesBySize(1:10), "`x` must be a matrix.")
})

test_that("Rate functions return ArraySpeciesBySize", {
    params <- NS_params_small

    expect_true(is.ArraySpeciesBySize(enc_small))
    expect_true(is.ArraySpeciesBySize(getFeedingLevel(params)))
    expect_true(is.ArraySpeciesBySize(getCriticalFeedingLevel(params)))
    expect_true(is.ArraySpeciesBySize(getEReproAndGrowth(params)))
    expect_true(is.ArraySpeciesBySize(getERepro(params)))
    expect_true(is.ArraySpeciesBySize(getEGrowth(params)))
    expect_true(is.ArraySpeciesBySize(getFMort(params)))
    expect_true(is.ArraySpeciesBySize(pred_mort_small))
    expect_true(is.ArraySpeciesBySize(getMort(params)))
    expect_true(is.ArraySpeciesBySize(getFlux(params)))
    expect_true(is.ArraySpeciesBySize(intake_max(params)))
    expect_true(is.ArraySpeciesBySize(metab(params)))
    expect_true(is.ArraySpeciesBySize(search_vol(params)))
    expect_true(is.ArraySpeciesBySize(ext_mort(params)))
    expect_true(is.ArraySpeciesBySize(ext_encounter(params)))
    expect_true(is.ArraySpeciesBySize(maturity(params)))
    expect_true(is.ArraySpeciesBySize(repro_prop(params)))
    expect_true(is.ArraySpeciesBySize(ext_diffusion(params)))
})

test_that("All rate functions have consistent dimnames", {
    params <- NS_params_small
    expected_dimnames <- dimnames(params@metab)

    expect_identical(dimnames(enc_small), expected_dimnames)
    expect_identical(dimnames(getFeedingLevel(params)), expected_dimnames)
    expect_identical(dimnames(getCriticalFeedingLevel(params)), expected_dimnames)
    expect_identical(dimnames(getEReproAndGrowth(params)), expected_dimnames)
    expect_identical(dimnames(getERepro(params)), expected_dimnames)
    expect_identical(dimnames(getEGrowth(params)), expected_dimnames)
    expect_identical(dimnames(getFMort(params)), expected_dimnames)
    expect_identical(dimnames(pred_mort_small), expected_dimnames)
    expect_identical(dimnames(getMort(params)), expected_dimnames)
    expect_identical(dimnames(getFlux(params)), expected_dimnames)
})

test_that("print.ArraySpeciesBySize works", {
    expect_output(print(enc_small), "Encounter rate")
    expect_output(print(enc_small), "species x")
    expect_output(print(enc_small), "g/year")
})

test_that("print.ArraySpeciesBySize truncates a wide size grid", {
    enc <- getEncounter(NS_params)
    # 12 species is below the truncation threshold, so all species still show
    expect_output(print(enc), "Cod")
    expect_output(print(enc), "Sprat")
    # but the 100 sizes must be truncated with an explanatory footer
    expect_output(print(enc), "showing .* of 100 sizes")
    expect_output(print(enc), "log-spaced")
    expect_output(print(enc), "as.data.frame")
})

test_that("print.ArraySpeciesBySize truncates a large number of species", {
    params_big <- suppressMessages(newTraitParams(no_sp = 60))
    enc_big <- getEncounter(params_big)
    expect_output(print(enc_big), "showing 8 of 60 species")
})

test_that("summary.ArraySpeciesBySize works", {
    s <- summary(enc_small)
    expect_s3_class(s, "summary.ArraySpeciesBySize")
    expect_identical(s$value_name, "Encounter rate")
    expect_identical(nrow(s$per_species), nrow(NS_params_small@species_params))
    expect_output(print(s), "Encounter rate")
})

test_that("summary.ArraySpeciesBySize covers the same sizes as plot()", {
    # The asymmetry this fixes: plot() has always dropped the values outside a
    # species' size range, and summary() reported them, so the two described
    # different arrays. The out-of-range values are also the extreme ones — the
    # encounter rate a fish far larger than its `w_max` would have — so it was
    # the summary's Min and Max that disagreed.
    s <- summary(enc_small)$per_species
    pd <- plot(enc_small, return_data = TRUE)
    for (sp in s$Species) {
        vals <- pd[[2]][pd$Species == sp]
        expect_equal(s$Max[s$Species == sp], max(vals))
        expect_equal(s$Min[s$Species == sp], min(vals))
    }

    # `all.sizes = TRUE` gives the whole grid back, which for a rate that grows
    # with size means a strictly larger maximum.
    all_sizes <- summary(enc_small, all.sizes = TRUE)$per_species
    expect_equal(all_sizes$Max, unname(apply(unclass(enc_small), 1, max)))
    # Never smaller, and strictly larger for a species whose `w_max` is below
    # the top of the grid, which is where a rate that grows with size peaks.
    expect_true(all(all_sizes$Max >= s$Max))
    expect_true(any(all_sizes$Max > s$Max))
})

test_that("summary.ArraySpeciesBySize reports NA for an empty species", {
    # Every size class masked out must give NA, not the -Inf/Inf that `min()`
    # and `max()` of an empty selection return, with a warning.
    x <- enc_small
    x[1, ] <- NA
    s <- expect_no_warning(summary(x, all.sizes = TRUE)$per_species)
    expect_true(is.na(s$Min[1]) && is.na(s$Mean[1]) && is.na(s$Max[1]))
    expect_false(anyNA(s$Min[-1]))
})

test_that("summary.ArraySpeciesBySize keeps species the model does not know", {
    # An array whose rows are not model species has no range to be outside of.
    x <- ArraySpeciesBySize(matrix(1:6, nrow = 2,
                                   dimnames = list(c("a", "b"), NULL)))
    expect_identical(summary(x)$per_species$Max, c(5L, 6L))
})

test_that("str.ArraySpeciesBySize works", {
    expect_output(str(enc_small), "ArraySpeciesBySize")
    expect_output(str(enc_small), "Encounter rate")
    expect_output(str(enc_small), "g/year")
    expect_output(str(enc_small), "params")

    out <- capture.output(str(enc_small))
    expect_false(any(grepl("intake_max", out)))
})

test_that("plot.ArraySpeciesBySize returns ggplot", {
    # params attribute is used for styling
    p <- plot(enc_small)
    expect_s3_class(p, "ggplot")

    # With species selection
    p2 <- plot(enc_small, species = c("Cod", "Herring"))
    expect_s3_class(p2, "ggplot")

    # return_data works
    df <- plot(enc_small, return_data = TRUE)
    expect_true(is.data.frame(df))
    expect_true(all(c("w", "Species") %in% names(df)))
    # 2nd column is the named value (e.g. "Encounter rate")
    expect_true(is.numeric(df[[2]]))
})

test_that("length-axis array plots only transform number densities", {
    density <- initialN(NS_params_small)
    density_w <- plot(density, size_axis = "w", return_data = TRUE)
    density_l <- plot(density, size_axis = "l", return_data = TRUE)
    sp_idx <- match(density_l$Species,
                    NS_params_small@species_params$species)
    jacobian <- NS_params_small@species_params$b[sp_idx] * density_w$w /
        density_l$l
    expect_equal(density_l[[2]], unname(density_w[[2]] * jacobian))
    expect_identical(plot(density, size_axis = "l")$scales$get_scales("y")$name,
                     "Number density [1/cm]")

    rate_w <- plot(enc_small, size_axis = "w", return_data = TRUE)
    rate_l <- plot(enc_small, size_axis = "l", return_data = TRUE)
    expect_equal(rate_l[[2]], rate_w[[2]])
})

test_that("ArraySpeciesBySize records the type of its values", {
    mat <- matrix(1:30, nrow = 3, ncol = 10,
                  dimnames = list(NS_params_small@species_params$species,
                                  signif(NS_params_small@w[1:10], 3)))

    # An explicit type is used as given
    for (type in array_types) {
        expect_identical(array_type(ArraySpeciesBySize(mat, type = type)), type)
    }
    expect_error(ArraySpeciesBySize(mat, type = "rate"), "must be one of")
    expect_error(ArraySpeciesBySize(mat, type = NA), "must be one of")
    expect_error(ArraySpeciesBySize(mat, type = c("density", "value")),
                 "must be one of")

    # Without one, a density is recognised from the other metadata, so that
    # arrays built before the attribute existed keep behaving as they did
    expect_identical(array_type(
        ArraySpeciesBySize(mat, value_name = "Number density")), "density")
    expect_identical(array_type(
        ArraySpeciesBySize(mat, value_name = "Resource capacity",
                           units = "1/g")), "density")
    expect_identical(array_type(
        ArraySpeciesBySize(mat, value_name = "Encounter rate",
                           units = "g/year")), "value")

    # An explicit type overrides that guess in both directions
    expect_identical(array_type(
        ArraySpeciesBySize(mat, value_name = "Number density",
                           type = "value")), "value")
    expect_identical(array_type(
        ArraySpeciesBySize(mat, value_name = "Flux", units = "1/year",
                           type = "density")), "density")

    # An array from before the attribute existed carries no attribute at all
    legacy <- ArraySpeciesBySize(mat, value_name = "Number density")
    attr(legacy, "type") <- NULL
    expect_identical(array_type(legacy), "density")

    # Subsetting preserves the type, arithmetic strips it with the class
    tagged <- ArraySpeciesBySize(mat, type = "proportion")
    expect_identical(array_type(tagged[1:2, ]), "proportion")
    expect_null(attr(tagged * 2, "type"))
})

test_that("only a density array is converted onto a length axis", {
    mat <- matrix(1, nrow = 3, ncol = 10,
                  dimnames = list(NS_params_small@species_params$species,
                                  signif(NS_params_small@w[1:10], 3)))
    # Mizer arrays live on the weight grid, so a density in one is always a
    # density with respect to weight
    expect_identical(array_density_wrt(ArraySpeciesBySize(mat,
                                                          type = "density")),
                     "w")
    for (type in c("value", "proportion")) {
        expect_identical(array_density_wrt(ArraySpeciesBySize(mat,
                                                              type = type)),
                         NA_character_)
    }
})

test_that("array_units restates the units on the plotted axis", {
    mat <- matrix(1, nrow = 3, ncol = 10,
                  dimnames = list(NS_params_small@species_params$species,
                                  signif(NS_params_small@w[1:10], 3)))
    density <- ArraySpeciesBySize(mat, value_name = "Number density",
                                  units = "1/g", type = "density")
    expect_identical(array_units(density, "w"), "1/g")
    expect_identical(array_units(density, "l"), "1/cm")

    # Compound units keep the rest of the unit string
    rate_density <- ArraySpeciesBySize(mat, units = "g^-1/year",
                                       type = "density")
    expect_identical(array_units(rate_density, "l"), "cm^-1/year")

    # Values that are not a density keep their units on either axis
    rate <- ArraySpeciesBySize(mat, units = "g/year")
    expect_identical(array_units(rate, "l"), "g/year")
    proportion <- ArraySpeciesBySize(mat, units = "", type = "proportion")
    expect_identical(array_units(proportion, "l"), "")
})

test_that("array plots use the declared type", {
    density <- initialN(NS_params_small)
    sp <- NS_params_small@species_params
    by_weight <- plot(density, size_axis = "w", return_data = TRUE)
    by_length <- plot(density, size_axis = "l", return_data = TRUE)
    sp_idx <- match(by_length$Species, sp$species)
    jacobian <- unname(sp$b[sp_idx]) * by_weight$w / by_length$l
    expect_equal(by_length[[2]], by_weight[[2]] * jacobian)

    # Declaring the same values not to be a density leaves them alone
    not_density <- density
    attr(not_density, "type") <- "value"
    expect_equal(plot(not_density, size_axis = "l", return_data = TRUE)[[2]],
                 by_weight[[2]])
})

test_that("plot() can express a density per logarithmic size", {
    density <- initialN(NS_params_small)
    sp <- NS_params_small@species_params
    by_weight <- plot(density, size_axis = "w", return_data = TRUE)

    # Per log weight is the density times the weight
    per_log <- plot(density, per_log_size = TRUE, return_data = TRUE)
    expect_identical(names(per_log)[[1]], "w")
    expect_equal(per_log[[2]], by_weight[[2]] * by_weight$w)

    # Per log length is the density per length times the length
    per_log_l <- plot(density, size_axis = "l", per_log_size = TRUE,
                      return_data = TRUE)
    by_length <- plot(density, size_axis = "l", return_data = TRUE)
    expect_identical(names(per_log_l)[[1]], "l")
    expect_equal(per_log_l[[2]], by_length[[2]] * by_length$l)

    # per_log_size = FALSE is the density itself, the same as not asking
    expect_equal(plot(density, per_log_size = FALSE, return_data = TRUE)[[2]],
                 by_weight[[2]])

    # The label says which quantity is shown, since the units no longer do
    expect_identical(plot(density, per_log_size = TRUE)$scales$
                         get_scales("y")$name,
                     "Number density in log weight")
    expect_identical(plot(density, size_axis = "l", per_log_size = TRUE)$scales$
                         get_scales("y")$name,
                     "Number density in log length")

    # Asking for it on something that is not a density is an error, where it
    # used to be swallowed by `...`
    expect_error(plot(enc_small, per_log_size = TRUE),
                 "only applies to an array that holds a density")
    expect_error(plot(density, per_log_size = "yes"), "not a flag")

    # No weight-length relationship is needed on a weight axis, so the total
    # survives there
    with_total <- plot(density, total = TRUE, per_log_size = TRUE,
                       return_data = TRUE)
    expect_true(any(with_total$Species == "Total"))
})

test_that("the total is shown on a length axis and matches the lines drawn", {
    density <- initialN(NS_params_small)

    on_w <- plot(density, total = TRUE, return_data = TRUE)
    on_l <- plot(density, total = TRUE, size_axis = "l", return_data = TRUE)
    expect_true("Total" %in% on_w$Species)
    expect_true("Total" %in% on_l$Species)

    # The species here share one weight-length relationship, so the total is
    # the plain sum of the converted lines at each length
    total_l <- on_l[on_l$Species == "Total", ]
    lines_l <- on_l[on_l$Species != "Total", ]
    summed <- vapply(total_l$l, function(len) {
        sum(lines_l[[2]][lines_l$l == len])
    }, numeric(1))
    expect_equal(total_l[[2]], summed)

    # The total is the total of everything in the array, so it depends on
    # neither the species selection nor the size trimming
    full <- unname(colSums(unclass(density)))
    expect_equal(on_w[[2]][on_w$Species == "Total"], full)
    for (args in list(list(all.sizes = TRUE), list(species = "Cod"),
                      list(background = FALSE))) {
        d <- do.call(plot, c(list(density, total = TRUE, return_data = TRUE),
                             args))
        expect_equal(d[[2]][d$Species == "Total"], full)
    }
})

test_that("the total on a length axis interpolates onto a common grid", {
    # Give two species their own weight-length relationships so that their
    # length grids no longer coincide
    params <- NS_params_small
    species_params(params)$a <- c(0.008, 0.012, 0.01)
    species_params(params)$b <- c(3.1, 2.9, 3)
    density <- initialN(params)

    on_w <- plot(density, total = TRUE, return_data = TRUE)
    on_l <- plot(density, total = TRUE, size_axis = "l", return_data = TRUE)
    # The weight axis still has one point per size, the length axis needs the
    # union of the three grids
    expect_identical(sum(on_w$Species == "Total"), length(params@w))
    expect_gt(sum(on_l$Species == "Total"), length(params@w))

    # At the largest length only one species reaches, so the total is exactly
    # that species there. The total covers the whole array, so compare it
    # against the untrimmed lines.
    total_l <- on_l[on_l$Species == "Total", ]
    all_l <- plot(density, all.sizes = TRUE, size_axis = "l",
                  return_data = TRUE)
    lines_l <- all_l[all_l$Species != "Total", ]
    at <- max(total_l$l)
    expect_identical(sum(lines_l$l == at), 1L)
    expect_equal(total_l[[2]][total_l$l == at], lines_l[[2]][lines_l$l == at])
})

test_that("a proportion is plotted against the whole of [0, 1]", {
    feeding_level <- getFeedingLevel(NS_params_small)
    expect_identical(array_type(feeding_level), "proportion")
    expect_identical(plot(feeding_level)$scales$get_scales("y")$limits,
                     c(0, 1))

    # A value above 1 is a real feature of the model and stays visible
    above_one <- ArraySpeciesBySize(unclass_rate(feeding_level) * 3,
                                    value_name = "Test proportion",
                                    type = "proportion",
                                    params = NS_params_small)
    limits <- plot(above_one)$scales$get_scales("y")$limits
    expect_identical(limits[[1]], 0)
    # The sizes outside a species' size range are not plotted, so compare
    # against the values that are
    expect_gte(limits[[2]], max(plot(above_one, return_data = TRUE)[[2]]))

    # An explicit limit always wins, and a log axis is left alone
    expect_identical(plot(feeding_level, ylim = c(0.1, 0.5))$scales$
                         get_scales("y")$limits,
                     c(0.1, 0.5))
    expect_identical(plot(feeding_level, ylim = c(NA, 0.5))$scales$
                         get_scales("y")$limits,
                     c(0, 0.5))
    expect_equal(plot(feeding_level, log_y = TRUE)$scales$get_scales("y")$limits,
                 c(NA_real_, NA_real_))

    # Values that are not a proportion are unaffected
    expect_true(all(is.na(plot(enc_small)$scales$get_scales("y")$limits)))
})

test_that("comparisons reject two arrays that hold different types of value", {
    density <- initialN(NS_params_small)
    as_value <- density
    attr(as_value, "type") <- "value"
    # The type decides the Jacobian and the y-axis range, so there is no one
    # pair of axes that can carry both.
    expect_error(plot2(density, as_value), "different types")
    expect_error(plotRelative(density, as_value), "different types")
})

test_that("comparisons still only warn about a differing name or units", {
    enc <- getEncounter(NS_params_small)
    renamed <- enc
    attr(renamed, "value_name") <- "Something else"
    expect_warning(plot2(enc, renamed), "value name")

    rescaled <- enc
    attr(rescaled, "units") <- "kg/year"
    expect_warning(plot2(enc, rescaled), "y units")
})

test_that("plot.ArraySpeciesBySize supports base plot log argument", {
    p_y <- plot(enc_small, log = "y")
    expect_identical(p_y$scales$get_scales("x")$trans$name, "identity")
    expect_identical(p_y$scales$get_scales("y")$trans$name, "log-10")

    p_xy <- plot(enc_small, log = "xy")
    expect_identical(p_xy$scales$get_scales("x")$trans$name, "log-10")
    expect_identical(p_xy$scales$get_scales("y")$trans$name, "log-10")

    p_none <- plot(enc_small, log = "")
    expect_identical(p_none$scales$get_scales("x")$trans$name, "identity")
    expect_identical(p_none$scales$get_scales("y")$trans$name, "identity")

    expect_error(plot(enc_small, log = "z"), "containing only")
})

test_that("plot2.ArraySpeciesBySize compares compatible arrays", {
    p <- plot2(enc_small, enc_small, name1 = "Original", name2 = "Changed",
               species = "Cod", total = TRUE, background = FALSE,
               wlim = c(1, NA), log = "xy")
    expect_s3_class(p, "ggplot")
    expect_identical(levels(p$data$Model), c("Original", "Changed"))
    expect_true(all(p$data$Species %in% c("Cod", "Total")))
    expect_true(all(p$data$w >= 1))
    expect_identical(p$scales$get_scales("x")$trans$name, "log-10")
    expect_identical(p$scales$get_scales("y")$trans$name, "log-10")

    p_none <- plot2(enc_small, enc_small, species = "Cod", log = "")
    expect_identical(p_none$scales$get_scales("x")$trans$name, "identity")
    expect_identical(p_none$scales$get_scales("y")$trans$name, "identity")

    expect_error(plot2(enc_small, getBiomass(NS_sim_small)), "Both objects must be")
    expect_error(plot2(enc_small, enc_small, log = "z"), "containing only")
})

test_that("plotRelative.ArraySpeciesBySize plots symmetric relative difference", {
    enc2 <- enc_small
    enc2[] <- unclass(enc_small) * 2

    p <- plotRelative(enc_small, enc2, species = "Cod", total = TRUE,
                      background = FALSE, wlim = c(1, NA))
    expect_s3_class(p, "ggplot")
    expect_true(all(p$data$Species %in% c("Cod", "Total")))
    expect_true(all(p$data$w >= 1))
    expect_true(all(abs(p$data$rel_diff - 2 / 3) < 1e-12))
    expect_identical(p$scales$get_scales("x")$trans$name, "log-10")

    p_linear <- plotRelative(enc_small, enc2, species = "Cod", log_x = FALSE)
    expect_identical(p_linear$scales$get_scales("x")$trans$name, "identity")

    expect_error(plotRelative(enc_small, getBiomass(NS_sim_small)), "Both objects must be")
})

test_that("addPlot.ArraySpeciesBySize adds lines to an existing ggplot", {
    p <- plot(enc_small, species = "Cod")
    p2 <- addPlot(p, enc_small, species = "Cod", linetype = "dashed", alpha = 0.5)

    expect_s3_class(p2, "ggplot")
    expect_equal(length(p2$layers), length(p$layers) + 1)
    expect_identical(p2$layers[[length(p2$layers)]]$aes_params$linetype,
                     "dashed")
    expect_identical(p2$layers[[length(p2$layers)]]$aes_params$alpha,
                     0.5)
    p2$layers[[1]]$aes_params$alpha <- 0.25
    expect_null(p$layers[[1]]$aes_params$alpha)
    expect_error(addPlot("not a plot", enc_small), "ggplot")
    expect_error(addPlot(plot(getBiomass(NS_sim_small)), enc_small), "x variable `w`")
    expect_warning(addPlot(p, pred_mort_small, species = "Cod"), "y units")
})

test_that("plot.ArraySpeciesBySize supports full size grid", {
    df <- plot(pred_rate_small, return_data = TRUE, all.sizes = TRUE)
    expect_equal(sort(unique(df$w)), NS_params_small@w_full)

    p <- plot(pred_rate_small)
    expect_s3_class(p, "ggplot")
})

test_that("plot.ArraySpeciesBySize errors for unknown size grid", {
    mat <- matrix(1, nrow = nrow(NS_params_small@initial_n), ncol = 3,
                  dimnames = list(sp = dimnames(NS_params_small@initial_n)$sp,
                                  w = as.character(1:3)))
    x <- ArraySpeciesBySize(mat, params = NS_params_small)

    expect_error(
        plot(x, return_data = TRUE),
        "Can not determine the size grid"
    )
})

test_that("ArraySpeciesBySize has interactive plotly methods", {
    expect_s3_class(plotHover(enc_small), "plotly")
})

test_that("as.data.frame.ArraySpeciesBySize works", {
    df <- as.data.frame(enc_small)
    expect_true(is.data.frame(df))
    expect_true(all(c("w", "value", "Species") %in% names(df)))
    expect_equal(nrow(df), prod(dim(enc_small)))

    df_pr <- as.data.frame(pred_rate_small)
    expect_equal(sort(unique(df_pr$w)), NS_params_small@w_full)
})

test_that("ArraySpeciesBySize subsetting preserves class for 2D", {
    # Subsetting rows keeps class
    sub <- enc_small[1:3, ]
    expect_true(is.ArraySpeciesBySize(sub))
    expect_identical(nrow(sub), 3L)

    # Subsetting to a single row with drop = TRUE returns vector
    sub1 <- enc_small[1, ]
    expect_false(is.ArraySpeciesBySize(sub1))
    expect_true(is.numeric(sub1))

    # Subsetting to a single row with drop = FALSE keeps matrix
    sub1_nodrop <- enc_small[1, , drop = FALSE]
    expect_true(is.ArraySpeciesBySize(sub1_nodrop))
})

test_that("ArraySpeciesBySize arithmetic strips class", {
    # Arithmetic should strip ArraySpeciesBySize and return a plain matrix
    double_enc <- enc_small * 2
    expect_false(is.ArraySpeciesBySize(double_enc))
    expect_true(is.matrix(double_enc))
    expect_equal(double_enc, unclass(enc_small) * 2, ignore_attr = TRUE)

    # Addition with a matrix should work
    mat <- matrix(1, nrow = nrow(enc_small), ncol = ncol(enc_small))
    result <- enc_small + mat
    expect_false(is.ArraySpeciesBySize(result))
    expect_true(is.matrix(result))

    # Comparison operators should work
    expect_true(is.logical(enc_small > 0))
})

test_that("ArraySpeciesBySize value_name attribute", {
    expect_identical(attr(enc_small, "value_name"), "Encounter rate")

    fl <- getFeedingLevel(NS_params_small)
    expect_identical(attr(fl, "value_name"), "Feeding level")

    g <- getEGrowth(NS_params_small)
    expect_identical(attr(g, "value_name"), "Growth rate")

    mort <- getMort(NS_params_small)
    expect_identical(attr(mort, "value_name"), "Total mortality")
})

# Second-order plotting placement (#382) -----------------------------------

test_that("producers tag the representation honestly", {
    # Bin-averaged sinks are tagged "average"; flux/point quantities "point".
    expect_identical(attr(getPredMort(NS_params_small), "representation"),
                     "average")
    expect_identical(attr(getMort(NS_params_small), "representation"),
                     "average")
    expect_identical(attr(getFMort(NS_params_small), "representation"),
                     "average")
    expect_identical(attr(ext_mort(NS_params_small), "representation"),
                     "average")
    expect_identical(attr(getEncounter(NS_params_small), "representation"),
                     "point")
    expect_identical(attr(getEGrowth(NS_params_small), "representation"),
                     "point")
})

test_that("bin_midpoints are the geometric bin centres", {
    p <- NS_params_small
    beta <- p@w_full[2] / p@w_full[1]
    expect_equal(bin_midpoints(p), p@w * sqrt(beta))
    expect_equal(bin_midpoints(p, w_full = TRUE), p@w_full * sqrt(beta))
})

test_that("the size axis honours representation only under second_order_w", {
    p <- NS_params_small
    # Default first-order model: averaged quantities still plot at the nodes,
    # so existing plots are unchanged.
    expect_equal(get_ArraySpeciesBySize_w(getPredMort(p)), p@w)
    expect_equal(get_ArraySpeciesBySize_w(getEncounter(p)), p@w)

    second_order_w(p) <- c(bin_average = TRUE)
    # Now averaged quantities move to the geometric bin centres ...
    expect_equal(get_ArraySpeciesBySize_w(getPredMort(p)), bin_midpoints(p))
    # ... while point-valued quantities stay on the nodes.
    expect_equal(get_ArraySpeciesBySize_w(getEncounter(p)), p@w)
})

test_that("as.data.frame x-values follow the representation tag", {
    p <- NS_params_small
    second_order_w(p) <- c(bin_average = TRUE)
    df_avg <- as.data.frame(getPredMort(p))
    df_pt  <- as.data.frame(getEncounter(p))
    expect_equal(sort(unique(df_avg$w)), sort(unname(bin_midpoints(p))))
    expect_equal(sort(unique(df_pt$w)),  sort(unname(p@w)))
})

test_that("subsetting preserves the representation tag", {
    p <- NS_params_small
    second_order_w(p) <- c(bin_average = TRUE)
    mort <- getPredMort(p)
    sub <- mort[1:2, ]
    expect_identical(attr(sub, "representation"), "average")
    expect_equal(get_ArraySpeciesBySize_w(sub), bin_midpoints(p))
})

# parsePlotLog ----

test_that("parsePlotLog passes the defaults through when `log` is NULL", {
    expect_identical(parsePlotLog(NULL), list(log_x = FALSE, log_y = FALSE))
    expect_identical(parsePlotLog(NULL, log_x = TRUE),
                     list(log_x = TRUE, log_y = FALSE))
    expect_identical(parsePlotLog(NULL, log_x = TRUE, log_y = TRUE),
                     list(log_x = TRUE, log_y = TRUE))
})

test_that("parsePlotLog treats a logical `log` as the legacy y-axis toggle", {
    expect_identical(parsePlotLog(TRUE), list(log_x = FALSE, log_y = TRUE))
    expect_identical(parsePlotLog(FALSE), list(log_x = FALSE, log_y = FALSE))
    # the legacy form overrides the defaults, including log_x
    expect_identical(parsePlotLog(TRUE, log_x = TRUE, log_y = FALSE),
                     list(log_x = FALSE, log_y = TRUE))
})

test_that("parsePlotLog selects axes from a character `log`", {
    expect_identical(parsePlotLog("x"), list(log_x = TRUE, log_y = FALSE))
    expect_identical(parsePlotLog("y"), list(log_x = FALSE, log_y = TRUE))
    expect_identical(parsePlotLog("xy"), list(log_x = TRUE, log_y = TRUE))
    expect_identical(parsePlotLog("yx"), list(log_x = TRUE, log_y = TRUE))
    # an empty string switches both axes off, whatever the defaults were
    expect_identical(parsePlotLog("", log_x = TRUE, log_y = TRUE),
                     list(log_x = FALSE, log_y = FALSE))
})

test_that("parsePlotLog rejects malformed `log` arguments", {
    msg <- "`log` must be a single logical value or a character string"
    expect_error(parsePlotLog(NA), msg)
    expect_error(parsePlotLog(c(TRUE, FALSE)), msg)
    expect_error(parsePlotLog(NA_character_), msg)
    expect_error(parsePlotLog(c("x", "y")), msg)
    expect_error(parsePlotLog("z"), msg)
    expect_error(parsePlotLog("xz"), msg)
    expect_error(parsePlotLog(1), msg)
})

# apply_wlim ----

test_that("apply_wlim restricts the data to the weight limits", {
    data <- data.frame(w = c(1, 2, 5, 10, 20), value = 1:5)
    expect_identical(apply_wlim(data, c(2, 10)), data[2:4, ])
    # the limits are inclusive
    expect_identical(apply_wlim(data, c(1, 20)), data)
})

test_that("apply_wlim leaves a side unrestricted when its limit is NA", {
    data <- data.frame(w = c(1, 2, 5, 10, 20), value = 1:5)
    expect_identical(apply_wlim(data, c(NA, NA)), data)
    expect_identical(apply_wlim(data, c(5, NA)), data[3:5, ])
    expect_identical(apply_wlim(data, c(NA, 5)), data[1:3, ])
})

test_that("apply_wlim can return no rows at all", {
    data <- data.frame(w = c(1, 2, 5), value = 1:3)
    expect_identical(nrow(apply_wlim(data, c(100, 200))), 0L)
})

# Comparison plots across two models ----

test_that("plot2 transforms each array with its own weight-length parameters", {
    params2 <- NS_params_small
    # Halving `a` doubles the length a given weight corresponds to.
    given_species_params(params2)$a <- species_params(NS_params_small)$a / 2
    enc1 <- getEncounter(NS_params_small)
    enc2 <- getEncounter(params2)

    dat1 <- ArraySpeciesBySize_plot_data(enc1, size_axis = "l")
    dat2 <- ArraySpeciesBySize_plot_data(enc2, size_axis = "l")
    # The second model really does put the same weights at other lengths.
    expect_false(isTRUE(all.equal(dat1$l, dat2$l)))
    b <- species_params(NS_params_small)$b[[1]]
    expect_equal(dat2$l, dat1$l * 2^(1 / b))

    p <- plot2(enc1, enc2, size_axis = "l")
    lengths <- p$data$l[p$data$Model == "Second"]
    expect_equal(sort(unique(lengths)), sort(unique(dat2$l)))
})

test_that("plotRelative interpolates when the two length grids differ", {
    params2 <- NS_params_small
    given_species_params(params2)$a <- species_params(NS_params_small)$a / 2
    enc1 <- getEncounter(NS_params_small)
    enc2 <- getEncounter(params2)

    # The two grids share almost no coordinate, so matching by equality would
    # leave essentially nothing to compare.
    rel <- plotRelative(enc1, enc2, size_axis = "l")
    expect_gt(nrow(rel$data), 10)
    expect_true(all(is.finite(rel$data$rel_diff)))
    # Every retained coordinate lies inside the range both series cover.
    dat1 <- ArraySpeciesBySize_plot_data(enc1, size_axis = "l")
    dat2 <- ArraySpeciesBySize_plot_data(enc2, size_axis = "l")
    for (sp in unique(rel$data$Species)) {
        shown <- rel$data$l[rel$data$Species == sp]
        l1 <- dat1$l[dat1$Species == sp]
        l2 <- dat2$l[dat2$Species == sp]
        expect_true(all(shown >= max(min(l1), min(l2))))
        expect_true(all(shown <= min(max(l1), max(l2))))
    }
})

test_that("plotRelative reproduces the exact values when the grids match", {
    enc1 <- getEncounter(NS_params_small)
    enc2 <- getEncounter(NS_params_small) * 2
    enc2 <- ArraySpeciesBySize(enc2, value_name = attr(enc1, "value_name"),
                               units = attr(enc1, "units"),
                               params = NS_params_small)
    rel <- plotRelative(enc1, enc2)
    # 2 (2N - N) / (N + 2N) = 2/3 everywhere.
    expect_equal(unique(round(rel$data$rel_diff, 12)), 2 / 3)
})

test_that("plotRelative applies the density Jacobian of each model", {
    params2 <- NS_params_small
    # A different exponent gives a different dw/dl, which no longer cancels.
    given_species_params(params2)$b <- species_params(NS_params_small)$b + 0.2
    n1 <- initialN(NS_params_small)
    n2 <- initialN(params2)
    rel_w <- plotRelative(n1, n2, size_axis = "w")
    rel_l <- plotRelative(n1, n2, size_axis = "l")
    expect_true(all(abs(rel_w$data$rel_diff) < 1e-12))
    expect_false(all(abs(rel_l$data$rel_diff) < 1e-12))
})

test_that("the comparison methods honour highlight", {
    enc <- getEncounter(NS_params_small)
    doubled <- ArraySpeciesBySize(unclass_rate(enc) * 2,
                                  value_name = attr(enc, "value_name"),
                                  units = attr(enc, "units"),
                                  params = NS_params_small)

    p2 <- plot2(enc, doubled, highlight = "Cod")
    expect_gt(drawn_linewidth(p2, "Cod", NS_params_small),
              drawn_linewidth(p2, "Sprat", NS_params_small))

    # The relative plot draws its zero line first, so the values are layer 2.
    pr <- plotRelative(enc, doubled, highlight = "Cod")
    expect_gt(drawn_linewidth(pr, "Cod", NS_params_small, layer = 2),
              drawn_linewidth(pr, "Sprat", NS_params_small, layer = 2))
})

test_that("the array plots all share one preparation", {
    enc <- getEncounter(NS_params_small)
    prepared <- ArraySpeciesBySize_plot_data(
        enc, species = "Cod", wlim = c(1, NA), total = TRUE, size_axis = "l")
    from_plot <- plot(enc, species = "Cod", wlim = c(1, NA), total = TRUE,
                      size_axis = "l", return_data = TRUE)
    expect_equal(from_plot, prepared)

    # The comparison of an array with itself is that same data, twice.
    p2 <- plot2(enc, enc, species = "Cod", wlim = c(1, NA), total = TRUE,
                size_axis = "l")
    expect_equal(sum(p2$data$Model == "First"), nrow(prepared))
    shown <- p2$data[p2$data$Model == "First", names(prepared)]
    # The renderer turns `Legend` into a factor to order the legend.
    shown$Legend <- as.character(shown$Legend)
    expect_equal(shown, prepared, ignore_attr = TRUE)
})
