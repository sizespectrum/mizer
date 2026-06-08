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
    expect_true(is.ArraySpeciesBySize(getMaxIntakeRate(params)))
    expect_true(is.ArraySpeciesBySize(getMetabolicRate(params)))
    expect_true(is.ArraySpeciesBySize(getSearchVolume(params)))
    expect_true(is.ArraySpeciesBySize(getExtMort(params)))
    expect_true(is.ArraySpeciesBySize(getExtEncounter(params)))
    expect_true(is.ArraySpeciesBySize(getMaturityProportion(params)))
    expect_true(is.ArraySpeciesBySize(getReproductionProportion(params)))
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

test_that("summary.ArraySpeciesBySize works", {
    s <- summary(enc_small)
    expect_s3_class(s, "summary.ArraySpeciesBySize")
    expect_identical(s$value_name, "Encounter rate")
    expect_identical(nrow(s$per_species), nrow(NS_params_small@species_params))
    expect_output(print(s), "Encounter rate")
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
