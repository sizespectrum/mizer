test_that("ArraySpeciesBySize constructor works", {
    mat <- matrix(1:120, nrow = 12, ncol = 10)
    rownames(mat) <- NS_params@species_params$species
    colnames(mat) <- signif(NS_params@w[1:10], 3)
    
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
    params <- NS_params

    expect_true(is.ArraySpeciesBySize(getEncounter(params)))
    expect_true(is.ArraySpeciesBySize(getFeedingLevel(params)))
    expect_true(is.ArraySpeciesBySize(getCriticalFeedingLevel(params)))
    expect_true(is.ArraySpeciesBySize(getEReproAndGrowth(params)))
    expect_true(is.ArraySpeciesBySize(getERepro(params)))
    expect_true(is.ArraySpeciesBySize(getEGrowth(params)))
    expect_true(is.ArraySpeciesBySize(getFMort(params)))
    expect_true(is.ArraySpeciesBySize(getPredMort(params)))
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
    params <- NS_params
    expected_dimnames <- dimnames(params@metab)
    
    expect_identical(dimnames(getEncounter(params)), expected_dimnames)
    expect_identical(dimnames(getFeedingLevel(params)), expected_dimnames)
    expect_identical(dimnames(getCriticalFeedingLevel(params)), expected_dimnames)
    expect_identical(dimnames(getEReproAndGrowth(params)), expected_dimnames)
    expect_identical(dimnames(getERepro(params)), expected_dimnames)
    expect_identical(dimnames(getEGrowth(params)), expected_dimnames)
    expect_identical(dimnames(getFMort(params)), expected_dimnames)
    expect_identical(dimnames(getPredMort(params)), expected_dimnames)
    expect_identical(dimnames(getMort(params)), expected_dimnames)
    expect_identical(dimnames(getFlux(params)), expected_dimnames)
})

test_that("print.ArraySpeciesBySize works", {
    enc <- getEncounter(NS_params)
    expect_output(print(enc), "Encounter rate")
    expect_output(print(enc), "species x")
    expect_output(print(enc), "g/year")
})

test_that("summary.ArraySpeciesBySize works", {
    enc <- getEncounter(NS_params)
    s <- summary(enc)
    expect_s3_class(s, "summary.ArraySpeciesBySize")
    expect_identical(s$value_name, "Encounter rate")
    expect_identical(nrow(s$per_species), nrow(NS_params@species_params))
    expect_output(print(s), "Encounter rate")
})

test_that("plot.ArraySpeciesBySize returns ggplot", {
    enc <- getEncounter(NS_params)

    # params attribute is used for styling
    p <- plot(enc)
    expect_s3_class(p, "ggplot")

    # With species selection
    p2 <- plot(enc, species = c("Cod", "Herring"))
    expect_s3_class(p2, "ggplot")

    # return_data works
    df <- plot(enc, return_data = TRUE)
    expect_true(is.data.frame(df))
    expect_true(all(c("w", "value", "Species") %in% names(df)))
})

test_that("addPlot.ArraySpeciesBySize adds lines to an existing ggplot", {
    enc <- getEncounter(NS_params)
    pred_mort <- getPredMort(NS_params)

    p <- plot(enc, species = "Cod")
    p2 <- addPlot(p, enc, species = "Cod", linetype = "dashed", alpha = 0.5)

    expect_s3_class(p2, "ggplot")
    expect_equal(length(p2$layers), length(p$layers) + 1)
    expect_identical(p2$layers[[length(p2$layers)]]$aes_params$linetype,
                     "dashed")
    expect_identical(p2$layers[[length(p2$layers)]]$aes_params$alpha,
                     0.5)
    expect_error(addPlot("not a plot", enc), "ggplot")
    expect_error(addPlot(plot(getBiomass(NS_sim)), enc), "x variable `w`")
    expect_warning(addPlot(p, pred_mort, species = "Cod"), "y units")
})

test_that("plot.ArraySpeciesBySize supports full size grid", {
    pred_rate <- getPredRate(NS_params)

    df <- plot(pred_rate, return_data = TRUE, all.sizes = TRUE)
    expect_equal(sort(unique(df$w)), NS_params@w_full)

    p <- plot(pred_rate)
    expect_s3_class(p, "ggplot")
})

test_that("plot.ArraySpeciesBySize errors for unknown size grid", {
    mat <- matrix(1, nrow = nrow(NS_params@initial_n), ncol = 3,
                  dimnames = list(sp = dimnames(NS_params@initial_n)$sp,
                                  w = as.character(1:3)))
    x <- ArraySpeciesBySize(mat, params = NS_params)

    expect_error(
        plot(x, return_data = TRUE),
        "Can not determine the size grid"
    )
})

test_that("ArraySpeciesBySize has interactive plotly methods", {
    enc <- getEncounter(NS_params)

    expect_s3_class(ggplotly(enc), "plotly")
})

test_that("as.data.frame.ArraySpeciesBySize works", {
    enc <- getEncounter(NS_params)
    df <- as.data.frame(enc)
    expect_true(is.data.frame(df))
    expect_true(all(c("w", "value", "Species") %in% names(df)))
    expect_equal(nrow(df), prod(dim(enc)))

    pred_rate <- getPredRate(NS_params)
    df_pr <- as.data.frame(pred_rate)
    expect_equal(sort(unique(df_pr$w)), NS_params@w_full)
})

test_that("ArraySpeciesBySize subsetting preserves class for 2D", {
    enc <- getEncounter(NS_params)
    
    # Subsetting rows keeps class
    sub <- enc[1:3, ]
    expect_true(is.ArraySpeciesBySize(sub))
    expect_identical(nrow(sub), 3L)
    
    # Subsetting to a single row with drop = TRUE returns vector
    sub1 <- enc[1, ]
    expect_false(is.ArraySpeciesBySize(sub1))
    expect_true(is.numeric(sub1))
    
    # Subsetting to a single row with drop = FALSE keeps matrix
    sub1_nodrop <- enc[1, , drop = FALSE]
    expect_true(is.ArraySpeciesBySize(sub1_nodrop))
})

test_that("ArraySpeciesBySize arithmetic strips class", {
    enc <- getEncounter(NS_params)
    
    # Arithmetic should strip ArraySpeciesBySize and return a plain matrix
    double_enc <- enc * 2
    expect_false(is.ArraySpeciesBySize(double_enc))
    expect_true(is.matrix(double_enc))
    expect_equal(double_enc, unclass(enc) * 2, ignore_attr = TRUE)
    
    # Addition with a matrix should work
    mat <- matrix(1, nrow = nrow(enc), ncol = ncol(enc))
    result <- enc + mat
    expect_false(is.ArraySpeciesBySize(result))
    expect_true(is.matrix(result))
    
    # Comparison operators should work
    expect_true(is.logical(enc > 0))
})

test_that("ArraySpeciesBySize value_name attribute", {
    enc <- getEncounter(NS_params)
    expect_identical(attr(enc, "value_name"), "Encounter rate")
    
    fl <- getFeedingLevel(NS_params)
    expect_identical(attr(fl, "value_name"), "Feeding level")
    
    g <- getEGrowth(NS_params)
    expect_identical(attr(g, "value_name"), "Growth rate")
    
    mort <- getMort(NS_params)
    expect_identical(attr(mort, "value_name"), "Total mortality")
})
