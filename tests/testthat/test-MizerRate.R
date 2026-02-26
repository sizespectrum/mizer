test_that("MizerRate constructor works", {
    mat <- matrix(1:120, nrow = 12, ncol = 10)
    rownames(mat) <- NS_params@species_params$species
    colnames(mat) <- signif(NS_params@w[1:10], 3)
    
    rate <- MizerRate(mat, rate_name = "Test rate", units = "g/year")
    
    expect_s3_class(rate, "MizerRate")
    expect_true(is.MizerRate(rate))
    expect_false(is.MizerRate(mat))
    expect_true(is.matrix(rate))
    expect_identical(dim(rate), dim(mat))
    expect_identical(attr(rate, "rate_name"), "Test rate")
    expect_identical(attr(rate, "units"), "g/year")
})

test_that("MizerRate constructor validates input", {
    expect_error(MizerRate(1:10), "`x` must be a matrix.")
})

test_that("Rate functions return MizerRate", {
    params <- NS_params
    
    expect_true(is.MizerRate(getEncounter(params)))
    expect_true(is.MizerRate(getFeedingLevel(params)))
    expect_true(is.MizerRate(getCriticalFeedingLevel(params)))
    expect_true(is.MizerRate(getEReproAndGrowth(params)))
    expect_true(is.MizerRate(getERepro(params)))
    expect_true(is.MizerRate(getEGrowth(params)))
    expect_true(is.MizerRate(getFMort(params)))
    expect_true(is.MizerRate(getPredMort(params)))
    expect_true(is.MizerRate(getMort(params)))
    expect_true(is.MizerRate(getFlux(params)))
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

test_that("print.MizerRate works", {
    enc <- getEncounter(NS_params)
    expect_output(print(enc), "Encounter rate")
    expect_output(print(enc), "species x")
    expect_output(print(enc), "g/year")
})

test_that("summary.MizerRate works", {
    enc <- getEncounter(NS_params)
    s <- summary(enc)
    expect_s3_class(s, "summary.MizerRate")
    expect_identical(s$rate_name, "Encounter rate")
    expect_identical(nrow(s$per_species), nrow(NS_params@species_params))
    expect_output(print(s), "Encounter rate")
})

test_that("plot.MizerRate returns ggplot", {
    enc <- getEncounter(NS_params)
    
    # With params for styling
    p <- plot(enc, NS_params)
    expect_s3_class(p, "ggplot")
    
    # Without params (basic plot)
    p2 <- plot(enc)
    expect_s3_class(p2, "ggplot")
    
    # With species selection
    p3 <- plot(enc, NS_params, species = c("Cod", "Herring"))
    expect_s3_class(p3, "ggplot")
    
    # return_data works
    df <- plot(enc, NS_params, return_data = TRUE)
    expect_true(is.data.frame(df))
    expect_true(all(c("w", "value", "Species") %in% names(df)))
})

test_that("as.data.frame.MizerRate works", {
    enc <- getEncounter(NS_params)
    df <- as.data.frame(enc)
    expect_true(is.data.frame(df))
    expect_true(all(c("w", "value", "Species") %in% names(df)))
    expect_equal(nrow(df), prod(dim(enc)))
})

test_that("MizerRate subsetting preserves class for 2D", {
    enc <- getEncounter(NS_params)
    
    # Subsetting rows keeps class
    sub <- enc[1:3, ]
    expect_true(is.MizerRate(sub))
    expect_identical(nrow(sub), 3L)
    
    # Subsetting to a single row with drop = TRUE returns vector
    sub1 <- enc[1, ]
    expect_false(is.MizerRate(sub1))
    expect_true(is.numeric(sub1))
    
    # Subsetting to a single row with drop = FALSE keeps matrix
    sub1_nodrop <- enc[1, , drop = FALSE]
    expect_true(is.MizerRate(sub1_nodrop))
})

test_that("MizerRate arithmetic strips class", {
    enc <- getEncounter(NS_params)
    
    # Arithmetic should strip MizerRate and return a plain matrix
    double_enc <- enc * 2
    expect_false(is.MizerRate(double_enc))
    expect_true(is.matrix(double_enc))
    expect_equal(double_enc, unclass(enc) * 2, ignore_attr = TRUE)
    
    # Addition with a matrix should work
    mat <- matrix(1, nrow = nrow(enc), ncol = ncol(enc))
    result <- enc + mat
    expect_false(is.MizerRate(result))
    expect_true(is.matrix(result))
    
    # Comparison operators should work
    expect_true(is.logical(enc > 0))
})

test_that("MizerRate rate_name attribute", {
    enc <- getEncounter(NS_params)
    expect_identical(attr(enc, "rate_name"), "Encounter rate")
    
    fl <- getFeedingLevel(NS_params)
    expect_identical(attr(fl, "rate_name"), "Feeding level")
    
    g <- getEGrowth(NS_params)
    expect_identical(attr(g, "rate_name"), "Growth rate")
    
    mort <- getMort(NS_params)
    expect_identical(attr(mort, "rate_name"), "Total mortality")
})
