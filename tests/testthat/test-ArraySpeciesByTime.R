test_that("ArraySpeciesByTime constructor works", {
    mat <- matrix(1:120, nrow = 10, ncol = 12)
    rownames(mat) <- as.character(1:10)
    colnames(mat) <- NS_params@species_params$species

    bio <- ArraySpeciesByTime(mat, value_name = "Biomass", units = "g")

    expect_s3_class(bio, "ArraySpeciesByTime")
    expect_true(is.ArraySpeciesByTime(bio))
    expect_false(is.ArraySpeciesByTime(mat))
    expect_true(is.matrix(bio))
    expect_identical(dim(bio), dim(mat))
    expect_identical(attr(bio, "value_name"), "Biomass")
    expect_identical(attr(bio, "units"), "g")
})

test_that("ArraySpeciesByTime constructor validates input", {
    expect_error(ArraySpeciesByTime(1:10), "`x` must be a matrix.")
})

test_that("Sim accessor functions return ArraySpeciesByTime", {
    expect_true(is.ArraySpeciesByTime(getBiomass(NS_sim)))
    expect_true(is.ArraySpeciesByTime(getSSB(NS_sim)))
    expect_true(is.ArraySpeciesByTime(getN(NS_sim)))
    expect_true(is.ArraySpeciesByTime(getYield(NS_sim)))
})

test_that("Sim accessor functions have correct dimnames", {
    bio <- getBiomass(NS_sim)
    expect_identical(colnames(bio), NS_params@species_params$species)
    expect_identical(nrow(bio), length(NS_sim@n[, 1, 1]))
})

test_that("print.ArraySpeciesByTime works", {
    bio <- getBiomass(NS_sim)
    expect_output(print(bio), "Biomass")
    expect_output(print(bio), "times x")
    expect_output(print(bio), "g")
})

test_that("summary.ArraySpeciesByTime works", {
    bio <- getBiomass(NS_sim)
    s <- summary(bio)
    expect_s3_class(s, "summary.ArraySpeciesByTime")
    expect_identical(s$value_name, "Biomass")
    expect_identical(nrow(s$per_species), nrow(NS_params@species_params))
    expect_output(print(s), "Biomass")
    expect_output(print(s), "times x")
})

test_that("plot.ArraySpeciesByTime returns ggplot", {
    bio <- getBiomass(NS_sim)

    p <- plot(bio)
    expect_s3_class(p, "ggplot")

    # With species selection
    p2 <- plot(bio, species = c("Cod", "Herring"))
    expect_s3_class(p2, "ggplot")

    # return_data works
    df <- plot(bio, return_data = TRUE)
    expect_true(is.data.frame(df))
    expect_true(all(c("Year", "Biomass", "Species") %in% names(df)))
})

test_that("ArraySpeciesByTime has interactive plotly methods", {
    bio <- getBiomass(NS_sim)

    expect_s3_class(ggplotly(bio), "plotly")
})

test_that("plot.ArraySpeciesByTime time filtering works", {
    bio <- getBiomass(NS_sim)
    t <- as.numeric(rownames(bio))

    mid <- median(t)
    df_start <- plot(bio, start_time = mid, return_data = TRUE)
    expect_true(all(df_start$Year >= mid))

    df_end <- plot(bio, end_time = mid, return_data = TRUE)
    expect_true(all(df_end$Year <= mid))

    expect_error(plot(bio, start_time = mid, end_time = mid - 1),
                 "start_time must be less than end_time")
})

test_that("plot.ArraySpeciesByTime total works", {
    bio <- getBiomass(NS_sim)
    df <- plot(bio, total = TRUE, return_data = TRUE)
    expect_true("Total" %in% df$Species)
})

test_that("plot.ArraySpeciesByTime errors on unknown species", {
    bio <- getBiomass(NS_sim)
    expect_error(plot(bio, species = "NotASpecies"),
                 "None of the selected species are in the array.")
})

test_that("as.data.frame.ArraySpeciesByTime works", {
    bio <- getBiomass(NS_sim)
    df <- as.data.frame(bio)
    expect_true(is.data.frame(df))
    expect_true(all(c("time", "value", "Species") %in% names(df)))
    expect_equal(nrow(df), prod(dim(bio)))
    expect_identical(sort(unique(df$Species)),
                     sort(NS_params@species_params$species))
})

test_that("ArraySpeciesByTime subsetting preserves class for 2D", {
    bio <- getBiomass(NS_sim)

    # Row subsetting keeps class
    sub <- bio[1:5, ]
    expect_true(is.ArraySpeciesByTime(sub))
    expect_identical(nrow(sub), 5L)

    # Column subsetting keeps class
    sub_col <- bio[, 1:3]
    expect_true(is.ArraySpeciesByTime(sub_col))
    expect_identical(ncol(sub_col), 3L)

    # Subsetting to single column with drop = TRUE returns vector
    sub1 <- bio[, 1]
    expect_false(is.ArraySpeciesByTime(sub1))
    expect_true(is.numeric(sub1))

    # Subsetting to single column with drop = FALSE keeps matrix
    sub1_nodrop <- bio[, 1, drop = FALSE]
    expect_true(is.ArraySpeciesByTime(sub1_nodrop))
})

test_that("ArraySpeciesByTime arithmetic strips class", {
    bio <- getBiomass(NS_sim)

    # Multiplication strips class
    doubled <- bio * 2
    expect_false(is.ArraySpeciesByTime(doubled))
    expect_true(is.matrix(doubled))
    expect_equal(doubled, unclass(bio) * 2, ignore_attr = TRUE)

    # Addition strips class
    mat <- matrix(1, nrow = nrow(bio), ncol = ncol(bio))
    result <- bio + mat
    expect_false(is.ArraySpeciesByTime(result))
    expect_true(is.matrix(result))

    # Comparison operators work
    expect_true(is.logical(bio > 0))
})

test_that("ArraySpeciesByTime value_name attributes are set correctly", {
    expect_identical(attr(getBiomass(NS_sim), "value_name"), "Biomass")
    expect_identical(attr(getSSB(NS_sim), "value_name"), "Spawning stock biomass")
    expect_identical(attr(getYield(NS_sim), "value_name"), "Yield rate")
    expect_identical(attr(getN(NS_sim), "value_name"), "Abundance")
})
