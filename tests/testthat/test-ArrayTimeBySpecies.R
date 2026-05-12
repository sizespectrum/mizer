test_that("ArrayTimeBySpecies constructor works", {
    mat <- matrix(1:120, nrow = 10, ncol = 12)
    rownames(mat) <- as.character(1:10)
    colnames(mat) <- NS_params@species_params$species

    bio <- ArrayTimeBySpecies(mat, value_name = "Biomass", units = "g")

    expect_s3_class(bio, "ArrayTimeBySpecies")
    expect_true(is.ArrayTimeBySpecies(bio))
    expect_false(is.ArrayTimeBySpecies(mat))
    expect_true(is.matrix(bio))
    expect_identical(dim(bio), dim(mat))
    expect_identical(attr(bio, "value_name"), "Biomass")
    expect_identical(attr(bio, "units"), "g")
})

test_that("ArrayTimeBySpecies constructor validates input", {
    expect_error(ArrayTimeBySpecies(1:10), "`x` must be a matrix.")
})

test_that("Sim accessor functions return ArrayTimeBySpecies", {
    expect_true(is.ArrayTimeBySpecies(getBiomass(NS_sim)))
    expect_true(is.ArrayTimeBySpecies(getSSB(NS_sim)))
    expect_true(is.ArrayTimeBySpecies(getN(NS_sim)))
    expect_true(is.ArrayTimeBySpecies(getYield(NS_sim)))
})

test_that("Sim accessor functions have correct dimnames", {
    bio <- getBiomass(NS_sim)
    expect_identical(colnames(bio), NS_params@species_params$species)
    expect_identical(nrow(bio), length(NS_sim@n[, 1, 1]))
})

test_that("print.ArrayTimeBySpecies works", {
    bio <- getBiomass(NS_sim)
    expect_output(print(bio), "Biomass")
    expect_output(print(bio), "times x")
    expect_output(print(bio), "g")
})

test_that("summary.ArrayTimeBySpecies works", {
    bio <- getBiomass(NS_sim)
    s <- summary(bio)
    expect_s3_class(s, "summary.ArrayTimeBySpecies")
    expect_identical(s$value_name, "Biomass")
    expect_identical(nrow(s$per_species), nrow(NS_params@species_params))
    expect_output(print(s), "Biomass")
    expect_output(print(s), "times x")
})

test_that("plot.ArrayTimeBySpecies returns ggplot", {
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

test_that("plot.ArrayTimeBySpecies supports base plot log argument", {
    bio <- getBiomass(NS_sim)

    p_xy <- plot(bio, log = "xy")
    expect_identical(p_xy$scales$get_scales("x")$trans$name, "log-10")
    expect_identical(p_xy$scales$get_scales("y")$trans$name, "log-10")

    p_none <- plot(bio, log = "")
    expect_identical(p_none$scales$get_scales("x")$trans$name, "identity")
    expect_identical(p_none$scales$get_scales("y")$trans$name, "identity")

    expect_error(plot(bio, log = TRUE), "`log` must be a character string")
})

test_that("plot2.ArrayTimeBySpecies compares compatible arrays", {
    bio <- getBiomass(NS_sim)
    years <- as.numeric(rownames(bio))

    p <- plot2(bio, bio, name1 = "Original", name2 = "Changed",
               species = "Cod", total = TRUE, start_time = years[2],
               end_time = years[5], log = "xy")
    expect_s3_class(p, "ggplot")
    expect_identical(levels(p$data$Model), c("Original", "Changed"))
    expect_true(all(p$data$Species %in% c("Cod", "Total")))
    expect_true(all(p$data$Year >= years[2]))
    expect_true(all(p$data$Year <= years[5]))
    expect_identical(p$scales$get_scales("x")$trans$name, "log-10")
    expect_identical(p$scales$get_scales("y")$trans$name, "log-10")

    p_none <- plot2(bio, bio, species = "Cod", log = "")
    expect_identical(p_none$scales$get_scales("x")$trans$name, "identity")
    expect_identical(p_none$scales$get_scales("y")$trans$name, "identity")

    warnings <- character()
    withCallingHandlers(
        plot2(bio, getYield(NS_sim), species = "Cod"),
        warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    expect_true(any(grepl("value name", warnings)))
    expect_true(any(grepl("y units", warnings)))
    expect_error(plot2(bio, getEncounter(NS_params)), "Both objects must be")
})

test_that("plotRelative.ArrayTimeBySpecies plots symmetric relative difference", {
    bio <- getBiomass(NS_sim)
    bio2 <- bio
    bio2[] <- unclass(bio) * 2
    years <- as.numeric(rownames(bio))

    p <- plotRelative(bio, bio2, species = "Cod", total = TRUE,
                      start_time = years[2], end_time = years[5])
    expect_s3_class(p, "ggplot")
    expect_true(all(p$data$Species %in% c("Cod", "Total")))
    expect_true(all(p$data$Year >= years[2]))
    expect_true(all(p$data$Year <= years[5]))
    expect_true(all(abs(p$data$rel_diff - 2 / 3) < 1e-12))
    expect_identical(p$scales$get_scales("x")$trans$name, "identity")

    p_log <- plotRelative(bio, bio2, species = "Cod", log_x = TRUE)
    expect_identical(p_log$scales$get_scales("x")$trans$name, "log-10")

    expect_error(plotRelative(bio, getEncounter(NS_params)), "Both objects must be")
})

test_that("addPlot.ArrayTimeBySpecies adds lines to an existing ggplot", {
    bio <- getBiomass(NS_sim)
    yield <- getYield(NS_sim)

    p <- plot(bio, species = "Cod")
    p2 <- addPlot(p, bio, species = "Cod", linetype = "dashed", alpha = 0.5)

    expect_s3_class(p2, "ggplot")
    expect_equal(length(p2$layers), length(p$layers) + 1)
    expect_identical(p2$layers[[length(p2$layers)]]$aes_params$linetype,
                     "dashed")
    expect_identical(p2$layers[[length(p2$layers)]]$aes_params$alpha,
                     0.5)
    p2$layers[[1]]$aes_params$alpha <- 0.25
    expect_null(p$layers[[1]]$aes_params$alpha)
    expect_error(addPlot("not a plot", bio), "ggplot")
    expect_error(addPlot(plot(getEncounter(NS_params)), bio), "x variable `Year`")

    warnings <- character()
    withCallingHandlers(
        addPlot(p, yield, species = "Cod"),
        warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    expect_true(any(grepl("y variable", warnings)))
    expect_true(any(grepl("y units", warnings)))
})

test_that("ArrayTimeBySpecies has interactive plotly methods", {
    bio <- getBiomass(NS_sim)

    expect_s3_class(ggplotly(bio), "plotly")
})

test_that("plot.ArrayTimeBySpecies time filtering works", {
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

test_that("plot.ArrayTimeBySpecies total works", {
    bio <- getBiomass(NS_sim)
    df <- plot(bio, total = TRUE, return_data = TRUE)
    expect_true("Total" %in% df$Species)
})

test_that("plot.ArrayTimeBySpecies errors on unknown species", {
    bio <- getBiomass(NS_sim)
    expect_error(plot(bio, species = "NotASpecies"),
                 "None of the selected species are in the array.")
})

test_that("as.data.frame.ArrayTimeBySpecies works", {
    bio <- getBiomass(NS_sim)
    df <- as.data.frame(bio)
    expect_true(is.data.frame(df))
    expect_true(all(c("time", "value", "Species") %in% names(df)))
    expect_equal(nrow(df), prod(dim(bio)))
    expect_identical(sort(unique(df$Species)),
                     sort(NS_params@species_params$species))
})

test_that("ArrayTimeBySpecies subsetting preserves class for 2D", {
    bio <- getBiomass(NS_sim)

    # Row subsetting keeps class
    sub <- bio[1:5, ]
    expect_true(is.ArrayTimeBySpecies(sub))
    expect_identical(nrow(sub), 5L)

    # Column subsetting keeps class
    sub_col <- bio[, 1:3]
    expect_true(is.ArrayTimeBySpecies(sub_col))
    expect_identical(ncol(sub_col), 3L)

    # Subsetting to single column with drop = TRUE returns vector
    sub1 <- bio[, 1]
    expect_false(is.ArrayTimeBySpecies(sub1))
    expect_true(is.numeric(sub1))

    # Subsetting to single column with drop = FALSE keeps matrix
    sub1_nodrop <- bio[, 1, drop = FALSE]
    expect_true(is.ArrayTimeBySpecies(sub1_nodrop))
})

test_that("ArrayTimeBySpecies arithmetic strips class", {
    bio <- getBiomass(NS_sim)

    # Multiplication strips class
    doubled <- bio * 2
    expect_false(is.ArrayTimeBySpecies(doubled))
    expect_true(is.matrix(doubled))
    expect_equal(doubled, unclass(bio) * 2, ignore_attr = TRUE)

    # Addition strips class
    mat <- matrix(1, nrow = nrow(bio), ncol = ncol(bio))
    result <- bio + mat
    expect_false(is.ArrayTimeBySpecies(result))
    expect_true(is.matrix(result))

    # Comparison operators work
    expect_true(is.logical(bio > 0))
})

test_that("ArrayTimeBySpecies value_name attributes are set correctly", {
    expect_identical(attr(getBiomass(NS_sim), "value_name"), "Biomass")
    expect_identical(attr(getSSB(NS_sim), "value_name"), "Spawning stock biomass")
    expect_identical(attr(getYield(NS_sim), "value_name"), "Yield rate")
    expect_identical(attr(getN(NS_sim), "value_name"), "Abundance")
})
