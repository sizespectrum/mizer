# Cache sim accessor results once to avoid recomputing them in every test block.
bio_small <- getBiomass(NS_sim_small)
ssb_small <- getSSB(NS_sim_small)
yield_small <- getYield(NS_sim_small)
abundance_small <- getN(NS_sim_small)

test_that("ArrayTimeBySpecies constructor works", {
    mat <- matrix(1:30, nrow = 10, ncol = 3)
    rownames(mat) <- as.character(1:10)
    colnames(mat) <- NS_params_small@species_params$species

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
    expect_true(is.ArrayTimeBySpecies(bio_small))
    expect_true(is.ArrayTimeBySpecies(ssb_small))
    expect_true(is.ArrayTimeBySpecies(abundance_small))
    expect_true(is.ArrayTimeBySpecies(yield_small))
})

test_that("Sim accessor functions have correct dimnames", {
    expect_identical(colnames(bio_small), NS_params_small@species_params$species)
    expect_identical(nrow(bio_small), length(NS_sim_small@n[, 1, 1]))
})

test_that("print.ArrayTimeBySpecies works", {
    expect_output(print(bio_small), "Biomass")
    expect_output(print(bio_small), "times x")
    expect_output(print(bio_small), "g")
})

test_that("print.ArrayTimeBySpecies truncates a long time series", {
    sim_long <- suppressMessages(
        project(NS_params_small, t_max = 60, dt = 0.5, t_save = 0.5,
               progress_bar = FALSE))
    bio <- getBiomass(sim_long)
    expect_identical(nrow(bio), 121L)
    out <- paste(capture.output(print(bio)), collapse = "\n")
    # the earliest and latest time steps should both be visible ...
    expect_match(out, "\\b0\\b")
    expect_match(out, "\\b60\\b")
    # ... via an evenly spaced sample, reported in a trailing footer
    expect_match(out, "showing 8 of 121 times")
    expect_match(out, "evenly spaced")
})

test_that("summary.ArrayTimeBySpecies works", {
    s <- summary(bio_small)
    expect_s3_class(s, "summary.ArrayTimeBySpecies")
    expect_identical(s$value_name, "Biomass")
    expect_identical(nrow(s$per_species), nrow(NS_params_small@species_params))
    expect_output(print(s), "Biomass")
    expect_output(print(s), "times x")
})

test_that("str.ArrayTimeBySpecies works", {
    expect_output(str(bio_small), "ArrayTimeBySpecies")
    expect_output(str(bio_small), "Biomass")
    expect_output(str(bio_small), "g")
    expect_output(str(bio_small), "params")

    out <- capture.output(str(bio_small))
    expect_false(any(grepl("intake_max", out)))
})

test_that("plot.ArrayTimeBySpecies returns ggplot", {
    p <- plot(bio_small)
    expect_s3_class(p, "ggplot")

    # With species selection
    p2 <- plot(bio_small, species = c("Cod", "Herring"))
    expect_s3_class(p2, "ggplot")

    # return_data works
    df <- plot(bio_small, return_data = TRUE)
    expect_true(is.data.frame(df))
    expect_true(all(c("Year", "Biomass", "Species") %in% names(df)))
})

test_that("plot.ArrayTimeBySpecies supports base plot log argument", {
    p_xy <- plot(bio_small, log = "xy")
    expect_identical(p_xy$scales$get_scales("x")$trans$name, "log-10")
    expect_identical(p_xy$scales$get_scales("y")$trans$name, "log-10")

    p_none <- plot(bio_small, log = "")
    expect_identical(p_none$scales$get_scales("x")$trans$name, "identity")
    expect_identical(p_none$scales$get_scales("y")$trans$name, "identity")

    expect_error(plot(bio_small, log = 1), "must be a single logical value or a character string")
})

test_that("plot2.ArrayTimeBySpecies compares compatible arrays", {
    years <- as.numeric(rownames(bio_small))

    p <- plot2(bio_small, bio_small, name1 = "Original", name2 = "Changed",
               species = "Cod", total = TRUE,
               tlim = c(years[2], years[4]), log = "xy")
    expect_s3_class(p, "ggplot")
    expect_identical(levels(p$data$Model), c("Original", "Changed"))
    expect_true(all(p$data$Species %in% c("Cod", "Total")))
    expect_true(all(p$data$Year >= years[2]))
    expect_true(all(p$data$Year <= years[4]))
    expect_identical(p$scales$get_scales("x")$trans$name, "log-10")
    expect_identical(p$scales$get_scales("y")$trans$name, "log-10")

    p_none <- plot2(bio_small, bio_small, species = "Cod", log = "")
    expect_identical(p_none$scales$get_scales("x")$trans$name, "identity")
    expect_identical(p_none$scales$get_scales("y")$trans$name, "identity")

    warnings <- character()
    withCallingHandlers(
        plot2(bio_small, yield_small, species = "Cod"),
        warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    expect_true(any(grepl("value name", warnings)))
    expect_true(any(grepl("y units", warnings)))
    expect_error(plot2(bio_small, getEncounter(NS_params_small)), "Both objects must be")
})

test_that("plotRelative.ArrayTimeBySpecies plots symmetric relative difference", {
    bio2 <- bio_small
    bio2[] <- unclass(bio_small) * 2
    years <- as.numeric(rownames(bio_small))

    p <- plotRelative(bio_small, bio2, species = "Cod", total = TRUE,
                      tlim = c(years[2], years[4]))
    expect_s3_class(p, "ggplot")
    expect_true(all(p$data$Species %in% c("Cod", "Total")))
    expect_true(all(p$data$Year >= years[2]))
    expect_true(all(p$data$Year <= years[4]))
    expect_true(all(abs(p$data$rel_diff - 2 / 3) < 1e-12))
    expect_identical(p$scales$get_scales("x")$trans$name, "identity")

    p_log <- plotRelative(bio_small, bio2, species = "Cod", log_x = TRUE)
    expect_identical(p_log$scales$get_scales("x")$trans$name, "log-10")

    expect_error(plotRelative(bio_small, getEncounter(NS_params_small)), "Both objects must be")
})

test_that("addPlot.ArrayTimeBySpecies adds lines to an existing ggplot", {
    p <- plot(bio_small, species = "Cod")
    p2 <- addPlot(p, bio_small, species = "Cod", linetype = "dashed", alpha = 0.5)

    expect_s3_class(p2, "ggplot")
    expect_equal(length(p2$layers), length(p$layers) + 1)
    expect_identical(p2$layers[[length(p2$layers)]]$aes_params$linetype,
                     "dashed")
    expect_identical(p2$layers[[length(p2$layers)]]$aes_params$alpha,
                     0.5)
    p2$layers[[1]]$aes_params$alpha <- 0.25
    expect_null(p$layers[[1]]$aes_params$alpha)
    expect_error(addPlot("not a plot", bio_small), "ggplot")
    expect_error(addPlot(plot(getEncounter(NS_params_small)), bio_small), "x variable `Year`")

    warnings <- character()
    withCallingHandlers(
        addPlot(p, yield_small, species = "Cod"),
        warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    expect_true(any(grepl("y variable", warnings)))
    expect_true(any(grepl("y units", warnings)))
})

test_that("ArrayTimeBySpecies has interactive plotly methods", {
    expect_s3_class(plotHover(bio_small), "plotly")
})

test_that("plot.ArrayTimeBySpecies time filtering works", {
    t <- as.numeric(rownames(bio_small))

    mid <- median(t)
    df_start <- plot(bio_small, tlim = c(mid, NA), return_data = TRUE)
    expect_true(all(df_start$Year >= mid))

    df_end <- plot(bio_small, tlim = c(NA, mid), return_data = TRUE)
    expect_true(all(df_end$Year <= mid))

    expect_error(plot(bio_small, tlim = c(mid, mid - 1)),
                 "tlim\\[1\\] must be less than tlim\\[2\\]")
})

test_that("plot.ArrayTimeBySpecies total works", {
    df <- plot(bio_small, total = TRUE, return_data = TRUE)
    expect_true("Total" %in% df$Species)
})

test_that("plot.ArrayTimeBySpecies errors on unknown species", {
    expect_error(plot(bio_small, species = "NotASpecies"),
                 "None of the selected species are in the array.")
})

test_that("as.data.frame.ArrayTimeBySpecies works", {
    df <- as.data.frame(bio_small)
    expect_true(is.data.frame(df))
    expect_true(all(c("time", "value", "Species") %in% names(df)))
    expect_equal(nrow(df), prod(dim(bio_small)))
    expect_identical(sort(unique(df$Species)),
                     sort(NS_params_small@species_params$species))
})

test_that("ArrayTimeBySpecies subsetting preserves class for 2D", {
    # Row subsetting keeps class
    sub <- bio_small[1:3, ]
    expect_true(is.ArrayTimeBySpecies(sub))
    expect_identical(nrow(sub), 3L)

    # Column subsetting keeps class
    sub_col <- bio_small[, 1:3]
    expect_true(is.ArrayTimeBySpecies(sub_col))
    expect_identical(ncol(sub_col), 3L)

    # Subsetting to single column with drop = TRUE returns vector
    sub1 <- bio_small[, 1]
    expect_false(is.ArrayTimeBySpecies(sub1))
    expect_true(is.numeric(sub1))

    # Subsetting to single column with drop = FALSE keeps matrix
    sub1_nodrop <- bio_small[, 1, drop = FALSE]
    expect_true(is.ArrayTimeBySpecies(sub1_nodrop))
})

test_that("ArrayTimeBySpecies arithmetic strips class", {
    # Multiplication strips class
    doubled <- bio_small * 2
    expect_false(is.ArrayTimeBySpecies(doubled))
    expect_true(is.matrix(doubled))
    expect_equal(doubled, unclass(bio_small) * 2, ignore_attr = TRUE)

    # Addition strips class
    mat <- matrix(1, nrow = nrow(bio_small), ncol = ncol(bio_small))
    result <- bio_small + mat
    expect_false(is.ArrayTimeBySpecies(result))
    expect_true(is.matrix(result))

    # Comparison operators work
    expect_true(is.logical(bio_small > 0))
})

test_that("ArrayTimeBySpecies value_name attributes are set correctly", {
    expect_identical(attr(bio_small, "value_name"), "Biomass")
    expect_identical(attr(ssb_small, "value_name"), "Spawning stock biomass")
    expect_identical(attr(yield_small, "value_name"), "Yield rate")
    expect_identical(attr(abundance_small, "value_name"), "Abundance")
})

# Background species, selection and non-positive values ----

test_that("a background species is drawn once, under the Background legend", {
    sim <- markBackground(NS_sim_small, species = "Sprat")
    dat <- plot(getBiomass(sim), return_data = TRUE)
    sprat <- dat[dat$Species == "Sprat", ]
    # One row per saved time, not two: it used to be selected under its own
    # name and then appended again as "Background".
    expect_equal(nrow(sprat), nrow(sim@n))
    expect_identical(unique(sprat$Legend), "Background")
    expect_setequal(unique(as.character(dat$Species)),
                    c("Sprat", "Herring", "Cod"))
})

test_that("background = FALSE removes the background species", {
    sim <- markBackground(NS_sim_small, species = "Sprat")
    dat <- plot(getBiomass(sim), background = FALSE, return_data = TRUE)
    expect_false("Sprat" %in% dat$Species)
    expect_false("Background" %in% dat$Legend)
})

test_that("a background species can be selected explicitly", {
    sim <- markBackground(NS_sim_small, species = "Sprat")
    dat <- plot(getBiomass(sim), species = "Sprat", return_data = TRUE)
    expect_identical(unique(as.character(dat$Species)), "Sprat")
    expect_identical(unique(dat$Legend), "Background")
    # And a selection that leaves it out really leaves it out.
    dat <- plot(getBiomass(sim), species = "Cod", return_data = TRUE)
    expect_identical(unique(as.character(dat$Species)), "Cod")
})

test_that("a background series with no retained values does not error", {
    sim <- markBackground(NS_sim_small, species = "Sprat")
    bio <- getBiomass(sim)
    bio[, "Sprat"] <- 0
    # Every Sprat value is now below the cutoff a logarithmic axis applies, so
    # the background group is empty. That used to abort the plot.
    dat <- plot(bio, return_data = TRUE)
    expect_false("Sprat" %in% dat$Species)
    expect_s3_class(plot(bio), "ggplot")
})

test_that("non-positive values survive a linear y axis", {
    signed <- getBiomass(NS_sim_small)
    signed[, "Cod"] <- 0
    signed[, "Sprat"] <- -signed[, "Sprat"]

    linear <- plot(signed, log_y = FALSE, return_data = TRUE)
    expect_true(all(linear$Biomass[linear$Species == "Cod"] == 0))
    expect_true(all(linear$Biomass[linear$Species == "Sprat"] < 0))

    # A logarithmic axis still drops them, since it cannot show them.
    logged <- plot(signed, log_y = TRUE, return_data = TRUE)
    expect_false("Cod" %in% logged$Species)
    expect_false("Sprat" %in% logged$Species)
})

test_that("plotRelative keeps values a logarithmic axis would have dropped", {
    first <- getBiomass(NS_sim_small)
    second <- first
    second[, "Cod"] <- 0
    rel <- plotRelative(first, second)
    # 2 (0 - N) / (N + 0) = -2 wherever the second model has nothing.
    cod <- rel$data$rel_diff[rel$data$Species == "Cod"]
    expect_equal(unique(cod), -2)
})

test_that("the total is the total over every species the array holds", {
    bio <- getBiomass(NS_sim_small)
    all_species <- plot(bio, total = TRUE, return_data = TRUE)
    one_species <- plot(bio, species = "Cod", total = TRUE, return_data = TRUE)
    expect_equal(all_species$Biomass[all_species$Species == "Total"],
                 one_species$Biomass[one_species$Species == "Total"])
})

test_that("the time comparison methods honour highlight", {
    bio <- getBiomass(NS_sim_small)
    doubled <- ArrayTimeBySpecies(unclass_time(bio) * 2,
                                  value_name = attr(bio, "value_name"),
                                  units = attr(bio, "units"),
                                  params = NS_params_small)
    p2 <- plot2(bio, doubled, highlight = "Cod")
    expect_gt(drawn_linewidth(p2, "Cod", NS_params_small),
              drawn_linewidth(p2, "Sprat", NS_params_small))

    pr <- plotRelative(bio, doubled, highlight = "Cod")
    expect_gt(drawn_linewidth(pr, "Cod", NS_params_small, layer = 2),
              drawn_linewidth(pr, "Sprat", NS_params_small, layer = 2))
})
