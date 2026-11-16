test_that("snapshot biomass plots consume their own measurement metadata", {
    rlang::local_bindings(
        snapshot_pool_biomass = function(...) {
            data.frame(series = c("small", "large"), value = c(2, 5),
                       label = c("Small pool", "Large pool"),
                       colour = c("#0072B2", "#D55E00"))
        }, .env = globalenv()
    )
    p <- setComponent(NS_params_small, "pool", initial_value = 1,
                      dynamics_fun = "constant_other",
                      biomass_fun = "snapshot_pool_biomass")
    biomass <- getBiomass(p, include = "component:pool")
    before <- p$linecolour
    expect_s3_class(biomass, "MizerBiomass")
    plotted <- plot(biomass, return_data = TRUE)
    expect_equal(plotted$Biomass, c(2, 5))
    expect_identical(plotted$Label, c("Small pool", "Large pool"))
    expect_identical(p$linecolour, before)
    expect_s3_class(plot(biomass), "ggplot")
    expect_error(plot(biomass, species = "absent"), "Unknown biomass")
    expect_equal(plot(biomass, total = TRUE, return_data = TRUE)$Biomass,
                 c(2, 5, 7))
    vdiffr::expect_doppelganger("Biomass pool measurements", plot(biomass))
})

test_that("biomass vectors preserve coverage when subset and strip it in ops", {
    biomass <- getBiomass(NS_params_small, include = c("species", "resource"))
    selected <- biomass[c(4, 1, 4)]
    expect_identical(attr(selected, "series")$id, names(selected))
    expect_equal(as.numeric(selected), as.numeric(biomass)[c(4, 1, 4)])
    expect_error(plot(selected, total = TRUE), "same physical source")
    plain <- biomass * 2
    expect_identical(names(plain), names(biomass))
    expect_null(attr(plain, "series"))
    expect_null(attr(plain, "params"))
    expect_false(inherits(plain, "MizerBiomass"))
    expect_identical(log(biomass), log(as.numeric(biomass)) |>
                         setNames(names(biomass)))
    expect_null(attr(log(biomass), "units"))
    expect_equal(as.data.frame(biomass)$value, as.numeric(biomass))
    expect_false(any(grepl("intake_max", capture.output(str(biomass)))))
})

test_that("unavailable snapshot measurements explain missing totals", {
    p <- setComponent(NS_params_small, "unknown", initial_value = 1,
                      dynamics_fun = "constant_other")
    biomass <- getBiomass(p, include = c("resource", "components"))
    expect_output(print(biomass), "Unavailable.*component:unknown")
    expect_error(plot(biomass, total = TRUE), "total is unavailable")
    expect_match(plot(biomass)$labels$caption, "No biomass callback")
    expect_false(tail(as.data.frame(biomass)$available, 1))
    expect_null(measurement_caption(biomass["resource"]))
    unavailable <- getBiomass(p, include = "components")
    expect_no_warning(ggplot2::ggplot_build(plot(unavailable)))
})
