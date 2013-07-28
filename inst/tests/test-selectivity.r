context("Selectivity functions")

test_that("knife-edge selectivity function is working properly",{
    data(species_params_gears)
    species_params_gears$sel_func <- "knife_edge"
    species_params_gears$knife_edge_size <- 1000
    species_params_gears$knife_edge_size[species_params_gears$gear == "industrial"] <- 4
    # Chop off l25, l50, a and b columns
    species_params_gears <- species_params_gears[,!(colnames(species_params_gears) %in% c("l25","l50","a","b"))]
    data(inter)
    params <- MizerParams(species_params_gears, inter)
    industrial_species <- as.character(species_params_gears$species[species_params_gears$gear == "industrial"])
    pelagic_species <- as.character(species_params_gears$species[species_params_gears$gear == "pelagic"])
    beam_trawl_species <- as.character(species_params_gears$species[species_params_gears$gear == "beam_trawl"])
    otter_trawl_species <- as.character(species_params_gears$species[species_params_gears$gear == "otter_trawl"])

    expect_that(all(params@selectivity["industrial",industrial_species,params@w >= 4] == 1), is_true())
    expect_that(all(params@selectivity["industrial",industrial_species,params@w < 4] == 0), is_true())
    expect_that(all(params@selectivity["pelagic",pelagic_species,params@w >= 1000] == 1), is_true())
    expect_that(all(params@selectivity["pelagic",pelagic_species,params@w < 1000] == 0), is_true())
    expect_that(all(params@selectivity["beam_trawl",beam_trawl_species,params@w >= 1000] == 1), is_true())
    expect_that(all(params@selectivity["beam_trawl",beam_trawl_species,params@w < 1000] == 0), is_true())
    expect_that(all(params@selectivity["otter_trawl",otter_trawl_species,params@w >= 1000] == 1), is_true())
    expect_that(all(params@selectivity["otter_trawl",otter_trawl_species,params@w < 1000] == 0), is_true())

    sim <- project(params,t_max=10)
    fm <- getFMortGear(sim)
    expect_that(all(fm[10,"industrial",industrial_species,sim@params@w >= 4] > 0), is_true())
    expect_that(all(fm[10,"industrial",industrial_species,sim@params@w < 4] == 0), is_true())
    expect_that(all(fm[10,"pelagic",pelagic_species,sim@params@w >= 1000] > 0), is_true())
    expect_that(all(fm[10,"pelagic",pelagic_species,sim@params@w < 1000] == 0), is_true())
    expect_that(all(fm[10,"beam_trawl",beam_trawl_species,sim@params@w >= 1000] > 0), is_true())
    expect_that(all(fm[10,"beam_trawl",beam_trawl_species,sim@params@w < 1000] == 0), is_true())
    expect_that(all(fm[10,"otter_trawl",otter_trawl_species,sim@params@w >= 1000] > 0), is_true())
    expect_that(all(fm[10,"otter_trawl",otter_trawl_species,sim@params@w < 1000] == 0), is_true())
})
