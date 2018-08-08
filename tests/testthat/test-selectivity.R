context("Selectivity functions")

test_that("knife-edge selectivity function is working properly",{
    data(NS_species_params_gears)
    NS_species_params_gears$sel_func <- "knife_edge"
    NS_species_params_gears$knife_edge_size <- 1000
    NS_species_params_gears$knife_edge_size[NS_species_params_gears$gear == "Industrial"] <- 500
    # Chop off l25, l50, a and b columns - the trawl selectivity
    NS_species_params_gears <- NS_species_params_gears[,!(colnames(NS_species_params_gears) %in% c("l25","l50","a","b"))]
    data(inter)
    params <- MizerParams(NS_species_params_gears, inter)
    industrial_species <- as.character(NS_species_params_gears$species[NS_species_params_gears$gear == "Industrial"])
    pelagic_species <- as.character(NS_species_params_gears$species[NS_species_params_gears$gear == "Pelagic"])
    beam_trawl_species <- as.character(NS_species_params_gears$species[NS_species_params_gears$gear == "Beam"])
    otter_trawl_species <- as.character(NS_species_params_gears$species[NS_species_params_gears$gear == "Otter"])

    expect_that(all(params@selectivity["Industrial",industrial_species,params@w < 500] == 0), is_true())
    expect_that(all(params@selectivity["Industrial",industrial_species,params@w >= 500] == 1), is_true())
    expect_that(all(params@selectivity["Pelagic",pelagic_species,params@w >= 1000] == 1), is_true())
    expect_that(all(params@selectivity["Pelagic",pelagic_species,params@w < 1000] == 0), is_true())
    expect_that(all(params@selectivity["Beam",beam_trawl_species,params@w >= 1000] == 1), is_true())
    expect_that(all(params@selectivity["Beam",beam_trawl_species,params@w < 1000] == 0), is_true())
    expect_that(all(params@selectivity["Otter",otter_trawl_species,params@w >= 1000] == 1), is_true())
    expect_that(all(params@selectivity["Otter",otter_trawl_species,params@w < 1000] == 0), is_true())

    sim <- project(params, t_max = 10, effort = 1)
    fm <- getFMortGear(sim)
    expect_that(all(fm[10,"Industrial",industrial_species,sim@params@w < 500] == 0), is_true())
    expect_that(all(fm[10,"Industrial",industrial_species,sim@params@w >= 500] > 0), is_true())
    expect_that(all(fm[10,"Pelagic",pelagic_species,sim@params@w >= 1000] > 0), is_true())
    expect_that(all(fm[10,"Pelagic",pelagic_species,sim@params@w < 1000] == 0), is_true())
    expect_that(all(fm[10,"Beam",beam_trawl_species,sim@params@w >= 1000] > 0), is_true())
    expect_that(all(fm[10,"Beam",beam_trawl_species,sim@params@w < 1000] == 0), is_true())
    expect_that(all(fm[10,"Otter",otter_trawl_species,sim@params@w >= 1000] > 0), is_true())
    expect_that(all(fm[10,"Otter",otter_trawl_species,sim@params@w < 1000] == 0), is_true())
})
