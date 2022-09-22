context("Selectivity functions")

species_params <- NS_species_params
w <- NS_params@w

# knife-edge ----
test_that("knife-edge selectivity function is working properly",{
    expect_length(knife_edge(w, 20, species_params = species_params[1, ]),
                  length(w))
    
    NS_species_params_gears$sel_func <- "knife_edge"
    NS_species_params_gears$knife_edge_size <- 1000
    NS_species_params_gears$knife_edge_size[NS_species_params_gears$gear == "Industrial"] <- 500
    # Chop off l25, l50, a and b columns - the trawl selectivity
    NS_species_params_gears <- NS_species_params_gears[,!(colnames(NS_species_params_gears) %in% c("l25","l50","a","b"))]
    params <- newMultispeciesParams(NS_species_params_gears, inter)
    industrial_species <- as.character(NS_species_params_gears$species[NS_species_params_gears$gear == "Industrial"])
    pelagic_species <- as.character(NS_species_params_gears$species[NS_species_params_gears$gear == "Pelagic"])
    beam_trawl_species <- as.character(NS_species_params_gears$species[NS_species_params_gears$gear == "Beam"])
    otter_trawl_species <- as.character(NS_species_params_gears$species[NS_species_params_gears$gear == "Otter"])
    
    expect_true(all(params@selectivity["Industrial",industrial_species,params@w < 500] == 0))
    expect_true(all(params@selectivity["Industrial",industrial_species,params@w >= 500] == 1))
    expect_true(all(params@selectivity["Pelagic",pelagic_species,params@w >= 1000] == 1))
    expect_true(all(params@selectivity["Pelagic",pelagic_species,params@w < 1000] == 0))
    expect_true(all(params@selectivity["Beam",beam_trawl_species,params@w >= 1000] == 1))
    expect_true(all(params@selectivity["Beam",beam_trawl_species,params@w < 1000] == 0))
    expect_true(all(params@selectivity["Otter",otter_trawl_species,params@w >= 1000] == 1))
    expect_true(all(params@selectivity["Otter",otter_trawl_species,params@w < 1000] == 0))
    
    sim <- project(params, t_max = 10, effort = 1)
    fm <- getFMortGear(sim)
    expect_true(all(fm[10,"Industrial",industrial_species,sim@params@w < 500] == 0))
    expect_true(all(fm[10,"Industrial",industrial_species,sim@params@w >= 500] > 0))
    expect_true(all(fm[10,"Pelagic",pelagic_species,sim@params@w >= 1000] > 0))
    expect_true(all(fm[10,"Pelagic",pelagic_species,sim@params@w < 1000] == 0))
    expect_true(all(fm[10,"Beam",beam_trawl_species,sim@params@w >= 1000] > 0))
    expect_true(all(fm[10,"Beam",beam_trawl_species,sim@params@w < 1000] == 0))
    expect_true(all(fm[10,"Otter",otter_trawl_species,sim@params@w >= 1000] > 0))
    expect_true(all(fm[10,"Otter",otter_trawl_species,sim@params@w < 1000] == 0))
})

# sigmoid_length ----
test_that("sigmoid_length works", {
    expect_error(sigmoid_length(w, 20, 30, species_params = species_params),
                 "The selectivity function needs the weight-length parameters ")
    species_params$a <- 0.5
    species_params$b <- 3
    expect_length(sigmoid_length(w, 20, 30, species_params = species_params[1, ]),
                  length(w))
})

# double_sigmoid_length ----
test_that("double_sigmoid_length works", {
    expect_error(double_sigmoid_length(w, 20, 30, 40, 50, 
                                       species_params = species_params),
                 "The selectivity function needs the weight-length parameters ")
    species_params$a <- 0.5
    species_params$b <- 3
    expect_length(double_sigmoid_length(w, 20, 30, 40, 50,
                                        species_params = species_params[1, ]),
                  length(w))
    expect_error(double_sigmoid_length(w, 20, 30, 40, 30,
                                        species_params = species_params[1, ]),
                 "l50_right not less than l25_right")
})

# sigmoid_weight ----
test_that("sigmoid_weight works", {
    expect_length(sigmoid_weight(w, sigmoidal_weight = 20, 
                                 sigmoidal_sigma = 2), length(w))
})
