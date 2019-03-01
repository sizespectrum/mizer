data(NS_species_params_gears)
data(inter)

params <- MizerParams(NS_species_params_gears, inter)
sim <- project(params)

params <- MizerParams(NS_species_params_gears, inter, no_w = 500)
sim1 <- project(params)

params <- MizerParams(NS_species_params_gears, inter, no_w = 1000)
sim2 <- project(params)

params <- MizerParams(NS_species_params_gears, inter, no_w = 5000)
sim3 <- project(params)

tmax <- 101
round((getBiomass(sim)[tmax,]/getBiomass(sim1)[tmax,]),3)
round((getBiomass(sim2)[tmax,]/getBiomass(sim3)[tmax,]),3)
