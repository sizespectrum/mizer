# This is useful for playing with the example from the app
# without running the app, which allows for easier bug fixing

#### set background ####
p_bg  <- setBackground(
    set_scaling_model(no_sp = 10, min_w_inf = 10, max_w_inf = 1e5,
                      min_egg = 1e-4, min_w_mat = 10^(0.4),
                      knife_edge_size = 100, kappa = 500)
)

######### add mullet ####
# some data from fishbase at 
# http://www.fishbase.org/summary/Mullus-barbatus+barbatus.html
# some parameter info is in this table
# https://www.dropbox.com/s/iqiydcrxqrx0k0w/paramsTable.jpg?dl=0
# length to weight conversion constants from 
# http://www.fishbase.org/popdyn/LWRelationshipList.php?ID=790&GenusName=Mullus&SpeciesName=barbatus+barbatus&fc=332
a_m <- 0.0085
b_m <- 3.11
# asymptotic length from
# http://www.fishbase.org/popdyn/PopGrowthList.php?ID=790&GenusName=Mullus&SpeciesName=barbatus+barbatus&fc=332
L_inf_m <- 24.3
# length at maturity from 
# http://www.fishbase.org/summary/Mullus-barbatus+barbatus.html
L_mat <- 11.1
species_params <- data.frame(
    species = "mullet",
    w_min = 0.002, # mizer's default egg weight, used in NS
    # w_inf = 251.94, #is the old value we used. Where is it from ? It differs to below
    w_inf = a_m*L_inf_m^b_m, # from fishbase
    # w_mat = 16.48, #is the old value we used. Where is it from ? It differs to below
    w_mat = a_m*L_mat^b_m, # from fishbase
    beta = 283, # = beta_gurnard from North sea. Silvia says gurnard is similar.
    sigma = 1.8, # = sigma_gurnard from North sea. Silvia says gurnard is similar.
    z0 = 0,
    alpha = 0.6, # unknown, mizer default=0.6
    erepro = 0.1, # unknown
    sel_func = "knife_edge", # not used but required
    knife_edge_size = 100, # we can choose
    gear = "knife_edge_gear",
    k = 0,
    r_max = 10^50,
    k_vb = 0.6,
    a = a_m,
    b = b_m
)
# k_vb is from 
# http://www.fishbase.org/popdyn/PopGrowthList.php?ID=790&GenusName=Mullus&SpeciesName=barbatus+barbatus&fc=332
p <- addSpecies(p_bg, species_params, SSB = 100)
plotSpectra(p)

############# add hake ####
# Merluccius merluccius  (European hake)
# http://www.fishbase.org/summary/Merluccius-merluccius.html
#! Currently hake and mullet are both using the same feeding level. What to do about it ?
# length to weight conversion: w=a*L^b, a = 0.0046, b = 3.12
# http://www.fishbase.org/popdyn/LWRelationshipList.php?ID=30&GenusName=Merluccius&SpeciesName=merluccius&fc=184
a <- 0.0046
b <- 3.12
# characteristic weights from
# http://www.fishbase.org/Reproduction/MaturityList.php?ID=30&GenusName=Merluccius&SpeciesName=merluccius&fc=184
# http://www.fishbase.org/graph/graphLengthFM01.php?RequestTimeout=50000&ID=30&genusname=Merluccius&speciesname=merluccius&fc=184&gm_lm=29.832069860776&gm_loo=81.220460002349
L_inf <- 81.2
L_mat <- 29.83
# Some information below is from Richard Law's document (RLD) at
# https://www.dropbox.com/s/g701wgcnhr12qpg/species%20%282%29.pdf?dl=0

species_params <- data.frame(
    species = "hake",
    w_min = 0.001, # mizer default
    w_inf = a*L_inf^b, # from fishbase
    w_mat = a*L_mat^b, # from fishbase
    beta = exp(2.4), #RLD and Blanchard thesis p 88
    sigma = 1.1, #RLD and Blanchard thesis p 88
    z0 = 0,
    alpha = 0.6, # unknown, using mizer default=0.6
    erepro = 0.1, # unknown
    sel_func = "knife_edge", # not used but required
    knife_edge_size = 100, # can choose
    gear = "knife_edge_gear",
    k = 0,
    r_max = 10^50, #why do I need r_max after combining before
    k_vb = 0.1, # from FB website below
    a = a,
    b = b
)
#k_vb <- 0.1 # from FB website below
# http://www.fishbase.org/popdyn/PopGrowthList.php?ID=30&GenusName=Merluccius&SpeciesName=merluccius&fc=184
p <- addSpecies(p, species_params, SSB = 200)
plotSpectra(p)

#### run simulation ####
s <- project(p, t_max = 15, effort = 0)

#### plot biomass ####
plotBiomass(s)

ggplot(p@species_params, aes(x = species, y = erepro)) + 
    geom_col() + geom_hline(yintercept = 1, color="red")
