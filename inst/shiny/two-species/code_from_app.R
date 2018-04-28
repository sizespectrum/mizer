# This is useful for playing with the example from the app
# without running the app, which allows for easier bug fixing

#### set background ####
p_bg  <- setBackground(
    set_scaling_model(no_sp = 10, no_w = 400,
                      min_w_inf = 10, max_w_inf = 1e5,
                      min_egg = 1e-4, min_w_mat = 10^(0.4),
                      knife_edge_size = 10^5, kappa = 5000,
                      lambda = 2.08, f0 = 0.6, h = 34)
)

rfac <- 1.01
effort <- 0.4
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
    species = "Mullet",
    w_min = 0.001, # mizer's default egg weight, used in NS
    # w_inf = 251.94, #is the old value we used. Where is it from ? It differs to below
    w_inf = a_m*L_inf_m^b_m, # from fishbase
    # w_mat = 16.48, #is the old value we used. Where is it from ? It differs to below
    w_mat = a_m*L_mat^b_m, # from fishbase
    beta = 283, # = beta_gurnard from North sea. Silvia says gurnard is similar.
    sigma = 1.8, # = sigma_gurnard from North sea. Silvia says gurnard is similar.
    z0 = 0,
    alpha = 0.6, # unknown, mizer default=0.6
    erepro = 0.1, # unknown
    sel_func = "sigmoid_length", # not used but required
    gear = "sigmoid_gear",
    l25 = 10,
    l50 = 16.6,
    k = 0,
    k_vb = 0.6,
    a = a_m,
    b = b_m,
    gamma = 0.0017,
    h = 50
)
# k_vb is from 
# http://www.fishbase.org/popdyn/PopGrowthList.php?ID=790&GenusName=Mullus&SpeciesName=barbatus+barbatus&fc=332
p <- addSpecies(p_bg, species_params, SSB = 1400, effort=effort, 
                rfac = rfac, iterate=FALSE)
# plotSpectra(p)

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
    species = "Hake",
    w_min = 0.001, # mizer default
    w_inf = a*L_inf^b, # from fishbase
    w_mat = a*L_mat^b, # from fishbase
    beta = exp(2.4), #RLD and Blanchard thesis p 88
    sigma = 1.1, #RLD and Blanchard thesis p 88
    z0 = 0,
    alpha = 0.6, # unknown, using mizer default=0.6
    erepro = 0.1, # unknown
    sel_func = "sigmoid_length", # not used but required
    gear = "sigmoid_gear",
    l25 = 10,
    l50 = 16.6,
    k = 0,
    r_max = 10^50, #why do I need r_max after combining before
    k_vb = 0.1, # from FB website below
    a = a,
    b = b,
    gamma = 0.003,
    h = 20
)
#k_vb <- 0.1 # from FB website below
# http://www.fishbase.org/popdyn/PopGrowthList.php?ID=30&GenusName=Merluccius&SpeciesName=merluccius&fc=184
p <- addSpecies(p, species_params, SSB = 600, effort=effort, 
                rfac = rfac, iterate=FALSE)
plotSpectra(p)

p@species_params$erepro <- 1000

#### run simulation ####
sim <- project(p, t_max = 50, t_save = 5, effort = effort)
plotBiomass(sim)
p@initial_n <- sim@n[11, , ]
p@initial_n_pp <- sim@n_pp[11, ]
p@species_params$r_max <- Inf

# Retune the values of erepro so that we get the correct level of
# recruitment without stock-recruitment relationship
mumu <- getZ(p, p@initial_n, p@initial_n_pp, effort = effort)
gg <- getEGrowth(p, p@initial_n, p@initial_n_pp)
rdi <- getRDI(p, p@initial_n, p@initial_n_pp)
# TODO: vectorise this
no_sp <- length(p@species_params$species)
erepro <- 1:no_sp  # set up vector of right dimension
for (i in (1:no_sp)) {
    gg0 <- gg[i, p@species_params$w_min_idx[i]]
    mumu0 <- mumu[i, p@species_params$w_min_idx[i]]
    DW <- p@dw[p@species_params$w_min_idx[i]]
    erepro[i] <- p@species_params$erepro[i] *
        (p@initial_n[i, p@species_params$w_min_idx[i]] *
             (gg0 + DW * mumu0)) / rdi[i]
}
p@species_params$erepro <- erepro
s1 <- project(p, t_max=15, effort = 0.4)
plotBiomass(s1)
yield_old <- getYield(s1)[11, ]

a <- p@species_params["Hake", "a"]
b <- p@species_params["Hake", "b"]
l50 <- 20.50
sd <- 0.331
l25 = l50 - log(3) * sd
p@species_params["Hake", "l50"] <- l50
p@species_params["Hake", "l25"] <- l25
p@selectivity["sigmoid_gear", "Hake", ] <- sigmoid_length(p@w, l25, l50, a, b)
s2 <- project(p, t_max=15, effort = 0.4)
plotBiomass(s2)
plotYield(s2)
y <- getYield(s2)
y[1, ] <- yield_old
ym <- reshape2::melt(y, varnames = c("Year", "Species"), 
                     value.name = "Yield")
ym <- subset(ym, ym$Yield > 0)
p <- ggplot(ym) + 
    geom_line(aes(x=Year, y=Yield, colour=Species, linetype=Species)) +
    scale_y_continuous(name="Yield [g/year]", limits = c(0, NA))
p

s <- project(p, t_max = 1, effort = 0.4)
t_max = 1; t_save=0.1; effort = 0.4; shiny_progress = NULL
dt = 0.1
initial_n <- p@initial_n; initial_n_pp <- p@initial_n_pp
object <- p

#### plot biomass ####

plot(p@w, s@n[1,11,], type="l", log="xy")
lines(p@w, s@n[12,11,], col="blue")
sum(p@w * p@dw * s@n[1,11, ] * p@psi[11, ])
sum(p@w * p@dw * s@n[2,11, ] * p@psi[11, ])
sum(p@w * p@dw * s@n[3,11, ] * p@psi[11, ])
plot(p@w, s@n[3,11, ] - s@n[2,11, ], type="l", log="x")
lines(p@w, s@n[2,11, ] - s@n[1,11, ], log="x", col="blue")
lines(p@w, s@n[11,11, ] - s@n[1,11, ], log="x", col="red")

ggplot(p@species_params, aes(x = species, y = erepro)) + 
    geom_col() + geom_hline(yintercept = 1, color="red")
