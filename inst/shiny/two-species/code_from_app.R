# This is useful for playing with the example from the app
# without running the app, which allows for easier bug fixing

#### set fishing parameters ####
effort <- 0.4
l50 <- c(16.6, 15.48, 20.5, 15.85)
names(l50) <- c("hake_old", "mullet_old", "hake_new", "mullet_new")
sd <- c(0.462, 2.10, 0.331, 2.05)
l25 = l50 - log(3) * sd

#### set background ####
p_bg  <- setBackground(
    set_scaling_model(no_sp = 10, no_w = 400,
                      min_w_inf = 10, max_w_inf = 1e5,
                      min_egg = 1e-4, min_w_mat = 10^(0.4),
                      knife_edge_size = Inf, kappa = 10000,
                      lambda = 2.08, f0 = 0.6, h = 34)
)

rfac <- 1.01

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
    l25 = l25["mullet_old"],
    l50 = l50["mullet_old"],
    k = 0,
    k_vb = 0.6,
    a = a_m,
    b = b_m,
    gamma = 0.0017,
    h = 50,
    linecolour = "red",
    linetype = "solid"
)
# k_vb is from 
# http://www.fishbase.org/popdyn/PopGrowthList.php?ID=790&GenusName=Mullus&SpeciesName=barbatus+barbatus&fc=332
p <- addSpecies(p_bg, species_params, SSB = 2800, effort=effort, 
                rfac = rfac)
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
    l25 = l25["hake_old"],
    l50 = l50["hake_old"],
    k = 0,
    k_vb = 0.1, # from FB website below
    a = a,
    b = b,
    gamma = 0.003,
    h = 20,
    linecolour = "blue",
    linetype = "solid"
)
#k_vb <- 0.1 # from FB website below
# http://www.fishbase.org/popdyn/PopGrowthList.php?ID=30&GenusName=Merluccius&SpeciesName=merluccius&fc=184
p_old <- addSpecies(p, species_params, SSB = 1200, effort=effort, 
                rfac = rfac, iterate=FALSE)
plotSpectra(p_old)

# p_old@species_params$erepro <- 1000

#### run simulation ####
sim_old <- project(p_old, t_max = 50, t_save = 5, effort = effort)
plotBiomass(sim_old)

# save(sim_old, file="hake_mullet.RData)
# sim_old <- readRDS(file="hake_mullet.RData")

#### run new simulation ####
p <- sim_old@params
no_sp <- length(p@species_params$species)
no_t <- dim(sim_old@n)[1]
p@initial_n <- sim_old@n[no_t, , ]
p@initial_n_pp <- sim_old@n_pp[no_t, ]
p@species_params$r_max <- Inf

# Retune the values of erepro so that we get the correct level of
# recruitment without stock-recruitment relationship
mumu <- getZ(p, p@initial_n, p@initial_n_pp, effort = effort)
gg <- getEGrowth(p, p@initial_n, p@initial_n_pp)
rdi <- getRDI(p, p@initial_n, p@initial_n_pp)
# TODO: vectorise this

for (i in (1:no_sp)) {
    gg0 <- gg[i, p@species_params$w_min_idx[i]]
    mumu0 <- mumu[i, p@species_params$w_min_idx[i]]
    DW <- p@dw[p@species_params$w_min_idx[i]]
    p@species_params$erepro[i] <- p@species_params$erepro[i] *
        (p@initial_n[i, p@species_params$w_min_idx[i]] *
             (gg0 + DW * mumu0)) / rdi[i]
}

# s1 <- project(p, t_max=15, effort = 0.4)
# plotBiomass(s1)
# getYield(s1)[dim(s1@n)[1], ]
# getSSB(s1)[dim(s1@n)[1], ]


a <- p@species_params["Hake", "a"]
b <- p@species_params["Hake", "b"]
p@species_params["Hake", "l50"] <- l50["hake_old"]
p@species_params["Hake", "l25"] <- l25["hake_old"]
p@selectivity["sigmoid_gear", "Hake", ] <- 
    sigmoid_length(p@w, l25["hake_old"], l50["hake_old"], a, b)
a <- p@species_params["Mullet", "a"]
b <- p@species_params["Mullet", "b"]
p@species_params["Mullet", "l50"] <- l50["mullet_old"]
p@species_params["Mullet", "l25"] <- l25["mullet_old"]
p@selectivity["sigmoid_gear", "Mullet", ] <- 
    sigmoid_length(p@w, l25["mullet_old"], l50["mullet_old"], a, b)
s2 <- project(p, t_max=15, t_save = 0.1, effort = 0.4)
plotBiomass(s2)

a <- p@species_params["Hake", "a"]
b <- p@species_params["Hake", "b"]
p@species_params["Hake", "l50"] <- l50["hake_new"]
p@species_params["Hake", "l25"] <- l25["hake_new"]
p@selectivity["sigmoid_gear", "Hake", ] <- 
    sigmoid_length(p@w, l25["hake_new"], l50["hake_new"], a, b)
a <- p@species_params["Mullet", "a"]
b <- p@species_params["Mullet", "b"]
p@species_params["Mullet", "l50"] <- l50["mullet_new"]
p@species_params["Mullet", "l25"] <- l25["mullet_new"]
p@selectivity["sigmoid_gear", "Mullet", ] <- 
    sigmoid_length(p@w, l25["mullet_new"], l50["mullet_new"], a, b)
s2 <- project(p, t_max=15, t_save = 0.1, effort = 0.4)
plotBiomass(s2)

# Plot yield
y <- getYield(s2)
y[1, ] <- yield_old
ym <- reshape2::melt(y, varnames = c("Year", "Species"), 
                     value.name = "Yield")
ym <- subset(ym, ym$Yield > 0)
ggplot(ym) + 
    geom_line(aes(x=Year, y=Yield, colour=Species, linetype=Species)) +
    scale_y_continuous(name="Yield [g/year]", limits = c(0, NA)) +
    scale_colour_manual(values = p@linecolour) +
    scale_linetype_manual(values = p@linetype)


# Plot changes in abundance
year <- 10
no_w <- length(p@w)
w_sel <- seq.int(1, no_w, by = floor(no_w/50))
w <- p@w[w_sel]
change <- s2@n[10*year+1, ,w_sel]/s2@n[1, ,w_sel] - 1
change_total <- colSums(s2@n[10*year+1, ,w_sel], na.rm = TRUE) /
                     colSums(s2@n[1, ,w_sel], na.rm = TRUE) - 1
ch <- rbind(change, "Total" = change_total)
names(dimnames(ch)) <- names(dimnames(change))
cf <- reshape2::melt(ch)
cf <- subset(cf, !is.nan(value))
names(cf)[1] <- "Species"
ggplot(cf, aes(x = w, y = value)) +
    geom_line(aes(colour = Species, linetype = Species)) +
    geom_hline(yintercept = 0) +
    scale_x_log10(name = "Size [g]", labels = prettyNum,
                  breaks = 10^(-3:4)) +
    scale_y_continuous(name = "Percentage change", limits = c(-0.50, 0.60),
                       labels = scales::percent, breaks = (-7:9)/10) +
    scale_colour_manual(values = p@linecolour) +
    scale_linetype_manual(values = p@linetype) +
    # geom_vline(xintercept = l50) +
    # geom_vline(xintercept = l50_old) +
    theme(text = element_text(size = 14))


# Selectivity curve
w_min_idx <- sum(p@w < 0.5)
w_max_idx <- which.min(p@w < 200)
w_sel <- seq(w_min_idx, w_max_idx, by = floor((w_max_idx-w_min_idx)/50))
w <- p@w[w_sel]
selectivity <- p@selectivity[2, , w_sel]
sf <- reshape2::melt(selectivity)
sf$Gear <- "T90 modfied"
selectivity_old <- sim@params@selectivity[2, , w_sel]
sf_old <- reshape2::melt(selectivity_old)
sf_old$Gear <- "Standard"
sf <- rbind(sf, sf_old)
names(sf)[1] <- "Species"
sf <- subset(sf, value > 0)
ggplot(sf, aes(x = w, y = value)) +
    geom_line(aes(colour = Species, linetype = Gear)) +
    scale_x_continuous(name = "Size [g]", labels = prettyNum)  +
    scale_colour_manual(values = p@linecolour)

a <- p@species_params$a
names(a) <- p@species_params$species
b <- p@species_params$b
names(b) <- p@species_params$species
ggplot(sf, aes(x = (w/a[Species])^(1/b[Species]), y = value)) +
    geom_line(aes(colour = Species, linetype = Gear)) +
    scale_x_continuous(name = "Length [cm]", labels = prettyNum)  +
    scale_colour_manual(values = p@linecolour)


# Catch profile
w_min_idx <- sum(p@w < 4)
w_max_idx <- which.min(p@w < 200)
w_sel <- seq(w_min_idx, w_max_idx, by = 1)
w <- p@w[w_sel]
catch <- p@selectivity[2, , w_sel] * s2@n[10*year+1, ,w_sel] * 0.4
catchf <- reshape2::melt(catch)
catchf$Gear <- "T90 modfied"
catch_old <- sim@params@selectivity[2, , w_sel] * sim@n[11, ,w_sel] * 0.4
catchf_old <- reshape2::melt(catch_old)
catchf_old$Gear <- "Standard"
catchf <- rbind(catchf, catchf_old)
names(catchf)[1] <- "Species"
catchf <- subset(catchf, value > 0)
ggplot(catchf, aes(x = w, y = value)) +
    geom_line(aes(colour = Species, linetype = Gear)) +
    scale_x_continuous(name = "Size [g]", labels = prettyNum)  +
    scale_y_continuous()
    scale_colour_manual(values = s2@params@linecolour)

ggplot(catchf, aes(x = w, y = value * w)) +
    geom_line(aes(colour = Species, linetype = Gear)) +
    scale_x_continuous(name = "Size [g]", labels = prettyNum)  +
    scale_y_continuous(name = "Yield distribution")
    scale_colour_manual(values = s2@params@linecolour)

a <- p@species_params$a
names(a) <- p@species_params$species
b <- p@species_params$b
names(b) <- p@species_params$species
ggplot(sf, aes(x = (w/a[Species])^(1/b[Species]), y = value)) +
    geom_line(aes(colour = Species, linetype = Gear)) +
    scale_x_continuous(name = "Length [cm]", labels = prettyNum)  +
    scale_colour_manual(values = s2@params@linecolour)


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
