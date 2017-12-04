
library(ggplot2)
library(grid)
library(methods)
library(plyr)
library(reshape2)
library(mizer)
params_data <- read.csv("./vignettes/NS_species_params.csv")

inter <- read.csv("./vignettes/inter.csv", row.names=1)
inter <- as(inter, "matrix")

####
params <- MizerParams(params_data, interaction = inter, no_w = 1000)

ptm <- proc.time() 
sim <- project(params, effort = 1, t_max = 10, dt = 0.1, t_save = 1)
Tsim <- (proc.time() - ptm)[[3]]

plot(sim)

#############

n<- sim@n[11, , ]
n_pp <- sim@n_pp[11, ]

L <- 10^(3)

############


ptm <- proc.time() 
for (re in (1:L)){
    preyy <- getPhiPrey(object = params, n = n, n_pp = n_pp)
}
TgetPhiPrey <- ((proc.time() - ptm)[[3]])/L

###############

ptm <- proc.time() 
for (re in (1:L)){
    fd <- getFeedingLevel(object = params, n = n, n_pp = n_pp, phi_prey = preyy)
}
TgetFeedingLevel <- ((proc.time() - ptm)[[3]])/L

###############

ptm <- proc.time() 
for (re in (1:L)){
    gp <- getPredRate(object = params, n = n, n_pp = n_pp, feeding_level = fd)
}
TgetPredRate <- ((proc.time() - ptm)[[3]])/L

###############

ptm <- proc.time() 
for (re in (1:L)){
    m22 <- getM2(object = params, n = n, n_pp = n_pp, pred_rate = gp)
}
TgetM2 <- ((proc.time() - ptm)[[3]])/L

###############

ptm <- proc.time() 
for (re in (1:L)){
    m22 <- getM2(object = params, n = n, n_pp = n_pp)
}
TgetM2all <- ((proc.time() - ptm)[[3]])/L

###############

c(Tsim,TgetPhiPrey,TgetFeedingLevel,TgetPredRate,TgetM2, TgetM2all)

Tsimf <- 8.29
TgetPhiPreyf <- 0.00515
TgetFeedingLevelf <- 0.00011
TgetPredRatef <- 0.01138
TgetM2f <- 0.00003
TgetM2allf <- 0.01538

Tsimm <- 113.50000
TgetPhiPreym <- 0.80291
TgetFeedingLevelm <- 0.00014
TgetPredRatem <- 0.22052
TgetM2m <- 0.00003
TgetM2allm <- 1.07239


