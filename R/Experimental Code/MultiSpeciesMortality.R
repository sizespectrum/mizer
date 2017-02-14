source('./R/MizerParams-class.R')
source('./R/MizerSim-class.R')
source('./R/project_methods.R')
source('./R/selectivity_funcs.R')
source('./R/summary_methods.R')
source('./R/wrapper_functions.R')
source('./R/plots.R')
source('./R/project.R')
library(ggplot2)
library(grid)
library(methods)
library(plyr)
library(reshape2)

#params <- set_community_model(z0 = 0.1, f0 = 0.7, alpha = 0.2, recruitment = 4e7,
#                              sigma = 2, no_w=100)

#params <- set_trait_model(no_sp = 2, min_w_inf = 10, max_w_inf = 1e5)
#sim <- project(params)

params_data <- read.csv("./vignettes/NS_species_params.csv")
inter <- read.csv("./vignettes/inter.csv", row.names=1)
inter <- as(inter, "matrix")
params <- MizerParams(params_data, interaction = inter)
sim <- project(params, effort = 1, t_max = 10, dt = 0.1, t_save = 1)
plot(sim)
object <- sim@params 

# extract n_pp and n from sim object 
nt <- dim(sim@n_pp)[1]  # index of last time step
n_pp <- sim@n_pp[nt, ]
n <- sim@n[nt, , ]
# sim@n[time,species, weight ]
# we need to get species index back even though there is only one species
dim(n) <- c(dim(object@interaction)[1], dim(sim@n)[3])

noSpecies <- dim(object@interaction)[1]
muVals <- matrix(0, nrow = noSpecies, ncol = length(params@w))

w <- params@w
x <- log(w)
x <- x - x[1]
dx <- x[2]-x[1]

for (i in 1:noSpecies){
  for (j in 1:noSpecies){
    Beta <- log(object@species_params$beta)[j]
    sigma <- object@species_params$sigma[j]
    
  
    
    Delta <- dx*round(min(2*sigma, Beta)/dx)
    Beta <- dx*round(Beta/dx)
    Delta <- Beta
    min_cannibal <- 1+floor((Beta-Delta)/dx)
    #min_cannibal <- 0
    
    P <- x[length(x)] + 2*Delta
    no_P <- 1+ceiling(P/dx)  # P/dx should already be integer 
    x_P <- (1:no_P)*dx#+Beta-Delta-dx
    
    phi <- rep(0, length(x_P))
    
    phi[abs(x_P+Beta-P)<Delta] <- exp(-(x_P[abs(x_P+Beta-P)<Delta] + Beta - P)^2/(2*sigma^2)) 
    
    feeding_level <- getFeedingLevel(object, n=n, n_pp=n_pp)[j,]
    f <- (1-feeding_level)*object@search_vol[j,]*n[j,]*w*object@interaction[j,i]
    f <- c(f[min_cannibal:length(x)], rep(0, length(phi)-length(x)))
    
    
    mortalityIntegral <- dx*Re(fft(fft(phi)*fft(f), inverse=TRUE)/no_P)
    
    mu <- c(mortalityIntegral[(no_P-min_cannibal):no_P], mortalityIntegral[1:(length(x)-min_cannibal-1)])
  
    muVals[i, ] <- muVals[i, ]+mu  
  }
   
}



m2 <- getM2(object, n, n_pp)

dim(m2)

plot(muVals[1, ])
lines(m2[1, ])
