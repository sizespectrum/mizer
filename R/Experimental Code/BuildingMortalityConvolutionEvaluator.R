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

params <- set_community_model(z0 = 0.1, f0 = 0.7, alpha = 0.2, recruitment = 4e7,
                              sigma = 2, no_w=100)
sim <- project(params)
object <- sim@params 

# extract n_pp and n from sim object 
nt <- dim(sim@n_pp)[1]  # index of last time step
n_pp <- sim@n_pp[nt, ]
n <- sim@n[nt, , ]
# we need to get species index back even though there is only one species
dim(n) <- c(1, length(n))

# Do the mizer calculation of m2
# from lines 211 and 213 of project_methods.R
n_total_in_size_bins <- sweep(n, 2, object@dw, '*', check.margin=FALSE) # N_i(w)dw
pred_rate <- sweep(object@pred_kernel,c(1,2),(1-feeding_level)*object@search_vol*n_total_in_size_bins,"*", check.margin=FALSE)
# from lines 293 and 297 of project_methods.R
idx_sp <- (length(object@w_full) - length(object@w) + 1):length(object@w_full)
m2 <- t(object@interaction) %*% colSums(aperm(pred_rate, c(2,1,3)),dims=1)[,idx_sp]

Beta <- log(object@species_params$beta)
sigma <- object@species_params$sigma
w <- params@w
x <- log(w)
x <- x - x[1]
dx <- x[2]-x[1]

Delta <- dx*round(min(2*sigma, Beta)/dx)

Beta <- dx*round(Beta/dx)
Delta <- Beta
min_cannibal <- 1+floor((Beta-Delta)/dx)
#min_cannibal <- 0

P <- x[length(x)] + 2*Delta
no_P <- 1+ceiling(P/dx)  # P/dx should already be integer 
x_P <- (1:no_P)*dx#+Beta-Delta-dx

feeding_level <- getFeedingLevel(object, n=n, n_pp=n_pp)

# pred_rate <- sweep(object@pred_kernel,c(1,2),(1-feeding_level)*object@search_vol*n_total_in_size_bins,"*", check.margin=FALSE)

fmat <- sweep(object@pred_kernel,c(1,2),(1-feeding_level)*object@search_vol*n,"*", check.margin=FALSE)

dim(object@pred_kernel)
length(f)


for(i in 1:dim(n)[1]){
  fmat[i, ] <- c(fmat[i,min_cannibal:length(x)], rep(0, no_P-length(x)+min_cannibal-1))
}

#min_cannibal

#dim(fmat)

phi <- rep(0, length(x_P))

phi[abs(x_P+Beta-P)<Delta] <- exp(-(x_P[abs(x_P+Beta-P)<Delta] + Beta - P)^2/(2*sigma^2)) 

f <- (1-feeding_level)*object@search_vol*n*w
f <- c(f[min_cannibal:length(x)], rep(0, no_P-length(x)+min_cannibal-1))
f <- c(f[min_cannibal:length(x)], rep(0, length(phi)-length(x)))





mortalityIntegral <- dx*Re(fft(fft(phi)*fft(f), inverse=TRUE)/no_P)

mu <- c(mortalityIntegral[(no_P-min_cannibal):no_P], mortalityIntegral[1:(length(x)-min_cannibal-1)])

plot(mu)
lines(m2[1,])


