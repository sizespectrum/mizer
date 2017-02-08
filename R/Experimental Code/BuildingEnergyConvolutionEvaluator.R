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
# sim@n[time,species, weight ]
# we need to get species index back even though there is only one species
dim(n) <- c(1, length(n))

mizerResult <- rowSums(sweep(object@pred_kernel,3,object@dw_full*object@w_full*n_pp,"*", check.margin=FALSE),dims=2)
#presumably the above over estimates, because it is based on a left-Riemman sum. The following uses right-Riemman sum, provided 
# that the last entry integrated via rowSums in `mizerResult` is negligable

Presweep <- array(dim=c(1,dim(object@pred_kernel)[2],(dim(object@pred_kernel)[3]-1)))
Presweep[1,,] <- object@pred_kernel[,,-1]
mizerUnderEstimate2 <- rowSums(sweep(Presweep,3,(object@dw_full[-(length(object@dw_full))])*(object@w_full[-1])*(n_pp[-1]),"*", check.margin=FALSE),dims=2)

######################
w0 <- object@w[1]

Beta <- log(object@species_params$beta)
sigma <- object@species_params$sigma
wFull <- object@w_full
xFull <- log(wFull)
dx <- xFull[2]-xFull[1]
s <- exp(-(xFull - xFull[1] - Beta)^2/(2*sigma^2))

f <- wFull*wFull*n_pp

#background energy via fft
energy <- dx*Re(fft(fft(s)*fft(f), inverse=TRUE)/length(s))
idx_sp <- (length(object@w_full) - length(object@w) + 1):length(object@w_full)

#compare
plot(log(object@w),energy[idx_sp])
lines(log(object@w), mizerResult)
lines(log(object@w), mizerUnderEstimate2)

######### adding energy from eating fish ###

fishEaten <- rep(0, length.out = length(wFull))
fishEaten[idx_sp] <- (object@interaction %*% n)[1, ]
f2 <- wFull*wFull*(n_pp + fishEaten)
fullEnergy <- dx*Re(fft(fft(s)*fft(f2), inverse=TRUE)/length(s))
plot(log(object@w),fullEnergy[idx_sp])

################################## compare with mizer result

n_eff_prey <- sweep(object@interaction %*% n, 2, object@w * object@dw, "*", check.margin=FALSE) 
idx_sp <- (length(object@w_full) - length(object@w) + 1):length(object@w_full)
mizerResultFish <- rowSums(sweep(object@pred_kernel[,,idx_sp,drop=FALSE],c(1,3),n_eff_prey,"*", check.margin=FALSE),dims=2)[1,]

lines(log(object@w),mizerResultFish+mizerResult)


