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

object@phiMortality

# extract n_pp and n from sim object 
nt <- dim(sim@n_pp)[1]  # index of last time step
n_pp <- sim@n_pp[nt, ]
n <- sim@n[nt, , ]
# sim@n[time,species, weight ]
# we need to get species index back even though there is only one species
dim(n) <- c(dim(object@interaction)[1], dim(sim@n)[3])

mizerResult <- rowSums(sweep(object@pred_kernel,3,object@dw_full*object@w_full*n_pp,"*", check.margin=FALSE),dims=2)
#presumably the above over estimates, because it is based on a left-Riemman sum. The following uses right-Riemman sum, provided 
# that the last entry integrated via rowSums in `mizerResult` is negligable

#Presweep <- array(dim=c(1,dim(object@pred_kernel)[2],(dim(object@pred_kernel)[3]-1)))
#Presweep[1,,] <- object@pred_kernel[,,-1]
#mizerUnderEstimate2 <- rowSums(sweep(Presweep,3,(object@dw_full[-(length(object@dw_full))])*(object@w_full[-1])*(n_pp[-1]),"*", check.margin=FALSE),dims=2)

######################
w0 <- object@w[1]

Beta <- log(object@species_params$beta)
sigma <- object@species_params$sigma
wFull <- object@w_full
xFull <- log(wFull)
xFull <- xFull - xFull[1]
dx <- xFull[2]-xFull[1]

idx_sp <- (length(object@w_full) - length(object@w) + 1):length(object@w_full)

s <- exp(-(xFull - Beta)^2/(2*sigma^2))

smat <- matrix(0, nrow = dim(n)[1], ncol=length(xFull))
for(i in 1:dim(n)[1]){
  smat[i, ] <- exp(-(xFull - Beta[i])^2/(2*sigma[i]^2))
}

f <- wFull*wFull*n_pp

#background energy via fft
energy <- dx*Re(fft(fft(s)*fft(f), inverse=TRUE)/length(s))

#compare
#plot(log(object@w),energy[idx_sp])
#lines(log(object@w), mizerResult)
#lines(log(object@w), mizerUnderEstimate2)

######### adding energy from eating fish ###

fullEnergyMat <- matrix(0, nrow = dim(n)[1], ncol=length(object@w))

for(i in 1:(dim(n)[1])) {
    fishEaten <- rep(0, length.out = length(wFull))
    fishEaten[idx_sp] <- (object@interaction %*% n)[i, ]
    f2 <- wFull*wFull*(n_pp + fishEaten)
    fullEnergy <- dx*Re(fft(fft(smat[i,])*fft(f2), inverse=TRUE)/length(smat[i,]))
    fullEnergyMat[i,] <- fullEnergy[idx_sp]
}

plot(log(object@w),fullEnergyMat[1,])

plot(log(object@w),fullEnergyMat[2,])

################################## compare with mizer result

gPP <- getPhiPrey(object, n, n_pp)

i <- 6
plot(log(object@w), gPP[i,])
lines(log(object@w), fullEnergyMat[i,])
