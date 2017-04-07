
source('./R/MizerParams-class.r')
source('./R/MizerSim-class.R')
source('./R/project_methods.R')
source('./R/selectivity_funcs.R')
source('./R/summary_methods.R')
source('./R/wrapper_functions.R')
source('./R/project_methods.R')
source('./R/plots.R')
source('./R/project.R')
library(ggplot2)
library(grid)
library(methods)
library(plyr)
library(reshape2)
params_data <- read.csv("./vignettes/NS_species_params.csv")
params_data$sel_func <- "sigmoid_length"
#params_data$l25 <- c(7.6, 9.8, 8.7, 10.1, 11.5, 19.8, 16.4, 19.8, 11.5,
#                     19.1, 13.2, 35.3)
params_data$l25 <- c(7.6, 9.8, 8.7, 10.1, 11.5, 27.0, 16.4, 19.8, 11.5,
                     19.1, 13.2, 35.3)
params_data$l50 <- c(8.1, 11.8, 12.2, 20.8, 17.0, 29.0, 25.8, 29.0, 17.0,
                     24.3, 22.9, 43.6)
params_data$a <- c(0.007, 0.001, 0.009, 0.002, 0.010, 0.006, 0.008, 0.004,
                   0.007, 0.005, 0.005, 0.007)
params_data$b <- c(3.014, 3.320, 2.941, 3.429, 2.986, 3.080, 3.019, 3.198,
                   3.101, 3.160, 3.173, 3.075)
inter <- read.csv("./vignettes/inter.csv", row.names=1)
inter <- as(inter, "matrix")
#### sing
nme <- rownames(inter)[6]
inter <- as(inter[6,6], "matrix")
rownames(inter) <- nme
colnames(inter) <-nme
####
params <- MizerParams(params_data[6,], interaction = inter, no_w = 200)
sim <- project(params, effort = 1, t_max = 80, dt = 0.1, t_save = 1)
plot(sim)
##################
####### put initializeation code here
######  Evaluate it using new smatMlong, or just directly

###### Write a sample method

####### put it through on a new branch


########################################## define smatMlong here in a way that can be inserted around line 635 of params code

Beta <- log(res@species_params$beta)            
# Here we use default rr[j] = beta[j] + 3*sigma[j]
rr <- Beta + 3*res@species_params$sigma 

wfull <- res@w_full
xfull <- log(wfull)
xfull <- xfull - xfull[1]
dx <- xfull[2]-xfull[1]

rr <- dx*ceiling(rr/dx)
# Determine period used
P <- max(xfull[length(xfull)] + rr)
# Determine number of x points used in period
no_P <- 1+ceiling(P/dx)  # P/dx should already be integer
noSpecies <- dim(res@interaction)[1]

# initially fill matrices with zeros
phiMortality <- matrix(0, nrow = noSpecies, ncol = no_P)
fphiMortality <- matrix(0, nrow = noSpecies, ncol = no_P)

for (j in 1:noSpecies){
    # Prepare local phi, which will equal j th row of matrix, used in loop
    phi <- rep(0, no_P)
    # Our phi is a periodic extension of the normal feeding kernel, so, for 0<=x<=P we use phi[x-P] as our
    # value of the period P extension of phi, since support(phi)=[-rr,0]
    phi[x_P-P >= -rr[j]] <- exp(-(Beta[j]-P+x_P[x_P-P >= -rr[j]])^2/(2*(res@species_params$sigma[j])^2))
    # This phi value is added to our output
    phiMortality[j, 1:length(phi)] <- phi
    # We also save the fft of this vector, so we don't have to use too many fft s in the time evolution
    fphiMortality[j, ] <- fft(phiMortality[j, ])
}

res@smatMlong <- phiMortality
res@fsmatMlong <- fphiMortality

####################################

########################## Next part of code determines get_pred_rate_fft[j,k] for all j in (1:no_sp) and all k in (1:length(wfull))
# here j is the species of the predator, and wfull[k] is the size of the prey, here the sum of interaction[j,i]*get_pred_rate_fft[j,k]
# over j death rate of a size k fish, provided k >= FishEggIndex
# also, the sum of get_pred_rate_fft[j,k], over j is getM2Background[k] == death rate of a generic size wfull[k] individual 
# due to predation when the interaction matrix is unitary, this also gives the anhilation rate for dead fish, 
# however for dead fish we may only wish to use these values in the fish size range 

object <- params
noSpecies <- dim(object@interaction)[1]

# Prepare a (noSpecies times length(w)) matrix, that will be used to ouput the mortality integral data
muVals <- matrix(0, nrow = noSpecies, ncol = length(object@w_full))

# Obtain weight vector wfull, and the corresponding vector xfull in log-space
wfull <- object@w_full
xfull <- log(wfull)
xfull <- xfull - xfull[1]
# Values of x are evenly spaced, with a difference of dx, starting with zero
dx <- xfull[2]-xfull[1]
# Obtain feeding level
feeding_level <- getFeedingLevel(object, n=n, n_pp=n_pp)

# no_P is the number of x points sampled over a period P (period P is used in spectral integration)
no_P <- length(object@smatMlong[1,])

LL <- length(object@w)
LLfull <- length(object@w_full)

#xlonger[1] == 0
#xlonger[no_P] == dx*(no_P-1) == P == xfull[LLfull] + ceiling(max((log(object@beta) + 3*object@sigma))/dx)*dx
# no_P == 1+ ceiling(P/dx)
#xlonger[(1:LLfull)] == xfull

xlonger <- (0:(no_P-1))*dx

FishEggIndex <- LLfull-LL+1

#x == xlonger[(FishEggIndex:LLfull)]

# Make a matrix to compute intermediate values, which are mortality integrals from predation via different 
# species, whose terms equal the integral within the sum of (3.12), discounting the interaction matrix theta
muIntermediate <- matrix(0, nrow = noSpecies, ncol = length(object@w_full))

########################################## Express full convolution integral ######


# Fill out this intermediate matrix
for (j in 1:noSpecies){
    # We express the intermediate values as a a convolution integral involving
    # two functions: q and fsmatM. Here q is all the integrand of (3.12), except the
    # feeding kernel and theta, and we sample it from 0 to P, but it is only 
    # non-zero from fishEggSize to X, where P = X + beta + 3*sigma
    q <- rep(0, no_P)
    q[(FishEggIndex: LLfull)] <- (1-feeding_level[j,])*object@search_vol[j,]*n[j,]*w
    # For convolution, we imagine f is a period P function, and sample it on [0,P]
    # We use fast fourier transforms to evalute this convolution integral
    mortalityIntegral <- dx*Re(fft((object@fsmatM[j,])*fft(q), inverse=TRUE)/no_P)
    # muIntermediate[j, ] measures how much i dies from being predated on by species j when object@interaction[j,i]=1
    muIntermediate[j, ] <- mortalityIntegral[1:length(object@w_full)]
}

OUTPUT@get_pred_rate_fft <- muIntermediate