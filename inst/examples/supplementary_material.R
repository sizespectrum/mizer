# Supplementary material for
# mizer: an R package for multispecies, trait-based and community size spectrum ecological modelling
# Finlay Scott, Julia L. Blanchard and Ken H. Andersen
# Methods in Ecology and Evolution

#----------------------------------------------------
# Preliminaries
#----------------------------------------------------

library(mizer)

# library(ggplot2)

rm(list=ls())

#----------------------------------------------------
# Spectra and Tropic cascades
#----------------------------------------------------

# We use the community model, the trait-based model and the multispecies model (parameterised for the
# North Sea) to explore the impacts of fishing.

# In the first row, for each model we plot:
#     the unfished biomass spectra of each species
#     the unfished total community biomass spectrum
#     the fished total community biomass spectrum
#     the unfished background resource spectrum
# In the second row, for each model we plot the relative abundances in the fished and unfished case
# (to simulate trophic cascades).

# Parameters for the community and trait-based models.
# Some of the parameters have the same value.
# See the package vignette for a description of the parameters.
h <- 20
sigma <- 1.3
beta <- 100
q <- 0.8 
n <- 3/4
lambda <- 2+q-n
f0 <- 0.6
kappa <- 5
r_pp <- 4
max_w_com <- 50e3
max_w_trait <- 100e3
alpha_com <- 0.17
alpha_trait <- 0.6
ks_com <- 0 # standard metabolism is turned off in the Community model
ks_trait <- 2.4
z0pre_com <- 3
z0pre_trait <- 2
z0_com <- z0pre_com * (max_w_com*0.9) ^(n-1)

# Make the Community model parameters object using the wrapper function: set_community_model()
# The value of rec_mult has been set so that N in the smallest size is approximately the same as the resource spectrum in the smallest size.
comm_params <- set_community_model(max_w=max_w_com, beta=beta, sigma=sigma, h=h, alpha=alpha_com, r_pp=r_pp, n=n, q=q, z0=z0_com, kappa=kappa, f0=f0, rec_mult=0.01)
# The time step of the simulations is 0.1, and the results are stored every time step (default values).
# The solver is not well behaved with Community model and it can be quite  unstable.
# We therefore run a long projection and average over time series.
comm_time_sim <- 800 
# First we project the community model through time without fishing (effort is set to 0).
comm_sim0 <- project(comm_params, t_max=comm_time_sim, effort=0, dt=0.1)
# Plot the results as a quick check that something sensible has happened
# Note that the plots with size on the x-axis are snap shots of the final time step.
# As the Community model is not well behaved, these can look a little strange. Averaging over time will solve this problem.
plot(comm_sim0) 
# Zoom in on the biomass near the end of the simulation
plotBiomass(comm_sim0, start_time=700)
# Project again, this time with some fishing.
# The default is a knife-edge selectivity with only species larger than 1000g getting caught.
comm_sim1 <- project(comm_params, t_max=comm_time_sim, effort=1.5, dt=0.1)
# Again plot the results
plot(comm_sim1)
plotBiomass(comm_sim1, start_time=700)
# We want the average abundances over the last x steps to remove  the 'wrinkles' caused by the solver.
# Abundances are stored in the 'n' slot which contains a three dimensional array (time by species by size).
# We can use the abundances to calculate the relative abundances between the fished and unfished case at size:
avg_steps <- 400
# Mean abundance of the background resource with no fishing
mean_npp_comm0 <- apply(comm_sim0@n_pp[(comm_time_sim-avg_steps+2):(comm_time_sim+1),],2,mean)
# Mean abundance of the community with no fishing
mean_n_comm0 <- apply(comm_sim0@n[(comm_time_sim-avg_steps+2):(comm_time_sim+1),,],2,mean)
# Mean abundance of the community with fishing
mean_n_comm1 <- apply(comm_sim1@n[(comm_time_sim-avg_steps+2):(comm_time_sim+1),,],2,mean)
# Calculate the relative abundance
comm_relative_abundance <- mean_n_comm1 / mean_n_comm0

## Plot biomass spectra
#plot(comm_params@w, mean_n_comm0*comm_params@w, log="xy", type = "l")
#lines(comm_params@w, mean_n_comm1*comm_params@w, col="blue")
## Plot abundance spectra
#par(mfrow=c(2,1))
#plot(comm_params@w, mean_n_comm0, log="xy", type = "l")
#lines(comm_params@w, mean_n_comm1, col="blue")
## Plot relative abundance
#plot(comm_params@w, comm_relative_abundance, log="xy", type="l", ylim=c(0.1,10))
#lines(x=c(1e-10,1e10),y=c(1,1), lty=2)

# Make the trait-based model parameters object using the wrapper function: set_trait_model()
# k0 has been set so that the resource spectrum and the community spectra form a continuum.
trait_params <- set_trait_model(max_w=max_w_trait, no_sp=10, min_w_inf=5, max_w_inf=max_w_trait*0.9, n=n, q=q, f0=f0, h=h, ks=ks_trait, beta=beta, sigma=sigma, alpha=alpha_trait, r_pp=r_pp, z0pre=z0pre_trait, kappa=kappa, k0=5000)
# Project through time. The trait-based model is better behaved than the Community model so we don't need to project for so long.
trait_time_sim <- 200
# Project without fishing
trait_sim0 <- project(trait_params, t_max=trait_time_sim, effort=0, dt=1)
# Quick check - looks stable
plot(trait_sim0)
# Turn on fishing
# The default is a knife-edge selectivity with only species larger than 1000g getting caught.
trait_sim1 <- project(trait_params, t_max=trait_time_sim, effort=0.6, dt=1)
plot(trait_sim1)
# Although the Trait-based model is stable, we still get the time averaged abundances anyway.
avg_steps <- 50
# Mean abundance of the background resource with no fishing
mean_npp_trait0 <- apply(trait_sim0@n_pp[(trait_time_sim-avg_steps+2):(trait_time_sim+1),],2,mean)
# Mean abundance of the community with no fishing
mean_n_trait0 <- apply(trait_sim0@n[(trait_time_sim-avg_steps+2):(trait_time_sim+1),,],c(2,3),mean)
# Mean abundance of the community with fishing
mean_n_trait1 <- apply(trait_sim1@n[(trait_time_sim-avg_steps+2):(trait_time_sim+1),,],c(2,3),mean)
# Sum over the community to get the total abundances
mean_total_n_trait0 <- apply(mean_n_trait0,2,sum)
mean_total_n_trait1 <- apply(mean_n_trait1,2,sum)
# Calculate the relative abundances
trait_relative_abundance <-  mean_total_n_trait1 / mean_total_n_trait0

# Make the multispecies model based on the North Sea parameterisation
# This has different parameters to the Community and Trait-based models
kappa_ns <- 9.27e10
q_ns <- 0.8
n_ns <- 2/3
z0pre <- 0.6 
lambda_ns <- 2+q_ns-n_ns
# Load the species parameters (they come with mizer)
data(NS_species_params)
data(inter)
# Set the fishing gear to be a knife-edge with only species larger than 1000g getting caught.
NS_species_params$knife_edge_size <- 1000
# Make the MizerParams object
ms_params <- MizerParams(NS_species_params, inter, max_w=1e6, kappa=kappa_ns, q=q_ns, n=n_ns, z0pre=z0pre)
# Project through time. The multispcies model is normally well behaved.
ms_time_sim <- 150
# Project forward without any fishing
ms_sim0 <- project(ms_params, t_max=ms_time_sim, effort=0)
# Check the results
plot(ms_sim0)
# Now project forward with fishing
ms_sim1 <- project(ms_params, t_max=ms_time_sim, effort=1.0)
plot(ms_sim1)
# Although the multispecies model is stable, we still get the time averaged abundances anyway.
avg_steps <- 50
# Mean abundance of the background resource with no fishing
mean_npp_ms0 <- apply(ms_sim0@n_pp[(ms_time_sim-avg_steps+2):(ms_time_sim+1),],2,mean)
# Mean abundance of the community with no fishing
mean_n_ms0 <- apply(ms_sim0@n[(ms_time_sim-avg_steps+2):(ms_time_sim+1),,],c(2,3),mean)
# Mean abundance of the community with fishing
mean_n_ms1 <- apply(ms_sim1@n[(ms_time_sim-avg_steps+2):(ms_time_sim+1),,],c(2,3),mean)
# Sum over the community to get the total abundances
mean_total_n_ms0 <- apply(mean_n_ms0,2,sum)
mean_total_n_ms1 <- apply(mean_n_ms1,2,sum)
# Calculate the relative total abundance
ms_relative_abundance <-  mean_total_n_ms1 / mean_total_n_ms0

#----------------------------------------------------------
# Make FIGURE 1 for the paper
#----------------------------------------------------------

# Biomass
mean_b_comm0 <- mean_n_comm0 * comm_sim0@params@w
mean_b_comm1 <- mean_n_comm1 * comm_sim1@params@w
mean_total_b_trait0 <- mean_total_n_trait0 * trait_sim0@params@w
mean_total_b_trait1 <- mean_total_n_trait1 * trait_sim1@params@w
mean_total_b_ms0 <- mean_total_n_ms0 * ms_sim0@params@w
mean_total_b_ms1 <- mean_total_n_ms1 * ms_sim1@params@w
# Plot biomass relative to resource carrying capacity at 1e-3g
comm_base_b <- (comm_sim0@params@w_full * comm_sim0@params@cc_pp)[31]
trait_base_b <- (trait_sim0@params@w_full * trait_sim0@params@cc_pp)[31]
ms_base_b <- (ms_sim0@params@w_full * ms_sim0@params@cc_pp)[31]

width <- 14
height <- 7
#png(filename = "trophic_cascades.png", width = width, height = height, units="cm", res = 1000, pointsize=8)
#tiff(filename = "trophic_cascades.tiff", width = width, height = height, units="cm", res=500, pointsize=8)
postscript(file = "trophic_cascades.eps", width = width/2.5, height = height/2.5, pointsize=8, horizontal = FALSE, onefile=FALSE, paper='special')

nf <- layout(matrix(1:6,2,3,byrow=TRUE), rep(width/3,3),c(3/6,3/6)*height,TRUE)
#layout.show(nf)
#par(mfrow=c(2,3))

ylim <- c(1e-10,10)
xlim <- c(1e-3,1e4)
cascade_ylim <- c(0.1,10)
fat_lwd <- 2
fished_lty <- 3
resource_colour <- "green"
resource_lwd <- 1

# Community
par(mar=c(1,5,2,1))
plot(x=comm_sim0@params@w, y= (mean_b_comm0) / (comm_base_b), log="xy", type="n", ylab="Biomass relative to carrying capacity", xlim=xlim, ylim=ylim, main = "(a)", xlab="")
# Resource
lines(x=comm_sim0@params@w_full, y= (mean_npp_comm0 * comm_sim0@params@w_full) / comm_base_b, col=resource_colour, lwd=resource_lwd)
# Unfished community spectrum relative to unfished biomass in size 1
lines(x=comm_sim0@params@w, y=(mean_b_comm0) / (comm_base_b), col="red", lwd=fat_lwd)
# Fished community spectrum relative to unfished biomass in size 1
lines(x=comm_sim0@params@w, y=(mean_b_comm1) / (comm_base_b), col="blue", lwd=fat_lwd, lty=fished_lty)

# Trait
plot(x=trait_sim0@params@w, y= (mean_total_b_trait0) / (trait_base_b), log="xy", type="n",  ylab="", xlim=xlim, ylim=ylim, main = "(b)", xlab="")
# Resource
lines(x=trait_sim0@params@w_full, y= (mean_npp_trait0 * trait_sim0@params@w_full) / trait_base_b, col=resource_colour, lwd=resource_lwd)
# Unfished community spectrum relative to unfished biomass in size 1
lines(x=trait_sim0@params@w, y=(mean_total_b_trait0) / (trait_base_b), col="red", lwd=fat_lwd)
# Fished community spectrum relative to unfished biomass in size 1
lines(x=trait_sim1@params@w, y=(mean_total_b_trait1) / (trait_base_b), col="blue", lwd=fat_lwd, lty=fished_lty)
# Unfished species relative to unfished biomass in size1
for (i in 1:10){
    lines(x=trait_sim0@params@w, y=(mean_n_trait0[i,]*trait_sim0@params@w) / (trait_base_b))
}

# Multispecies
# Biomass
plot(x=ms_sim0@params@w, y= (mean_total_b_ms0) / ms_base_b, log="xy", type="n",  ylab="", xlim=xlim, ylim=ylim, main = "(c)", xlab="")
# Resource
lines(x=ms_sim0@params@w_full, y=(mean_npp_ms0 * ms_sim0@params@w_full) / ms_base_b, col=resource_colour, lwd=resource_lwd)
# Unfished community spectrum relative to unfished biomass in size 1
lines(x=ms_sim0@params@w, y=(mean_total_b_ms0) / ms_base_b, col="red", lwd=fat_lwd)
# Fished community spectrum relative to unfished biomass in size 1
lines(x=ms_sim1@params@w, y=(mean_total_b_ms1) / ms_base_b, col="blue", lwd=fat_lwd, lty=3)
# Unfished species relative to unfished biomass in size1
for (i in 1:12){
    lines(x=ms_sim0@params@w, y=(mean_n_ms0[i,]*ms_sim0@params@w) / (ms_base_b))
}

# Cascades
# Community
par(mar=c(5,5,5,1))
plot(x=comm_sim0@params@w, y=comm_relative_abundance, log="xy", type="n", ylab="Relative abundance", xlim=xlim, ylim=cascade_ylim, main = "(d)", xlab="Body mass (g)")
lines(x=comm_sim0@params@w, y=comm_relative_abundance)
lines(x=c(min(comm_sim0@params@w),max(comm_sim0@params@w)), y=c(1,1),lty=2)
lines(x=c(1000,1000),y=c(1e-20,1e20),lty=3)

# Trait
plot(x=trait_sim0@params@w, y=trait_relative_abundance, log="xy", type="n", ylab="", xlim=xlim, ylim=cascade_ylim, main = "(e)", xlab="Body mass (g)")
lines(x=trait_sim0@params@w, y=trait_relative_abundance)
lines(x=c(min(trait_sim0@params@w),max(trait_sim0@params@w)), y=c(1,1),lty=2)
lines(x=c(1000,1000),y=c(1e-20,1e20),lty=3)

# Multispecies
plot(x=ms_sim0@params@w, y=ms_relative_abundance, log="xy", type="n", xlab = "Body mass (g)", ylab="", xlim=xlim, ylim=cascade_ylim, main = "(f)")
lines(x=ms_sim0@params@w, y=ms_relative_abundance)
lines(x=c(min(ms_sim0@params@w),max(ms_sim0@params@w)), y=c(1,1),lty=2)
lines(x=c(1000,1000),y=c(1e-20,1e20),lty=3)

dev.off()


#----------------------------------------------------
# Code for FIGURE 2
#----------------------------------------------------

# Here we illustrate the performance of the model with fishing effort
# changing over time.
# We use the multispecies North Sea model with multiple gears.
# Initially, only the 

# Set up MS model. To have multiple gears. Run to equib.
# Project forward with time changing efforts.
# e.g. initially no demersal fishery, only pelagic.
# Then turn on demersal fishing
# Calc indicators: slope, MMW etc through time.
# Multi panel plot through time: fishing effort of gears, biomass of each species, indicators etc.

data(NS_species_params)
data(inter)

# Add an extra 'gears' column to the data set to specify the name of the fishing gear for each species.
# We have 4 gears: Industrial, Pelagic, Otter and Beam
NS_species_params$gear <- c("Industrial", "Industrial", "Industrial", "Pelagic", "Beam", "Otter", "Beam", "Otter", "Beam", "Otter", "Otter", "Otter")
# We also set the selectivity parameters - we still use knife edge 
# selectivty but the position of the 'edge' changes depending on the gear
NS_species_params$knife_edge_size <- NA
NS_species_params[NS_species_params$gear == "Industrial", "knife_edge_size"] <- 10
NS_species_params[NS_species_params$gear == "Pelagic", "knife_edge_size"] <- 20
NS_species_params[NS_species_params$gear == "Otter", "knife_edge_size"] <- 25 
NS_species_params[NS_species_params$gear == "Beam", "knife_edge_size"] <- 30
# Check the gear columns has been set correctly
NS_species_params
# Make the parameter object
ms_params <- MizerParams(NS_species_params, inter)

# Set up fishing effort for each gear through time as two dimensional array 
# Run for 100 years to get equilibrium
# No industrial fishing
time_to_equib <- 100


# No fishing to start with
equib_effort <- array(0, dim=c(time_to_equib, 4), dimnames=list(time=1:time_to_equib, gear=c("Industrial", "Pelagic", "Beam", "Otter")))
#equib_effort <- array(NA, dim=c(time_to_equib, 4), dimnames=list(time=1:time_to_equib, gear=c("Industrial", "Pelagic", "Beam", "Otter")))
#equib_effort[,"Industrial"] <- 0
#equib_effort[,"Pelagic"] <- 0.2
#equib_effort[,"Otter"] <- 0.2
#equib_effort[,"Beam"] <- 0.2
# Project to equilibrium
ms_equib <- project(ms_params, effort=equib_effort)
# Plot everything
plot(ms_equib)
# Pull out the equilibrium population abundances - we will use them as the initial abundances for the projection
n_equib <- ms_equib@n[time_to_equib+1,,]
n_pp_equib <- ms_equib@n_pp[time_to_equib+1,]

# Set up fishing history
# After 10 years a pelagic fishery starts (increases from 0 to 1 over 10 yrs)
# After 20 years a beam fishery starts (increases from 0 to 0.5 over 10 yrs)
# After 30 years an otter fishery starts (increases from 0 to 0.5 over 10 yrs)
# After 40 years an industrial fishery starts (increases from 0 to 0.5 over 10 yrs)

project_time <- 100
fishing_effort <- array(0, dim=c(project_time, 4), dimnames=list(time=1:project_time, gear=c("Industrial", "Pelagic", "Beam", "Otter")))
fishing_effort[,"Pelagic"] <- c(rep(0,10),seq(from = 0, to = 1, length = 10), rep(1,80))
fishing_effort[,"Beam"] <- c(rep(0,30),seq(from = 0, to = 1, length = 10), rep(1,60))
fishing_effort[,"Otter"] <- c(rep(0,50),seq(from = 0, to = 1, length = 10), rep(1,40))
fishing_effort[,"Industrial"] <- c(rep(0,70),seq(from = 0, to = 1, length = 10), rep(1,20))
# Have a quick look
plot(x = 1:project_time, y = seq(from=0,to=1,length=project_time), type="n", xlab = "Years", ylab="Fishing effort")
for (i in 1:4){
    lines(x=1:project_time, y = fishing_effort[,i], lty=i)
}
legend(x="bottomright",legend=c("Industrial", "Pelagic", "Beam", "Otter"), lty=1:4)
# What happens
ms_sim <- project(ms_params, effort=fishing_effort, initial_n = n_equib, initial_n_pp = n_pp_equib)
plotBiomass(ms_sim)
# And Yield?
yield <- getYieldGear(ms_sim)
# Melt this down for easy plotting with ggplot2
library(reshape2)
yield_melt <- melt(yield)
ggplot(yield_melt) + geom_line(aes(x=time, y=value, colour = sp)) + scale_y_log10()
# But remove 0 rows
yield_melt <- yield_melt[yield_melt$value >0,]
library(ggplot2)
ggplot(yield_melt) + geom_line(aes(x=time, y=value, linetype = sp, colour = gear)) + scale_y_log10()

fishing_effort[,"Pelagic"] <- c(rep(0.2,30), rep(0.5,70))
fishing_effort[,"Otter"] <- 0.2
fishing_effort[,"Beam"] <- c(rep(0.2,50),rep(0.6,50))



#fishing_effort <- array(NA, dim=c(project_time, 4), dimnames=list(time=1:project_time, gear=c("Industrial", "Pelagic", "Beam", "Otter")))
#fishing_effort[,"Industrial"] <- c(rep(0,10),rep(0.8,90))
#fishing_effort[,"Pelagic"] <- c(rep(0.2,30), rep(0.5,70))
#fishing_effort[,"Otter"] <- 0.2
#fishing_effort[,"Beam"] <- c(rep(0.2,50),rep(0.6,50))
ms_sim <- project(ms_params, effort=fishing_effort, initial_n = n_equib, initial_n_pp = n_pp_equib)

# Plots
# Get indicators - slope, mean max weight, mean weight, lfi, yield, ssb throught time
# size range
min_w <- 10
max_w <- 5000
threshold_w <- 100

ssb <- getSSB(ms_sim)
yield <- getYield(ms_sim)
lfi <- getProportionOfLargeFish(ms_sim, min_w=min_w, max_w=max_w, threshold_w=threshold_w)
mw <- getMeanWeight(ms_sim, min_w=min_w, max_w=max_w)
mmw <- getMeanMaxWeight(ms_sim, min_w=min_w, max_w=max_w, measure="biomass")
slope <- getCommunitySlope(ms_sim, min_w=min_w, max_w=max_w)[,"slope"]

#png(filename="time_series.png", width = 7, height = 14, units="cm",res=800)
#par(mfcol =c(7,1))
nf <- layout(matrix(1:7,7,1,byrow=TRUE), rep(7,7),c(1.2,rep(1,5),1.8),TRUE)
#layout.show(nf)
# Effort of gears
par(mar=c(0,4,0.5,1))
plot(x = 1:project_time, y=1:project_time, type="n", ylim=c(0,1), xlab="", ylab="Effort", xaxt="n")
for (i in 1:4){
    lines(x = 1:project_time, y=ms_sim@effort[,i], col=i)
}
# ssb
par(mar=c(0,4,0,1))
plot(x = 1:project_time, y=1:project_time, type="n", ylim=c(min(ssb),max(ssb)),log="y", ylab="SSB", xlab="", xaxt="n")
for (i in 1:12){
    lines(x = 1:project_time, y=ssb[2:(project_time+1),i])
}
# yield
par(mar=c(0,4,0,1))
plot(x = 1:project_time, y=1:project_time, type="n", ylim=c(min(yield[11:100,]),max(yield[11:100,])),log="y", ylab="Yield", xlab="", xaxt="n")
for (i in 1:12){
    lines(x = 1:project_time, y=yield[1:project_time,i])
}
# lfi
par(mar=c(0,4,0,1))
plot(x = 1:project_time, y=lfi[2:(project_time+1)], type="l", ylab="LFI", xlab="", xaxt="n")
# mw
par(mar=c(0,4,0,1))
plot(x = 1:project_time, y=mw[2:(project_time+1)], type="l", ylab="Mean weight", xlab="", xaxt="n")
# mmw
par(mar=c(0,4,0,1))
plot(x = 1:project_time, y=mmw[2:(project_time+1)], type="l", ylab="Mean max weight", xlab="", xaxt="n")
# slope
par(mar=c(5,4,0,1))
plot(x = 1:project_time, y=slope[2:(project_time+1)], type="l", ylab="Slope", xlab="Years")
#dev.off()







#----------------------------------------------------
# COMMUNITY MODEL EXAMPLE
#----------------------------------------------------

# Create a MizerParams object using the default values
comm_params <- set_community_model()
# Project through time with and without fishing effort
# The time step of the simulations is 1, and the results are
# stored every time step (default values).
comm_sim0 <- project(comm_params, t_max=50, effort=0)
comm_sim1 <- project(comm_params, t_max=50, effort=1)
# We want the abundances in the final time step
# Abundances are stored in the 'n' slot
# which contains a three dimensional array (time by species by size).
# We can use the abundances to calculate the relative abundances at size:
comm_relative_abundance <- comm_sim1@n[51,,] / comm_sim0@n[51,,]

# Plot FIGURE 1 - a trophic cascade
plot(x=comm_sim0@params@w, y=comm_relative_abundance, log="x", type="n", xlab = "Size (g)", ylab="Relative abundance")
lines(x=comm_sim0@params@w, y=comm_relative_abundance)
lines(x=c(min(comm_sim0@params@w),max(comm_sim0@params@w)), y=c(1,1),lty=2)

#----------------------------------------------------
# TRAIT BASED MODEL EXAMPLE
#----------------------------------------------------

# Set up a trait-based model with 10 species, with asymptotic sizes ranging from 10 g to 100 kg. All the other parameters have default values. 
trait_params <- set_trait_model(no_sp=10, min_w_inf=10, max_w_inf=1e5)
# Project forward for 100 time steps by which time the community has reached equilibrium:
# Without fishing
trait_sim0 <- project(trait_params, t_max=100, effort=0)
# With fishing
trait_sim1 <- project(trait_params, t_max=100, effort=0.75)
# Plot the spectra FIGURE 2
plotSpectra(trait_sim1, biomass=FALSE) 
# For the paper - some tweaks to the plot
p <- plotSpectra(trait_sim1, biomass=FALSE)
p <- p + theme_bw() + theme(legend.position="none")
ggsave("figure2.png", plot = p, width = 7, height = 7, units = "cm", scale = 2)
# 

# Calculate the slopes for individuals between 10 g and 10 kg
slope0 <- getCommunitySlope(trait_sim0, min_w=10, max_w=10e3)
slope1 <- getCommunitySlope(trait_sim1, min_w=10, max_w=10e3)

#----------------------------------------------------
# MULTISPECIES EXAMPLE
#----------------------------------------------------

# Load the multispecies data set and interaction matrix
data(NS_species_params)
data(inter)

# Add an extra 'gears' column to the data set to specify the gear name for each species
NS_species_params$gear <- c("Industrial", "Industrial", "Industrial", "Pelagic", "Beam", "Otter", "Beam", "Otter", "Beam", "Otter", "Otter", "Otter")
# Check the gear columns has been set correctly
NS_species_params
# Make the MizerParams object
ms_params <- MizerParams(NS_species_params, inter)
# Set up fishing effort for each gear through time as a two dimensional array 
# This will be used to perform a 40 year projection, with the first 20 years as transients
fishing_effort <- array(NA, dim=c(40, 4), dimnames=list(time=1:40, gear=c("Industrial", "Pelagic", "Beam", "Otter")))
# Set the fishing effort of each of the gears through time
fishing_effort[,"Industrial"] <- c(rep(0.5,20), seq(from=0.5, to=2, length=20))
fishing_effort[,"Pelagic"] <- c(rep(1,20), seq(from=1, to=2, length=20))
fishing_effort[,"Beam"] <- c(rep(2,20), seq(from=2, to=0.5, length=20))
fishing_effort[,"Otter"] <- c(rep(1,20), seq(from=2, to=0.5, length=20))
# And project
ms_sim <- project(ms_params, effort=fishing_effort)
# Plot everything
plot(ms_sim)

# Plot for FIGURE 3
# Pull out detail from summary methods
bm <- getBiomass(ms_sim)
fl <- getFeedingLevel(ms_sim)[41,,]
m2 <- getM2(ms_sim)[41,,]
fm <- getFMort(ms_sim)[40,,]
# Helpful variables for the plot
sim_time <- 21:40
max_bm <- max(bm[sim_time,])
min_bm <- min(bm[sim_time,])
size <- ms_sim@params@w
# Plot formatting variables
cols <- rep(c("black","blue","red","green"),each = 3)
ltys <- rep(c(1,3,4),4)
# Set up the plot layout
layout_matrix <- matrix(c(1,2,3,4), byrow=TRUE, ncol = 1)
nf <- layout(layout_matrix, widths = 1, heights = c(1.3,1,0.8,1.2),respect = FALSE)
# Then figure out  margins: B, L, T, R
par(mar=c(4,5,2,2))
plot(x = sim_time, y = rep(1, length(sim_time)), type="n", ylim=c(min_bm,max_bm), log = "y", ylab = "Biomass", xlab = "Time (yr)")
for (i in 1:dim(bm)[2]){
    lines(x = sim_time, y = bm[sim_time+1,i], col = cols[i], lty = ltys[i])
}
par(mar=c(0,5,2,2))
plot(x=size, y = rep(1,length(size)), type="n", ylim = c(min(fl), max(fl)), log="x", ylab = "Feeding level", xlab = "", xaxt="n")
    for (i in 1:dim(fl)[1]){
        lines(x = size, y = fl[i,], col = cols[i], lty = ltys[i])
    }
par(mar=c(0,5,0,2))
plot(x=size, y = rep(1,length(size)), type="n", ylim = c(min(m2), max(m2)), log="x", ylab = "M2", xlab = "", xaxt="n")
    for (i in 1:dim(m2)[1]){
        lines(x = size, y = m2[i,], col = cols[i], lty = ltys[i])
    }
par(mar=c(5,5,0,2))
plot(x=size, y = rep(1,length(size)), type="n", ylim = c(min(fm), max(fm)), log="x", ylab = "Fishing mortality", xlab = "Size (g)")
    for (i in 1:dim(fm)[1]){
        lines(x = size, y = fm[i,], col = cols[i], lty = ltys[i])
    }
legend(x="left",legend = dimnames(bm)$sp, col = cols, lty=ltys, ncol = 2, bty="n")







