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
rm(list=ls()) # start afresh
library(mizer)
data(NS_species_params)
data(inter)

# Add an extra 'gears' column to the data set to specify the name of the fishing gear for each species.
# We have 4 gears: Industrial, Pelagic, Otter and Beam
NS_species_params$gear <- c("Industrial", "Industrial", "Industrial", "Pelagic", "Beam", "Otter", "Beam", "Otter", "Beam", "Otter", "Otter", "Otter")
# We also set the selectivity parameters - we still use knife edge selectivty but the position of the 'edge' changes depending on the gear.
# Size is specified as mass. Assume a L-W relationship 
# w = a * l ^ b with a = 0.01 and b = 3
NS_species_params$knife_edge_size <- NA
NS_species_params[NS_species_params$gear == "Industrial", "knife_edge_size"] <- 10 # l = 10
NS_species_params[NS_species_params$gear == "Pelagic", "knife_edge_size"] <- 80 # l = 20
NS_species_params[NS_species_params$gear == "Otter", "knife_edge_size"] <- 270 # l = 30 
NS_species_params[NS_species_params$gear == "Beam", "knife_edge_size"] <- 155 # l = 25 
# It would be more sophisticated to base the selectivity on working group assumptions and also change
# the selectivity ogive. But for simplicity here we just use knife-edge.
# Check the gear columns has been set correctly
NS_species_params
# Make the parameter object
NS_params <- MizerParams(NS_species_params, inter)
# Out of interest we can see that the fishing catchability has been correctly set
# i.e. gears that do not catch a species have a catchability of 0
NS_params@catchability

# We want our simulation to start at the 'unfished' state so to start with we run for 100 years with no fishing effort to get stocks at equilibrium.
time_to_equib <- 100
# Project to equilibrium (with default timestep of 0.1)
NS_equib <- project(NS_params, effort=0, t_max = time_to_equib)
# Plot everything - check we are at equilibrium
plot(NS_equib)
# Pull out the equilibrium population abundances - we will use them as the initial abundances for the projection
n_equib <- NS_equib@n[time_to_equib+1,,]
n_pp_equib <- NS_equib@n_pp[time_to_equib+1,]

# Set up fishing effort for each gear.
# Here, we want the fishing effort of each gear to change through time as additional gears start operating.
# We therefore set up the fishing effort as a two dimensional array: time by gear
# Set up fishing history:
# After 10 years a pelagic fishery starts (effort increases from 0 to 1 over 10 yrs)
# After 30 years a beam fishery starts (increases from 0 to 0.75 over 10 yrs)
# After 50 years an otter fishery starts (increases from 0 to 0.9 over 10 yrs)
# After 70 years an industrial fishery starts (increases from 0 to 1.5 over 10 yrs)
gear_names <- c("Pelagic","Beam","Otter","Industrial")
gear_lty <- 1:4
names(gear_lty) <- gear_names
project_time <- 100
fishing_effort <- array(0, dim=c(project_time, 4), dimnames=list(time=1:project_time, gear=gear_names))
fishing_effort[,"Pelagic"] <- c(rep(0,10),seq(from = 0, to = 1, length = 10), rep(1,80))
fishing_effort[,"Beam"] <- c(rep(0,30),seq(from = 0, to = 0.75, length = 10), rep(0.75,60))
fishing_effort[,"Otter"] <- c(rep(0,50),seq(from = 0, to = 0.9, length = 10), rep(0.9,40))
fishing_effort[,"Industrial"] <- c(rep(0,70),seq(from = 0, to = 1.5, length = 10), rep(1.5,20))
# Have a quick look at the fishing effort by gear
plot(x = 1:project_time, y = seq(from=0,to=1,length=project_time), type="n", xlab = "Years", ylab="Fishing effort", ylim=c(0,1.5))
for (i in gear_names){
    lines(x=1:project_time, y = fishing_effort[,i], lty=gear_lty[i])
}
legend(x="bottomright",legend=gear_names, lty=gear_lty)
# Run the simulation, passing in the effort array and the initial population abundances
NS_sim <- project(NS_params, effort=fishing_effort, initial_n = n_equib, initial_n_pp = n_pp_equib)
# Plot everything - what happened?
plot(NS_sim)
# A closer look at the biomass through time - can we see when the different gears start?
plotBiomass(NS_sim)

## And Yield?
#yield <- getYieldGear(NS_sim)
## Do some tidying up as not all gears catch all species
## To help with the plot set yields of non-caught species to NA
#catchability <- NS_sim@params@catchability
#catchability[catchability == 0] <- NA
#yield <- sweep(yield, c(2,3), catchability, "*")
## Melt this down for easy plotting with ggplot2
#library(reshape2)
#library(ggplot2)
# yield_melt <- melt(yield)
## Do some tidying up because not all gears caught all species
## But remove 0 rows
#ggplot(yield_melt) + geom_line(aes(x=time, y=value, linetype = sp, colour = gear)) + scale_y_log10()

# Plot against time
# a) fishing effort by gear
# b) SSB - species and gear
# c) Yield - species and gear
# d) Maximum weight
# e) and Mean maximum weight
# f) Slope
# g) LFI

# Size range for calculating community metrics
min_w <- 10
max_w <- 5000
threshold_w <- 100 # for large fish indicator - what is a large fish?
# Get the data to be plotted
ssb <- getSSB(NS_sim)
# rescale ssb to be relative to the unfished ssb
rescale_ssb <- sweep(ssb,2,ssb[1,],"/")

yield <- getYieldGear(NS_sim)
# rescale yield relative to the maximum yield over the time series
max_yield <- apply(yield,c(2,3),max)
rescale_yield <- sweep(yield,c(2,3), max_yield, "/")
# Do some tidying up of yield as not all gears catch all species
# To help with the plot set yields of non-caught species to NA
yield[yield==0] <- NA
rescale_yield[rescale_yield==0] <- NA

# yield[yield==0] <- NA
#catchability <- NS_sim@params@catchability
#catchability[catchability == 0] <- NA
#yield <- sweep(yield, c(2,3), catchability, "*")
# Community metrics 
lfi <- getProportionOfLargeFish(NS_sim, min_w=min_w, max_w=max_w, threshold_w=threshold_w)
mw <- getMeanWeight(NS_sim, min_w=min_w, max_w=max_w)
mmw <- getMeanMaxWeight(NS_sim, min_w=min_w, max_w=max_w, measure="biomass")
slope <- getCommunitySlope(NS_sim, min_w=min_w, max_w=max_w)[,"slope"]

# Scale these relative to the unfished community
rescale_lfi <- lfi / lfi[1]
rescale_mw <- mw / mw[1]
rescale_mmw <- mmw / mmw[1]

# Set colours and lty for the species
# Pelagic Beam Otter Industrial  Catching 12 species
# Each species has same lty as the gear that catches it
# Each species within a gear has a different colour
species_names <- as.character(NS_params@species_params$species)
species_lty <- rep(NA,12)
names(species_lty) <- species_names
#cols <- 1:5
cols <- c("black","blue","magenta","green","red")
species_col <- rep(NA,12)
names(species_col) <- species_names
for (i in gear_names){
    gear_idx <- (NS_params@species_params$gear == i)
    species_col[NS_params@species_params$species[gear_idx]] <- cols[1:(sum(gear_idx))]
    species_lty[NS_params@species_params$species[gear_idx]] <- gear_lty[i]
}


#------------------------------------------------------------------
# FIGURE 2 

# Figure size
width <- 7
height <- 20

# Function to add fishing effort lines to the time plots
add_effort_lines <- function(){
    lwd <- 0.5
    lines(x=c(11,11), y=c(-1e20,1e20), lty=gear_lty[1], lwd=lwd)
    lines(x=c(31,31), y=c(-1e20,1e20), lty=gear_lty[2], lwd=lwd)
    lines(x=c(51,51), y=c(-1e20,1e20), lty=gear_lty[3], lwd=lwd)
    lines(x=c(71,71), y=c(-1e20,1e20), lty=gear_lty[4], lwd=lwd)
}

png(filename="Figure212.png", width = width, height = height, units="cm",res=800, pointsize=8)

#rel_heights <- c(1,2,2,2,2,0.5) # rel heights of panels
#rel_heights <- c(1,0.5,0.5,0.5,0.5,2,2,2,0.5) # rel heights of panels
#rel_heights <- c(1,rep(0.5,8),0.5,1,2,0.5) # rel heights of panels
rel_heights <- c(0.7,rep(0.5,8),0.5,0.8,1) # rel heights of panels

heights = (height / sum(rel_heights)) * rel_heights
nf <- layout(matrix(1:length(rel_heights),length(rel_heights),1,byrow=TRUE), widths = width, heights=heights,TRUE)
right_margin <- 4
left_margin <- 4
# Other plotting parameters
legend_txt_cex <- 0.7
leg_line <-  0.3
seg_len <- 2.5
leg_bty <- "n"
leg_box_lwd <- 0

# (a) Effort of gears
par(mar=c(0,left_margin,0.5,right_margin))
plot(x = 1:project_time, y=1:project_time, type="n", ylim=c(0,max(fishing_effort)), xlab="", ylab="Effort", xaxt="n")
for (i in gear_names){
    lines(x = 1:project_time, y=NS_sim@effort[,i], lty=gear_lty[i])
}
add_effort_lines()
legend(x="bottomright", legend = gear_names, lty=gear_lty, cex=legend_txt_cex, pt.lwd=leg_line, seg.len=seg_len, bty=leg_bty, box.lwd=leg_box_lwd)
# legend(x="topleft", legend = gear_names, lty=gear_lty, cex=legend_txt_cex, pt.lwd=leg_line, seg.len=seg_len, bty=leg_bty, box.lwd=leg_box_lwd, ncol=2)

# (b) yield
#par(mar=c(0,left_margin,0,right_margin))
#yield_min <- min(yield[yield>0], na.rm=TRUE)
#yield_max <- max(yield[yield>0], na.rm=TRUE)
#rescale_yield_min <- min(rescale_yield[rescale_yield>0], na.rm=TRUE)
#rescale_yield_max <- max(rescale_yield[rescale_yield>0], na.rm=TRUE)
#plot(x = 1:project_time, y=1:project_time, type="n", ylim=c(rescale_yield_min, rescale_yield_max), ylab="Yield", xlab="", xaxt="n")
#for (gear in gear_names){
#    for (i in species_names){
#        lines(x = 1:project_time, y=rescale_yield[1:project_time,gear,i], col=species_col[i], lty=species_lty[i])
#    }
#}
#add_effort_lines()

# (b) alternative yield - each gear has a separate panel
rescale_yield_min <- min(rescale_yield[rescale_yield>0], na.rm=TRUE)
rescale_yield_max <- max(rescale_yield[rescale_yield>0], na.rm=TRUE)
for (gear in gear_names){
    species_in_gear <- NS_params@species_params$species[NS_params@species_params$gear==gear]
    par(mar=c(0,left_margin,0,right_margin))
    plot(x = 1:project_time, y=1:project_time, type="n", ylim=c(rescale_yield_min, rescale_yield_max), ylab="", xlab="", xaxt="n", yaxt="n")
    axis(4)
    mtext(gear, side=2, line=1, cex=0.7)
    if (gear == gear_names[2]){
        mtext("Relative Yield", side=4, line=3, cex=0.7, adj=-3)
    }
    add_effort_lines()
    for (i in species_in_gear){
        lines(x = 1:project_time, y=rescale_yield[1:project_time,gear,i], col=species_col[i], lty=species_lty[i])
    }
    legend(x="bottomright", legend=species_in_gear, lty=species_lty[species_in_gear], col=species_col[species_in_gear], cex = legend_txt_cex, ncol=1, pt.lwd=leg_line, seg.len = seg_len, bty=leg_bty, box.lwd=leg_box_lwd)
}

# (c) Relative SSB
#par(mar=c(0,left_margin,0,right_margin))
#plot(x = 1:project_time, y=1:project_time, type="n", ylim=c(min(rescale_ssb),max(rescale_ssb)), ylab="Relative SSB", xlab="", xaxt="n")
#for (i in species_names){
#    lines(x = 1:project_time, y=rescale_ssb[2:(project_time+1),i], col=species_col[i], lty=species_lty[i])
#}
#add_effort_lines()
#legend(x="topleft", legend=species_names, lty=species_lty, col = species_col, bty='o', cex = legend_txt_cex, ncol=2, pt.lwd=leg_line, seg.len = seg_len)
# (c) alternative - each gear has a separate panel
#par(mfrow=c(4,1))
for (gear in gear_names){
    species_in_gear <- NS_params@species_params$species[NS_params@species_params$gear==gear]
    par(mar=c(0,left_margin,0,right_margin))
    plot(x = 1:project_time, y=1:project_time, type="n", ylim=c(min(rescale_ssb),max(rescale_ssb)), ylab="", xlab="", xaxt="n")
    if (gear == gear_names[2]){
        mtext("Relative SSB", side=2, line=3, cex=0.7, adj=-3)
    }
    mtext(gear, side=4, line=1, cex=0.7)
    add_effort_lines()
        for (i in species_in_gear){
            lines(x = 1:project_time, y=rescale_ssb[2:(project_time+1),i], col=species_col[i], lty=species_lty[i])
        }
}

# (d) Slope
par(mar=c(0,left_margin,0,right_margin))
ylim <- range(slope)
plot(x = 1:project_time, y=1:project_time, type="n", ylab="", xlab="Years", ylim=ylim, yaxt="n", xaxt="n")
axis(4)
mtext("Community slope", side=4, line=3, cex=0.7)
add_effort_lines()
lines(x = 1:project_time, y = slope[2:(project_time+1)])

# (e) LFI, MW, MMW, 
#slope_scale <- 1 / max(abs(slope[2:(project_time+1)]))
#slope_plot <- slope[2:(project_time+1)]
#slope_adj <- abs(min(slope_plot))
#slope_plot <- slope_plot + slope_adj
#slope_scale <- 1 / max(slope_plot)
par(mar=c(4,left_margin,0,right_margin))
ylim <- c(0,max(rescale_lfi,rescale_mw, rescale_mmw, slope))
plot(x = 1:project_time, y=1:project_time, type="n", ylab="Relative metrics", xlab="Years", ylim=ylim)
lines(x = 1:project_time, y = rescale_lfi[2:(project_time+1)], col=1)
lines(x = 1:project_time, y = rescale_mw[2:(project_time+1)], col=2)
lines(x = 1:project_time, y = rescale_mmw[2:(project_time+1)], col=3)
#lines(x = 1:project_time, y = slope[2:(project_time+1)], col=4)
# lines(x = 1:project_time, y = slope_plot * slope_scale, col=4)
add_effort_lines()
legend(x="bottomright", legend = c("LFI", "MW", "MMW"), lty=1, col=c(1,2,3), cex = legend_txt_cex, ncol=1, pt.lwd=leg_line, seg.len = seg_len, bty=leg_bty, box.lwd=leg_box_lwd)

# (f) Spectra at various points in time
# Size distribution at start and end - five lines
# times of something happening
# 20, 40, 60, 80
b0 <- apply(NS_sim@n[2,,],2,sum) * NS_sim@params@w
b1 <- apply(NS_sim@n[21,,],2,sum) * NS_sim@params@w
b2 <- apply(NS_sim@n[41,,],2,sum) * NS_sim@params@w
b3 <- apply(NS_sim@n[61,,],2,sum) * NS_sim@params@w
b4 <- apply(NS_sim@n[81,,],2,sum) * NS_sim@params@w
xlim <- c(1,5e4)
ylim <- c(5e5,max(c(b0,b1,b2,b3,b4)))
par(mar=c(4,left_margin,1,right_margin))
plot(x=NS_sim@params@w, y = NS_sim@params@w, type="n", ylab="Total biomass", xlab = "Size (g)", log="xy", ylim = ylim, xlim=xlim)
lines(x=NS_sim@params@w, y = b0, col=1)
lines(x=NS_sim@params@w, y = b1, col=2)
lines(x=NS_sim@params@w, y = b2, col=3)
lines(x=NS_sim@params@w, y = b3, col=4)
lines(x=NS_sim@params@w, y = b4, col=6)
legend(x="bottomleft", legend=c("Unfished", "Year 20: Pelagic", "Year 40: Beam", "Year 60: Otter", "Year 80: Industrial"), lty=1, col=c(1,2,3,4,6), cex = legend_txt_cex, pt.lwd=leg_line, seg.len = seg_len, bty=leg_bty, box.lwd=leg_box_lwd)
dev.off()

