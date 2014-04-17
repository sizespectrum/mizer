# Supplementary material for
# mizer: an R package for multispecies, trait-based and community size spectrum ecological modelling
# Finlay Scott, Julia L. Blanchard and Ken H. Andersen
# Methods in Ecology and Evolution

#----------------------------------------------------
# Preliminaries
#----------------------------------------------------

# Loading the library and clearing out workspace

library(mizer)
rm(list=ls())

#----------------------------------------------------
# Example 1: Spectra and Tropic cascades
#----------------------------------------------------

# We use the three types of model in mizer:
#     the community model
#     the trait-based model
#     and the multispecies model (parameterised for the North Sea)
# We explore the impacts of fishing by projecting each model forward in time with and without fishing, then calculate the change in abundances.
# The results are used to create Figure 1 in the paper
# For each model we plot:
# In the first row:
#     the unfished biomass spectra of each species
#     the unfished total community biomass spectrum
#     the fished total community biomass spectrum
#     the unfished background resource spectrum
# In the second row:
#     the relative abundances in the fished and unfished case (to simulate trophic cascades).

# Make the Community model parameters object using the wrapper function: set_community_model()
# (See the package vignette for a description of the parameters.)
comm_params <- set_community_model(max_w=50e3, # spectra parameters
                                   beta=100, sigma=1.3, # size selection parameters
                                   h=20, alpha=0.17,  n=0.75, q=0.8, z0=0.2, f0=0.6,
                                   r_pp=4, kappa=5, # resource parameters
                                   rec_mult=0.01)
# The value of rec_mult has been set so that N in the smallest size is approximately the same as the resource spectrum in the smallest size (found through trial and error).
# Project the Community model to equilibrium
# The solver is not well behaved with the Community model and it can be quite unstable.
# We therefore need ro run a long projection and average over time series.
comm_time_sim <- 800 
# First we project the community model through time without fishing (effort is set to 0).
# The time step of the simulations is 0.1, and the results are stored every time step (default values).
comm_sim0 <- project(comm_params, t_max=comm_time_sim, effort=0, dt=0.1)
# Plot the results as a quick check (plot() has been overloaded for MizerSim objects)
plot(comm_sim0) 
# Note that the plots with size on the x-axis are snap shots of the final time step.
# As the Community model is not well behaved, these can look a little strange.
# Averaging over time will solve this problem.
# Zoom in on the biomass near the end of the simulation
plotBiomass(comm_sim0, start_time=700)
# Project again, this time with some fishing.
# The default selectivity is a knife-edge selectivity with only species larger than 1000g getting caught (see vignette for examples of changing the selectivity of the fishing gear).
comm_sim1 <- project(comm_params, t_max=comm_time_sim, effort=1.5, dt=0.1)
# Again, plot the results
plot(comm_sim1)
plotBiomass(comm_sim1, start_time=700)
# Solver is still unstable. We will take the mean later on.

# Make the trait-based model parameters object using the wrapper function set_trait_model()
# Use similar parameters as the Community model
trait_params <- set_trait_model(max_w=100e3, no_sp=10, min_w_inf=5, max_w_inf=100e3*0.9, # spectra parameters
                                beta=100, sigma=1.3, # size selection parameters
                                h=20, alpha=0.6, n=0.75, q=0.8, z0pre=2, f0=0.6, ks=2.4,
                                r_pp=4, kappa=5, # resource parameters
                                k0=5000) # recruitment parameters
# k0 has been set so that the resource spectrum and the community spectra form a continuum (found through trial and error).
# Project through time to equilibrium. The trait-based model is better behaved than the Community model so we don't need to project for so long.
trait_time_sim <- 200
# Project without fishing (time step is 1 as it is better behaved)
trait_sim0 <- project(trait_params, t_max=trait_time_sim, effort=0, dt=1)
# Plotting shows that a stable looking equilibrium has been reached 
plot(trait_sim0)
plotBiomass(trait_sim0, start_time=150)
# Turn on fishing
# The default is a knife-edge selectivity with only species larger than 1000g getting caught.
trait_sim1 <- project(trait_params, t_max=trait_time_sim, effort=0.6, dt=1)
plot(trait_sim1)
# Looks different to the unfished case
plotBiomass(trait_sim1, start_time=150)

# Make the multispecies model based on the North Sea parameterisation
# This has different parameters to the Community and Trait-based models 
# Load the species parameters (they come with mizer) and the interaction matrix (see Blanchard et al., 2014)
data(NS_species_params)
data(inter)
# Add a column to the species parameters so that the fishing gear is a knife-edge with only species larger than 1000g getting caught.
NS_species_params$knife_edge_size <- 1000
# Make the MizerParams object using the constructor
ms_params <- MizerParams(NS_species_params, inter, # the species parameters and interaction matrix
                         max_w=1e6, kappa=9.27e10, q=0.8, n=2/3, z0pre=0.6)
# There are lots of note explaining that many parameters have been set to the default values
# Project through time to equilibrium. The multispcies model is normally well behaved.
ms_time_sim <- 150
# Project forward without any fishing
ms_sim0 <- project(ms_params, t_max=ms_time_sim, effort=0)
# Check the results
plot(ms_sim0)
# Now project forward with fishing
ms_sim1 <- project(ms_params, t_max=ms_time_sim, effort=0.75)
plot(ms_sim1)

# Process the model results for the Figure 1
# For the Community model we want the average abundances over the last x steps to remove the 'wrinkles' caused by the solver.
# Abundances are stored in the 'n' slot of the MizerSim object. Abundances are stored as a three dimensional array (time by species by size).
# For example:
dim(comm_sim1@n) # Only 1 species in the Community model
dimnames(comm_sim1@n)
# The weights of individuals in the community can be found in the @params@w slot
comm_sim1@params@w # Just the community
comm_sim1@params@w_full # Full community with background spectrum

# We want to calculate average abundances (biomass) for the fished and unfished simulations.
# Then we want to look at the relative abundances.
# We can use the abundances to calculate the relative abundances between the fished and unfished case at size:
no_avg_steps <- 400 # no of steps to average over
avg_steps_indices <- (comm_time_sim-no_avg_steps+2):(comm_time_sim+1) # the range to average over (note that the first step is the inital abundances so we ignore that step).
# Mean abundance of the background resource with no fishing
comm0_npp <- apply(comm_sim0@n_pp[avg_steps_indices,],2,mean) * comm_sim0@params@w_full
# Mean abundance of the community with no fishing
comm0_n <- apply(comm_sim0@n[avg_steps_indices,,],2,mean) * comm_sim0@params@w
# Mean abundance of the community with fishing
comm1_n <- apply(comm_sim1@n[avg_steps_indices,,],2,mean) * comm_sim1@params@w
# Calculate the relative total abundance
comm_relative_abundance <- comm1_n / comm0_n

# The Trait-based model is stable so there is no need to average the abundances
# Abundance of the background resource with no fishing in the final time step (the first time step in the initial abundance)
trait0_npp <- trait_sim0@n_pp[trait_time_sim+1,] * trait_sim0@params@w_full
# Abundance of the community with no fishing in the final time step
trait0_n <- sweep(trait_sim0@n[trait_time_sim+1,,],2,trait_sim0@params@w,"*")
# Abundance of the community with fishing in the final time step
trait1_n <- sweep(trait_sim1@n[trait_time_sim+1,,],2,trait_sim1@params@w,"*")
# Sum over the community to get the total abundances
trait0_n_total <- apply(trait0_n,2,sum)
trait1_n_total <- apply(trait1_n,2,sum)
# Calculate the relative total abundances
trait_relative_abundance <-  trait1_n_total / trait0_n_total

# The multispecies model is stable so there is no need to average the abundances
# Abundance of the background resource with no fishing in the final time step (the first time step in the initial abundance)
ms0_npp <- ms_sim0@n_pp[ms_time_sim+1,] * ms_sim0@params@w_full
# Abundance of the community with no fishing in the final time step
ms0_n <- sweep(ms_sim0@n[ms_time_sim+1,,],2,ms_sim0@params@w,"*")
# Abundance of the community with fishing in the final time step
ms1_n <- sweep(ms_sim1@n[ms_time_sim+1,,],2, ms_sim1@params@w,"*")
# Sum over the community to get the total abundances
ms0_n_total <- apply(ms0_n,2,sum)
ms1_n_total <- apply(ms1_n,2,sum)
# Calculate the relative total abundances
ms_relative_abundance <-  ms1_n_total / ms0_n_total

#save.image(file="fig1_workspace.RData")

#----------------------------------------------------------
# Make FIGURE 1 for the paper
#----------------------------------------------------------

# Plot biomass relative to resource carrying capacity in the smallest size of community spectrum
# The resource carrying capacity is found in the @params@cc_pp slot
# It's in the same position for each model as each model has the same size range for the background spectrum
cc_index <- which(comm_sim0@params@w_full == comm_sim0@params@w[1])
comm_base <- (comm_sim0@params@w_full * comm_sim0@params@cc_pp)[cc_index]
trait_base <- (trait_sim0@params@w_full * trait_sim0@params@cc_pp)[cc_index]
ms_base <- (ms_sim0@params@w_full * ms_sim0@params@cc_pp)[cc_index]

# Figure dimensions for printing
width <- 14
height <- 7
#postscript(file = "Figure1.eps", width = width/2.5, height = height/2.5, pointsize=8, horizontal = FALSE, onefile=FALSE, paper='special')
png(filename="Figure1.png", width = width, height = height, units="cm",res=800, pointsize=8)

# Figure layput
nf <- layout(matrix(1:6,2,3,byrow=TRUE), rep(width/3,3),c(3/6,3/6)*height,TRUE)

ylim <- c(1e-10,10)
xlim <- c(1e-3,1e4)
cascade_ylim <- c(0.1,10)
fat_lwd <- 2
fished_lty <- 3
resource_colour <- "green"
resource_lwd <- 1

# Community
par(mar=c(1,5,2,1))
plot(x=comm_sim0@params@w, y= comm0_n / comm_base, log="xy", type="n", ylab="Biomass relative to carrying capacity", xlim=xlim, ylim=ylim, main = "(a)", xlab="")
# Resource
lines(x=comm_sim0@params@w_full, y= comm0_npp / comm_base, col=resource_colour, lwd=resource_lwd)
# Unfished community spectrum relative to unfished biomass in size 1
lines(x=comm_sim0@params@w, y=comm0_n / comm_base, col="red", lwd=fat_lwd)
# Fished community spectrum relative to unfished biomass in size 1
lines(x=comm_sim0@params@w, y=comm1_n / comm_base, col="blue", lwd=fat_lwd, lty=fished_lty)

# Trait
plot(x=trait_sim0@params@w, y= trait0_n_total / trait_base, log="xy", type="n",  ylab="", xlim=xlim, ylim=ylim, main = "(b)", xlab="")
# Resource
lines(x=trait_sim0@params@w_full, y= trait0_npp / trait_base, col=resource_colour, lwd=resource_lwd)
# Unfished community spectrum relative to unfished biomass in size 1
lines(x=trait_sim0@params@w, y=trait0_n_total / trait_base, col="red", lwd=fat_lwd)
# Fished community spectrum relative to unfished biomass in size 1
lines(x=trait_sim1@params@w, y=trait1_n_total / trait_base, col="blue", lwd=fat_lwd, lty=fished_lty)
# Unfished species relative to unfished biomass in size1
for (i in 1:10){
    lines(x=trait_sim0@params@w, y=trait0_n[i,] / trait_base)
}

# Multispecies
# Biomass
plot(x=ms_sim0@params@w, y= ms0_n_total / ms_base, log="xy", type="n",  ylab="", xlim=xlim, ylim=ylim, main = "(c)", xlab="")
# Resource
lines(x=ms_sim0@params@w_full, y=ms0_npp / ms_base, col=resource_colour, lwd=resource_lwd)
# Unfished community spectrum relative to unfished biomass in size 1
lines(x=ms_sim0@params@w, y=ms0_n_total / ms_base, col="red", lwd=fat_lwd)
# Fished community spectrum relative to unfished biomass in size 1
lines(x=ms_sim1@params@w, y=ms1_n_total / ms_base, col="blue", lwd=fat_lwd, lty=3)
# Unfished species relative to unfished biomass in size1
for (i in 1:12){
    lines(x=ms_sim0@params@w, y=ms0_n[i,] / ms_base)
}

# Cascades
# Community
par(mar=c(5,5,5,1))

plot(x=comm_sim0@params@w, y=comm_relative_abundance, log="xy", type="n", ylab="Relative abundance", xlim=xlim, ylim=cascade_ylim, main = "", xlab="Body mass (g)", yaxt="n")
axis(side = 2, at=c(0.2,1,5))
lines(x=comm_sim0@params@w, y=comm_relative_abundance)
lines(x=c(min(comm_sim0@params@w),max(comm_sim0@params@w)), y=c(1,1),lty=2)
lines(x=c(1000,1000),y=c(1e-20,1e20),lty=3)

# Trait
plot(x=trait_sim0@params@w, y=trait_relative_abundance, log="xy", type="n", ylab="", xlim=xlim, ylim=cascade_ylim, main = "", xlab="Body mass (g)", yaxt="n")
axis(side = 2, at=c(0.2,1,5))
lines(x=trait_sim0@params@w, y=trait_relative_abundance)
lines(x=c(min(trait_sim0@params@w),max(trait_sim0@params@w)), y=c(1,1),lty=2)
lines(x=c(1000,1000),y=c(1e-20,1e20),lty=3)

# Multispecies
plot(x=ms_sim0@params@w, y=ms_relative_abundance, log="xy", type="n", xlab = "Body mass (g)", ylab="", xlim=xlim, ylim=cascade_ylim, main = "", yaxt="n")
axis(side = 2, at=c(0.2,1,5))
lines(x=ms_sim0@params@w, y=ms_relative_abundance)
lines(x=c(min(ms_sim0@params@w),max(ms_sim0@params@w)), y=c(1,1),lty=2)
lines(x=c(1000,1000),y=c(1e-20,1e20),lty=3)

dev.off()


#----------------------------------------------------------------------
# Example 2: Changing effort time with the multispecies North Sea model
#----------------------------------------------------------------------

# Here we illustrate the performance of the multispecies model with fishing effort changing over time.
# We use the multispecies North Sea model with multiple gears

# Start afresh
rm(list=ls()) 
library(mizer)
data(NS_species_params)
data(inter)

# In this example we will set up 4 gears that catch different groups of species
# To do this we add an extra 'gears' column to the data set and specify the name of the fishing gear for each species.
# Our species are:
NS_species_params$species
# We have 4 gears: Industrial, Pelagic, Otter and Beam
# The species will be caught by:
NS_species_params$gear <- c("Industrial", "Industrial", "Pelagic", "Pelagic", "Beam", "Otter", "Beam", "Otter", "Beam", "Otter", "Otter", "Otter")
# We set the selectivity parameters for each gear.
# We will still use a knife edge selectivity pattern but but the position of the 'edge' size depends on the gear.
# Size is specified as mass. To set the size as mass we assume a Length-Weight relationship 
# w = a * l ^ b with a = 0.01 and b = 3
# Add the knife_edge_size column
NS_species_params$knife_edge_size <- NA
# Set the selection size for each gear (alternatively, the selection size could be set by species rather than being the same within a gear)
NS_species_params[NS_species_params$gear == "Industrial", "knife_edge_size"] <- 10 # length = 10
NS_species_params[NS_species_params$gear == "Pelagic", "knife_edge_size"] <- 80 # length = 20
NS_species_params[NS_species_params$gear == "Otter", "knife_edge_size"] <- 270 # length = 30 
NS_species_params[NS_species_params$gear == "Beam", "knife_edge_size"] <- 155 # length = 25 
# Check the gear columns has been set correctly
NS_species_params
# Now make the parameter object using default values for everything else
NS_params <- MizerParams(NS_species_params, inter)
# Out of interest we can see that the fishing catchability has been correctly set by gear and species (i.e. what is catching who)
NS_params@catchability

# Our initial abundances in the simulation will be the equilibrium unfished abundances
time_to_equib <- 100
# Project to equilibrium
NS_equib <- project(NS_params, t_max=time_to_equib)
plot(NS_equib)
# Pull out the equilibrium population abundances - we will use them as the initial abundances for the projection
n_equib <- NS_equib@n[time_to_equib+1,,]
n_pp_equib <- NS_equib@n_pp[time_to_equib+1,]

# Now we set up an arrat of fishing effort through time.
# The effort of each gear will be different
# After 10 years a pelagic fishery starts (increases from 0 to 1 over 10 yrs)
# After 30 years a beam fishery starts (increases from 0 to 0.75 over 10 yrs)
# After 50 years an otter fishery starts (increases from 0 to 0.9 over 10 yrs)
# After 70 years an industrial fishery starts (increases from 0 to 1.5 over 10 yrs)

# Pelagic, Beam, Otter, Industrial
gear_names <- c("Pelagic","Beam","Otter","Industrial")
project_time <- 100
fishing_effort <- array(0, dim=c(project_time, 4), dimnames=list(time=1:project_time, gear=gear_names))
fishing_effort[,"Pelagic"] <- c(rep(0,10),seq(from = 0, to = 1, length = 10), rep(1,80))
fishing_effort[,"Beam"] <- c(rep(0,30),seq(from = 0, to = 0.75, length = 10), rep(0.75,60))
fishing_effort[,"Otter"] <- c(rep(0,50),seq(from = 0, to = 0.9, length = 10), rep(0.9,40))
fishing_effort[,"Industrial"] <- c(rep(0,70),seq(from = 0, to = 1.5, length = 10), rep(1.5,20))

# Have a quick look at what the fishing efforts will do through time
plot(x = 1:project_time, y = seq(from=0,to=1,length=project_time), type="n", xlab = "Years", ylab="Fishing effort", ylim=c(0,1.5))
for (i in 1:4){
    lines(x=1:project_time, y = fishing_effort[,i], lty=i)
}
legend(x="bottomright",legend=c("Pelagic", "Beam", "Otter", "Industrial"), lty=1:4)

# Run the simulation
NS_sim <- project(NS_params, effort=fishing_effort, initial_n = n_equib, initial_n_pp = n_pp_equib)
# Have a look
plot(NS_sim)
plotBiomass(NS_sim)

# Calculate the things we're going to plot
ssb <- getSSB(NS_sim)
# Rescale ssb to be relative to the unfished ssb
rescale_ssb <- sweep(ssb,2,ssb[1,],"/")
yield <- getYieldGear(NS_sim)
# Rescale yield relative to the maximum yield over the time series
max_yield <- apply(yield,c(2,3),max)
rescale_yield <- sweep(yield,c(2,3), max_yield, "/")
# Do some tidying up of yield as not all gears catch all species
# To help with the plot set yields of non-caught species to NA
yield[yield==0] <- NA
rescale_yield[rescale_yield==0] <- NA

# Community indicators
# Size range for indicators
min_w <- 10
max_w <- 5000
threshold_w <- 100
lfi <- getProportionOfLargeFish(NS_sim, min_w=min_w, max_w=max_w, threshold_w=threshold_w)
mw <- getMeanWeight(NS_sim, min_w=min_w, max_w=max_w)
mmw <- getMeanMaxWeight(NS_sim, min_w=min_w, max_w=max_w, measure="biomass")
slope <- getCommunitySlope(NS_sim, min_w=min_w, max_w=max_w)[,"slope"]
# Scale these relative to the unfished community
rescale_lfi <- lfi / lfi[1]
rescale_mw <- mw / mw[1]
rescale_mmw <- mmw / mmw[1]

#save.image(file="fig2_workspace.RData")

#----------------------------------------------------------
# Make FIGURE 2 for the paper
#----------------------------------------------------------

# Set colours and lty for the species
# Pelagic Beam Otter Industrial gears catching 12 species
# Each gear has an lty
# Each species has same lty as the gear that catches it
# Each species within a gear has a different colour
gear_lty <- 1:4
names(gear_lty) <- gear_names
species_names <- as.character(NS_params@species_params$species)
species_lty <- rep(NA,12)
names(species_lty) <- species_names
cols <- c("black","blue","magenta","green","red")
species_col <- rep(NA,12)
names(species_col) <- species_names
for (i in gear_names){
    gear_idx <- (NS_params@species_params$gear == i)
    species_col[NS_params@species_params$species[gear_idx]] <- cols[1:(sum(gear_idx))]
    #species_lty[NS_params@species_params$species[gear_idx]] <- gear_lty[i]
    species_lty[NS_params@species_params$species[gear_idx]] <- 1
}

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

#postscript(file = "Figure2.eps", width = width/2.5, height = height/2.5, pointsize=8, horizontal = FALSE, onefile=FALSE, paper='special')
#png(filename="Figure2.png", width = width, height = height, units="cm",res=800, pointsize=8)
png(filename="Figure2.png", width = width, height = height, units="cm",res=1000, pointsize=8)

# Figuring out figure panel heights - bit fiddly
rel_heights <- c(0.7,rep(0.5,8),0.5,0.8,1) # rel heights of panels
heights = (height / sum(rel_heights)) * rel_heights
nf <- layout(matrix(1:length(rel_heights),length(rel_heights),1,byrow=TRUE), widths = width, heights=heights,TRUE)
right_margin <- 4
left_margin <- 4.5
# Other plotting parameters
legend_txt_cex <- 0.7
leg_line <- 0.3
seg_len <- 2.5
leg_bty <- "n"
leg_box_lwd <- 0

# (a) Effort of gears
par(mar=c(0,left_margin,0.5,right_margin))
plot(x = 1:project_time, y=1:project_time, type="n", ylim=c(0,max(fishing_effort)), xlab="", ylab=expression(Effort~(y^{-1})), xaxt="n")
text(x=5,y=1.4,labels="(a)")

for (i in gear_names){
    lines(x = 1:project_time, y=NS_sim@effort[,i], lty=gear_lty[i])
}
add_effort_lines()
legend(x="bottomright", legend = gear_names, lty=gear_lty, cex=legend_txt_cex, pt.lwd=leg_line, seg.len=seg_len, bty=leg_bty, box.lwd=leg_box_lwd)

# (b) yield - each gear has a separate panel
rescale_yield_min <- min(rescale_yield[rescale_yield>0], na.rm=TRUE)
rescale_yield_max <- max(rescale_yield[rescale_yield>0], na.rm=TRUE)
main_labels <- c("(b)","(c)","(d)","(e)")
names(main_labels) <- gear_names
for (gear in gear_names){
    species_in_gear <- NS_params@species_params$species[NS_params@species_params$gear==gear]
    par(mar=c(0,left_margin,0,right_margin))
    plot(x = 1:project_time, y=1:project_time, type="n", ylim=c(rescale_yield_min, rescale_yield_max), ylab="", xlab="", xaxt="n", yaxt="n")
    axis(4)
    text(x=5,y=0.9,labels=main_labels[gear])
    mtext(gear, side=2, line=1, cex=0.6)
    if (gear == gear_names[2]){
        mtext("Relative Yield", side=4, line=3, cex=0.6, adj=-3)
    }
    add_effort_lines()
    for (i in species_in_gear){
        lines(x = 1:project_time, y=rescale_yield[1:project_time,gear,i], col=species_col[i], lty=species_lty[i])
    }
    legend(x="bottomright", legend=species_in_gear, lty=species_lty[species_in_gear], col=species_col[species_in_gear], cex = legend_txt_cex, ncol=1, pt.lwd=leg_line, seg.len = seg_len, bty=leg_bty, box.lwd=leg_box_lwd)
}

# (c) SSB - each gear has a separate panel
main_labels <- c("(f)","(g)","(h)","(i)")
names(main_labels) <- gear_names
for (gear in gear_names){
    species_in_gear <- NS_params@species_params$species[NS_params@species_params$gear==gear]
    par(mar=c(0,left_margin,0,right_margin))
    plot(x = 1:project_time, y=1:project_time, type="n", ylim=c(min(rescale_ssb),max(rescale_ssb)), ylab="", xlab="", xaxt="n")
    if (gear == gear_names[2]){
        mtext("Relative SSB", side=2, line=3, cex=0.6, adj=-3)
    }
    text(x=5,y=1.75,labels=main_labels[gear])
    mtext(gear, side=4, line=1, cex=0.6)
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
text(x=5,y=-1.5,labels="(j)")
mtext("Community slope", side=4, line=3, cex=0.6)
add_effort_lines()
lines(x = 1:project_time, y = slope[2:(project_time+1)])

# (e) LFI, MW, MMW,
par(mar=c(4,left_margin,0,right_margin))
ylim <- c(0,max(rescale_lfi,rescale_mw, rescale_mmw, slope))
plot(x = 1:project_time, y=1:project_time, type="n", ylab="Relative metrics", xlab="Years", ylim=ylim)
text(x=5,y=2.5,labels="(k)")
lines(x = 1:project_time, y = rescale_lfi[2:(project_time+1)], col=1)
lines(x = 1:project_time, y = rescale_mw[2:(project_time+1)], col=2)
lines(x = 1:project_time, y = rescale_mmw[2:(project_time+1)], col=3)
add_effort_lines()
legend(x="bottomright", legend = c("LFI", "MW", "MMW"), lty=1, col=c(1,2,3), cex = legend_txt_cex, ncol=1, pt.lwd=leg_line, seg.len = seg_len, bty=leg_bty, box.lwd=leg_box_lwd)

# (f) Spectra at various points in time
# Size distribution at start and end - five lines
# times of something happening
# 20, 40, 60, 80
biomass_time <- sweep(apply(NS_sim@n[c(2,21,41,61,81),,],c(1,3),sum),2,NS_sim@params@w,"*")
xlim <- c(1,5e4)
ylim <- c(5e5,max(biomass_time))/1000
par(mar=c(4,left_margin,1,right_margin))
plot(x=NS_sim@params@w, y = NS_sim@params@w, type="n", ylab="Total biomass (kg)", xlab = "Size (g)", log="xy", ylim = ylim, xlim=xlim)
text(x=2,y=1e9,labels="(l)")
cols <- c(1,2,3,4,6)
for (i in 1:5){
    lines(x=NS_sim@params@w, y = biomass_time[i,]/1000, col=cols[i])
}
legend(x="bottomleft", legend=c("Unfished", "Year 20: Pelagic", "Year 40: Beam", "Year 60: Otter", "Year 80: Industrial"), lty=1, col=cols, cex = legend_txt_cex, pt.lwd=leg_line, seg.len = seg_len, bty=leg_bty, box.lwd=leg_box_lwd)
dev.off()
# Done
#-----------------------------------------------------------------------------
