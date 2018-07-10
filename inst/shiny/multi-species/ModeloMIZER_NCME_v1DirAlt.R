# MULTISPECIES SIZE SPECTRUM MODEL / Setting up a model for the North Chilean Marine Ecosystem in MIZER
# ----------------------------------------------------------------------------------------------------------------------------------------------------
# The script is developed following the vignette MIZER v0.2 and MIZER v1.0 (Scott, Blanchard and Andersen, 2014) and improvement by G. Delius,
# R. Southwell and R. Law (April 2018).
# -----------------------------------------------------------------------------------------------------------------------------------------------------

# install.packages("mizer")
# help(package='mizer') # This command gives a technical summary of the package, including the available functions
# help(mizer)           # The second command gives a brief introduction to mizer.
# help(project)         # The third gives the documentation page for the method project
# help(plot, package="mizer") # Example To select the help page for the appropriate plotting method
# method ? getFeedingLevel  or class ? MizerParams  # TOther ways of accessing package documentation

library(mizer)
packageVersion('mizer')
# vignette("mizer_vignette") # Here we just called the vignette of MIZER package
#setwd("/Users/Mlla/WORK_MLLA/POSTDOC_PUC/Objetivo3/NCME_MIZER")
rm(list=ls())

#----------------------------------------------
# 1) Load the parameter input data
#----------------------------------------------
# Notes from vignettes for MIZER v0.2 and MIZER v1.0
# Species-specific parameters should be store in a single data.frame species (row) x paras (col). Best used a *.csv file.
# MizerParams class is used for storing model parameters.
# The MizerParams class  stores: life-history parameters (Winf,Wmat),
# Size-based biological parameters: search volume, V(w)
# Stock recruitmeny functions, Growth and dynamics of the resource spectrum
# Fishing gear pars: selectivity and catchability.

# --- Species Parameter file --- 
# I load only the minimun of parameters suggested.
# Input parameters as in the vignette

params_data<-read.csv("speciesNCME_edited2.csv")[,c(2,7,8,13,14,10,29)]

# --- Species Interaction Matrix ---
# The interaction matrix describes the interaction of each pair of species in the model.
# The values are between 0 (species do not overlap and therefore do not interact with each other) 
# to 1 (species overlap perfectly). 
# If all the values in the interaction matrix are set to 1 then predator-prey interactions are determined entirely by size-preference.
# The interaction matrix must be of type array or matrix. Better to read it from *.csv file.
# It should be noted that the order of species in the interaction matrix has to be the same as the order in the species parameters data.frame.
# This is a matrix of 1s, implying that all species fully interact with each other, i.e. the species are spread
# homogeneously across the model area
# it has to be a class matrix object

inter=as(read.csv('theta.csv',row.names=1),"matrix")

# --- Other parameters that are default ---
# Pars that have default values are: k=0 (ctivilty coefficient), alpha=0.6 (assimitaion efficecny),
# erepro=1 (reproductive efficiency, ~ steepness), w_min= recruitment size (default smallest size of the com SS),
# sel_func=knife_edge, with wmin (selectivity function), gear, default names is the name of the species, catchability=1

# Pars that are essential but if they are not provided (see notes in my notebook), values for
# them are estimated using the values in the other columns. These are: h (maximun intake), gamma (volumetric search rate),
# ks coeficient for standard metabolism, z0 mortality ub,i.  

#------------------------------------------------------------------------------------------
# 2) Setting up model parameters -- MizerParams () (does not store pars that change with time)
#------------------------------------------------------------------------------------------
# Determining a value for the kappa argument can be dificult and may need to be estimated through
# some kind calibration process. The default value kappa is for the North Sea model (kappa=1e11).

# ---------------------------
# Setting the carrying capacity, and creating the MIZER objects params to run simulation
# Kappa would be estimated along with the values r_max, calibration process.
# ---------------------------
species=read.csv("speciesNCME_edited2.csv")

wfmin <- min(species$Wegg)                                # min mass (g) of a fish (min egg size)
wfmax <- max(species$w_inf)                               # max mass (g) of a fish (largest Winf in species community; Julia has 2*max(species$Winf))
xpmin  <- -23					                                    # smallest plankton prey (10^-10 g) Chilean data
xpmax  <- log(0.1)					                              # largest plankton prey	Chilean data
nx     <- 295


params <- MizerParams(params_data, interaction=inter, kappa = 3e10, min_w=wfmin,max_w=wfmax, min_w_pp=exp(xpmin), 
                      w_pp_cutoff=exp(xpmax))

# --- Adding the differents gears ---
# Here the fisheries are identified regarding each species
# NF = No Fishery
# P  = Pelagic
# Each selectivity function has a range of arguments. Values for these arguments mustbe included as columns in the species parameters data.frame

params_data_gears <- params_data
params_data_gears$gear <-read.csv("speciesNCME_edited2.csv")[,18]
params_gears <- MizerParams(params_data_gears, interaction = inter, kappa = 3e10)

# To explore what is in params
params_gears@species_params  # We can see that this object contains the original species data.frame (with w inf and so on), plus any default
                       # values that may not have been included in the original data.frame.
summary(params_gears)

#----------------------------------------------
# 3) Running the simulations  --- project()
#----------------------------------------------
# Project function is used to run simulations. This method takes a MizerParams object
# a project forward through time, starting from an initial populations.
# and with pre-determined fishing effort pattern.
# It requires MizerParams object, Fishing effort, Initial Population and time argument (sim time step, length of the sim,
# how often record the output)

# arguments of project()
# ---------------------------
# object A MizerParams object
# effort The effort of each fishing gear through time: a) Project without an effort argument, b) Project with a constant effort, c) Project with time varying effor
# t_max The maximum time the projection runs for. The default value is 100. However, this argument is not needed if an array is used for the effort argument,
#        in which case this argument is ignored
# dt Time step of the solver. The default value is 0.1.
# t_save The frequency with which the output is stored. The default value is 1. Must be an integer multiple of dt.
# initial_n  The initial populations of the species
# initial_n_pp  The initial population of the background spectrum. It should be a numeric vector of the same length as the w_full slot of the MizerParams argument.
# shiny_progress A shiny progress object used to update shiny progress bar
# Fishing effort
# Initial Population
# Time arguments
# help('project')
# time arguments: dt (time step, default is 0.1), t_max (default=100).
# t_save: this sets how frecuently project() stores the state of the model in the resulting MizerSim object. Defaul value is 1

# ---------------------------
# Setting the Initial abundance (E=0)
# ---------------------------

#sim<- project(params, effort = 0, t_max = 10, dt = 0.05, t_save = 0.1)
