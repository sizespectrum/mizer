data(NS_species_params_gears)
data(inter)
params <- MizerParams(NS_species_params_gears, inter)
no_gear <- dim(params@catchability)[1]
no_sp <- dim(params@catchability)[2]
no_w <- length(params@w)
no_w_full <- length(params@w_full)
sim <- project(params, effort=1, t_max=20, dt = 0.5, t_save = 0.5)
sim2 <- project(params, effort=2, t_max=20, dt = 0.5, t_save = 0.5)



# 
# setMethod('getSSBFrame', signature(sim='MizerSim')
#           
#           setMethod('plotBiomass', signature(sim='MizerSim')
plotBiomass(sim)
#                     setMethod('plotYield', signature(sim='MizerSim', sim2='missing'),
plotYield(sim)
#                               
#                               setMethod('plotYield', signature(sim='MizerSim', sim2='MizerSim'),
plotYield(sim,sim2)

#                                         setMethod('plotYieldGear', signature(sim='MizerSim'),
plotYieldGear(sim)
#error here
#                                                   setMethod('plotSpectra', signature(object='MizerSim'),
plotSpectra(sim)
#                                                             setMethod('plotSpectra', signature(object='MizerParams'),
plotSpectra(sim@params)
#error here
#                                                                       setMethod('plotFeedingLevel', signature(sim='MizerSim'),
plotFeedingLevel(sim)
#                                                                                 setMethod('plotM2', signature(sim='MizerSim'),
plotM2(sim)
#                                                                                           setMethod('plotFMort', signature(sim='MizerSim'),
plotFMort(sim)
#                                                                                                     plot
plot(sim)
#                                                                                                     setMethod('plotGrowthCurves', signature(object = 'MizerSim'),
plotGrowthCurves(sim)
#                                                                                                               setMethod('plotGrowthCurves', signature(object = 'MizerParams'),
#                                                                                                                         
#      
plotGrowthCurves(sim@params)


# what to do about plotYieldGear ?
# properly implement plotSpectra
# test getSSBFrame
# reindent and clean up plots.R
