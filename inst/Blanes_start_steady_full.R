species_params <- read.csv("vignettes/Blanes_species_params_full.csv", sep = ";")

#species_params$r_max[] <- Inf


library(readr)
theta <- read_csv("vignettes/Blanes_theta_full3.csv", 
                               col_names = FALSE)
theta <- t(as.matrix(theta))

## Choose B ----

B <- c(detritus = 1, carrion = 1)

## Choose rho ----
my_rho  <- 10^20
no_sp <- dim(species_params)[1]
rho <- array(0, dim = c(no_sp,2))
#rho[1,2] <- my_rho
rho[,1] <- species_params$detQ*10^20
rho[,2] <- species_params$scavQ*10^20


resource_dynamics <- list(detritus = detritus_dynamics,
                            carrion = carrion_dynamics)
resource_params <- list("detritus_external" = 0, "detritus_proportion" = 0.5, 
                          "carrion_external" = 0)

params <- multispeciesParams(species_params, rho = rho,
                             resource_dynamics = resource_dynamics,
                             resource_params = resource_params, 
                             interaction = theta, interaction_p = rep(0,no_sp))
params@initial_B <- B

params <- steady(params, t_max = 220, t_per = 220)


sim <- project(params, t_max = 220)
plotBiomass(sim)

# detritus and carrion not yet implemented.

GB <- getBiomass(sim)
dim(GB)
GB[221,]
which.min(GB[221,])
plot(sim)

# igraph
#graph_from_adjacency_matrix(adjmatrix, mode = c("directed", "undirected",
#                                                "max", "min", "upper", "lower", "plus"), weighted = NULL, diag = TRUE,
#                            add.colnames = NULL, add.rownames = NA)

library(igraph)
plot.igraph(graph_from_adjacency_matrix(theta), arrow.size = 0.01, arrow.width = 0.01)
################################################################

params@initial_n

my.mat = matrix(c(0,0,1,0),nrow = 2, ncol = 2)
my.mat
plot.igraph(graph_from_adjacency_matrix(my.mat), arrow.size = 0.01, arrow.width = 0.01)


#enc <- getEncounter(params,params@initial_n
#,params@initial_n_pp
#,B = B)

params@interaction[19,]
