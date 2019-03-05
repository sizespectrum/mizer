data(NS_species_params_gears)
data(inter)

# Get a feel for size of rounding errors ####
params <- MizerParams(NS_species_params_gears, inter)
sim <- project(params)

params <- MizerParams(NS_species_params_gears, inter, no_w = 1000)
sim2 <- project(params)

tmax <- 101
round((getBiomass(sim)[tmax,]/getBiomass(sim2)[tmax,]),3)

# Compare results for available energy between old and fft code ####
params <- MizerParams(NS_species_params_gears, inter, no_w = 30)
params2 <- params
params2@ft_pred_kernel_e <- array()
afft <- getAvailEnergy(params, params@initial_n, params@initial_n_pp)
a <- getAvailEnergy(params2, params@initial_n, params@initial_n_pp)
all.equal(as.vector(afft),as.vector(a), tolerance = 0)
# The discrepancy is due to the kernels that extend to below the smallest plankton
i <- 4
plot(params@pred_kernel[i ,1, ])
all.equal(as.vector(afft[i, ]), as.vector(a[i, ]), tolerance = 0)
# This is fixed by extending the plankton spectrum
params <- MizerParams(NS_species_params_gears, inter, no_w = 30,
                      min_w_pp = 1e-13)  # small min_w_pp required
params2 <- params
params2@ft_pred_kernel_e <- array()
afft <- getAvailEnergy(params, params@initial_n, params@initial_n_pp)
a <- getAvailEnergy(params2, params@initial_n, params@initial_n_pp)
all.equal(as.vector(afft),as.vector(a), tolerance = 0)

# Compare results for predation rate between old and fft code ####
params2@ft_pred_kernel_p <- array()
f <- getFeedingLevel(params, params@initial_n, params@initial_n_pp)
pfft <- getPredRate(params, params@initial_n, params@initial_n_pp,
                    feeding_level = f)
p <- getPredRate(params2, params@initial_n, params@initial_n_pp,
                 feeding_level = f)
all.equal(as.vector(pfft),as.vector(p), tolerance = 0)
# The biggest discrepancy comes from species with a feeding kernel that
# extends to close to the predator's size
i <- 3
plot(params@pred_kernel[i, 1, ])
all.equal(as.vector(pfft[i, ]), as.vector(p[i, ]), tolerance = 0)
plot((pfft[i, ] - p[i, ])/p[i, ])
# We see that in the fft method there is a nonzero predation rate from this 
# species at a size larger than the largest individual
w_inf <- params@species_params$w_inf[i]
pfft[i, params@w_full > w_inf]
p[i, params@w_full > w_inf]


plot(params@initial_n[i, ], log = "y")
