# Load Mariella's species parameters
library(readr)
params_data <- read_csv("speciesNCME_Mariella.csv",
                        col_types = cols(
                            X1 = col_integer(),
                            species = col_character(),
                            Linf = col_double(),
                            Lmat = col_double(),
                            a2 = col_double(),
                            b2 = col_double(),
                            w_inf = col_double(),
                            w_mat = col_double(),
                            Wegg = col_double(),
                            k_vb = col_double(),
                            log10_beta = col_double(),
                            log10_sigma = col_double(),
                            beta = col_double(),
                            sigma = col_double(),
                            a1_maturity = col_double(),
                            b1_maturity = col_double(),
                            `Ni0 (egg/m3)` = col_double(),
                            `Fishery Gear` = col_character(),
                            `w_start_Fishery (cm)` = col_double(),
                            `w_L50%` = col_double(),
                            Fmort = col_double(),
                            obsCatch_tonnes = col_double(),
                            TotBiom_tonnes = col_double(),
                            w_TB_start = col_double(),
                            SSB_tonnes = col_double(),
                            w_SSB_start_g = col_double(),
                            eRepro = col_double(),
                            N0 = col_double(),
                            to = col_double(),
                            us_prey = col_double(),
                            ls_prey = col_double(),
                            to_corr = col_double(),
                            r_max_guess = col_double()
                        ))

# Set some general parameters
lambda <- 2.12  # Exponent of community power law
effort <- 1

# What Mariella calls sigma in this file is the exponential
# of the standard deviation in the log of the observed predator/prey ratio
# observed in stomach data. We convert that to mizer's sigma:
params_data$sigma <- log(params_data$sigma)

# What Mariella calls beta in this file is the exponential
# of the mean of the log of the observed predator/prey ratio 
# observed in stomach data. We convert that to mizer's beta:
params_data$beta <- params_data$beta *
    exp(-params_data$sigma^2 * (lambda - 5/3))


# Fishing selectivity parameters are missing from data frame
params_data$l25 <- c(1.0e+29,     1.9e+00,     4.0e+00,     5.0e+00,     2.9e+00,     3.2e+01,     4.9e+00,     4.9e+01 )
params_data$l50 <-  c(1.1e+29,     2.0e+00,     5.0e+00,     6.0e+00,     3.0e+00,     3.6e+01,     8.0e+00,     5.1e+01) 

# Remove the mesopelagics from the dataframe because we model them via the
# background species
# params_data <- subset(params_data, params_data$species != "Mesopelagic")

# Spread the background species over the entire range up to the largest fish
min_egg <- min(params_data$Wegg)
max_w_inf <- max(params_data$w_inf)

# We need to give even the smallest individuals a full range of planktonic prey
min_w_pp = min(params_data$Wegg /
                   (params_data$beta * exp(3 * params_data$sigma)))

p <- setBackground(
    set_scaling_model(min_w_pp = 1e-9,
        no_sp = 10, no_w = 400, min_w_inf = 2, max_w_inf = max_w_inf,
        min_egg = min_egg, min_w_mat = 2 / 10^0.8, 
        lambda = lambda, beta = 500, sigma = 2,
        knife_edge_size = Inf)
)
plotSpectra(p)


gears <- "knife_edge"
no_sp <- dim(params_data)[1]
for (i in (1:no_sp)) {
    if (params_data$species[i] == "Anchovy") {
        gear <- "sigmoid_gear_Anchovy"
    } else {
        gear <- "sigmoid_gear"
    }
    gears <- union(gears, gear)
    species_params <- data.frame(
        species = params_data$species[i],
        w_min = params_data$Wegg[i],
        w_inf = params_data$w_inf[i],
        w_mat = params_data$w_mat[i],
        beta = params_data$beta[i],
        sigma = params_data$sigma[i],
        z0 = 0,
        alpha = 0.6,
        erepro = 0.1, # unknown, determined later
        sel_func = "sigmoid_length",
        gear = gear,
        l25 = params_data$l25[i],
        l50 = params_data$l50[i],
        k = 0,
        k_vb = params_data$k_vb[i],
        t0 = params_data$to[i],
        a = params_data$a2[i],
        b = params_data$b2[i],
        catchability = params_data$Fmort[i],
        catch_observed = params_data$obsCatch_tonnes[i] * 1e-6,
        biomass_observed = params_data$TotBiom_tonnes[i] * 1e-6,
        biomass_cutoff = params_data$w_TB_start[i]
    )

    p <- addSpecies(p, species_params, effort = effort, rfac = Inf)
}
plotSpectra(p)

# Run to steady state
p <- steady(p, effort = effort, t_max = 500,  tol = 1e-3)
plotSpectra(p)

# Remove background species
p <- removeSpecies(p, is.na(p@A))
p <- steady(p, effort = effort, t_max = 500,  tol = 1e-3)

# rescale kappa so that Anchovy has desired abundance
sp <- 2
biomass <- cumsum(p@initial_n[sp, ] * p@w * p@dw)
cutoff_idx <- which.min(p@w <= p@species_params$biomass_cutoff[sp])
scale <- p@species_params$biomass_observed[sp] / 
    (max(biomass) - biomass[cutoff_idx])
p@initial_n <- p@initial_n * scale
p@initial_n_pp <- p@initial_n_pp * scale
p@cc_pp <- p@cc_pp * scale
# To keep the same per-capity behaviour, we have to scale down the
# search volume
p@species_params$gamma <- p@species_params$gamma / scale
p@search_vol <- p@search_vol / scale
p@kappa <- p@kappa * scale

# Check that this worked
sim <- project(p, t_max = 15, t_save = 0.1, effort = effort)
plotBiomass(sim)

# Save params object
saveRDS(p, file = "humboldt_params.rds")

###############

p <- readRDS("humboldt_params.rds")
class(p)

effort <- 1.4
sim <- project(p, t_max = 15, t_save = 0.1, effort = effort)
plot(sim)

p@species_params$erepro
plotGrowthCurves(p, species = "Sardine")

# Now change one of the parameters
sp <- "Sardine"

species_params <- p@species_params

species_params[sp, "gamma"] <- 
    species_params[sp, "gamma"] * 2
# species_params["JMackerel", "w_mat"] <-
#     species_params["JMackerel", "w_mat"] * 0.9
# species_params["JMackerel", "w_min"] <-
#     species_params["JMackerel", "w_min"] * 10

# effort <- c(knife_edge_gear = 0, sigmoid_gear = 1.4, sigmoid_gear_Anchovy = 1.1)

pc <- MizerParams(
    species_params,
    p = p@p,
    n = p@n,
    q = p@q,
    lambda = p@lambda,
    f0 = p@f0,
    kappa = p@kappa,
    min_w = min(p@w),
    max_w = max(p@w),
    no_w = length(p@w),
    min_w_pp = min(p@w_full),
    w_pp_cutoff = max(p@w_full),
    r_pp = (p@rr_pp / (p@w_full ^ (p@p - 1)))[1]
)
pc@linetype <- p@linetype
pc@linecolour <- p@linecolour
pc@A <- p@A
pc@sc <- p@sc
pc@cc_pp <- p@cc_pp
pc@mu_b <- p@mu_b

pc@initial_n <- p@initial_n
pc@initial_n_pp <- p@initial_n_pp

# The spectrum for the changed species is calculated with new
# parameters but in the context of the original community
# Compute death rate for changed species
mumu <- getMort(pc, p@initial_n, p@initial_n_pp, effort = effort)[sp, ]
# compute growth rate for changed species
gg <- getEGrowth(pc, p@initial_n, p@initial_n_pp)[sp, ]
# Compute solution for changed species
w_inf_idx <- sum(pc@w < pc@species_params[sp, "w_inf"])
idx <- pc@w_min_idx[sp]:(w_inf_idx - 1)
if (any(gg[idx] == 0)) {
    stop("Can not compute steady state due to zero growth rates")
}
pc@initial_n[sp, ] <- 0
pc@initial_n[sp, pc@w_min_idx[sp]:w_inf_idx] <- 
    c(1, cumprod(gg[idx] / ((gg + mumu * pc@dw)[idx + 1])))
if (any(is.infinite(pc@initial_n))) {
    stop("Candidate steady state holds infinities")
}
if (any(is.na(pc@initial_n) || is.nan(pc@initial_n))) {
    stop("Candidate steady state holds none numeric values")
}
plotSpectra(pc)

# Run to steady state
pcs <- steady(pc, effort = effort, t_max = 20, tol = 1e-2)

pcs@species_params$erepro
plotSpectra(pcs)


# Testing changing a general parameter

p_old <- p

p <- setBackground(
    set_scaling_model(
        #min_w_pp = input$min_w_pp, no_sp = input$no_bg_sp, no_w = 400,
        # min_w_pp = 1e-12, no_sp = input$no_bg_sp, no_w = 400,
        # min_w_inf = 2, max_w_inf = 6e5,
        # min_egg = 1e-4, min_w_mat = 2 / 10^0.6,
        # lambda = input$lambda, knife_edge_size = Inf,
        # f0 = input$f0, h = input$h_bkgd, r_pp = 10^input$log_r_pp
        min_w_pp = 1e-12,
        no_sp = 10, no_w = 400, min_w_inf = 2, max_w_inf = 6e5,
        min_egg = 1e-4, min_w_mat = 2 / 10^0.6, 
        lambda = 2.12,
        knife_edge_size = Inf,
        
    )
)
# Loop over all foreground species and add them one-by-one to the new
# background
effort <- 0
names(effort) <- "knife_edge_gear"
all_efforts <- c(0, 1.4, 1.1)
names(all_efforts) <- c("knife_edge_gear", "sigmoid_gear", "sigmoid_gear_Anchovy")
no_sp <- length(p_old@A)
for (sp in (1:no_sp)[!is.na(p_old@A)]) {
    if (!(p_old@species_params[sp, "gear"] %in% names(effort))) {
        effort <- c(effort, all_efforts[p_old@species_params[sp, "gear"]])
    }
    p <- addSpecies(p, p_old@species_params[sp, ],
                    effort = effort,
                    rfac = Inf)
}

# Run to steady state
p <- steady(p, effort = effort, 
            t_max = 100, tol = 1e-2)

