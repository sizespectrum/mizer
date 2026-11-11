---
title: "Dynamic multi-species size-spectrum modelling in R: steady-state calibration, advection–diffusion numerics, stability analysis and modular extensions with mizer"
shorttitle: "Dynamic size-spectrum modelling with mizer"
author:
  - name: Gustav W. Delius
    orcid: 0000-0003-4092-8228
    affiliation: "1,*"
  - name: Finlay Scott
    orcid: 0000-0001-9950-9023
    affiliation: "2"
  - name: Julia L. Blanchard
    orcid: 0000-0003-0532-4824
    affiliation: "3,4"
  - name: Ken H. Andersen
    orcid: 0000-0002-8478-3430
    affiliation: "5"
affiliation:
  - id: "1"
    institution: "Department of Mathematics, University of York, York, UK"
  - id: "2"
    institution: "Oceanic Fisheries Programme, The Pacific Community (SPC), Nouméa, New Caledonia"
  - id: "3"
    institution: "Institute for Marine and Antarctic Studies, University of Tasmania, Hobart, TAS, Australia"
  - id: "4"
    institution: "Centre for Marine Socioecology, University of Tasmania, Hobart, TAS, Australia"
  - id: "5"
    institution: "National Institute of Aquatic Resources (DTU Aqua), Technical University of Denmark, Kongens Lyngby, Denmark"
correspondence: "*Correspondence: Gustav W. Delius (gustav.delius@york.ac.uk)"
running_headline: "Dynamic size-spectrum modelling with mizer"
bibliography: references.bib
---

## Abstract

1. Individual body size governs metabolic demand, prey selection, somatic growth, fecundity and vulnerability to both natural predators and fishing gears across aquatic ecosystems, where fish routinely grow across four to six orders of magnitude in mass from egg to adult. Dynamic multi-species size-spectrum models (MSSMs) capture these ontogenetic shifts mechanistically, yet bridging theoretical size-spectrum partial differential equations (PDEs) and empirical ecosystem management has been hindered by strong nonlinear feedbacks during calibration, artificial numerical diffusion in standard discretisations, and a lack of tools for diagnosing dynamical stability or composing ecosystem extensions.
2. Here we present the mathematical formulation, calibration methodology, numerical algorithms, stability diagnostics and modular software architecture of `mizer` (version 3), an open-source R package (with C++ computational kernels) for dynamic size-spectrum modelling. `mizer` generalises consumer dynamics to coupled McKendrick–von Foerster advection–diffusion PDEs incorporating both predation-induced jump-growth variance and external growth variability, accelerated by $O(N_w \log N_w)$ fast Fourier transform convolutions. A steady-state-first calibration workflow decouples trophic equilibrium tuning (`tuneSteadyState()`) and exact abundance-scaling symmetries (`calibrateBiomass()`) from stock–recruitment resilience (`reproduction_level()`). Flux-limited second-order spatial reconstructions (`second_order_w`) and L-stable two-stage TR-BDF2 time integration eliminate numerical cohort smearing and damping, while algebraic Newton root-finding (`solver = "newton"`) and continuous-time Jacobian eigensolvers (`getStability()`) characterise equilibria regardless of dynamical stability.
3. Applied to pulsed-cohort benchmarks and a 12-species North Sea food web (`NS_params`, $1{,}171$ coupled state variables), `mizer` reproduces observed species biomasses, ontogenetic diet shifts and growth trajectories at machine-precision steady state. Stability and bifurcation analyses (`scanModel()`) reveal a supercritical Hopf bifurcation in the North Sea food web at relative fishing effort $E^* = 1.065$, above which multi-species cohort resonances generate limit cycles whose squared relative amplitude follows canonical normal-form scaling ($R^2 = 0.999$)—oscillations that conventional first-order semi-implicit solvers either erase near equilibrium or understate four-fold in amplitude.
4. Combined with a four-tier S3 extension architecture supporting coupled benthic, detrital, thermal and bioenergetic modules, `mizer` provides a unified, reproducible methodology linking ecological theory, fisheries assessment and global climate-impact ensembles.

**Keywords:** allometric scaling, ecosystem-based fisheries management, food-web dynamics, Hopf bifurcation, marine ecosystems, McKendrick–von Foerster equation, trophic interactions

---

## 1 Introduction

Across marine and freshwater ecosystems, trophic structure is organised at least as strongly by individual body mass as by taxonomic identity [@sheldon1972; @andersen2019]. Unlike terrestrial vertebrates, most teleost fishes spawn eggs of roughly $0.5\text{–}2\text{ mg}$ regardless of whether adults mature at a few grams or several hundred kilograms [@andersen2006; @andersen2016]. As a growing fish traverses four to six orders of magnitude in body mass, it moves continuously through trophic guilds—beginning life as planktivorous prey alongside forage species, competing with benthic or pelagic mesopredators as a juvenile, and eventually preying upon both competitors and smaller conspecifics as an apex piscivore [@hartvig2011]. Anthropogenic pressures act directly on this size structure: size-selective fishing truncates adult biomass spectra, releases intermediate size classes from predation mortality, alters somatic growth rates via prey competition, and reshapes community resilience to warming and deoxygenation [@blanchard2014; @jacobsen2014; @lindmark2022].

Quantitative ecosystem modelling frameworks have historically approached this complexity from opposite poles. Static trophic mass-balance and biomass-pool models such as Ecopath with Ecosim [@christensen2004] and Rpath [@lucey2020] resolve dozens of functional groups from routine diet and catch data, while spatially explicit end-to-end platforms such as Atlantis [@audzijonyte2019] and individual-based models such as OSMOSE [@shin2004] incorporate detailed biogeochemistry, life-history stages and spatial movement. However, discretising populations into coarse age stanzas or unstructured biomass pools obscures continuous ontogenetic feedbacks between prey availability, individual growth, maturation timing and size-dependent predation mortality, while demanding hundreds to thousands of empirical parameters that are difficult to constrain outside data-rich shelf seas [@spence2016; @pennino2026]. Conversely, classical size-spectrum partial differential equations (PDEs) derive community structure parsimoniously from individual bioenergetics and allometric scaling laws [@andersen2006; @datta2010; @hartvig2011], yet early implementations remained largely restricted to idealised power-law communities or single-species theoretical analyses.

The R package `mizer` was originally introduced over a decade ago as an application note describing its initial prototype (version 0.2), which united community, trait-based and multi-species size-spectrum formulations under a forward-projection interface [@scott2014]. Since that initial release, `mizer` has been adopted across more than 40 empirical and theoretical studies worldwide [@pennino2026], spanning mixed-fishery harvest trade-offs [@blanchard2014; @dejuan2023; @wo2024], trophic cascades, early-warning recovery signals and nutrient-yield management [@szuwalski2017; @clements2019; @robinson2022], freshwater food webs [@benoit2022], benthic–pelagic coupling [@audzijonyte2023], and regional-to-global climate-change ensembles within the Fisheries and Marine Ecosystem Model Intercomparison Project (Fish-MIP) [@tittensor2018; @lotze2019; @reum2020; @ortegacisneros2025].

This transition from theoretical exploration to routine empirical calibration and ensemble projection exposed four fundamental methodological challenges that could not be addressed within the architecture of early size-spectrum software:

1. **Ill-conditioned forward calibration.** Projecting a newly parameterised multi-species PDE forward in time with active density-dependent stock–recruitment kinetics rarely converges to observed community structure: imperfect initial feeding levels distort somatic growth and spawning biomass, triggering recruitment collapses or competitive exclusions before trophic parameters can be adjusted.
2. **Individual growth dispersion and artificial numerical diffusion.** Purely advective McKendrick–von Foerster formulations treat somatic growth as deterministic, even though discrete prey captures and individual heterogeneity disperse cohorts across body mass via an advection–diffusion process [@datta2010]. At the same time, standard first-order upwind finite-difference schemes inject artificial numerical diffusivity of magnitude $d_{\text{num}} \approx g(w)\,w\ln\beta + g(w)^2\Delta t$, which overwhelms biological growth dispersion on practical grids and damps genuine cohort oscillations.
3. **Hidden dynamical instabilities and solver artefacts.** Size-structured predator–prey feedbacks readily generate oscillatory cohort resonances and supercritical or subcritical Hopf bifurcations [@datta2011; @canales2020]. Time-marching alone cannot locate unstable equilibria, and coarse first-order implicit time steps can make dynamically unstable ecosystems appear spuriously stable.
4. **Monolithic extensibility.** Incorporating ecosystem-specific mechanisms—such as detrital carrion pools, multiple size-structured resources, seasonal reproduction or physiological thermal performance—previously required forking core package code, fragmenting reproducibility across studies.

Here we present the mathematical, numerical and architectural framework of `mizer` (version 3), representing a complete methodological transformation of dynamic size-spectrum modelling in R. Section 2 formulates the generalised multi-species advection–diffusion food-web equations, multi-gear fishing selectivity and parsimonious allometric parameterisation. Section 3 introduces a steady-state-first calibration workflow that exploits exact scaling symmetries and decouples trophic equilibrium fitting from dynamical stock–recruitment resilience. Section 4 derives the conservative finite-volume discretisation, single-quadrature fast Fourier transform (FFT) convolutions, flux-limited second-order spatial reconstruction and L-stable TR-BDF2 time integrator, alongside smoothness conditions for state-dependent rate functions. Section 5 presents direct algebraic Newton root-finding, continuous-time Jacobian stability diagnostics and nonlinear bifurcation continuation, illustrated on a calibrated 12-species North Sea model. Finally, Section 6 details `mizer`'s state-bearing S3 architecture, four-tier extension system and machine-readable agent skills, and Section 7 situates `mizer` within the broader landscape of ecological modelling tools.

---

## 2 Mathematical Formulation of the Generalised Size-Spectrum Model

### 2.1 Advection–diffusion consumer dynamics and boundary conditions

`mizer` represents an aquatic community through the time-dependent size-spectrum number densities $N_i(w, t)$ ($\text{numbers}\,\text{g}^{-1}$, or area/volume-normalised equivalents $\text{m}^{-2}\,\text{g}^{-1}$ or $\text{m}^{-3}\,\text{g}^{-1}$) of $i = 1, \dots, S$ consumer species as functions of individual body mass $w$ ($\text{g}$) and time $t$ ($\text{yr}$), coupled to a size-structured resource spectrum $N_R(w, t)$ and optional additional ecosystem state variables $N_{\text{other}}(t)$. Rather than assuming purely deterministic growth trajectories, `mizer` governs each consumer spectrum via a conservative McKendrick–von Foerster advection–diffusion PDE (a forward Kolmogorov / Fokker–Planck equation with mortality sink) derived from stochastic jump-growth predator–prey dynamics [@datta2010]:

$$
\frac{\partial N_i(w, t)}{\partial t} + \frac{\partial J_i(w, t)}{\partial w} = -\mu_i(w, t)\,N_i(w, t), \qquad w_{\min,i} < w < w_{\max,i},
$$

where $\mu_i(w, t)$ ($\text{yr}^{-1}$) is total per-capita mortality and $J_i(w, t)$ ($\text{yr}^{-1}$) is the net flux of individuals across mass threshold $w$:

$$
J_i(w, t) = g_i(w, t)\,N_i(w, t) - \frac{1}{2}\frac{\partial}{\partial w}\bigl(d_i(w, t)\,N_i(w, t)\bigr).
$$

Here $g_i(w, t)$ ($\text{g}\,\text{yr}^{-1}$) is the mean somatic growth rate (advective drift) and $d_i(w, t)$ ($\text{g}^2\,\text{yr}^{-1}$) is the somatic growth diffusivity. Notice that the derivative $\partial_w$ acts on the product $d_i N_i$, as follows directly from a second-order Kramers–Moyal expansion of discrete feeding jumps [@datta2010].

The total diffusivity $d_i(w, t) = d_{p,i}(w, t) + d_{\text{ext},i}(w)$ comprises two additive biological contributions. First, when predation diffusion is enabled (`use_predation_diffusion(params) <- TRUE`), `mizerDiffusion()` computes the endogenous variance of mass increments arising from discrete captures of prey of mass $w_p$, assimilated with efficiency $\alpha_i$ and reduced to somatic retention by the reproductive investment fraction $\psi_i(w)$:

$$
d_{p,i}(w, t) = \alpha_i^2\bigl(1 - \psi_i(w)\bigr)^2 \bigl(1 - f_i(w, t)\bigr)\,\gamma_i(w) \int \phi_i(w / w_p)\,N_c^{(i)}(w_p, t)\,w_p^2\,\mathrm{d}w_p,
$$

where $\gamma_i(w)$ is search volume, $\phi_i(w/w_p)$ is the predation kernel, $f_i(w, t) \in [0, 1)$ is the Type II feeding level, and $N_c^{(i)}(w_p, t)$ is the effective prey abundance weighted by species interaction strengths. Under a power-law prey spectrum $N_c(w_p) \propto w_p^{-\lambda}$, $d_{p,i}(w)$ scales allometrically as $w^{n + 3 - \lambda}$ (yielding $w^{7/4}$ for canonical exponents $n = 3/4$ and $\lambda = 2 + q - n \approx 2.05$). Second, `setExtDiffusion()` supplies an external diffusivity $d_{\text{ext},i}(w)$ capturing non-modelled individual growth variability—such as fine-scale prey patchiness, behavioural foraging differences or genetic variation—defaulting to the power law $d_{\text{ext},i}(w) = D_{\text{ext},i}\,w^{n+1}$ parameterised by species column `D_ext`.

At the species-specific egg mass $w_{\min,i}$, recruits enter the spectrum via a Robin flux boundary condition equating the inward boundary flux to the realised density-dependent recruitment rate $R_{\text{dd},i}(t)$ ($\text{numbers}\,\text{yr}^{-1}$):

$$
J_i(w_{\min,i}, t) = R_{\text{dd},i}(t), \qquad N_i(w, t) = 0 \quad \text{for } w < w_{\min,i}.
$$

At the upper grid cutoff $w_{\max,i}$, a zero-flux condition $J_i(w_{\max,i}, t) = 0$ prevents artificial leakage of biomass. Crucially, `mizer` distinguishes between a species' physiological asymptotic mass $w_{\infty,i}$ (`w_inf`, where deterministic somatic growth ceases as reproductive allocation reaches unity at `w_repro_max`) and its computational upper boundary $w_{\max,i}$ (`w_max`, defaulting to $1.5\,w_{\infty,i}$). Providing headroom $w_{\max,i} > w_{\infty,i}$ accommodates individuals that grow past $w_{\infty,i}$ through positive growth diffusion without artificial pile-up against a reflecting wall.

### 2.2 Size-selective predation, bioenergetics and reproduction

Ecological coupling among species emerges from size-resolved encounter and mortality integrals. Consumer $i$ of mass $w$ searches water at volumetric rate $\gamma_i(w) = \gamma_i\,w^q$ and selects prey of mass $w_p$ according to a dimensionless predation kernel $\phi_i(w / w_p)$, most commonly the lognormal kernel (`lognormal_pred_kernel()`, parameterised by preferred predator–prey mass ratio $\beta_i$ and logarithmic width $\sigma_i$), or alternatives including `box_pred_kernel()`, `power_law_pred_kernel()` and `truncated_lognormal_pred_kernel()`. Total encountered food rate $E_i(w, t)$ ($\text{g}\,\text{yr}^{-1}$) sums resource prey, consumer prey weighted by dimensionless spatial/behavioural overlap $\theta_{ij} \in [0, 1]$ (`interaction_matrix(params)`), fixed external food `ext_encounter(params)` and optional dynamical component contributions $E_{\text{other},i}$:

$$
E_i(w, t) = \gamma_i(w) \int \phi_i(w / w_p) \left[ \theta_{iR}\,N_R(w_p, t) + \sum_{j=1}^{S} \theta_{ij}\,N_j(w_p, t) \right] w_p\,\mathrm{d}w_p + E_{\text{ext},i}(w) + E_{\text{other},i}(w, t).
$$

Consumption saturates via a Holling Type II functional response with maximum intake rate $I_{\max,i}(w) = h_i\,w^n$, yielding the dimensionless feeding level

$$
f_i(w, t) = \frac{E_i(w, t)}{E_i(w, t) + I_{\max,i}(w)}.
$$

Assimilated intake $\alpha_i\,f_i(w, t)\,I_{\max,i}(w)$ first pays basal and activity metabolic maintenance $k_{\text{met},i}(w) = k_{s,i}\,w^p + k_i\,w$, leaving net production rate

$$
E_{r,i}(w, t) = \max\left(0,\; \alpha_i\,f_i(w, t)\,I_{\max,i}(w) - k_{\text{met},i}(w)\right).
$$

Net production is partitioned between somatic growth $g_i(w, t) = \bigl(1 - \psi_i(w)\bigr) E_{r,i}(w, t)$ and reproductive energy $E_{\text{repro},i}(w, t) = \psi_i(w)\,E_{r,i}(w, t)$ via a smooth maturation ogive unified with life-history allocation:

$$
\psi_i(w) = \left[ 1 + \left(\frac{w}{w_{\text{mat},i}}\right)^{-U_i} \right]^{-1} \left(\frac{w}{w_{\text{repro\_max},i}}\right)^{1 - n},
$$

clamped at $\psi_i(w) = 1$ for $w \ge w_{\text{repro\_max},i}$ (where $w_{\text{repro\_max},i} = w_{\infty,i}$ by default) and $\psi_i(w) = 0$ below `w_min_idx`. Under constant feeding level $f_i(w) = f_0$ and metabolic exponent $p = n$, this allocation rule generates a juvenile-to-adult growth curve that closely matches a von Bertalanffy trajectory approaching $w_{\infty,i}$ [@hartvig2011; @andersen2019].

Reproductive output integrates adult investment across the population to give density-independent egg production $R_{\text{di},i}(t)$ ($\text{numbers}\,\text{yr}^{-1}$),

$$
R_{\text{di},i}(t) = \frac{\epsilon_i}{2\,w_{\min,i}} \int_{w_{\min,i}}^{w_{\max,i}} N_i(w, t)\,E_{\text{repro},i}(w, t)\,\mathrm{d}w,
$$

where the factor $1/2$ accounts for female spawning biomass and $\epsilon_i \in [0, 1]$ (`erepro`) represents reproductive efficiency encompassing gonad conversion and early larval survival prior to $w_{\min,i}$. A physiological stock–recruitment function (`RDD`, defaulting to `BevertonHoltRDD()`) maps $R_{\text{di},i}$ to realised recruitment:

$$
R_{\text{dd},i}(t) = \frac{R_{\text{di},i}(t)}{1 + R_{\text{di},i}(t) / R_{\max,i}},
$$

representing finite nursery habitat, egg cannibalism or unmodelled larval density dependence [@andersen2019; @canales2020].

Conservation of mass requires that every gram of consumer prey ingested inflicts a corresponding per-capita predation mortality $\mu_{p,j}(w_p, t)$ on prey species $j$:

$$
\mu_{p,j}(w_p, t) = \sum_{i=1}^{S} \theta_{ji} \int \phi_i(w / w_p)\,\bigl(1 - f_i(w, t)\bigr)\,\gamma_i(w)\,N_i(w, t)\,\mathrm{d}w,
$$

yielding total consumer mortality $\mu_i(w, t) = \mu_{p,i}(w, t) + \mu_{\text{ext},i}(w) + F_i(w, t) + \mu_{\text{other},i}(w, t)$, where $\mu_{\text{ext},i}(w)$ (`ext_mort()`, defaulting to allometric background mortality $z_{0,i}\,w^{n-1}$) captures senescence and non-modelled predators, and $F_i(w, t)$ is fishing mortality.

### 2.3 Resource spectra, multi-gear fishing selectivity and allometric parameterisation

Lower-trophic planktonic prey are represented on an extended size grid $[w_{\min,R}, w_{\text{cutoff},R}]$ by the resource density $N_R(w, t)$, governed by default via semi-chemostat renewal (`resource_semichemostat()`, with `resource_logistic()` and `resource_constant()` available as alternatives):

$$
\frac{\partial N_R(w, t)}{\partial t} = r_R(w)\bigl[c_R(w) - N_R(w, t)\bigr] - \mu_R(w, t)\,N_R(w, t),
$$

with intrinsic regeneration rate $r_R(w) = r_{p}\,w^{n-1}$, carrying capacity $c_R(w) = \kappa\,w^{-\lambda}$ up to $w_{\text{cutoff},R}$, and consumer predation mortality $\mu_R(w_p, t)$ defined analogously to $\mu_{p,j}(w_p, t)$. Because semi-chemostat renewal decouples resource size classes from one another and represents continuous influx from unmodelled nutrient–phytoplankton pools, it avoids oscillatory extinction cascades at the base of the food web [@cuesta2018].

Harvesting is parameterised through a dedicated `gear_params` data frame mapping multiple fleets or gears $g$ onto target species $i$. Fishing mortality decomposes cleanly into time-varying fleet effort $E_g(t)$, species–gear catchability $Q_{g,i}$ and size-based gear selectivity $S_{g,i}(w) \in [0, 1]$:

$$
F_i(w, t) = \sum_{g} E_g(t)\,Q_{g,i}\,S_{g,i}(w), \qquad Y_{g,i}(t) = \int_{w_{\min,i}}^{w_{\max,i}} E_g(t)\,Q_{g,i}\,S_{g,i}(w)\,N_i(w, t)\,w\,\mathrm{d}w,
$$

with built-in selectivity functions covering trawls, gillnets and dome-shaped hooks (`sigmoid_length()`, `double_sigmoid_length()`, `knife_edge()` and `sigmoid_weight()`), converting length thresholds ($L_{25}, L_{50}$) to mass via $w = a_i l^{b_i}$.

Building a multi-species model via `newMultispeciesParams()` requires as few as two empirical inputs per species—`species` and asymptotic mass `w_inf` (or asymptotic length `l_inf` with $a_i, b_i$)—while seamlessly incorporating richer life-history, diet and fishery data where available. Physiological exponents default to metabolic scaling values ($n = 3/4$, $p = 3/4$, $q = 0.8$, $\lambda = 2 + q - n = 2.05$, $\alpha_i = 0.6$, $\beta_i = 30$, $\sigma_i = 2$) [@hartvig2011; @andersen2016]. When maximum intake $h_i$ is omitted, `mizer` derives it from the von Bertalanffy parameter `k_vb` (or `age_mat` and `w_mat`) so juveniles at reference feeding level $f_0 = 0.6$ realise observed somatic growth, and computes the search coefficient $\gamma_i$ (`get_gamma_default()`) needed to achieve $f_0$ on the reference resource spectrum above critical feeding level $f_c = 0.2$. Separating user-supplied columns (`given_species_params(params)`) from derived defaults (`calculated_species_params(params)`) ensures that updating a fundamental trait via `species_params(params) <-` recomputes downstream defaults and rate arrays without overwriting empirical inputs.

---

## 3 Steady-State-First Calibration and Resilience Workflow

### 3.1 Decoupling trophic steady states from stock–recruitment dynamics

Calibrating a nonlinear multi-species PDE to empirical survey observations poses a classic inverse problem. In early size-spectrum workflows [@scott2014], practitioners guessed stock–recruitment ceilings $R_{\max,i}$, resource parameters $\kappa$ and search coefficients $\gamma_i$, and projected the full coupled PDE forward in time hoping abundances would settle near observed values. Because predation mortality $\mu_{p,i}(w)$, feeding level $f_i(w)$ and density-dependent recruitment $R_{\text{dd},i}$ form tight positive and negative feedback loops, slight initial parameter mismatches routinely cause species to starve, overgrow or go extinct during transient spin-up, confounding trophic calibration with arbitrary recruitment assumptions.

`mizer` v3 resolves this bottleneck through a **steady-state-first** calibration paradigm built around two distinct equilibrium operators with complementary mathematical contracts: `tuneSteadyState()` and `findSteadyState()` (Table 1).

During model construction and calibration, `tuneSteadyState(params)` solves for a consumer steady state $N_i^*(w)$ while **holding each species' recruitment flux $R_{\text{dd},i}^*$ and the resource spectrum $N_R^*(w)$ fixed** at their current values (derived from the initial abundances at $w_{\min,i}$). With $R_{\text{dd},i}$ and $N_R$ held constant, species interact solely through adult-on-juvenile predation and prey competition, avoiding stock–recruitment collapse. Once the consumer size distributions $N_i^*(w)$ converge so that $\partial_t N_i^* = 0$ across all size bins, `tuneSteadyState()` inverts the boundary equations algebraically:

1. From the converged spectrum $N_i^*(w)$ and realised growth/reproductive rate $E_{\text{repro},i}^*(w)$, `getRequiredRDD()` and `getRDI()` determine the reproductive efficiency $\epsilon_i$ (`erepro`) and/or recruitment parameter (`R_max`) needed so the Beverton–Holt curve reproduces $R_{\text{dd},i}^*$ endogenously at steady state (preserving the user's chosen `reproduction_level` by default via `preserve = "reproduction_level"`).
2. Simultaneously, solving $\partial_t N_R^* = 0$ updates the resource carrying capacity array $c_R(w) = N_R^*(w)\bigl[1 + \mu_R^*(w)/r_R(w)\bigr]$ so that the held resource spectrum $N_R^*(w)$ is itself an exact equilibrium under consumer predation $\mu_R^*(w)$ (or `resource_constant` if desired).

By contrast, `findSteadyState(params)` **modifies no biological parameters whatsoever**: it allows consumer spectra $N_i$, births $R_{\text{dd},i}$ and resource $N_R$ to co-evolve freely to the true dynamical equilibrium of the already-parameterised system, serving as the tool for counterfactual scenario analysis (e.g. evaluating long-term community reorganisation under altered fishing effort). Both functions attach a structured `"convergence"` attribute recording the attractor class, residual norm and termination criterion, verified independently by `isSteady()` and `getSteadyResidual()`.

### 3.2 Exploiting global scaling symmetry and iterative species matching

Once an initial steady state is established via `tuneSteadyState()`, `mizer` exploits an exact **global scaling symmetry** of the size-spectrum equations to match overall ecosystem scale in a single step. Multiplying consumer densities $N_i(w)$, resource density $N_R(w)$, resource capacity $\kappa$ (`cc_pp`) and maximum recruitment $R_{\max,i}$ uniformly by a positive scalar $c > 0$, while simultaneously dividing volumetric search rates $\gamma_i(w)$ by $c$, leaves all encounter rates $E_i(w)$, feeding levels $f_i(w)$, somatic growth rates $g_i(w)$, predation mortalities $\mu_{p,i}(w)$ and reproductive efficiencies $\epsilon_i$ **strictly invariant**. `calibrateBiomass(params)` (and `calibrateNumber(params)`) computes the exact ratio $c = \sum_i B_{\text{obs},i} / \sum_i B_{\text{mod},i}$ over all species with empirical survey targets `biomass_observed` (integrated above species-specific survey catchability cutoffs `biomass_cutoff`, such as trawl selectivity thresholds) and rescales the model instantaneously without disturbing the steady state (Figure 1a).

Matching individual species' relative abundances and growth rates then breaks global symmetry because altering a single species' abundance changes the prey field and predation mortality experienced by every other guild. `mizer` provides two targeted adjustment operators that are alternated with `tuneSteadyState()` in a rapidly converging fixed-point loop:

- `matchBiomasses(params)` (or `matchNumbers(params)`) scales each species' spectrum $N_i(w)$ proportionally by $B_{\text{obs},i} / B_{\text{mod},i}$, adjusting its boundary recruitment flux accordingly.
- `matchGrowth(params)` rescales each species' physiological foraging and metabolic time-scale rates (`h`, `gamma`, `ks` and `k` in unison, preserving feeding level $f_i(w)$ and critical feeding threshold $f_c$) so that emergent size-at-age trajectories integrated from $g_i^*(w)$ pass cleanly through empirical maturation targets (`w_mat`, `age_mat` or `k_vb`).

Because each adjustment slightly perturbs interspecific predation mortalities $\mu_{p,j}(w)$, repeating `params |> matchGrowth() |> matchBiomasses() |> tuneSteadyState()` five to ten times drives total relative biomass error to machine precision while keeping growth curves, ontogenetic diets and trophic interactions mutually consistent (Figure 1b–d). Observed fleet catches (`gear_params(params)$yield_observed`) provide an independent validation check via `plotYieldObservedVsModel()` and `plotFMort()`.

| Workflow stage | Primary `mizer` functions | Parameters updated | Mathematical invariant preserved |
| :--- | :--- | :--- | :--- |
| **1. Model assembly** | `newMultispeciesParams()`, `gear_params()<-` | Allometric defaults (`h`, `gamma`, `ks`, `z0`), size grids, gear selectivities | Separation of `given_species_params` vs `calculated_species_params` |
| **2. Initial trophic equilibrium** | `tuneSteadyState()` | `erepro` / `R_max`, `cc_pp` (`N_i^*(w)` solved) | Held recruitment $R_{\text{dd},i}^*$ and resource $N_R^*(w)$ become exact steady states |
| **3. Ecosystem-scale calibration** | `calibrateBiomass()` / `calibrateNumber()` | Uniform scaling $c$ applied to $N_i, N_R, \kappa, R_{\max,i}$ and $c^{-1}$ to $\gamma_i$ | Exact steady state, $f_i(w)$, $g_i(w)$, $\mu_{p,i}(w)$ and `erepro` unchanged |
| **4. Species-level iterative matching** | `matchGrowth()`, `matchBiomasses()` + `tuneSteadyState()` | Species-specific `h, gamma, ks, k`, spectra $N_i^*(w)$, `erepro`, `cc_pp` | Empirical `biomass_observed` (above `biomass_cutoff`) and maturation schedule matched |
| **5. Resilience specification** | `reproduction_level()<-` / `setBevertonHolt()` | Pairs $(\epsilon_i, R_{\max,i})$ along constant-$R_{\text{dd},i}^*$ hyperbola | Calibrated trophic steady state $N_i^*(w)$ remains byte-for-byte identical |
| **6. Scenario & stability analysis** | `findSteadyState()`, `getStability()`, `scanModel()`, `project()` | State trajectories $N_i(w, t)$ only (no hidden parameter re-tuning) | Full endogenous consumer–resource–reproduction feedbacks active |

: **Table 1.** Structured steady-state-first workflow for parameterising, calibrating and analysing an empirical multi-species size-spectrum model in `mizer`.

### 3.3 Orthogonal specification of reproductive resilience

A pivotal insight formalised in `mizer` is that **steady-state abundance observations alone cannot identify stock–recruitment curvature**, yet the two can be parameterised orthogonally [@andersen2019]. Along a Beverton–Holt curve at steady state, any pair $(\epsilon_i, R_{\max,i})$ satisfying

$$
R_{\text{dd},i}^* = \frac{\epsilon_i\,\tilde{R}_{\text{di},i}^*}{1 + (\epsilon_i\,\tilde{R}_{\text{di},i}^*) / R_{\max,i}}, \qquad \text{where } \tilde{R}_{\text{di},i}^* = R_{\text{di},i}^* / \epsilon_i,
$$

produces the **exact same** recruitment flux $R_{\text{dd},i}^*$ and therefore the exact same steady-state ecosystem $N_i^*(w)$. `mizer` re-parameterises this one-parameter family via the dimensionless **reproduction level** $\rho_i \in [0, 1)$, accessible via `reproduction_level(params)`:

$$
\rho_i = \frac{R_{\text{dd},i}^*}{R_{\max,i}}, \qquad \text{so that} \qquad R_{\max,i} = \frac{R_{\text{dd},i}^*}{\rho_i}, \quad \epsilon_i = \frac{R_{\text{dd},i}^*}{\tilde{R}_{\text{di},i}^*\,(1 - \rho_i)}.
$$

When $\rho_i \to 0$ ($R_{\max,i} \to \infty$), recruitment is linear in spawning output ($R_{\text{dd},i} = R_{\text{di},i}$), yielding high sensitivity to fishing mortality; as $\rho_i \to 1$, recruitment saturates strongly near $R_{\max,i}$, conferring steep compensatory resilience to harvesting (subject to physiological feasibility $\epsilon_i < 1$). Setting `reproduction_level(params) <- rho` (or `setBevertonHolt(params, erepro = ...)` / `R_max = ...`) moves along this iso-recruitment hyperbola **without altering the calibrated steady state**, enabling clean fitting of $\rho_i$ against historical fishing-mortality yield curves (`plotYieldVsF()`) or time-series catch data independently of static biomass calibration.

![Steady-state calibration and emergent trophic ecology of the 12-species North Sea model (`NS_params`). **(a)** Modelled versus observed survey biomasses after single-step global scaling with `calibrateBiomass()` (dashed line: 1:1 agreement; open circle: species lacking survey biomass target). **(b)** Exact agreement between modelled and observed species biomasses achieved by alternating `matchBiomasses()` and `tuneSteadyState()`. **(c)** Resulting steady-state Sheldon biomass density spectra ($N_i(w)\,w$, $\text{g}$) across body mass for all 12 fish species and the planktonic resource spectrum (`plotSpectra()`). **(d)** Emergent ontogenetic diet shift of Atlantic cod (*Gadus morhua*, `plotDiet()`), transitioning continuously from planktonic resource in larvae ($<0.1\text{ g}$) through forage fish (sprat, sandeel, Norway pout) at intermediate sizes to herring, gadoids and conspecific cannibalism in large adults ($>10^3\text{ g}$).](figure1_northsea.png)

---

## 4 Conservative Finite-Volume Discretisation and Second-Order Solvers

### 4.1 Geometric weight grid, variable co-location and single-quadrature FFT convolutions

Consumer body mass in multi-species communities spans six to eight orders of magnitude ($10^{-3}\text{ g}$ to $10^5\text{ g}$), extended downward to $w_{\min,R} \approx 10^{-12}\text{ g}$ for unicellular plankton. `mizer` discretises mass onto a geometric (logarithmically uniform) grid where consecutive bin boundaries satisfy a constant ratio $\beta = w_{j+1} / w_j > 1$, so bin widths $\Delta w_j = w_{j+1} - w_j = (\beta - 1)w_j$ grow proportionally to body mass.

To guarantee strict conservation of abundance across arbitrary grid resolutions and nonzero diffusion, `mizer` v3 formulates spatial discretisation via a **cell-centred finite-volume scheme** rather than naive point collocation. Integrating the transport PDE over cell $j = [w_j, w_{j+1}]$ and dividing by $\Delta w_j$ yields the exact semi-discrete balance law for the cell-averaged density $N_{i,j}(t) = \Delta w_j^{-1} \int_{w_j}^{w_{j+1}} N_i(w, t)\,\mathrm{d}w$:

$$
\frac{\mathrm{d} N_{i,j}}{\mathrm{d} t} + \frac{J_{i,j+1} - J_{i,j}}{\Delta w_j} = -\mu_{i,j}\,N_{i,j}.
$$

This balance law enforces a strict structural distinction between two classes of discrete variables:

- **Interface quantities** evaluated at left cell boundaries $w_j$: the advective somatic growth velocity $g_{i,j} = g_i(w_j, t)$ and total interface flux $J_{i,j}$.
- **Cell-averaged quantities** co-located with $N_{i,j}$ over $[w_j, w_{j+1}]$ (associated with geometric cell centres $w_j^c = \sqrt{w_j w_{j+1}}$, the midpoint in $\ln w$): per-capita mortality $\mu_{i,j}$ and growth diffusivity $d_{i,j}$, whose interface diffusive flux takes the centred gradient form $J_{i,j}^{\text{diff}} = -\frac{1}{2}(d_{i,j} N_{i,j} - d_{i,j-1} N_{i,j-1}) / \Delta w_{j-1}^c$, with $\Delta w_{j-1}^c = w_j^c - w_{j-1}^c$.

Evaluating encounter rates $E_i(w)$, predation mortalities $\mu_{p,j}(w_p)$ and predation diffusivities $d_{p,i}(w)$ directly across $N_{\text{full}}$ plankton-to-whale size bins would require $O(S\,N_{\text{full}}^2)$ operations per time step, dominating runtime. On the geometric grid $w_k = w_0 \beta^k$, however, any scale-invariant predation kernel $\phi_i(w / w_p) = \phi_i(\beta^{k - m})$ becomes a discrete **linear convolution** in the logarithmic index $k - m$. `mizer` pre-computes the discrete Fourier transforms of the reversed and forward kernels at parameter setup (`ft_pred_kernel_e`, `ft_pred_kernel_p` and `ft_pred_kernel_d`) and evaluates all three trophic integrals in $O(N_{\text{full}} \log N_{\text{full}})$ time via FFTs (`mvfft()`).

Combining finite-volume bin averages with convolution integrals requires care: if both a size-dependent factor (such as prey mass $w_p$ inside $E_i$) and the quadrature rule bin-average the same power law, the integral acquires a silent uniform bias of $(1 + \beta)/2$ ($\approx 10\%$ on a 100-bin grid where $\beta \approx 1.19$). `mizer` enforces a **single-quadrature invariant** across all integrators and indicators (`sizeIntegral()`, `getBiomass()`, `getSSB()`, `getYield()`, `encounter_kernel()`): each bin integral is performed in one canonical place, averaging pure power laws analytically over $[w_j, w_{j+1}]$ and composite rates trapezoidally when `second_order_w(params)$bin_average` is enabled.

### 4.2 Eliminating artificial numerical diffusion: van Leer flux reconstruction

In the historical first-order scheme (`second_order_w = FALSE`, retained as default for exact backwards compatibility with legacy models), the advective interface density $J_{i,j}^{\text{adv}} = g_{i,j} N_i(w_j)$ uses pure upwind reconstruction $N_i(w_j) \approx N_{i,j-1}$ paired with semi-implicit Euler stepping. Taylor expansion reveals that this discretisation solves a modified PDE carrying an **artificial numerical diffusivity**:

$$
d_{\text{num}}^{\text{Euler}}(w) \approx g(w)\,w\ln\beta + g(w)^2\,\Delta t.
$$

On a standard 100-bin consumer grid ($\ln\beta \approx 0.17$), $g(w)\,w\ln\beta$ routinely exceeds biological predation diffusion $d_{p,i}(w)$ by an order of magnitude, smearing pulsed annual cohorts within two to three years even when biological diffusion vanishes and damping oscillatory cohort waves.

Setting `second_order_w(params) <- TRUE` eliminates leading spatial numerical diffusion by reconstructing the interface density at $w_j$ (the logarithmic midpoint between $w_{j-1}^c$ and $w_j^c$) with flux-limited weight $\chi_{i,j}$:

$$
J_{i,j}^{\text{adv}} = g_{i,j}\left[ N_{i,j-1} + \frac{1}{2}\,\chi_{i,j}\,(N_{i,j} - N_{i,j-1}) \right].
$$

Setting $\chi_{i,j} \equiv 1$ (`flux = "centred"`) yields unlimited second-order centred advection, whereas the default `flux = "van_leer"` applies the total-variation-diminishing (TVD) flux limiter of @vanleer1974:

$$
\chi(r) = \frac{r + |r|}{1 + |r|}, \qquad r_{i,j} = \frac{N_{i,j-1} - N_{i,j-2}}{N_{i,j} - N_{i,j-1}},
$$

forced to $\chi_{i,j} = 0$ at the recruitment boundary ($j \le j_{\min,i} + 2$). Where the spectrum is smooth ($r_{i,j} \approx 1$), $\chi \to 1$ recovers second-order accuracy free of upwind diffusion; at cohort extrema or steep fronts, $\chi \to 0$ reverts smoothly to monotone upwind transport (Figure 2a). Evaluating $\chi_{i,j}$ from start-of-step (or midpoint) densities keeps the transport system linear in $N^{t+1}$ for the $O(N_w)$ C++ Thomas tridiagonal solver (`project_simple.cpp`), while an exponential moving-average relaxation $\chi^{(t)} = \frac{1}{2}\chi(N) + \frac{1}{2}\chi^{(t-1)}$ prevents limiter chatter during steady-state convergence.

### 4.3 L-stable second-order time integration and smoothness of custom rates

Advancing stiff advection–diffusion–reaction grids accurately in time requires an integrator that avoids explicit diffusion step restrictions without introducing artificial damping or boundary ringing. `project()` offers three time-stepping schemes via `method =`:

1. `"euler"` (default): first-order semi-implicit step evaluating ecological rates at $t$, solving $(I - \Delta t\,L^t) N^{t+1} = N^t + \Delta t\,q^t$ implicitly via the C++ tridiagonal solver, and updating the resource via its exact fixed-mortality semi-chemostat solution.
2. `"predictor_corrector"`: an Euler predictor $\hat{N}^{t+1}$ forms midpoint rates $r^{t+1/2} = \frac{1}{2}(r^t + \hat{r}^{t+1})$ for a Crank–Nicolson consumer step and midpoint-mortality resource update. Because Crank–Nicolson is A-stable rather than L-stable ($\lim_{z \to -\infty} |R_{\text{CN}}(z)| = 1$), stiff boundaries or coarse steps ($\Delta t \gtrsim 0.4\text{ yr}$) excite oscillatory ringing (Figure 2b).
3. `"tr_bdf2"` (aliased as `"second_order"`): a two-stage composite Runge–Kutta / backward-differentiation integrator [@bank1985; @hosea1996] driven by midpoint ecological rates. Stage 1 advances via a trapezoidal sub-step to $t + \gamma\Delta t$, and Stage 2 completes the step via BDF2:

$$
\left(I - \frac{\gamma}{2}\Delta t\,L\right) N^{t+\gamma} = \left(I + \frac{\gamma}{2}\Delta t\,L\right) N^t + \gamma\Delta t\,q,
$$

$$
\bigl(I - \alpha\Delta t\,L\bigr) N^{t+1} = \frac{1}{\gamma(2 - \gamma)} N^{t+\gamma} - \frac{(1 - \gamma)^2}{\gamma(2 - \gamma)} N^t + \alpha\Delta t\,q.
$$

Choosing $\gamma = 2 - \sqrt{2} \approx 0.5858$ makes $\alpha = (1 - \gamma)/(2 - \gamma) = \gamma/2$, so **both implicit stages share the identical tridiagonal matrix** $I - \alpha\Delta t\,L$. Moreover, TR-BDF2 is **L-stable** ($\lim_{z \to -\infty} |R_{\text{TR-BDF2}}(z)| = 0$), achieving clean $O(\Delta t^2)$ global convergence across fine time steps while damping stiff high-frequency boundary transients without ringing at coarse $\Delta t$ (Figure 2b).

Regardless of time-stepper order, formal convergence and Jacobian differentiability require ecological rate functions to vary **continuously** with state variables $N$. When users implement custom management or ecological rules via `setRateFunction()`—such as an abrupt fishery closure `if (SSB < B_lim) effort <- 0`—crossing the discontinuity mid-step degrades every integrator to $O(\Delta t)$ chattering whose apparent plateau shifts with $\Delta t$, stalls Newton root-finding, and corrupts finite-difference Jacobians (Figure 2c, left). Regularising step switches into continuous piecewise-linear ramps—such as an ICES hockey-stick harvest control rule ramping effort smoothly between $B_{\lim}$ and $B_{\text{trigger}}$—restores full $\Delta t$-independent convergence and well-posed stability diagnostics (Figure 2c, right).

![Advection–diffusion cohort propagation, second-order time-stepping convergence and regularisation of state-dependent rate switches. **(a)** Five pulsed annual cohorts in a single-species model (`no_w = 400`, `second_order_w = TRUE`, `method = "tr_bdf2"`) at $t = 5\text{ yr}$, showing sharp cohort preservation under pure deterministic advection (`No diffusion`, solid curve) versus natural cohort broadening under endogenous jump-growth variance (`Predation diffusion`, dashed curve; `plotSpectra2()`). **(b)** Log–log scaling of relative $L_2$ error in the final consumer spectrum against time step $\Delta t$, comparing `"euler"`, `"predictor_corrector"` and L-stable `"tr_bdf2"` (`method = "second_order"`), which maintains $O(\Delta t^2)$ accuracy at small $\Delta t$ without large-step instability. **(c)** Spawning stock biomass ($SSB / SSB_0$) of North Sea cod under a discontinuous threshold closure (`Hard switch`, left) exhibiting $\Delta t$-dependent chattering above $B_{\lim}$, versus step-invariant convergence under a continuous linear ramp (`Ramped rule`, right).](figure2_numerics.png)

---

## 5 Equilibrium Root-Finding, Linear Stability and Bifurcation Analysis

### 5.1 Algebraic Newton root-finding on stable and unstable equilibria

Size-structured food webs do not always settle onto stationary fixed points: size-selective predation, cannibalism and reproductive time lags frequently destabilise equilibria via Hopf bifurcations, generating persistent multi-year cohort cycles [@datta2011; @canales2020]. Whenever an equilibrium is dynamically unstable, forward time-projection (`solver = "project"`) diverges away from the fixed point onto a limit cycle or alternative state, making it impossible to inspect steady-state trophic rates or compute local eigenvalues along a parameter branch.

Both `tuneSteadyState()` and `findSteadyState()` therefore provide a direct algebraic root-finder (`solver = "newton"`). At any candidate state, setting $N_{i,j}^{t+1} = N_{i,j}^t = N_{i,j}^*$ in the semi-discrete finite-volume balance equations eliminates $\Delta t$ identically, yielding a time-step-independent tridiagonal transport system $\tilde{A}_{i,j} N_{i,j-1}^* + \tilde{B}_{i,j} N_{i,j}^* + \tilde{C}_{i,j} N_{i,j+1}^* = \tilde{S}_{i,j}$ coupled nonlinearly through the ecological rate integrals $F(N^*) = 0$. In `findSteadyState(..., solver = "newton")`, `mizer` solves the full coupled nonlinear algebraic system for active consumer and resource log-densities using damped Newton–Broyden iteration (`nleqslv`). Because algebraic root-finding depends only on local nonsingularity of the Jacobian $\partial F / \partial N$ and not on the signs of its eigenvalues, `solver = "newton"` converges quadratically to equilibria on **both** the stable and unstable sides of a bifurcation, typically driving `getSteadyResidual()` to $10^{-12}\text{–}10^{-14}$ in seconds.

### 5.2 Continuous-time Jacobian eigenvalues versus discrete solver damping

Once a fixed point $N^*$ is located, `getStability(params)` characterises its local linear stability by forming the Jacobian matrix $J = \partial F / \partial X \big|_{X^*}$ of the **continuous-time** semi-discretised ODE system $\mathrm{d}X/\mathrm{d}t = F(X)$ across all active consumer and resource bins ($n_{\text{active}}$ state variables). Because abundances span more than 15 orders of magnitude from picoplankton to apex gadoids, uniform additive finite-difference perturbations would underflow in large fish bins or overwhelm plankton bins; `getStability()` therefore perturbs each active bin multiplicatively ($\delta X_k = h\,X_k^*$ with default relative step $h = 10^{-6}$) and rescales columns to obtain scale-invariant derivatives completely independent of any projection time step $\Delta t$.

From the full eigenspectrum $\{\lambda_m\}$ of $J$, `getStability()` reports:

- `max_real_part` $= \max_m \operatorname{Re}(\lambda_m)$ ($\text{yr}^{-1}$), determining continuous-time asymptotic stability ($\max_m \operatorname{Re}(\lambda_m) < 0$);
- `leading_oscillatory_eigenvalue` $\lambda_1 = \sigma_1 + \mathrm{i}\omega_1$, together with its implied `dominant_period` $T = 2\pi / |\omega_1|$ ($\text{yr}$) and complex spatial–taxonomic eigenvector.

Crucially, `mizer` pairs `getStability()` with `getDiscreteStability(params, dt, method)`, which differentiates the **discrete one-step time-marching map** $X^{t+1} = \Phi_{\Delta t}(X^t)$ directly and computes its spectral radius $\rho(\partial \Phi_{\Delta t} / \partial X)$. Comparing the two exposes a major methodological hazard in size-spectrum simulation: for a lightly unstable continuous-time mode $\lambda = \sigma + \mathrm{i}\omega$ with $0 < \sigma \ll |\omega|$, backward Euler maps $\lambda$ to the amplification factor $(1 - \Delta t\,\lambda)^{-1}$, whose modulus satisfies

$$
\left|\frac{1}{1 - \sigma\Delta t - \mathrm{i}\omega\Delta t}\right| < 1 \quad \Longleftrightarrow \quad (1 - \sigma\Delta t)^2 + \omega^2\Delta t^2 > 1 \quad \Longleftrightarrow \quad \Delta t > \frac{2\sigma}{\sigma^2 + \omega^2}.
$$

Whenever the time step exceeds $2\sigma / (\sigma^2 + \omega^2)$, **first-order implicit time-stepping flips the sign of the net growth rate**, artificially stabilising a genuinely unstable ecosystem.

### 5.3 Case study: Hopf bifurcation and limit cycles in the North Sea food web

We illustrate the combined power of Newton root-finding, Jacobian eigenanalysis and L-stable second-order integration on the calibrated 12-species North Sea model (`NS_params`, $n_{\text{active}} = 1{,}171$ active size classes across sprat, sandeel, Norway pout, herring, dab, whiting, sole, gurnard, plaice, haddock, cod and saithe plus plankton resource). At historical baseline fishing effort ($E = 1.0$), the steady state is linearly stable ($\max\operatorname{Re}(\lambda) = -0.0174\text{ yr}^{-1}$). Raising fishing effort across all fleets using `findSteadyState(..., solver = "newton")` and evaluating `getStability()` at $E \in \{1.00, 1.05, 1.10, 1.15\}$ yields leading real parts $\{-0.0174, -0.00354, +0.00862, +0.0194\}\,\text{yr}^{-1}$, pinpointing a **Hopf bifurcation** at critical effort $E^* = 1.065$.

At elevated fishing effort $E = 1.5$, `findSteadyState(params, effort = 1.5, solver = "newton")` still locates the exact unstable equilibrium, where `getStability()` identifies a conjugate pair of unstable oscillatory eigenvalues $\lambda_{1,2} = 0.0710 \pm 1.27\mathrm{i}\,\text{yr}^{-1}$ (`stable = FALSE`), predicting growing multi-species oscillations with linear period $T = 2\pi / 1.27 = 4.93\text{ yr}$. Passing this eigenmode to `getOscillationModeSim(params_f15, amplitude = 0.1)` reconstructs a standard `MizerSim` object over one linear period without running a nonlinear simulation, revealing out-of-phase predator–prey cohort resonances led by saithe, herring, haddock and cod (Figure 3a).

Projecting nonlinear dynamics at $E = 1.5$ with `projectUntilSettled(..., method = "tr_bdf2")` confirms that perturbations grow away from the unstable fixed point and settle onto a self-sustained multi-species **limit cycle** (`attractor = "limit_cycle"`, nonlinear period $5.4\text{ yr}$, relative amplitude $0.658$; Figure 3b). Scanning fishing effort continuously across $E \in [1.00, 1.50]$ via `scanModel(..., set_func = scanEffort(), method = "tr_bdf2", style = "envelope")` traces the full bifurcation diagram (Figure 3c): a stationary single-valued curve below $E^*$ opens smoothly into a symmetric oscillation envelope above $E^*$, the hallmark of a **supercritical** Hopf bifurcation free of hysteresis. Moreover, regressing the squared relative cycle amplitude $A^2 = \bigl[(B_{\max} - B_{\min})/\bar{B}\bigr]^2$ against fishing effort yields a linear normal-form fit with $R^2 = 0.999$ across the entire scanned domain ($E \in [1.10, 1.50]$; Figure 3d), extrapolating to zero amplitude at $E = 1.088$ in close agreement with the independent Jacobian root $E^* = 1.065$.

Repeating the exact same $E = 1.5$ runs with the default first-order `"euler"` time step demonstrates both discrete failure modes predicted in Section 5.2 (for $\sigma_1 = 0.071\text{ yr}^{-1}$ and $\omega_1 = 1.27\text{ rad yr}^{-1}$, Euler damps the mode whenever $\Delta t > 2\sigma_1/(\sigma_1^2 + \omega_1^2) = 0.088\text{ yr}$, which includes the standard $\Delta t = 0.1\text{ yr}$):

1. Starting from a $5\%$ perturbation of the unstable steady state, `"euler"` **decays back to the unstable fixed point**, reporting `attractor = "fixed_point"` and `isSteady() == TRUE`, whereas `"tr_bdf2"` grows at the Jacobian rate into the $5.2\text{-yr}$, $64\%$-amplitude limit cycle.
2. Starting far from equilibrium (the $E = 1.0$ state), `"euler"` retains a residual oscillation (`period = 5.1 yr`) but **understates its amplitude four-fold** ($0.166$ versus $0.658$).

Thus, combining `solver = "newton"`, `getStability()` and `method = "tr_bdf2"` prevents numerical damping from masquerading as ecological stability.

![Linear stability, supercritical Hopf bifurcation and nonlinear multi-species limit cycles in the calibrated 12-species North Sea model (`NS_params`). **(a)** Biomass oscillations over one $4.93\text{-yr}$ linear period reconstructed directly from the leading unstable Jacobian eigenvector ($\lambda_1 = 0.0710 + 1.27\mathrm{i}\,\text{yr}^{-1}$) at `effort = 1.5` via `getOscillationModeSim()`. **(b)** Full 30-year nonlinear multi-species limit cycle (`period = 5.4 yr`, relative amplitude $0.658$) obtained via `projectUntilSettled(..., method = "tr_bdf2")`. **(c)** Supercritical Hopf bifurcation envelope (minimum to maximum saithe biomass on the settled attractor) across fishing effort $E \in [1.0, 1.5]$ computed with `scanModel()`. **(d)** Verification of canonical supercritical Hopf normal-form scaling, showing linear growth of squared relative cycle amplitude $A^2$ with fishing effort ($R^2 = 0.999$).](figure3_stability.png)

---

## 6 Composable S3 Architecture, Ecosystem Extensions and Reproducible Workflows

### 6.1 State-bearing `MizerParams` objects, provenance and self-describing arrays

To make complex ecosystem models transparent, inspectable and composable with tidyverse and S3 workflows in R, `mizer` v3 transitioned its core `MizerParams` and `MizerSim` structures from rigid S4 classes to lightweight, validated **S3 list objects** while retaining full `@` and `@<-` accessor compatibility alongside standard `$` and `[[]]` indexing (`validParams()` upgrades legacy saved objects automatically). Crucially, a `MizerParams` object is **self-contained and state-bearing**: `initialN(params)`, `initialNResource(params)` and `initialNOther(params)` hold the model's current ecological state (typically its calibrated steady state), so every diagnostic, summary and plotting generic—`getBiomass()`, `getSSB()`, `getDiet()`, `getTrophicLevel()`, `plotSpectra()`, `plotGrowthCurves()`—operates identically on a static `MizerParams` object or across time on a `MizerSim` trajectory.

Diagnostic and summary outputs are returned as lightweight S3 array wrappers (`ArraySpeciesBySize`, `ArrayTimeBySpecies`, `ArrayTimeBySpeciesBySize`, `ArrayResourceBySize`) carrying unit, dimension and model metadata, enabling instant `plot()`, `summary()` and `as.data.frame()` dispatch without passing auxiliary `MizerParams` arguments. Parameter edits propagate deterministically downward across three architectural layers: modifying a scalar species trait via `species_params(params) <-` records provenance in `given_species_params` and automatically rebuilds dependent size-resolved rate arrays (`search_vol`, `intake_max`, `metab`, `psi`, `ext_mort`, `catchability`, `selectivity`), **unless** a rate array has been explicitly customised via its setter (e.g. `setSearchVolume(params, search_vol = custom_array)`), which stamps a `comment` **freeze attribute** protecting the empirical array from accidental overwrite and notifies the user via `info_level`-controlled structured signals.

### 6.2 Four-tier extension hierarchy and composable S3 subclass dispatch

Rather than forcing modellers to fork `mizer` when adding new ecological processes, `mizer` provides a **four-tier extension architecture** ordered from zero-code array injection to fully coupled dynamical ecosystem compartments (Table 2):

1. **Fixed external rate arrays (Tier 1).** Static size-dependent prey fields, apex-predator mortalities or individual growth variances from external observations or coupled models are injected directly via `ext_encounter(params) <-`, `ext_mort(params) <-` and `setExtDiffusion()`.
2. **Stateless state-dependent rate hooks (Tier 2).** Processes that depend dynamically on current abundances $(N, N_R, t)$ without introducing new state variables—such as starvation mortality when metabolic demand exceeds intake, or opportunistic scavenging—register free-standing functions directly in `other_mort(params)` and `other_encounter(params)`.
3. **Custom core rate functions (Tier 3).** Any of `mizer`'s fourteen internal rate slots (`Encounter`, `FeedingLevel`, `EReproAndGrowth`, `ERepro`, `EGrowth`, `PredRate`, `PredMort`, `FMort`, `Mort`, `RDI`, `RDD`, `ResourceMort`, `Diffusion` and `Rates`) can be replaced by name via `setRateFunction()`, which tests return dimensions immediately upon registration and reads custom parameters from `other_params(params)` (e.g. seasonal reproduction, ramped harvest control rules or Arrhenius thermal scaling).
4. **Coupled dynamical ecosystem components (Tier 4).** `setComponent()` adds arbitrary new state variables $N_{\text{other}}(t)$—scalars, vectors or size spectra such as benthic detritus pools fed by consumer faeces, carrion and fishery discards—integrated alongside consumers at every time step and coupled bidirectionally into `Encounter` and `Mort`.

When packaged for reuse, extensions prepend their S3 class name onto `MizerParams` (e.g. `c("mizerShelf", "MizerParams")`), chain S3 generics (`validParams()`, `setParams()`, `plotDiet()`, `getBiomass()`) cooperatively via `NextMethod()` so **multiple extensions compose on the same model**, and register package provenance via `recordExtension()`. Saving and restoring models via `saveParams()` / `readParams()` (and `setMetadata()`) automatically verifies, loads and upgrades required extension packages across sessions.

| Extension tier | Mechanism / API | State & coupling | Representative packages & ecological scope |
| :--- | :--- | :--- | :--- |
| **Tier 1: External arrays** | `ext_encounter()`, `ext_mort()`, `setExtDiffusion()` | Time-invariant size arrays added to encounter, mortality or growth variance | Coupled Earth-system / apex-predator forcing; empirical cohort dispersion (`D_ext`) |
| **Tier 2: Stateless hooks** | `other_mort()`, `other_encounter()` | State-dependent functions of $(N, N_R, t)$ without extra state variables | `mizerStarvation` (condition- or mass-loss starvation mortality when $E_{r,i} < 0$) |
| **Tier 3: Rate replacement** | `setRateFunction()`, `other_params()` | Replaces any of 14 core rate slots with registration-time shape checks | `therMizer` [@woodworthjefcoats2019; @reum2020] (ocean warming & aerobic scope); seasonal reproduction [@datta2016]; HCRs |
| **Tier 4: Dynamical pools & S3 subclasses** | `setComponent()`, `NextMethod()`, `recordExtension()` | Coupled ODE/PDE state pools $N_{\text{other}}(t)$ + composable S3 generics | `mizerShelf` [@audzijonyte2023] (benthic detritus, carrion & sea-floor habitats); `mizerMR` (multiple size-structured resources); `mizerEcopath` (Ecopath mass-balance import & catch calibration); `mizerBioenergetic` |

: **Table 2.** Four-tier modular extension hierarchy in `mizer` and representative packages across the `sizespectrum` ecosystem.

### 6.3 Single-source documentation and machine-readable AI agent skills

Lowering the technical barrier to mechanistic ecosystem modelling requires operational guidance that prevents subtle methodological errors (such as confusing `tuneSteadyState()` with `findSteadyState()` or introducing discontinuous custom rates). In `mizer` v3, all eleven workflow guides—covering model construction, parameter modification, fishing gears, steady-state calibration, simulation, stability analysis and extension authoring—are maintained as structured **single-source `SKILL.md` modules** shipped inside `inst/skills/` alongside a curated API index (`inst/llms.txt`). Build scripts automatically compile these modules into human-facing website vignettes while preserving `<!-- agent-only -->` diagnostic tables that allow large-language-model coding assistants (configured via `mizerAgents`) to inspect the exact installed `mizer` API and guide researchers step-by-step through rigorous model development.

---

## 7 Discussion

Over the twelve years since its initial version 0.2 prototype [@scott2014], `mizer` has evolved from an exploratory size-spectrum simulator into a comprehensive mathematical, numerical and software platform for empirical ecosystem science [@pennino2026]. Within the broader landscape of aquatic ecosystem models, `mizer` occupies a distinctive methodological niche that complements both aggregate trophic networks and high-dimensional end-to-end simulators. Compared with biomass-pool mass-balance frameworks such as Ecopath with Ecosim [@christensen2004] and Rpath [@lucey2020], `mizer` resolves continuous intra-specific size structure, emergent ontogenetic diet shifts, size-dependent cannibalism and growth-mediated life-history feedbacks from first principles—and through the `mizerEcopath` extension can now be initialised directly from existing Ecopath mass-balance tables. Relative to spatially explicit end-to-end frameworks such as Atlantis [@audzijonyte2019], individual-based OSMOSE [@shin2004] or multi-area statistical assessment toolboxes such as gadget3 [@lentin2022], `mizer`'s allometric defaults and global scaling symmetries reduce parameter dimensionality by one to two orders of magnitude, enabling rapid ensemble uncertainty quantification [@spence2016; @spence2021] and standardised global-to-regional climate intercomparisons in Fish-MIP [@tittensor2018; @lotze2019; @ortegacisneros2025]. Within the family of size- and stage-structured packages, `mizer` generalises discrete length-based marine fish models such as LeMaRns [@spence2020lemarns] and StrathSPACE [@speirs2016] to continuous weight-based advection–diffusion PDEs with arbitrary plankton/benthic resources and modular ecosystem pools, while complementing general physiologically structured bifurcation packages such as PSPManalysis [@deroos2021] with fast FFT trophic convolutions, empirical multi-gear fisheries calibration and forward transient projections.

The methodological advances presented here carry direct ecological implications beyond software convenience:

1. **Separating trophic calibration from reproductive resilience** (Section 3) resolves a long-standing identifiability tension in ecosystem-based fisheries management: static survey biomasses constrain trophic fluxes and relative species scale, whereas dynamic responses to fishing pressure depend independently on `reproduction_level()`, which can now be varied or fitted systematically without invalidating the calibrated food web.
2. **Accounting for biological growth diffusion while removing numerical diffusion** (Section 4) enables size-spectrum models to represent real cohort broadening (`D_ext` and `mizerDiffusion()`) rather than an uncontrolled grid artefact proportional to $\ln\beta$ and $\Delta t$.
3. **Pairing Newton root-finding with Jacobian stability analysis and L-stable TR-BDF2 integration** (Section 5) demonstrates that realistic, empirically calibrated multi-species fish communities—not merely idealised single-spectrum power laws—can cross supercritical Hopf bifurcations under elevated fishing pressure ($E^* = 1.065$ in `NS_params`), and shows that standard first-order implicit solvers can silently mask or severely attenuate these ecological oscillations.

Several practical limitations and frontiers merit emphasis. First, core `mizer` represents a spatially well-mixed or implicitly overlap-weighted (`interaction_matrix`) food web; systems governed by ontogenetic habitat migrations, marine protected areas or sharp oceanographic fronts require coupling multiple `mizer` patches or spatially resolved extensions. Second, while allometric defaults allow rapid model assembly from minimal trait data, empirical projections of mixed-fishery yields remain sensitive to gear selectivity curves (`gear_params`), diet overlap matrices ($\theta_{ij}$) and stock–recruitment resilience ($\rho_i$), underscoring the importance of multi-model ensembles and Bayesian uncertainty propagation [@spence2016; @spence2021]. Third, empirical estimation of external growth diffusivity `D_ext` from otolith back-calculation or mark–recapture size dispersion represents a promising empirical frontier now unlocked by `mizer`'s conservative advection–diffusion solver. Ongoing developments across the `sizespectrum` ecosystem—including automated catch-and-effort likelihood fitting, seasonal bioenergetics (`mizerSeasonal`), multi-resource coupling (`mizerMR`) and Ecopath interoperability (`mizerEcopath`)—build directly upon the mathematical and S3 extension foundations established here.

---

## Author Contributions

**G.W.D.**, **F.S.**, **J.L.B.** and **K.H.A.** conceived the `mizer` modelling framework and software architecture; **F.S.**, **K.H.A.** and **J.L.B.** developed the initial package implementation and multi-species calibration case studies; **G.W.D.** led the design and implementation of the `mizer` v2–v3 architecture, steady-state calibration workflow, advection–diffusion finite-volume and TR-BDF2 numerical solvers, Newton/Jacobian stability diagnostics and modular S3 extension system, with scientific and methodological contributions from **K.H.A.**, **J.L.B.** and **F.S.**; **G.W.D.** led the writing of the manuscript, and all authors contributed critically to revisions and approved the final version for publication.

## Acknowledgements

We gratefully acknowledge the contributions of the global `mizer` developer and user community—including Romain Forestier, Asta Audzijonyte, Michael A. Spence, Phoebe A. Woodworth-Jefcoats, Jonathan C. P. Reum, Max Lindmark, Camilla Novaglio, Samik Datta, Richard Law and Ryan Heneghan—whose empirical applications, bug reports, extension packages and workshops shaped the methods reported here. Development of `mizer` has been supported by the European Commission Horizon 2020 Research and Innovation Programme (MINOUW, Grant Agreement No. 634495), the European Union Horizon Europe Programme (NECCTON, Grant Agreement No. 101081273), the Australian Research Council Discovery Project (`DP170104240`), and UK Research and Innovation under the Sustainable Management of Marine Resources (SMMR) programme (Pyramids of Life, Grant No. `NE/V01708X/1`).

*Declaration of generative AI and AI-assisted technologies:* During recent software development of `mizer` v3 (including code refactoring, unit-test authoring and documentation generation) and the preparation of an initial structural draft and bibliography verification for this manuscript, the authors used large-language-model coding and writing assistants (Anthropic Claude Code and Google DeepMind Antigravity). All mathematical derivations, numerical implementations, empirical analyses, interpretations and manuscript text were critically reviewed, edited and verified by the human authors, who take full responsibility for the integrity and accuracy of the published work.

## Conflict of Interest Statement

The authors declare no conflict of interest.

## Data Availability Statement

`mizer` is free and open-source software licensed under GPL-3. Stable releases are distributed via the Comprehensive R Archive Network (<https://CRAN.R-project.org/package=mizer>), active development, issue tracking and full source code (including the `NS_params` North Sea dataset, pre-computed stability benchmarks and vignette scripts reproducing all figures) are hosted on GitHub (<https://github.com/sizespectrum/mizer>), and documentation is available at <https://sizespectrum.org/mizer/>. Upon acceptance, the exact release of `mizer` and reproduction scripts accompanying this article will be permanently archived on Zenodo with a persistent DOI (`10.5281/zenodo.XXXXXXX`).

## References
