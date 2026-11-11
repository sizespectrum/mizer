---
title: 'mizer: Dynamic multi-species size-spectrum modelling in R'
tags:
  - R
  - marine ecology
  - fisheries science
  - size spectrum
  - food web
  - ecosystem modelling
authors:
  - name: Gustav W. Delius
    orcid: 0000-0003-4092-8228
    corresponding: true
    affiliation: 1
  - name: Finlay Scott
    orcid: 0000-0001-9950-9023
    affiliation: 2
  - name: Julia L. Blanchard
    orcid: 0000-0003-0532-4824
    affiliation: 3
  - name: Ken H. Andersen
    orcid: 0000-0002-8478-3430
    affiliation: 4
affiliations:
  - name: Department of Mathematics, University of York, York, United Kingdom
    index: 1
    ror: 04m01e293
  - name: Oceanic Fisheries Programme, The Pacific Community, Nouméa, New Caledonia
    index: 2
    ror: 05ewdm369
  - name: Institute for Marine and Antarctic Studies and Centre for Marine Socioecology, University of Tasmania, Hobart, Australia
    index: 3
    ror: 01nfmeh72
  - name: National Institute of Aquatic Resources, Technical University of Denmark, Kongens Lyngby, Denmark
    index: 4
    ror: 04qtj9h94
date: 11 November 2026
bibliography: paper.bib
---

# Summary

In aquatic ecosystems, body size governs nearly every aspect of an organism's ecology: larger animals generally eat smaller ones, metabolic demands scale predictably with mass, and individual fish grow across five to seven orders of magnitude in weight from egg to adult. As a fish grows, both its prey and its predators shift continuously, so selective fishing or environmental change on any species or size class ripples across the entire food web.

`mizer` is an open-source R package for building, calibrating, simulating, and analysing dynamic size-structured models of aquatic communities. From species life-history traits and empirical observations such as survey biomasses, commercial catches, or growth curves, `mizer` simulates how individuals feed on smaller organisms, allocate assimilated energy between somatic growth and reproduction, and experience natural and fishing mortality through time. Ecologists, fisheries scientists, and students use `mizer` to explore mixed-fisheries management strategies, climate-warming scenarios, and biodiversity trade-offs via deterministic projections, stability analyses, and static or interactive visualisations (\autoref{fig:northsea}).

![Representative `mizer` outputs for the bundled 12-species North Sea model (`NS_params`). **(a)** Simulated biomass density across ten orders of magnitude in body weight, coupling the planktonic resource spectrum (bright green) to interacting fish spectra. **(b)** Emergent ontogenetic diet shift of Atlantic cod (*Gadus morhua*), progressing from plankton at larval sizes to forage fish, juvenile groundfish, and conspecifics at adult sizes.\label{fig:northsea}](figure1.png){ width=100% }

# Statement of need

Ecosystem-based fisheries management and climate-impact projection require models that mechanistically link individual bioenergetics, size-dependent predation ("life-history omnivory"), and multi-fleet size selectivity [@blanchard2014; @andersen2016; @andersen2019]. Solving the coupled integro-partial differential equations of multi-species size-spectrum theory [@datta2010; @hartvig2011] from scratch demands specialised numerical schemes on logarithmic grids, strict mass conservation across trophic transfers, and practical workflows for calibrating uncertain parameters against sparse empirical data.

`mizer` provides marine and freshwater ecologists, quantitative fisheries scientists, and theoretical biologists with a unified R framework spanning idealised community and trait-based models through to single-species and calibrated multi-species food webs. Whereas an early prototype (`mizer` v0.2) was introduced over a decade ago as an application note [@scott2014] offering only uncalibrated forward simulation of first-order advection equations, `mizer` v3 represents a generational redesign: it centres on a steady-state-first calibration workflow, generalises transport to advection--diffusion equations for stochastic individual growth [@datta2010], provides optional second-order conservative finite-volume quadratures and L-stable time-stepping, adds Newton equilibrium root-finding and Jacobian eigenvalue/limit-cycle stability analysis, and exposes a composable S3 extension architecture. This paper documents the modern `mizer` v3 software and establishes its standard citation.

# State of the field

Existing aquatic food-web and structured-population tools occupy distinct trade-offs among biological resolution, parameterisation effort, and speed:

- **Biomass-pool and end-to-end frameworks**: **Ecopath with Ecosim** [@christensen2004] and `Rpath` [@lucey2020] model trophic networks across species or coarse stanza pools for rapid mass-balance synthesis, whereas **Atlantis** [@audzijonyte2019] couples three-dimensional oceanography, biogeochemistry, ecology, and socio-economics in C at the cost of heavy parameterisation. Neither resolves continuous body size as the primary dynamical coordinate driving emergent growth and mortality.
- **Individual-based and statistical stock-ecosystem toolboxes**: **OSMOSE** [@shin2004] simulates size-opportunistic predation via spatial super-individuals, trading deterministic steady states and speed for stochastic replication, while **Gadget** / `gadget3` [@lentin2022] targets statistical likelihood fitting for single or small multi-stock assessments rather than whole-community size spectra.
- **Structured-population R packages**: `LeMaRns` [@spence2020lemarns] and **FishSUMS** [@speirs2016] implement discrete length- or weight-class shelf-sea models with prescribed or simplified growth, whereas `PSPManalysis` [@deroos2021] performs bifurcation analysis of general physiologically structured models for low-dimensional systems rather than multi-fleet fisheries calibration.

`mizer` was built as a dedicated framework rather than extending these packages because dynamic size-spectrum theory [@hartvig2011; @andersen2019] occupies an essential middle ground: allometric scaling keeps parameter counts parsimonious across dozens of species, while size-resolved bioenergetics let growth rates, maturation timing, diet shifts, and predation mortality respond endogenously to fishing and environmental forcing.

# Software design

`mizer` balances an accessible R calibration workflow, fast numerical integration across ten or more orders of magnitude in body mass, and modular extensibility.

## State-bearing parameter objects and steady-state calibration

Both model parameters and initial size-resolved abundances live inside a single S3 object (`MizerParams`), paired with `MizerSim` for trajectories and self-describing array classes (`ArraySpeciesBySize`, `ArrayTimeBySpecies`) that dispatch directly to `plot()` and interactive `plotly` helpers. Migrating from rigid S4 classes to S3 while preserving `@` accessors maintained backwards compatibility with existing scripts while allowing extension packages to subclass `MizerParams` via `NextMethod()`.

To combine allometric theory with incomplete data, `mizer` separates user-supplied traits (`given_species_params`) from defaults derived via metabolic scaling [@andersen2019], and marks manually set rate arrays as frozen so later parameter edits do not overwrite custom biology. Crucially, `mizer` decouples trophic calibration from recruitment density dependence: users first match observed biomasses, catches, and growth curves at a steady state under constant reproduction, and subsequently configure Beverton--Holt or custom reproduction rules governing stock resilience and maximum sustainable yield without perturbing that equilibrium.

## Mathematical formulation and hybrid R/FFT/C++ numerical engine

Consumer number densities $N_i(w, t)$ across body weight $w$ satisfy the conservative McKendrick--von Foerster advection--diffusion equation [@datta2010]:
\begin{equation}\label{eq:mvf}
\frac{\partial N_i}{\partial t} + \frac{\partial}{\partial w}\left(g_i(w,t)\,N_i - \frac{1}{2}\frac{\partial \bigl(d_i(w,t)\,N_i\bigr)}{\partial w}\right) = -\mu_i(w,t)\,N_i,
\end{equation}
coupled to size-structured background resources and recruitment boundary fluxes at egg weight $w_{\min,i}$, with emergent somatic growth $g_i$, growth-diffusion $d_i$, and mortality $\mu_i$.

On a geometric grid $w_{j+1}/w_j = \beta$, feeding kernels depending on the predator--prey mass ratio $w/w_p$ turn encounter, predation-mortality, and predation-diffusion integrals into discrete convolutions. Evaluating these via Fast Fourier Transforms lowers per-step complexity from $O(N_w^2)$ to $O(N_w \log N_w)$. Because diffusion imposes a stiff $\Delta t \sim \Delta w^2$ restriction on explicit schemes, `mizer` steps densities with a semi-implicit finite-volume discretisation solved in C++ (`Rcpp`) using an $O(N_w)$ tridiagonal Thomas solver that guarantees non-negativity.

To preserve reproducibility of published models while supporting high-accuracy work, `second_order_w` retains historical first-order upwind fluxes and node-sampled rates by default, while offering flux-limited Van Leer advection and exact or composite bin-averaged quadratures alongside an L-stable second-order `tr_bdf2` integrator. A Newton--Raphson solver (`steady(solver = "newton")`) also finds unstable equilibria, supporting Jacobian eigenvalue diagnostics, limit-cycle analysis, and parameter scans (`scanModel()`).

## Composable extensions and dual-purpose documentation

Four orthogonal mechanisms keep core `mizer` lean yet extensible: additive external encounter, mortality, and diffusion arrays; custom rate-function dispatch (`setRateFunction()`); coupled state components such as detritus or benthic pools (`setComponent()`); and S3 extension registration (`recordExtension()`). Workflow guides in `inst/skills/` serve as single-source files that compile into website vignettes while doubling (alongside `inst/llms.txt`) as structured domain instructions for AI coding assistants.

# Research impact statement

Continuously maintained on CRAN since 2014, `mizer` is core community infrastructure for size-based ecosystem modelling. An independent systematic review documented over 40 peer-reviewed `mizer` studies across marine and freshwater ecosystems [@pennino2026], including calibrated models of the North Sea [@blanchard2014; @spence2021], Baltic Sea [@lindmark2022], Mediterranean Sea [@dejuan2023], Eastern Bering Sea [@reum2020], China Seas [@szuwalski2017; @wo2024], central North Pacific [@woodworthjefcoats2019], temperate reefs [@audzijonyte2023], and freshwater lakes [@benoit2022]. These studies have shaped research on balanced harvesting [@jacobsen2014], trophic cascades [@szuwalski2017], micronutrient-sensitive fisheries [@robinson2022], recovery early-warning signals [@clements2019], and benthic--pelagic climate coupling [@audzijonyte2023].

Internationally, `mizer` is a core regional model in the Fisheries and Marine Ecosystem Model Intercomparison Project (Fish-MIP), informing IPCC and IPBES biodiversity and climate assessments [@tittensor2018; @ortegacisneros2025]. Its extension architecture supports more than ten companion packages—including `mizerExperimental`, `mizerMR`, `mizerShelf`, `therMizer`, `mizerEcopath`, `mizerStarvation`, and `mizerShiny`—alongside curated models in `mizerExamples` and a three-part online course (`mizerCourse`) used in postgraduate teaching and international workshops.

# AI usage disclosure

Generative AI coding assistants (including Anthropic Claude, OpenAI Codex, and Google Gemini via Antigravity) assisted `mizer` v3.x development through code refactoring, unit-test authoring, numerical verification scripts, and drafting documentation in `inst/skills/`. For this manuscript, Google Gemini (via Antigravity) was used to verify JOSS formatting rules, compile bibliographic metadata, compose \autoref{fig:northsea} from package reference plots, and generate an initial draft. All architectural decisions, mathematical formulations, numerical schemes, and scientific claims were conceived, reviewed, edited, and validated by the human authors alongside `mizer`'s automated test suite.

# Acknowledgements

We thank Richard Southwell, Samik Datta, Romain Forestier, Asta Audzijonyte, Max Lindmark, Phoebe Woodworth-Jefcoats, Jonathan Reum, Camilla Novaglio, Michael Spence, Axel Rossberg, Richard Law, and the wider `mizer` community for code, mathematical, and modelling contributions. Development was supported by the European Union's Horizon 2020 and Horizon Europe programmes (Grant Agreements No. 634495, MINOUW; No. 101081273, NECCTON), the Australian Research Council (DP170104240), and UK Research and Innovation's Sustainable Management of UK Marine Resources programme (Pyramids of Life, NE/V01708X/1).

# References
