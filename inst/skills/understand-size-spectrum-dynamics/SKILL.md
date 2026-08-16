---
name: understand-size-spectrum-dynamics
description: >-
  Understand the dynamics of mizer size-spectrum models, individual physiology,
  emergent growth and mortality, and ecosystem feedbacks. Use whenever reasoning
  about why a species collapses, explodes, or oscillates, how growth, mortality
  and feeding levels emerge from size-based predation, how energy is allocated
  between somatic growth and reproduction, why the steady-state spectrum has a
  specific slope (the McKendrick–von Foerster equation and Sheldon spectrum), or
  how reproduction level (R_max) and trophic cascades govern stability and
  indirect species interactions.
---

# Understanding size-spectrum dynamics

Mizer models fish communities structured by **individual body size** rather than
age or fixed trophic levels. Individual growth, predation mortality, and
reproduction are not prescribed static functions — they emerge dynamically from
predator–prey encounters across the size spectrum.

To diagnose model behaviour, tune steady states, or understand simulation
outcomes, you must reason from size-structured physiological and ecological
principles rather than classical age-based or biomass-pool intuition.

## The size spectrum as state variable

* **Number density $N_i(w)$**: The state variable for species $i$ is the density
  of individuals per unit weight.
  The actual number of fish in weight interval $[w, w + dw]$ is $N_i(w)\,dw$.
* **Biomass density**: $w N_i(w)$ is the biomass density. Plotted on a
  logarithmic size grid ($\log w$), the biomass in a bin of constant logarithmic
  width $d(\ln w) = dw / w$ is $w N_i(w)\,d(\ln w)$.
* **Sheldon's size spectrum**: Across an entire aquatic ecosystem, total biomass
  density is approximately uniform across logarithmic size octaves (the Sheldon
  hypothesis), which means total abundance scales roughly as $N(w) \propto w^{-\lambda}$
  with $\lambda \approx 2.05$.
* **Community vs species spectra**: While the combined community spectrum is a
  steady downward power law, each individual species forms a **dome-shaped**
  spectrum spanning from its egg size $w_{\min}$ up to a maximum size $w_{\max}$.

## The dynamic equation (McKendrick–von Foerster)

The time evolution of each species spectrum is governed by the continuity
equation with growth, diffusion, and mortality:

$$
\frac{\partial N_i(w,t)}{\partial t} + \frac{\partial J_i(w,t)}{\partial w} = -\mu_i(w,t) N_i(w,t)
$$

where $J_i(w,t)$ is the rate at which individuals of species $i$ grow from sizes
smaller than $w$ to sizes larger than $w$. This is known as the **flux** and is
given in terms of the growth and diffusion rates:

$$
J_i(w,t) = g_i(w,t) N_i(w,t) - \frac{1}{2} \frac{\partial \big(D_i(w,t) N_i(w,t)\big)}{\partial w},
$$
where $g_i(w)$ is somatic growth rate ($\text{g}\cdot\text{year}^{-1}$), and $\mu_i(w)$ is total mortality
rate ($\text{year}^{-1}$). The diffusion rate $D_i(w)$ 
 ($\text{g}^2\cdot\text{year}^{-1}$) captures the stochasticity in individual growth.

The boundary condition at the egg size $w_{\min}$ is

$$
J_i(w_{\min}, t) = R_i(t),
$$
where $R_i$ is the recruitment flux of eggs
($\text{individuals}\cdot\text{year}^{-1}$).

### Key dynamical properties

* **The Traffic Jam Effect**: The flux of fish through size is $J_i(w)$.
  Where growth rate slows down (or advective flux decreases), individuals pile up, increasing local density
  $N_i(w)$ (analogous to traffic slowing down on a motorway). Conversely, where
  growth accelerates, density drops. Growth diffusion acts to smooth out sharp peaks
  and spread cohorts over time.
* **Spectrum slope**: At steady state, without growth deceleration and with negligible diffusion:
  $$\frac{d \ln N_i}{d \ln w} \approx -\frac{\mu_i(w)}{g_i(w)/w} - 1$$
  The ratio of mortality $\mu(w)$ to mass-specific growth $g(w)/w$ dictates how
  steeply abundance drops with size. High mortality or stunted growth produces a
  steep drop (few surviving to large size).

## Individual energetics and physiology

Energy flows through individuals in sequential stages:

### Step 1 — Encounter & feeding level
A predator of size $w$ searches a volume $\gamma_i w^q$ and encounters prey of size $w_p$
according to its size preference kernel $\phi(w/w_p)$ (typically lognormal parameterised
by preferred ratio $\beta$ and width $\sigma$) and species interaction matrix $\theta_{ij}$.

Encountered food $E_i(w)$ is converted to a dimensionless **feeding level** $f_i(w) \in [0, 1]$
via a Holling type II functional response:

$$
f_i(w) = \frac{E_i(w)}{E_i(w) + h_i w^n}
$$

* $f_i(w) = 0$: complete starvation.
* $f_i(w) = 1$: absolute satiation (feeding rate limited only by maximum intake rate $h_i w^n$).
* Default mizer calibration targets a baseline feeding level $f_0 \approx 0.6$.

### Step 2 — Net available energy
Assimilated intake minus standard metabolic costs ($ks_{i} w^p$):

$$
E_{\text{net},i}(w) = \max\big(0, \; \alpha f_i(w) h_i w^n - k_{\text{metab},i} w^p\big)
$$

If $E_{\text{net}} \le 0$, the fish has zero surplus energy for growth or reproduction.

### Step 3 — Energy allocation (growth vs reproduction)
* **Immature ($w < w_{\text{mat}}$)**: 100% of net energy goes to somatic growth $g_i(w) = E_{\text{net},i}(w)$.
* **Mature ($w \ge w_{\text{mat}}$)**: A proportion $\psi_i(w)$ is allocated to
  reproduction, and the remainder $(1 - \psi_i(w))$ to somatic growth:
  $$g_i(w) = (1 - \psi_i(w)) E_{\text{net},i}(w)$$
  $\psi_i(w)$ rises smoothly from 0 near $w_{\text{mat}}$ to 1 at $w_{\text{repro\_max}}$,
  causing somatic growth to halt as fish reach $w_{\text{repro\_max}}$ (defaulting to $w_{\text{inf}}$).


## Emergent predation mortality

Predation mortality on prey of size $w$ is the sum of consumption rates by all
larger predators $w'$ for which $w$ is suitable prey:

$$
\mu_{p,i}(w) = \sum_j \int \theta_{ji} \, (1 - f_j(w')) \, \gamma_j (w')^q \, \phi(w'/w) \, N_j(w') \, dw'
$$

* Predation mortality $\mu_p(w)$ is **highest for larvae and juveniles** and drops
  rapidly as fish grow into size refuges.
* Satiated predators ($f \to 1$) exert less predation mortality per capita than
  hungry predators ($f \ll 1$).
* Total mortality is made up of predation mortality, fishing mortality and
  external mortality (mortality from sources not explicitly modelled).

## Stock-recruitment and density dependence

1. **Egg production**: Total rate of egg production from mature females:
   $$R_{p,i} = \frac{\epsilon}{2 w_{\min}} \int \psi_i(w) E_{\text{net},i}(w) N_i(w) \, dw$$
2. **Density-dependent recruitment**: Realized recruitment $R_i$ at egg size $w_{\min}$
   is mediated by Beverton–Holt density dependence:
   $$R_i = \frac{R_{\max,i} R_{p,i}}{R_{\max,i} + R_{p,i}}$$
3. **Reproduction level $r_i = R_i / R_{\max,i} \in [0, 1]$**:
   * **Low $r_i < 0.5$ (Strong density dependence)**: Realized recruitment is well below
     capacity. If fishing reduces adult biomass, density compensation buffers the
     loss ($R_i$ drops much less than $R_{p,i}$). Highly stable fixed points.
   * **High $r_i > 0.85$ (Weak density dependence)**: Realized recruitment is close
     to $R_{\max,i}$. Small changes in adult biomass do not affect recruitment at
     first, but populations are susceptible to **cohort blooms and persistent limit cycles**.

## Size-structured trophic feedbacks

Because trophic roles change with size, size-spectrum models exhibit indirect
feedbacks absent from classical food-web models:

* **Predation release vs starvation bottlenecks**: Fishing large predators reduces
  mortality $\mu_p$ on small fish. The resulting surge in small fish can overconsume
  the background plankton resource, causing feeding levels $f(w)$ to collapse and
  stunting the growth of the predator's own larvae.
* **Cultivation–depensation**: Large predators feed on forage fish, but forage fish
  compete with or prey upon the larvae of those same predators. Depleting adult
  predators can allow forage fish to suppress predator recruitment.
* **Travelling cohort waves**: A pulse of recruits consumes food as it grows along
  the size axis, creating a moving depression in prey abundance followed by a
  wake of reduced mortality.

## Diagnostic reasoning heuristics

When investigating model behaviour, check the emergent physiological and
dynamical properties before changing structural parameters:

| Symptom | Dynamical cause | What to inspect in R |
|---|---|---|
| **Species collapses during `project()`** | Starving larvae ($f(w_{\min}) \approx 0$), intense juvenile predation, or insufficient egg production ($R_p < R_{\text{crit}}$). | `plotFeedingLevel(params)`, `plotDiet(params)`, `getPredMort(params)` |
| **Biomass oscillates in regular cycles** | Reproduction level $r_i$ too close to 1 ($R_{\max}$ saturated), narrow predation kernel ($\sigma$ too small), or sharp knife-edge fishing. | `reproduction_level(params)`, `getStability(params)` (see the `analyse-stability` skill) |
| **Growth slows before $w_{\text{mat}}$** | Food limitation in intermediate sizes (resource spectrum depleted or $w_{\text{pp\_cutoff}}$ too low). | `plotGrowthCurves(params)`, `plotSpectra(params, power = 2)` |
| **Tuning one species drops another** | Direct predation overlap or intense resource competition in shared juvenile size bins. | `interaction_matrix(params)`, `plotDiet(params)` |
| **Species insensitive to fishing** | Strong compensatory reproduction ($r_i \ll 0.2$) or infinite plankton buffer. | `reproduction_level(params)` |
