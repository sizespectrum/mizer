---
name: understand-size-spectrum-dynamics
description: >-
  Understand how mizer models behave: which quantities you impose and which the
  model produces for itself, the food and predation feedback loops that couple
  species, what sets the slope of the steady-state spectrum and the timescale of
  its dynamics, and the two distinct kinds of density dependence. Use whenever
  reasoning about why a species collapses, explodes, or oscillates, why growth or
  mortality is not what you asked for, why changing one species moved another,
  why a model is insensitive to fishing, or what the Sheldon spectrum, feeding
  level, reproduction level (R_max) and trophic cascades mean in a mizer model.
---

# Understanding size-spectrum dynamics

Mizer models fish communities structured by **individual body size** rather than
age or fixed trophic levels. To diagnose model behaviour you have to reason from
size-structured principles, not from age-based or biomass-pool intuition. 

## What you impose vs. what the model produces

The single most useful habit is knowing which side of the line a quantity sits
on: what you set, and what the model works out for itself. This is the
distinction that explains most surprises.

**What you impose** is never an outcome. It is:

* **The predation process** — search volume $\gamma_i w^q$, the size-preference
  kernel (`beta`, `sigma`), the interaction matrix $\theta_{ij}$. Together these
  give a *rate of predation events per predator–prey size pair*: how fast a
  predator of size $w$ would eat prey of size $w_p$ if it met them.
* **Individual physiology** — maximum intake $h_i w^n$, metabolic cost
  $k_{s,i} w^p$, assimilation efficiency $\alpha$, and the split between growth
  and reproduction (`w_mat`, `w_repro_max`).
* **Everything outside predation among the modelled fish** — the resource
  (`kappa`, `lambda`, `w_pp_cutoff`, resource rate), external mortality `z0`,
  fishing effort and selectivity, `erepro` and `R_max`.

**What the model produces** is everything you would actually want to observe:
the feeding level $f_i(w)$, the growth rate $g_i(w)$ (including whether a fish
reaches `w_mat` at all), the predation mortality $\mu_p(w)$, the abundance
spectra $N_i(w)$ and their slope, the realised recruitment $R_i$ and reproduction
level $r_i$ — and the steady state itself.

### Growth and mortality are one process seen from two ends

The most important consequence, and the one most often missed: **growth and
predation mortality are not two mechanisms. They are the same predation events,
counted twice.**

When a predator of size $w$ eats a prey of size $w_p$, one event has two
bookkeeping consequences: the predator gains mass $w_p$, and the prey dies. So
mizer takes the *one* pairwise predation rate above and integrates it two
different ways:

| Emergent quantity | Sum the same predation rate over… | Weighted by | So it responds to |
|---|---|---|---|
| **Encounter → feeding level → growth** of species $i$ | all **prey** small enough for $i$ to eat | prey mass $w_p$ (mass gained) | the abundance of $i$'s **prey** — bottom-up |
| **Predation mortality** on species $i$ | all **predators** big enough to eat $i$ | nothing (deaths counted, not mass) | the abundance of $i$'s **predators** — top-down |

Same kernel $\phi$, same search volume $\gamma$, same interaction matrix
$\theta$; even the satiation factor is shared, since a satiated predator both
stops gaining mass and stops killing. The only differences are which variable is
integrated out — prey size or predator size — and whether you weight by the
prey's mass or just count the death.

This is why **you cannot make a species grow faster without making its prey die
faster**, and why the food and predation feedback loops below are not
independent dials. It also means a single diagnosis often covers both symptoms:
if growth is wrong, the mortality that same predation imposes is wrong too, and
the fix is in the shared parameters, not in one or the other.

Further consequences worth internalising:

* **You cannot set a growth rate.** `matchGrowth()` adjusts physiological
  coefficients until the *emergent* growth matches observation; it does not
  write a growth curve into the model. Same for mortality. See the
  `calibrate-model` skill.
* **Everything is coupled through abundance.** Changing one species' parameters
  changes its abundance, which changes the food available to its predators and
  the mortality on its prey. There are no isolated species.
* **A steady state is a fixed point of that coupling**, not a property you can
  set species by species. This is why `steady()` exists.

## The size spectrum as state variable

* **Number density $N_i(w)$**: the state variable for species $i$ is the density
  of individuals per unit weight; the number of fish in $[w, w + dw]$ is
  $N_i(w)\,dw$.
* **Biomass density**: $w N_i(w)$. On a logarithmic size grid, the biomass in a
  bin of constant logarithmic width is $w N_i(w)\,d(\ln w)$. Plotting choices
  matter: `plotSpectra(params, power = 2)` shows biomass per log size bin, which
  is the view in which the Sheldon spectrum is flat.
* **Sheldon's size spectrum**: across an aquatic ecosystem, biomass density is
  roughly uniform across logarithmic size octaves, so abundance scales as
  $N(w) \propto w^{-\lambda}$ with $\lambda \approx 2.05$. In mizer this is an
  *input*, not a prediction: `lambda` sets the resource spectrum the fish feed
  against. See "The community slope" below.
* **Community vs species spectra** — keep these apart, they behave differently
  and have different slopes: the *combined* community spectrum is a steady
  downward power law, but each *species* forms a **dome** spanning from its egg
  size `w_min` up to `w_max`. A species spectrum that looks like a
  power law over its whole range is a species that is not experiencing the
  maturity-driven growth slowdown — usually a sign something is wrong.

## Bookkeeping in a size bracket

All the dynamics come from one piece of accounting. Pick a size bracket
$[w, w + dw]$ and ask what changes the number of fish in it. Individuals do not
jump around the size axis: they enter only by growing in through the bottom
edge, and leave only by growing out through the top edge or by dying. So

> **rate of change = what grows in − what grows out − what dies**

That is the whole content of the McKendrick–von Foerster equation. The quantity
doing the growing-in and growing-out is the **flux** $J_i(w)$: the rate
(individuals per year) at which fish of species $i$ grow past size $w$. It is
just density times speed, $g_i(w) N_i(w)$ — how many fish are sitting at that
size, times how fast they are moving through it. The dying term is
$\mu_i(w) N_i(w)$.

Three things follow, and they cover most of what you need:

* **The traffic jam effect.** If more fish grow *in* at the bottom of a bracket
  than grow *out* at the top, the bracket fills up. That happens wherever growth
  **slows** — exactly like traffic on a motorway, where density rises where cars
  slow down, not where more cars are added. Where growth accelerates, density
  drops. **A bump in a spectrum is usually a growth slowdown, not a recruitment
  event**, and the dome shape of a species spectrum is mostly the growth
  slowdown at `w_repro_max` piling fish up before mortality removes them.
* **Nothing enters except at the egg size.** There is no immigration term, so
  the only way into the spectrum is through the bottom boundary, where the flux
  is set to the recruitment rate $R_i(t)$. Everything at every larger size got
  there by growing.
* **Steady state is not "nothing happening".** At a fixed point the three terms
  cancel at every size while fish are still growing and dying at full rate. A
  steady state can be maintained by fast growth against heavy mortality or by
  slow growth against light mortality — the balance is what is fixed, not the
  rates. This is why the *ratio* of mortality to growth, not either one alone,
  sets the shape (see below).

With growth diffusion switched on (off by default) individuals of the same size
grow at slightly different speeds, which adds a spreading term to the flux; the
bookkeeping is unchanged.

## The two feedback loops

Close the abundance dependence on the two integrals above — let the prey and
predator abundances themselves respond — and each becomes a feedback loop.
Nearly every dynamical surprise in mizer comes from one of the two. Both run
*through* body size, which is what makes them different from food-web models.

**The food loop (bottom-up, stabilising for the consumer).** More consumers →
prey depleted → encountered food $E_i(w)$ falls → feeding level $f_i(w)$ falls →
growth slows → fish spend longer at small sizes where mortality is high →
fewer consumers. Diagnose with `plotFeedingLevel()` and `resource_level()`.

**The predation loop (top-down, and the source of cascades).** More predators →
predation mortality $\mu_p$ on their prey sizes rises → prey abundance falls →
which relieves predation on whatever the prey were eating. Diagnose with
`plotPredMort()` and `plotDiet()`.

The loops act on **different sizes of the same species at the same time**. A fish
is prey early and predator late, so a change can travel up the spectrum through
one loop and back down through the other. That round trip is what produces
trophic cascades and, when the delay around it matches a generation time,
oscillations.

Two knobs decide how tightly the loops are coupled:

* **`beta`** (preferred predator/prey mass ratio, default **30**): a predator
  eats prey about 30× lighter — roughly 1.5 orders of magnitude down the
  spectrum. This sets how far apart in size two species have to be before one
  eats the other.
* **`sigma`** (kernel width in natural log, default **2**): wide. At $\pm 1$
  s.d. the kernel spans mass ratios from about 4 to 200, and it is cut off at
  ratio 1 (nothing eats prey larger than itself). A **narrow** kernel couples
  narrow size bands tightly and is a classic cause of limit cycles; a wide one
  averages over many prey sizes and damps oscillations.

## Individual energetics

Energy flows through an individual in three stages.

**1 — Encounter and feeding level.** A predator of size $w$ searches volume
$\gamma_i w^q$ and encounters prey according to the predation kernel and the interaction
matrix. Encountered food is turned into a dimensionless **feeding level**
$f_i(w) \in [0,1]$ by a Holling type II response — $f = 0$ is starvation, $f = 1$
is complete satiation. Mizer targets $f_0 \approx 0.6$ by default.

Feeding level is the first thing to look at in almost every diagnosis, because
it is dimensionless and has an absolute scale: you can tell at a glance whether
a number is wrong. $f < 0.2$ anywhere in the juvenile range means the species is
starving; $f > 0.9$ everywhere means the species is satiated and effectively
decoupled from its prey — its growth will not respond to food very much, and
neither will it exert much predation mortality (see below).

**2 — Net available energy.** Assimilated intake ($\alpha$, default 0.6) minus
standard metabolism, floored at zero. If it hits zero the fish has no surplus
for growth or reproduction and simply stops.

**3 — Allocation.** Below `w_mat` all net energy goes to somatic growth. Above
it, a fraction $\psi_i(w)$ goes to reproduction and the rest to growth, with
$\psi$ rising to 1 at `w_repro_max`, so growth halts there. This is what turns
the species spectrum into a dome: growth decelerating towards `w_repro_max` piles
individuals up (see the traffic-jam effect below) and then mortality removes
them.

## What sets the slope

Two different slopes get called "the slope", and they are not the same object.
The **species** slope is the one that follows from the bookkeeping, and it is
what this equation gives. The **community** slope, $-\lambda$, is a property of
the sum over all species and the resource; it is discussed after.

### The species slope

This is the bookkeeping above, made quantitative. For one species $i$, setting
"in − out − deaths" to zero and solving gives the local log-log slope of *that
species'* spectrum, ignoring diffusion:

$$\frac{d \ln N_i}{d \ln w} = -\frac{\mu_i(w)}{g_i(w)/w} - \frac{d \ln g_i}{d \ln w}$$

Read it as: **the slope is minus the ratio of mortality to mass-specific
growth**, minus how fast growth itself is accelerating. With allometric rates
$g \propto w^n$ the second term is just $n$ (the `n` species parameter;
`newMultispeciesParams()` sets 2/3 by default).

This one ratio explains most species-spectrum shapes. High mortality or stunted
growth ⇒ steep drop, few fish survive to large size. Fast growth relative to
mortality ⇒ shallow spectrum. Note that nothing here mentions $\lambda$: a
species' slope is set by its own mortality and growth, whatever the community
around it is doing.

### The community slope, and where $\lambda$ actually comes from

$\lambda$ is **not derived by the model**. It is a parameter you set (`lambda`,
default 2.05): the exponent of the resource spectrum the fish feed against. What
mizer does is not to predict $\lambda$ but to keep the model **self-consistent**
with it, in four steps:

1. **Assume** the prey a fish sees is distributed as $N(w) \propto w^{-\lambda}$.
2. **Choose $q$ to match it.** Mizer defaults `q` to $\lambda - 2 + n$. Against a
   $w^{-\lambda}$ background, a predator searching at rate $\gamma w^q$
   encounters food scaling as $w^{q + 2 - \lambda}$; that default is exactly the
   $q$ making it $w^n$ — the same power as maximum intake $h w^n$. Numerator and
   denominator of the feeding level then scale together, so **$f$ comes out
   independent of size**. (Verify on any fresh model: [`getFeedingLevel()`](../reference/getFeedingLevel.html) is flat
   at `f0` across the juvenile range.)
3. **Each species' slope therefore goes constant.** This is the equation above.
   With $f$ size-independent, $g \propto w^n$ and $\mu \propto w^{n-1}$, so
   $\mu/(g/w)$ is constant — a constant slope, i.e. a straight line on log-log
   axes over each species' juvenile range.
4. **The community reproduces the assumption.** Sum those species spectra —
   overlapping domes staggered by `w_max` — together with the resource, and you
   get back approximately $w^{-\lambda}$. The loop closes.

Step 4 is exact only in the idealised scaling model, where many species are
spread evenly over a wide range of `w_max`. In a model with a handful of real
species it is approximate: the resource carries the small-size end, and the
community spectrum is bumpy where species domes overlap unevenly. Do not expect
a fitted community slope to come out at exactly `lambda`.

The practical value of this is diagnostic. **The consistency is a chain, and
breaking any link shows up in the feeding level first.** Change `lambda` without
letting `q` follow, set `q` by hand, or deplete the resource away from its
assumed power law, and $f$ acquires a trend in size; species spectra then curve
instead of running straight, and the community spectrum stops being a power law.
So when a spectrum looks curved, check `plotFeedingLevel()` for a size trend
before touching anything structural.

## What sets the timescales

Size sets speed, and this is why age-based intuition misleads.

Mass-specific growth $g/w \propto w^{n-1} \approx w^{-1/3}$: a fish 1000× heavier
grows 10× more slowly in relative terms. Mortality scales the same way
($\mu \propto w^{n-1}$; mizer's default external mortality is
$z_0 = 0.6\,w_{\inf}^{\,n-1}$). So **small fish live fast and die young, large
fish live slow** — and the whole spectrum runs on a gradient of timescales rather
than a single one.

Time to grow from egg to `w_max` scales as $w_{\max}^{\,1-n} \approx
w_{\max}^{1/3}$. A species 1000× larger takes roughly 10× longer to mature.
Practical consequences:

* **Oscillation periods track generation time**, so a limit cycle in a large
  slow species has a long period and needs a correspondingly long `project()`
  run to even be visible.
* **Equilibration is set by the slowest species.** A `steady()` run that looks
  converged for small species may be nowhere near it for large ones.
* **Small species and juvenile size classes respond first** to any perturbation;
  the large-fish response arrives a generation later. A transient that looks
  like a cascade may just be different species running at different speeds.

## Density dependence comes in two kinds

This distinction decides how a model responds to fishing, and confusing the two
is a common source of wrong conclusions.

**1 — Imposed, via `R_max`.** Egg production $R_p$ is converted to realised
recruitment by a Beverton–Holt relation, summarised by the **reproduction
level** $r_i = R_i/R_{\max,i} \in [0,1]$ (`reproduction_level()`).

* **Low $r_i$ (< 0.5), strong density dependence**: recruitment is held well
  below capacity. Fishing that halves adult biomass barely moves recruitment —
  the stock is buffered. Very stable, and *insensitive to fishing* in a way that
  is imposed by assumption rather than earned by the dynamics.
* **High $r_i$ (> 0.85), weak density dependence**: recruitment tracks egg
  production almost proportionally. Susceptible to cohort blooms and persistent
  limit cycles.

**2 — Emergent, via the two feedback loops.** Even with $R_{\max} = \infty$
(mizer's default) a species is still regulated: more juveniles deplete the shared
resource (food loop) and attract predation (predation loop). This regulation is
*produced* by the model, is size-structured, and spills over onto other species.

Why it matters: both give a stable model, but they respond to fishing completely
differently. Imposed density dependence is species-local and unmoved by
ecosystem change; emergent density dependence weakens when you fish down the
predators supplying it. **A model tuned to stability with a high `R_max` is a
different ecosystem from one that is stable on its own.** If a species seems
inert to fishing, check `reproduction_level()` first — you may be looking at an
assumption, not a result.

## Size-structured trophic feedbacks

Because trophic role changes with size, mizer shows feedbacks absent from
classical food-web models:

* **Predation release with a starvation bottleneck**: fishing large predators
  lowers $\mu_p$ on small fish; the surge in small fish overconsumes the
  resource, feeding levels collapse, and the growth of the predator's *own*
  larvae is stunted. The cascade comes back around and bites the species that
  started it.
* **Cultivation–depensation**: large predators eat forage fish, but forage fish
  compete with — or prey on — the larvae of those same predators. Deplete the
  adults and the forage fish can suppress predator recruitment, holding the
  predator down even after fishing stops.
* **Travelling cohort waves**: a pulse of recruits eats its way up the size axis,
  leaving a moving depression in prey abundance and a wake of reduced mortality
  behind it. Note that the default numerical scheme damps these — see the
  `run-simulation` skill on numerical diffusion before concluding a wave is
  physically absent.

## Diagnostic reasoning heuristics

Check emergent properties before changing structural parameters.

| Symptom | Likely cause | What to inspect |
|---|---|---|
| **Species collapses during `project()`** | Starving larvae, intense juvenile predation, or too little egg production | `plotFeedingLevel()`, `plotDiet()`, `getPredMort()` |
| **Biomass oscillates in regular cycles** | Reproduction level near 1, narrow kernel (`sigma` too small), or knife-edge fishing | `reproduction_level()`, `getStability()` (see the `analyse-stability` skill) |
| **Growth slows before `w_mat`** | Food limitation at intermediate sizes; resource depleted or `w_pp_cutoff` too low | `plotGrowthCurves()`, `resource_level()`, `plotSpectra()` |
| **Growth is not what `matchGrowth()` asked for** | Growth is emergent — the food to support it is not there | `plotFeedingLevel()`, then the `calibrate-model` skill |
| **Tuning one species drops another** | Predation overlap or resource competition in shared juvenile size bins | `interaction_matrix()`, `plotDiet()` |
| **Species insensitive to fishing** | Strong imposed density dependence ($r_i$ low), or an effectively infinite resource | `reproduction_level()`, `resource_level()` |
| **Species spectrum curves instead of running straight** | `lambda`, `q` and `n` no longer mutually consistent, so feeding level is size-dependent | [`plotFeedingLevel()`](../reference/plotFeedingLevel.html) — is it flat in size? |
| **Community slope isn't `lambda`** | Expected only in the idealised scaling model; with few species the domes don't sum to a clean power law | [`getCommunitySlope()`](../reference/getCommunitySlope.html), [`plotSpectra()`](../reference/plotSpectra.html) — compare against species domes, not against `lambda` |
| **Species won't stay at steady state** | Fixed point is dynamically unstable, not a numerical failure | `getSteadyResidual()`, `steadyNewton()` |

<!-- agent-only -->
### Diagnostic procedure

Work outwards from the individual to the ecosystem; most problems resolve at
step 2 or 3.

1. **Is it emergent or imposed?** Decide which side of the table at the top of
   this skill the misbehaving quantity sits on. If emergent, the parameter you
   want to change is not the quantity you want to fix.
2. **Feeding level first.** `plotFeedingLevel(params)`. It is dimensionless with
   an absolute scale, so it is the fastest read. Flat near 0.6 is healthy; near 0
   is starvation; near 1 is decoupling; sloping in size means the
   `lambda`/`q`/`n` consistency is broken.
3. **Then growth.** `plotGrowthCurves(params)`. Does the fish reach `w_mat` in a
   plausible time? Growth failure downstream of a feeding-level problem is not a
   separate bug.
4. **Then mortality.** `plotPredMort(params)` and `getMort(params)`. Split
   predation, fishing and external contributions before blaming any one of them.
5. **Then reproduction.** `reproduction_level(params)`. Anything above ~0.85
   makes oscillations likely; anything below ~0.2 makes the species inert to
   fishing by assumption.
6. **Only then structural parameters** — `interaction_matrix()`, `beta`,
   `sigma`, resource settings.

Before concluding that dynamics are absent, check whether the default
first-order scheme is damping them (`second_order_w`; see the `run-simulation`
skill). Before concluding a state is unstable, check it is actually a fixed
point (`getSteadyResidual()`).
<!-- /agent-only -->

## Reference formulas

The equations behind the reasoning above, for when you need the exact form.

**Dynamics** — the McKendrick–von Foerster equation: the "in − out − deaths"
balance above, written as a continuity equation in size.

$$\frac{\partial N_i(w,t)}{\partial t} + \frac{\partial J_i(w,t)}{\partial w} = -\mu_i(w,t) N_i(w,t)$$

The derivative $\partial J_i/\partial w$ is the net loss by growth: flux out at
the top of the bracket minus flux in at the bottom. The right-hand side is the
loss by death. Growth flux $J_i = g_i N_i$, boundary condition
$J_i(w_{\min}, t) = R_i(t)$ at the egg size, $g_i$ in
$\text{g}\cdot\text{year}^{-1}$ and $\mu_i$ in $\text{year}^{-1}$. With growth
diffusion switched on (`D_ext > 0`, **off by default**) the flux gains a term
$-\tfrac{1}{2}\,\partial(D_i N_i)/\partial w$, which smooths sharp peaks and
spreads cohorts.

**The two integrals of the same predation rate.** Encountered food — integrate
over **prey** size $w_p$, weighting by the mass $w_p$ gained:

$$E_i(w) = \gamma_i(w) \int \Big(\theta_{ip} N_R(w_p) + \sum_j \theta_{ij} N_j(w_p)\Big)\, \phi_i(w/w_p)\, w_p\, dw_p$$

Predation mortality — integrate the same thing over **predator** size $w'$,
counting deaths rather than mass:

$$\mu_{p,i}(w) = \sum_j \int \theta_{ji}\,\big(1 - f_j(w')\big)\,\gamma_j (w')^q\,\phi_j(w'/w)\,N_j(w')\,dw'$$

The shared factors $\gamma$, $\phi$ and $\theta$ are what makes these two sides
of one coin. The $(1 - f_j)$ factor is worth noting: **satiated predators kill
less**, and a species whose predators are all near $f = 1$ experiences less
mortality than their abundance alone suggests. It applies to the growth side
too — realised intake is $f_i h_i w^n = (1 - f_i) E_i$ — so satiation throttles
both integrals together.

**Feeding level** from encountered food $E_i$ and maximum intake $h_i w^n$:

$$f_i(w) = \frac{E_i(w)}{E_i(w) + h_i w^n}$$

**Net energy** after assimilation ($\alpha$) and metabolism ($k_s w^p$, plus an
activity term $k w$ that defaults to zero):

$$E_{\text{net},i}(w) = \max\big(0,\; \alpha_i f_i(w) h_i w^n - k_{s,i} w^p\big)$$

**Growth**, after the reproductive allocation $\psi_i(w)$:

$$g_i(w) = \big(1 - \psi_i(w)\big) E_{\text{net},i}(w)$$

**Predation kernel** (lognormal, the default), for mass ratio
$\text{ppmr} = w/w_p \ge 1$ and zero below:

$$\phi_i(\text{ppmr}) = \exp\!\left[\frac{-\big(\ln(\text{ppmr}/\beta_i)\big)^2}{2\sigma_i^2}\right]$$

**Egg production** and Beverton–Holt density dependence:

$$R_{p,i} = \frac{\epsilon_i}{2 w_{\min}} \int \psi_i(w)\, E_{\text{net},i}(w)\, N_i(w)\, dw,
\qquad
R_i = \frac{R_{\max,i} R_{p,i}}{R_{\max,i} + R_{p,i}}$$

**Defaults worth remembering**: `beta` 30, `sigma` 2, `lambda` 2.05,
`w_pp_cutoff` 10 g, `alpha` 0.6, `f0` 0.6, `n` 2/3, `p` 0.7,
`q` $= \lambda - 2 + n$, `erepro` 1, `R_max` $\infty$, `w_min` 0.001 g,
`w_mat` $= w_{\inf}/4$, `z0` $= 0.6\,w_{\inf}^{\,n-1}$.

For the complete mathematical formulation see the
[model description](model_description.html) article.
