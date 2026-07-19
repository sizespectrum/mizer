# Calculation of Default Parameter Values

## Introduction & Philosophy

A key design feature of the `mizer` package is that a multi-species
size-spectrum model can be set up using only a small amount of
information. This is possible because `mizer` uses **allometric scaling
relations** and **size-based feeding rules** to choose sensible default
values for parameters that are not explicitly provided by the user.

For example, to set up a model using the
[`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.md)
function, the user only *needs* to provide the von Bertalanffy
asymptotic sizes (\\w\_{\infty}\\) for the species in the community.
Almost all other physiological, ecological, and mortality parameters can
be automatically filled in using defaults. (See the section on
[maximum-size parameters](#maximum-size-parameters) below for the
distinction between \\w\_{\infty}\\, \\w\_{\text{repro_max}}\\ and
\\w\_{\text{max}}\\.)

This vignette details the philosophy, the relationships, and the
mathematical derivations behind the default parameter values calculated
by `mizer`. For a description of the model equations themselves, see
[The General Mizer Size-spectrum
Model](https://sizespectrum.org/mizer/articles/model_description.md).

------------------------------------------------------------------------

## How species parameters are stored and used

Because most species parameters can be filled in with defaults, `mizer`
needs to keep track of which values *you* supplied and which it
calculated for you. It also needs a rule for what should happen when you
later change a value. This section explains the machinery — the two
species-parameter data frames, the accessor and setter functions, the
role of
[`setParams()`](https://sizespectrum.org/mizer/reference/setParams.md),
and the “comment” mechanism on the size-dependent rate arrays — so that
the rest of this vignette makes sense.

### Two data frames: given vs. full

Every `MizerParams` object stores **two** species-parameter data frames:

- **`given_species_params`** holds only the values that were supplied
  explicitly (by you, or by a wrapper function such as
  [`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.md)).
  Parameters you never mentioned simply do not appear here, or appear as
  `NA`.
- **`species_params`** is the *full* table that the model actually uses.
  It contains your given values **plus** all the defaults and derived
  values that `mizer` filled in.

You access them with three functions:

``` r

given_species_params(params)       # only what was supplied explicitly
calculated_species_params(params)  # only the values mizer filled in
species_params(params)             # the complete table (given + calculated)
```

[`calculated_species_params()`](https://sizespectrum.org/mizer/reference/species_params.md)
is just
[`species_params()`](https://sizespectrum.org/mizer/reference/species_params.md)
with the given entries blanked out, so the given and calculated frames
together reconstruct the full table.

Here is a minimal model set up from nothing but a name and an asymptotic
size:

``` r

params <- newMultispeciesParams(
    data.frame(species = "Cod", w_inf = 5000)
)
# The columns that count as "given":
names(given_species_params(params))
```

    [1] "species" "w_inf"   "n"       "p"      

``` r

# Everything else was calculated, for example:
intersect(c("w_mat", "w_min", "h", "gamma", "ks"),
          names(calculated_species_params(params)))
```

    [1] "w_mat" "w_min" "h"     "gamma" "ks"   

Notice that `n` and `p` appear among the *given* parameters even though
we did not put them in the data frame. This is deliberate: some species
parameters can also be set through arguments to the constructor
functions (here the `n` and `p` arguments of
[`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.md)),
and any species parameter set from such an argument is recorded as given
— **even when the argument was left at its default value in the function
signature**. In other words, “given” means “fixed by the way you called
the constructor”, not only “typed into the species parameter data
frame”. The practical consequence is that these values count as explicit
input and so are not silently changed by later recalculations.

(Constructor arguments that do *not* directly correspond to a species
parameter behave differently. For example `z0pre` and `z0exp` are used
to *compute* the external mortality parameter `z0`, which therefore
appears among the calculated parameters rather than the given ones.)

### Changing parameters: `species_params<-()` vs. `given_species_params<-()`

There are correspondingly two ways to change species parameters:

- **`species_params<-()`** (often preferred in scripts) intelligently
  detects which values you have changed, records those changes as
  *given* values (protecting them from being overwritten by defaults),
  and then silently triggers a recalculation of everything that depends
  on them.
- **`given_species_params<-()`** is very useful in interactive use. Like
  `species_params<-()` it triggers a recalculation of other parameters,
  but importantly it will warn you when you set a parameter whose effect
  is overridden by another value you have already given.

For example, in a model where the age at maturity is known, `mizer`
derives the maximum intake coefficient `h` from it (see the [“Dependency
Chain”](#the-dependency-chain) below). Changing the size at maturity via
`species_params<-()` then propagates to `h` (and onward to `gamma` and
`ks`) automatically:

``` r

growing <- newMultispeciesParams(
    data.frame(species = "Cod", w_inf = 5000, age_mat = 5)
)
species_params(growing)$h
```

         Cod
    26.68043 

``` r

# using species_params<-() detects the change and protects it
species_params(growing)$w_mat <- 1000
species_params(growing)$h           # recalculated
```

         Cod
    26.68043 

After this, `w_mat` has moved from the calculated frame into the given
frame of `growing`, because `species_params<-()` detected that you
supplied it explicitly.

A useful consequence of the given/calculated split is that supplying a
value *switches off* the corresponding default calculation. The
[“Dependency Chain”](#the-dependency-chain) below shows, for instance,
that `age_mat` is only used to calculate a default for `h`. If you give
`h` directly, `age_mat` is ignored; if you give `gamma`, the `f0` value
is ignored; and so on.

### `setParams()` and the size-dependent rate arrays

The species parameters are not used directly during a simulation.
Instead they are used to compute the **size-dependent rate arrays**
(search volume, maximum intake rate, metabolic rate, external mortality,
reproduction, …) that are stored in their own slots of the `MizerParams`
object. The function that does this translation is
[`setParams()`](https://sizespectrum.org/mizer/reference/setParams.md),
which simply calls each of the individual setter functions
([`setSearchVolume()`](https://sizespectrum.org/mizer/reference/setSearchVolume.md),
[`setMaxIntakeRate()`](https://sizespectrum.org/mizer/reference/setMaxIntakeRate.md),
[`setMetabolicRate()`](https://sizespectrum.org/mizer/reference/setMetabolicRate.md),
and so on) in the correct order. Both `given_species_params<-()` and
`species_params<-()` call
[`setParams()`](https://sizespectrum.org/mizer/reference/setParams.md)
for you, which is why changing a species parameter automatically updates
the rates.

### Overriding a rate directly: the “comment” mechanism

Sometimes you do not want a rate to follow the allometric default at all
— you have your own array of values. Every setter function accepts the
rate array directly, e.g.

``` r

sv <- getSearchVolume(params)
sv[] <- 2 * sv                      # your own custom values
params <- setSearchVolume(params, search_vol = sv)
```

When you supply an array this way, `mizer` attaches a
[`comment()`](https://rdrr.io/r/base/comment.html) to the stored slot
(the text `"set manually"` unless you set your own comment):

``` r

comment(params@search_vol)
```

    [1] "set manually"

This comment is a **flag that protects your values**. While a rate array
carries a comment,
[`setParams()`](https://sizespectrum.org/mizer/reference/setParams.md)
will *not* recalculate it from the species parameters, and the species
parameters that would normally feed into it are effectively ignored for
that rate. So after the override above, changing `gamma` no longer
affects the search volume:

``` r

given_species_params(params)$gamma <- 1e-10
species_params(params)$gamma        # the species parameter did change ...
```

      Cod
    1e-10 

``` r

identical(comment(params@search_vol), "set manually")  # ... but the rate is frozen
```

    [1] TRUE

To go back to the automatically calculated values, call the setter with
`reset = TRUE`, which removes the comment and recomputes the rate from
the species parameters:

``` r

params <- setSearchVolume(params, reset = TRUE)
comment(params@search_vol)          # NULL again
```

    NULL

In summary:

- Provide a *species parameter* (via `given_species_params<-()`) when
  you want `mizer` to keep deriving the rate for you, using your value
  as input.
- Provide a *rate array* (via the relevant setter) when you want to
  override the rate completely; this sets a comment that freezes the
  rate until you `reset` it.

------------------------------------------------------------------------

## Maximum-size parameters

`mizer` has three species parameters that all describe a “maximum size”
but play distinct roles. Keeping them separate avoids a lot of
confusion.

- **`w_inf`** is the von Bertalanffy asymptotic size of an *average*
  individual. It is the **required** maximum-size parameter and the one
  from which the other two are derived by default. It is used to
  constrain the growth rate and to calculate the age at maturity.
- **`w_repro_max`** is the size at which a *typical mature* individual
  invests all of its available energy into reproduction (see
  \[setReproduction()\]). It is **not** a hard ceiling on size: not all
  individuals are mature at `w_repro_max`, and diffusion in the growth
  process lets some individuals grow beyond it. It defaults to `w_inf`
  but you can set it independently — there is no expectation that it is
  larger or smaller than `w_inf`, because estimates of the von
  Bertalanffy parameters are often unreliable.
- **`w_max`** is purely a *computational* boundary. It sets the upper
  end of the size grid and the range of the plots; it does not by itself
  stop growth. It defaults to `1.5 * w_inf`.

A typical (but not enforced) ordering is therefore \\w\_{\text{mat}} \<
w\_{\infty} \approx w\_{\text{repro_max}} \le w\_{\text{max}}\\.

For backwards compatibility, if you supply `w_max` (or `w_repro_max`)
but not `w_inf`, then `w_inf` is taken from those columns and an
informational message is issued, so models and scripts written for
earlier versions of `mizer` continue to work unchanged.

These parameters have been handled differently over the history of
`mizer`:

- In **mizer 1**, `w_inf` was the only maximum-size parameter. It was
  both the size at which all growth stopped and the size beyond which no
  fish existed, and was set to the von Bertalanffy asymptotic size.
- In **mizer 2**, `w_max` was introduced to acknowledge that real fish
  can be larger than the von Bertalanffy `w_inf`. Existing models were
  upgraded by renaming `w_inf` to `w_max`, and the use of von
  Bertalanffy curves was discouraged. Later in the mizer 2 series,
  `w_repro_max` was added (the size at which mature fish invest all
  energy into reproduction) once diffusion made it possible for fish to
  exceed the size at which the growth rate reaches zero, and `w_max`
  became just a computational boundary.
- In **mizer 3**, with diffusion fully integrated, von Bertalanffy
  growth is again the standard way to constrain growth, and `w_inf` is
  once more the first-class maximum-size parameter. Fish are allowed to
  exist beyond the von Bertalanffy asymptotic size because diffusion can
  take them there.

------------------------------------------------------------------------

## Species Parameter Defaults Reference Table

When creating a `MizerParams` object,
[`validSpeciesParams()`](https://sizespectrum.org/mizer/reference/validSpeciesParams.md)
(and subsequent setup functions like
[`setSearchVolume()`](https://sizespectrum.org/mizer/reference/setSearchVolume.md),
[`setMaxIntakeRate()`](https://sizespectrum.org/mizer/reference/setMaxIntakeRate.md),
[`setMetabolicRate()`](https://sizespectrum.org/mizer/reference/setMetabolicRate.md),
[`setExtMort()`](https://sizespectrum.org/mizer/reference/setExtMort.md),
and
[`setReproduction()`](https://sizespectrum.org/mizer/reference/setReproduction.md))
check the species parameter data frame and fill in missing parameters.

### Where defaults live

Each default has a single home: the rate-setting function that uses the
parameter.
[`setMetabolicRate()`](https://sizespectrum.org/mizer/reference/setMetabolicRate.md)
defaults `p` and `k` because it is the only function that reads them;
[`setExtMort()`](https://sizespectrum.org/mizer/reference/setExtMort.md)
defaults `z0`, `z_ext` and `d`; and so on. This matters because the rate
setters are also called directly, without going through
[`setParams()`](https://sizespectrum.org/mizer/reference/setParams.md),
so a setter has to be able to supply the defaults it needs on its own.

Only parameters that no single rate setter owns are defaulted centrally
in
[`validSpeciesParams()`](https://sizespectrum.org/mizer/reference/validSpeciesParams.md).
That happens for five reasons:

- **several setters read it** — `n` is used by
  [`setExtDiffusion()`](https://sizespectrum.org/mizer/reference/setExtDiffusion.md),
  [`setExtEncounter()`](https://sizespectrum.org/mizer/reference/setExtEncounter.md),
  [`setExtMort()`](https://sizespectrum.org/mizer/reference/setExtMort.md),
  [`setMaxIntakeRate()`](https://sizespectrum.org/mizer/reference/setMaxIntakeRate.md),
  [`setReproduction()`](https://sizespectrum.org/mizer/reference/setReproduction.md)
  and
  [`setSearchVolume()`](https://sizespectrum.org/mizer/reference/setSearchVolume.md);
- **only the projection reads it** — no setter touches `alpha`; it is
  used when calculating growth and reproduction;
- **the size grid is built from it** — `w_min` and `w_max` are needed
  before any rate setter runs;
- **the length-weight conversion needs it** — `a` and `b` are used by
  [`validSpeciesParams()`](https://sizespectrum.org/mizer/reference/validSpeciesParams.md)
  itself to convert lengths such as `l_mat` into weights;
- **it is only used for reporting** — `is_background`.

A consequence is that
[`validSpeciesParams()`](https://sizespectrum.org/mizer/reference/validSpeciesParams.md)
applied to a bare data frame returns only the central defaults. The
setter-owned columns appear once the model is built, because
[`setParams()`](https://sizespectrum.org/mizer/reference/setParams.md)
calls every rate-setting function.

Below is a reference of the default values and the functions responsible
for calculating them:

| Parameter | Description | Default Value / Formula | Function |
|:---|:---|:---|:---|
| `w_max` | Computational upper size boundary | `1.5 * w_inf` | [`validSpeciesParams()`](https://sizespectrum.org/mizer/reference/validSpeciesParams.md) |
| `w_repro_max` | Size at which a typical mature fish invests all energy in reproduction | `w_inf` | [`validSpeciesParams()`](https://sizespectrum.org/mizer/reference/validSpeciesParams.md) |
| `w_mat` | Size at maturity | `w_inf / 4` | [`validSpeciesParams()`](https://sizespectrum.org/mizer/reference/validSpeciesParams.md) |
| `w_min` | Egg size / birth size | `0.001` g | [`validSpeciesParams()`](https://sizespectrum.org/mizer/reference/validSpeciesParams.md) |
| `alpha` | Assimilation efficiency | `0.6` | [`validSpeciesParams()`](https://sizespectrum.org/mizer/reference/validSpeciesParams.md) |
| `n` | Allometric scaling exponent of intake | `3/4` | [`validSpeciesParams()`](https://sizespectrum.org/mizer/reference/validSpeciesParams.md) |
| `a` | Length-weight conversion coefficient | `0.01` | [`validSpeciesParams()`](https://sizespectrum.org/mizer/reference/validSpeciesParams.md) |
| `b` | Length-weight conversion exponent | `3` | [`validSpeciesParams()`](https://sizespectrum.org/mizer/reference/validSpeciesParams.md) |
| `is_background` | Flag indicating background species | `FALSE` | [`validSpeciesParams()`](https://sizespectrum.org/mizer/reference/validSpeciesParams.md) |
| `h` | Maximum intake rate coefficient | Derived from growth/age at maturity (or `30`) | [`get_h_default()`](https://sizespectrum.org/mizer/reference/get_h_default.md) |
| `gamma` | Search volume coefficient | Derived from target feeding level \\f_0\\ | [`get_gamma_default()`](https://sizespectrum.org/mizer/reference/get_gamma_default.md) |
| `q` | Exponent of the search volume | `lambda - 2 + n` | [`setSearchVolume()`](https://sizespectrum.org/mizer/reference/setSearchVolume.md) |
| `f0` | Target feeding level | `0.6` (or derived from `gamma` if given) | [`get_f0_default()`](https://sizespectrum.org/mizer/reference/get_f0_default.md) |
| `fc` | Critical feeding level | `0.2` | [`get_ks_default()`](https://sizespectrum.org/mizer/reference/get_ks_default.md) |
| `ks` | Standard metabolic rate coefficient | Derived from critical feeding level \\f_c\\ | [`get_ks_default()`](https://sizespectrum.org/mizer/reference/get_ks_default.md) |
| `p` | Allometric scaling exponent of metabolism | `n` | [`setMetabolicRate()`](https://sizespectrum.org/mizer/reference/setMetabolicRate.md) |
| `k` | Activity/movement cost coefficient | `0` | [`setMetabolicRate()`](https://sizespectrum.org/mizer/reference/setMetabolicRate.md) |
| `z0` | Constant external mortality rate | `z0pre * w_inf^z0exp` | [`setExtMort()`](https://sizespectrum.org/mizer/reference/setExtMort.md) |
| `z_ext` | Size-dependent external mortality | `0` | [`setExtMort()`](https://sizespectrum.org/mizer/reference/setExtMort.md) |
| `d` | Exponent of size-dependent external mortality | `n - 1` | [`setExtMort()`](https://sizespectrum.org/mizer/reference/setExtMort.md) |
| `E_ext` | External encounter rate | `0` | [`setExtEncounter()`](https://sizespectrum.org/mizer/reference/setExtEncounter.md) |
| `D_ext` | External diffusion rate | `0` | [`setExtDiffusion()`](https://sizespectrum.org/mizer/reference/setExtDiffusion.md) |
| `interaction_resource` | Species interaction strength with resource | `1` | [`setInteraction()`](https://sizespectrum.org/mizer/reference/setInteraction.md) |
| `erepro` | Reproduction efficiency | `1` | [`setReproduction()`](https://sizespectrum.org/mizer/reference/setReproduction.md) |
| `m` | Exponent of the investment into reproduction | `1` | [`setReproduction()`](https://sizespectrum.org/mizer/reference/setReproduction.md) |
| `w_mat25` | Size at which 25% of individuals are mature | `w_mat / 3^(1/10)` | [`setReproduction()`](https://sizespectrum.org/mizer/reference/setReproduction.md) |
| `beta` | Preferred predator/prey mass ratio | `30` | [`default_pred_kernel_params()`](https://sizespectrum.org/mizer/reference/default_pred_kernel_params.md) |
| `sigma` | Width of predation kernel | `2` | [`default_pred_kernel_params()`](https://sizespectrum.org/mizer/reference/default_pred_kernel_params.md) |

------------------------------------------------------------------------

## The Dependency Chain

Calculating the default physiological parameters follows a sequential
dependency chain. One parameter value determines the next, starting from
basic size parameters and ending with standard metabolic rates:

    [w_min, w_mat, age_mat] (Given or from VB parameters)
             │
             ▼
     ┌───────────────┐
     │       h       │  (Maximum intake rate coefficient)
     └───────┬───────┘
             ├──────────────────────────┐
             ▼                          ▼
     ┌───────────────┐          ┌───────────────┐
     │     gamma     │          │      ks       │  (Metabolic rate coefficient)
     └───────────────┘          └───────────────┘
    (Search volume coefficient)

1.  **Age at maturity** (\\\text{age}\_{\text{mat}}\\) is either
    supplied directly, or derived from von Bertalanffy parameters via
    [`age_mat_vB()`](https://sizespectrum.org/mizer/reference/age_mat_vB.md).
2.  **Maximum intake rate** (\\h\\) is calculated from
    \\\text{age}\_{\text{mat}}\\, egg size \\w\_{\text{min}}\\, maturity
    size \\w\_{\text{mat}}\\, and the assimilation parameters.
3.  **Search volume** (\\\gamma\\) and **metabolism** (\\k_s\\) are then
    both computed using the calculated value of \\h\\.

------------------------------------------------------------------------

## Mathematical Derivations

### 1. Maximum Intake Rate Coefficient (\\h\\)

The maximum intake rate coefficient \\h_i\\ for species \\i\\ is
calculated by
[`get_h_default()`](https://sizespectrum.org/mizer/reference/get_h_default.md).
The derivation is based on the idea that an average individual must grow
from egg size \\w\_{\text{min}}\\ to maturity size \\w\_{\text{mat}}\\
in a given time \\\text{age}\_{\text{mat}}\\ when feeding at a constant
feeding level \\f_0\\ (which defaults to \\0.6\\).

By definition, the rate of somatic growth is: \\\frac{dw}{dt} = g(w)\\

Prior to maturation, individuals invest nothing in reproduction
(\\\psi(w) = 0\\), so growth is equal to the energy available for growth
and reproduction \\E_r(w)\\: \\g(w) = \max(0, \alpha f(w) h w^n -
\text{metab}(w))\\

Under the default assumptions: \* The feeding level is constant: \\f(w)
= f_0\\. \* The activity cost is zero: \\k = 0\\, so \\\text{metab}(w) =
k_s w^p\\. \* The metabolic exponent equals the intake exponent: \\p =
n\\.

Substituting these in gives: \\\frac{dw}{dt} = \alpha f_0 h w^n - k_s
w^n = (\alpha f_0 h - k_s) w^n\\

We also define the critical feeding level \\f_c\\ as the feeding level
at which assimilation exactly balances metabolic costs, i.e., \\g(w) =
0\\. This implies: \\\alpha f_c h w^n = k_s w^n \implies k_s = \alpha
f_c h\\

Using this expression for \\k_s\\ in the growth equation yields:
\\\frac{dw}{dt} = (\alpha f_0 h - \alpha f_c h) w^n = \alpha h (f_0 -
f_c) w^n\\

We can solve this differential equation by separating variables and
integrating from birth (\\t = 0\\, \\w = w\_{\text{min}}\\) to maturity
(\\t = \text{age}\_{\text{mat}}\\, \\w = w\_{\text{mat}}\\):
\\\int\_{w\_{\text{min}}}^{w\_{\text{mat}}} w^{-n} \\ dw = \alpha h
(f_0 - f_c) \int\_{0}^{\text{age}\_{\text{mat}}} dt\\

Evaluating the integrals: \\\frac{w\_{\text{mat}}^{1-n} -
w\_{\text{min}}^{1-n}}{1-n} = \alpha h (f_0 - f_c)
\text{age}\_{\text{mat}}\\

Solving for \\h\\ gives the formula implemented in
[`get_h_default()`](https://sizespectrum.org/mizer/reference/get_h_default.md):
\\h = \frac{w\_{\text{mat}}^{1-n} -
w\_{\text{min}}^{1-n}}{\text{age}\_{\text{mat}} (1-n) \alpha (f_0 -
f_c)}\\

If the age at maturity is not known, `mizer` attempts to estimate it
from von Bertalanffy parameters using: \\\text{age}\_{\text{mat}} = -
\frac{\ln\left(1 -
\left(\frac{w\_{\text{mat}}}{w\_{\infty}}\right)^{1/b}\right)}{k\_{\text{vb}}} +
t_0\\ If no growth information is provided at all, `mizer` defaults to
\\h = 30\\.

------------------------------------------------------------------------

### 2. Search Volume Coefficient (\\\gamma\\)

The search volume coefficient \\\gamma_i\\ is calculated by
[`get_gamma_default()`](https://sizespectrum.org/mizer/reference/get_gamma_default.md).
It is selected such that if the prey abundance were described by a
power-law resource spectrum: \\N_R(w_p) = \kappa w_p^{-\lambda}\\ then
the resulting encounter rate would lead to the target feeding level
\\f_0\\ (default \\0.6\\).

The encounter rate \\E_i(w)\\ of predator species \\i\\ at size \\w\\
is: \\E_i(w) = \gamma_i w^q \int \theta\_{iR} N_R(w_p) \phi_i(w, w_p)
w_p \\ dw_p\\

Substituting the power-law resource spectrum: \\E_i(w) = \gamma_i w^q
\theta\_{iR} \kappa \int w_p^{1-\lambda} \phi_i(w, w_p) \\ dw_p\\

If the predation kernel depends only on the predator/prey mass ratio \\x
= w/w_p\\, then we can substitute \\w_p = w/x\\ and \\dw_p = -w/x^2 \\
dx\\: \\\int w_p^{1-\lambda} \phi_i(w/w_p) \\ dw_p = w^{2-\lambda}
\int_0^\infty x^{\lambda - 2} \phi_i(x) \\ dx\\

Thus, the encounter rate scales with size as: \\E_i(w) = \gamma_i w^{q +
2 - \lambda} \theta\_{iR} \kappa \int_0^\infty x^{\lambda - 2} \phi_i(x)
\\ dx\\

Let us define the size-independent component of the encounter rate as:
\\A\_{\text{avail}, i} = \theta\_{iR} \kappa \int_0^\infty x^{\lambda -
2} \phi_i(x) \\ dx\\ so that \\E_i(w) = \gamma_i w^{q + 2 - \lambda}
A\_{\text{avail}, i}\\.

The feeding level is given by: \\f_i(w) = \frac{E_i(w)}{E_i(w) + h_i
w^n}\\

For the feeding level to be independent of size and equal to \\f_0\\, we
must require the size scaling exponents of encounter and maximum intake
to match: \\q + 2 - \lambda = n \implies q = n + \lambda - 2\\ When this
holds, the size \\w\\ cancels out: \\f_0 = \frac{\gamma_i
A\_{\text{avail}, i}}{\gamma_i A\_{\text{avail}, i} + h_i}\\

Solving this equation for \\\gamma_i\\ yields: \\\gamma_i =
\frac{h_i}{A\_{\text{avail}, i}} \frac{f_0}{1 - f_0}\\

In the code, \\A\_{\text{avail}, i}\\ is calculated numerically by
setting \\\gamma_i = 1\\ and evaluating the encounter rate at the
maximum size grid point divided by the appropriate power of \\w\\:
\\\text{avail_energy}\_i =
\frac{E\_{i}\|\_{\gamma_i=1}(w\_{\text{max}})}{w\_{\text{max}}^{q + 2 -
\lambda}}\\ Then, the default \\\gamma\\ is: \\\gamma_i =
\frac{h_i}{\text{avail_energy}\_i} \frac{f_0}{1 - f_0}\\

------------------------------------------------------------------------

### 3. Standard Metabolic Rate Coefficient (\\ks\\)

The standard metabolic rate coefficient \\k\_{s,i}\\ is calculated by
[`get_ks_default()`](https://sizespectrum.org/mizer/reference/get_ks_default.md).
It is chosen so that the critical feeding level \\f_c\\ (default
\\0.2\\) is exactly the feeding level required for assimilation to
balance metabolism at maturity size \\w\_{\text{mat}}\\.

At the critical feeding level \\f_c\\, the assimilation rate equals
metabolic losses. At size \\w\_{\text{mat}}\\: \\\alpha_i f_c h_i
w\_{\text{mat}}^n = \text{metab}\_i(w\_{\text{mat}})\\

Using the default metabolic cost equation \\\text{metab}\_i(w) =
k\_{s,i} w^p\\ (with activity \\k = 0\\): \\\alpha_i f_c h_i
w\_{\text{mat}}^n = k\_{s,i} w\_{\text{mat}}^p\\

Solving for \\k\_{s,i}\\ gives the formula: \\k\_{s,i} = f_c \alpha_i
h_i w\_{\text{mat}}^{n - p}\\

If \\n = p\\, this simplifies to: \\k\_{s,i} = f_c \alpha_i h_i\\

------------------------------------------------------------------------

### 4. Target Feeding Level (\\f_0\\)

For models where the search volume coefficient \\\gamma_i\\ is provided
by the user but the target feeding level \\f_0\\ is not,
[`get_f0_default()`](https://sizespectrum.org/mizer/reference/get_f0_default.md)
calculates the inverse of the search volume default.

Using the relationship derived above: \\f\_{0,i} = \frac{\gamma_i
A\_{\text{avail}, i}}{\gamma_i A\_{\text{avail}, i} + h_i}\\ Since the
encounter rate calculated with the given \\\gamma_i\\ is \\E_i(w) =
\gamma_i w^{q+2-\lambda} A\_{\text{avail}, i}\\, we evaluate:
\\\text{avail_energy}\_i =
\frac{E_i(w\_{\text{max}})}{w\_{\text{max}}^{q+2-\lambda}}\\ And then
calculate: \\f\_{0,i} = \frac{1}{h_i / \text{avail_energy}\_i + 1}\\

------------------------------------------------------------------------

## Other Defaults

- **External Mortality (\\z_0\\)**: Background mortality that is not due
  to predation or fishing. If not specified, it scales allometrically
  with the asymptotic size: \\z\_{0, i} = z\_{0\text{pre}} w\_{\infty,
  i}^{z\_{0\text{exp}}}\\ where \\z\_{0\text{pre}}\\ defaults to \\0.6\\
  and \\z\_{0\text{exp}}\\ defaults to \\n - 1\\ (typically \\-1/4\\).
- **Predation Kernel**: By default, `mizer` uses a lognormal predation
  kernel \\\phi_i(w, w_p)\\ with preferred predator-prey mass ratio
  \\\beta = 30\\ and width \\\sigma = 2\\.
- **Maturity Ogive**: The proportion of individuals mature at size \\w\\
  is a sigmoidal function: \\\text{maturity}(w) = \frac{1}{1 +
  (w/w\_{\text{mat}})^{-10}}\\

------------------------------------------------------------------------

## Defaults Editions

Sometimes improved methods for choosing defaults are added to `mizer`.
To prevent older code from changing behaviour unexpectedly when the
package is updated, `mizer` uses a versioning system controlled by the
[`defaults_edition()`](https://sizespectrum.org/mizer/reference/defaults_edition.md)
function.

- **Edition 1** (Current default): Historically chosen default
  relations.
- **Edition 2** (Available): Catchability is set to \\0.3\\ (instead of
  \\1\\), initial fishing effort is \\1\\ (instead of \\0\\), initial
  abundance is calculated using the steady-state solver instead of a
  simple power law, and `psi` is determined strictly by the maturity
  ogive and the reproductive investment proportion.

------------------------------------------------------------------------

## Worked Example

The following code illustrates how `mizer` automatically calculates
defaults when setting up a new model with only the species name and
asymptotic size:

``` r

# Create species parameter data frame with only w_inf specified
simple_species <- data.frame(
  species = c("Species A", "Species B"),
  w_inf = c(1000, 8000)
)

# Initialize parameters
params <- newMultispeciesParams(simple_species)
```

    For species where no growth information is available the parameter h has been set to h = 30.
    Because the age at maturity is not known, I need to fall back to using
    von Bertalanffy parameters, where available, and this is not reliable.
    No ks column so calculating from critical feeding level.
    Using z0 = z0pre * w_inf ^ z0exp for missing z0 values.
    Using f0, h, lambda, kappa and the predation kernel to calculate gamma.

``` r

# Inspect the automatically filled species parameters
species_params(params)[, c("species", "w_mat", "w_min", "h", "gamma", "ks")]
```

    An object of class "species_params" containing parameters for 2 species:
       species w_mat  h       ks w_min        gamma
     Species A   250 30 2.994823 0.001 7.214272e-11
     Species B  2000 30 2.794269 0.001 7.214272e-11
