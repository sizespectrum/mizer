# Where a species parameter default belongs

Use this skill when adding a new species parameter default, or moving or removing
an existing one. The instinct to tidy all the defaults into one central place is
wrong here, and acting on it silently breaks things.

## Rule of thumb

A default lives with the **rate setter that reads the parameter**. Only
parameters that no single rate setter owns are defaulted centrally, in
`species_params.data.frame()` in `R/species_params.R`.

## Decide by consumption, not by tidiness

Before assuming central is cleaner, find out who reads the column:

```
grep -rnE 'species_params\$<name>|species_params\[\["<name>"\]\]' R/
```

(Single quotes matter: in double quotes the shell eats the backslash and `$`
becomes an end-of-line anchor, so the search silently finds nothing.)

- **Exactly one `setX()` reads it** → that setter owns the default.
  (`p`, `k` → `setMetabolicRate()`; `z0`, `z_ext`, `d` → `setExtMort()`;
  `E_ext` → `setExtEncounter()`; `interaction_resource` → `setInteraction()`;
  `q`, `gamma` → `setSearchVolume()`; `erepro`, `m`, `w_mat25`, `R_max` →
  `setReproduction()`.)
- **No single setter owns it** → central. This happens for five reasons: several
  setters read it (`n`), only the projection reads it (`alpha`), the size grid is
  built from it before any setter runs (`w_min`, `w_max`), the table's own
  length-weight conversion needs it (`a`, `b`), or it is reporting-only
  (`is_background`).

The current list of central parameters and the five reasons live in the
**"Where defaults live"** section of `vignettes/default_parameters.qmd`. Keep
that table current instead of duplicating it here.

A default that needs the whole model (the `w` grid, `resource_params$lambda`, the
predation kernel) or the setter's own arguments **cannot** be central:
`species_params.data.frame()` receives only a data frame, and runs inside
`emptyParams()` before the model exists.

## Traps

### A setter-local default is not redundancy

The rate setters are public API and get called directly, without `setParams()`
and hence without `validParams()`/`validSpeciesParams()`.
`check_and_convert_species_params()` does not re-apply defaults either, so
`params@species_params$z_ext <- NULL` leaves a hole that only `setExtMort()`
fills. The setter's default is the function guarding its own precondition.

Tests pin this on purpose — `test-setExtMort.R` NULLs `z_ext` and `d`,
`test-setMetabolicRate.R` NAs `p`, and each requires the setter to restore it.

Deleting such a default tends to fail **silently** rather than loudly:
`params@species_params$z_ext <- NULL` makes `NULL != 0` give `logical(0)`, so
`any()` is `FALSE` and the power-law mortality term is dropped with no error.

### Two homes go divergent, not merely untidy

`species_params.data.frame()` set `p = n` while `setMetabolicRate()` set
`p = 3/4`. The central one ran first and won, so the disagreement stayed latent —
until the central one was removed, at which point the stale `3/4` would silently
have become the effective default. **If you remove one of two homes, check what
the other one says before trusting it.**

### A setter-owned default has a second registration

`setter_owned_species_params()` in `R/species_params.R` lists every column that a
rate setter derives or defaults. `species_params<-()` uses it to decide which
columns to drop when it rebuilds the table from `given_species_params`, so that
`setParams()` can derive them afresh. **Add your parameter to that list under its
owning setter, or the value calculated last time is carried over and stops
responding to the parameters it is derived from.** Nothing fails loudly if you
forget: the model just quietly stops propagating that change.

Only put a column there if a setter really refills it on the plain
`setX(params)` path — check that the default is set *before* the setter's early
returns for a supplied or commented rate array, not after.

### And a third, if the parameter feeds a rate array

`frozen_rate_params()` in `R/info_signals.R` maps each rate array that can be
frozen to the species parameters that feed it. `species_params<-()` uses it to
warn the user when a change cannot take effect because the array it feeds was
set by hand (#489). A missing entry only costs a warning, so the list is
deliberately not exhaustive — but **never list a parameter that does not in
fact feed the array**, because that warns about a change that did take effect.

`.claude/skills/info-signals.md` has the rest of that rule, and covers reporting
in general.

### The central defaults are an exported contract

`validSpeciesParams()` is exported and its roxygen lists exactly which defaults
it fills in. Change `R/validSpeciesParams.R`, the vignette table and `NEWS.md`
together.

## `given_species_params`

It holds what the user supplied, **plus the defaults of any function argument
that sets a species parameter**, even when the user did not override the
argument. `newMultispeciesParams()` applies its `n` and `p` arguments to the data
frame before `emptyParams()` snapshots it (see `R/newMultispeciesParams.R`),
which is why `n` and `p` appear in `given_species_params` of a model whose user
supplied neither. That is intended, not a bug.

Defaults that are **not** function arguments must stay out; they belong in
`species_params` only.

This matters because `given_species_params<-` discards the whole species
parameter table and rebuilds it with `validSpeciesParams(value)` followed by
`setParams()`. Anything neither present in `given_species_params` nor
re-supplied by a setter is **lost on that round trip**. `w_min` is a known
violation of the rule — see issue #460.

For the full user-facing conceptual model of parameter tiers, `given_species_params` vs `species_params`, and downward propagation, consult `inst/skills/change-parameters/SKILL.md` and `inst/skills/build-model/SKILL.md`.

## Verifying a change to defaults

Moving a default between homes should be numerically inert on the standard
construction path, because `setParams()` calls every rate setter. Prove it with
a golden-model diff rather than assuming:

1. Before touching anything, build a model (e.g.
   `newMultispeciesParams(NS_species_params_gears, inter)`) and save its
   `species_params` columns and rate arrays (`search_vol`, `intake_max`, `metab`,
   `mu_b`, `psi`, `maturity`, `ext_encounter`, `ext_diffusion`, `initial_n`,
   `interaction`).
2. Make the change, rebuild, and assert no column was lost and no array changed.
   Only genuinely new columns may differ.

Watch for defaults that reference a column that may not exist: with no maximum
size given, `species_params(data.frame(species = "a"))` has no `w_mat`, so any
default derived from `w_mat` must sit inside the `if ("w_inf" %in% names(sp))`
block.
