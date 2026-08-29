## Upgrading from mizer 3.1 to 3.2

### `species_params<-()` now detects and protects changes

Previously, modifying species parameters via `species_params<-()` updated the values in the model but bypassed `given_species_params()`. This meant that your changes were not protected, and any subsequent recalculation of defaults (for example, by a call to `given_species_params<-()`) would overwrite your custom values. Furthermore, changing a parameter like `w_inf` via `species_params<-()` did not automatically trigger a recalculation of downstream parameters like `w_mat` or `w_max`.

Now, `species_params<-()` intelligently diffs the new data frame against the old one to detect exactly which parameters you have changed. It automatically records those changed parameters in `given_species_params`, protecting them from future overwrites, and immediately recalculates any downstream defaults based on your changes. 

**How this affects existing code:**

1. If your existing code used `species_params<-()` to update a core parameter like `w_inf` and you expected `w_mat` or `w_max` to remain frozen at their old values, you will now see them automatically recalculate. If you wish to freeze downstream parameters, you must provide their frozen values explicitly in the same update.

2. If your code computes custom parameters and saves them via `species_params<-()`, those parameters will now be preserved and survive future recalculations.


### Setting resource parameters

Two related changes affect how you modify the resource size spectrum. Together
they make the resource scalars behave like the species parameters: a scalar is
an input, and the size-dependent arrays are computed from it.

#### Assigning to `resource_params()` now updates the resource arrays

Previously, assigning to `resource_params()` — or to one of its components, such
as `resource_params(params)$kappa <- ...` — only stored the new scalar values.
The size-dependent carrying capacity (`cc_pp`) and replenishment rate (`rr_pp`)
were left unchanged until you next called `setResource()`.

Now these assignments immediately rebuild the arrays from the scalars, exactly as
`species_params()<-` rebuilds the species rates:

- `kappa`, `lambda` and `w_pp_cutoff` rebuild the carrying capacity;
- `r_pp` and `n` rebuild the replenishment rate.

Arrays that you have set by hand are left untouched (see *Frozen arrays* below).

If your code changed a resource scalar and then called `setResource()` to apply
it, nothing breaks — you can drop the now-redundant `setResource()` call. If you
changed a resource scalar and relied on the arrays *not* changing until later,
review that code.

#### Assigning to `resource_params()` does not balance the resource

*Balancing* means adjusting the rate and capacity together so that the resource
replenishes at exactly the rate at which it is consumed, keeping it at its steady
state. Assigning to `resource_params()` rebuilds the arrays from the scalars but
does **not** balance, so the resource steady state generally shifts.

Balancing is now solely a feature of `setResource()`. To change a resource
coefficient *and* keep the resource balanced, call `setResource()` rather than
assigning to `resource_params()`:

```r
# Rebuild the capacity from a new coefficient and rebalance the rate,
# so the steady state is preserved:
params <- setResource(params, resource_capacity = new_kappa)

# Likewise, set a new rate coefficient and rebalance the capacity:
params <- setResource(params, resource_rate = new_r_pp)
```

#### The resource setters gained a `balance` argument

`resource_rate<-`, `resource_capacity<-`, `resource_level<-` and
`resource_dynamics<-` still balance by default (unchanged behaviour), but they
now accept a `balance` argument so you can switch balancing off:

```r
# Set the capacity but leave the rate untouched (do not rebalance):
resource_capacity(params, balance = FALSE) <- my_capacity
```

#### Frozen arrays are protected from incidental balancing

When you set the size dependence of the resource capacity or the resource rate
by hand (by assigning a full vector rather than a scalar), mizer marks it "set
manually" — it is *frozen* and will not be recomputed from the resource
parameters. Previously, an operation that re-balanced the resource *without*
being given a replacement rate or capacity — for example changing only
`resource_dynamics`, or calling `setResource()` with neither a rate nor a
capacity — would silently overwrite such a frozen array. It is now kept, and a
warning is issued instead.

To deliberately recompute a frozen array from the resource parameters, pass
`reset = TRUE` to `setResource()`.

### The `species_params` data frame is now an S3 subclass

The `species_params` data frame now has class `c("species_params",
"data.frame")` (and `gear_params` similarly). It behaves like an ordinary data
frame, but subsetting and subassignment go through class-preserving S3 methods
and can trigger reactive re-validation and conversions (for example filling in a
weight from a length). Code that relied on `class(species_params(params))` being
exactly `"data.frame"`, or that stripped attributes with the assumption of a
plain data frame, may need adjusting. When you need a plain frame, coerce
explicitly with `as.data.frame()`.

### Accessing a column with `$` now returns a named vector

Extracting a single column from a `species_params` or `gear_params` object with
`$` now returns a vector named by species (or by `"species, gear"` for
`gear_params`):

```r
species_params(params)$w_mat
#>   Sprat  Herring      Cod
#>    ...      ...      ...
```

The values are unchanged, but the names are new. This is convenient for
identifying entries, but code that compared such a vector with `identical()` to
an unnamed vector, or that used it as-is where names matter (for example as
row/column names elsewhere), may behave differently. Strip the names with
`unname()` if you need the old behaviour. The `species` column itself is
returned unnamed.

### Setting `sel_func` adds the required argument columns

Assigning a selectivity function name to a `gear_params` object now
automatically adds the argument columns that the function needs (as `NA`),
ready to be filled in:

```r
gp$sel_func <- "sigmoid_length"
# gp now has l25 and l50 columns, both NA
```

Previously these columns had to be added by hand. Code that checks which columns
are present in `gear_params`, or that expected setting `sel_func` to leave the
column set unchanged, will now see the extra columns (#431).

### Passing a data frame to `species_params()` / `given_species_params()` now validates it

Calling `species_params()` or `given_species_params()` on a plain data frame now
runs the same validation and defaults that `validSpeciesParams()` and
`validGivenSpeciesParams()` apply, rather than only checking for misspellings and
converting lengths to weights. `species_params(df)` fills in the default columns
(`w_max`, `alpha`, `n`, `p`, `interaction_resource`, `z_ext`, and the rest), and
`given_species_params(df)` applies the consistency corrections (for example
clamping `w_mat` below `w_inf`), derives `w_inf` from `w_max`/`w_repro_max` when
it is absent, and now stops if the frame has duplicate species rows. Models built
or modified through `newMultispeciesParams()`, `setParams()` and the
`species_params()<-` / `given_species_params()<-` setters are unaffected, because
those already ran this validation. Only code that called the two accessors
directly on a bare data frame will see the extra columns and stricter checks
(#432).

### Printing of mizer array objects shows the values

`print()` on the array objects returned by the rate getters (`ArraySpeciesBySize`,
`ArrayTimeBySpecies`, `ArrayResourceBySize`, `ArrayTimeByResourceBySize` and
`ArrayTimeBySpeciesBySize`, as returned by `getEncounter()`, `getBiomass()`,
`getFMort()`, `NResource()` and similar) now truncates the output instead of
flooding the console with all the array entries. If your code or reports relied
on the old printed format, use `as.data.frame()` to go back to the full output.

### Upper boundary condition at `w_max`

The size-spectrum solver now holds the abundance at zero above each species'
maximum size `w_max`. Without diffusion this happens automatically and results
are unchanged. With diffusion switched on this change stops a small amount of
density leaking to sizes above `w_max`, so results there change slightly. See
`vignette("numerical_details")`.
