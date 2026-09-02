# Save and restore mizer objects

`saveParams()` saves a MizerParams object to a file. This can then be
restored with `readParams()`. `saveSim()` and `readSim()` provide the
same lifecycle for MizerSim objects.

## Usage

``` r
saveParams(params, file)

readParams(file, install_extensions = FALSE)

saveSim(sim, file)

readSim(file, install_extensions = FALSE)
```

## Arguments

- params:

  A MizerParams object

- file:

  The name of the file or a connection where the object is saved to or
  read from.

- install_extensions:

  Logical. Should `readParams()` or `readSim()` attempt to install
  missing extension packages before loading them?

- sim:

  A MizerSim object

## Value

`saveParams()` and `saveSim()` return NULL invisibly. `readParams()`
returns a MizerParams object. `readSim()` returns a MizerSim object.

## Details

While these functions ultimately use
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) and
[`readRDS()`](https://rdrr.io/r/base/readRDS.html), they do extra work
to make the saved file more robust and more portable, so you should
always prefer them over calling
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html)/[`readRDS()`](https://rdrr.io/r/base/readRDS.html)
directly on a mizer object.

The complete S3 class vector is stored with the object, including any
extension classes.

## What `saveParams()` and `saveSim()` do beyond [`saveRDS()`](https://rdrr.io/r/base/readRDS.html)

- They **validate** the object before writing it, so a corrupted or
  inconsistent object is caught at save time rather than when you next
  try to use it.

- They **check that the required extension packages are installed** and
  stop with an informative error if they are not, so you do not save a
  file that you would be unable to read back.

- They **warn if the model relies on custom functions** (custom rate,
  dynamics, selectivity or predation-kernel functions that are not
  provided by mizer or a registered extension package). Such functions
  are not stored in the file, so to share the model you also need to
  share an R script or R Markdown file defining them.

Before saving a model you may want to set its metadata with
[`setMetadata()`](https://sizespectrum.org/mizer/reference/setMetadata.md).

## What `readParams()` and `readSim()` do beyond [`readRDS()`](https://rdrr.io/r/base/readRDS.html)

- They **upgrade** an object saved by an older version of mizer to the
  current structure (see
  [`upgradeParams()`](https://sizespectrum.org/mizer/reference/upgradeParams.md)),
  so that models saved years ago still load correctly.

- They **load the extension packages** that the model needs and,
  optionally, install any that are missing (see `install_extensions`),
  before revalidating the object.

- `readParams()` **reconciles the species parameters**, see
  [`reconcileSpeciesParams()`](https://sizespectrum.org/mizer/reference/reconcileSpeciesParams.md).
  A model may hold species parameter values that were written straight
  into the `species_params` slot and that mizer would therefore undo at
  the next recalculation. Those values are recorded as given species
  parameters instead, so that they survive.

## See also

The [guide to using mizer extension
packages](https://sizespectrum.org/mizer/articles/guide-use-extension-packages.html)

## Examples

``` r
# Save params to a temporary file and read them back
tmp <- tempfile(fileext = ".rds")
saveParams(NS_params, file = tmp)
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
params <- readParams(tmp)
identical(params, NS_params)
#> [1] FALSE

# Save and read back a simulation
tmp2 <- tempfile(fileext = ".rds")
saveSim(NS_sim, file = tmp2)
sim <- readSim(tmp2)
identical(sim, NS_sim)
#> [1] TRUE
```
