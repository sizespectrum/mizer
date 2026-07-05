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
  missing extension packages before registering the saved extension
  chain?

- sim:

  A MizerSim object

## Value

`saveParams()` and `saveSim()` return NULL invisibly. `readParams()`
returns a MizerParams object. `readSim()` returns a MizerSim object.

## Details

Issues a warning if the model you are saving relies on some custom
functions. Before saving a model you may want to set its metadata with
[`setMetadata()`](https://sizespectrum.org/mizer/reference/setMetadata.md).

## See also

"Using mizer extension packages":
[`vignette("using-extension-packages", package = "mizer")`](https://sizespectrum.org/mizer/articles/using-extension-packages.md)

## Examples

``` r
# Save params to a temporary file and read them back
tmp <- tempfile(fileext = ".rds")
saveParams(NS_params, file = tmp)
params <- readParams(tmp)
identical(params, NS_params)
#> [1] TRUE

# Save and read back a simulation
tmp2 <- tempfile(fileext = ".rds")
saveSim(NS_sim, file = tmp2)
sim <- readSim(tmp2)
identical(sim, NS_sim)
#> [1] TRUE
```
