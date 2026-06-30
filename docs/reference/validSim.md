# Validate MizerSim object and upgrade if necessary

Checks that the given MizerSim object is valid and upgrades it if
necessary. It also validates the embedded
[`MizerParams-class()`](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
object with
[`validParams()`](https://sizespectrum.org/mizer/reference/validParams.md).
If any entries of the consumer abundance array `sim@n` are non-finite, a
warning is issued and the simulation is truncated at the last time step
where `sim@n` is still finite.

## Usage

``` r
validSim(sim)
```

## Arguments

- sim:

  The MizerSim object to validate

## Value

A valid MizerSim object

## Details

Occasionally, during the development of new features for mizer, the
[MizerSim](https://sizespectrum.org/mizer/reference/MizerSim-class.md)
class or the
[MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
class gains extra slots. MizerSim objects created in older versions of
mizer are then no longer valid in the new version because of the missing
slots. You need to upgrade them with this function.

This function adds the missing slots and fills them with default values.
It also calls
[`validParams()`](https://sizespectrum.org/mizer/reference/validParams.md)
to upgrade the MizerParams object inside the MizerSim object. Any object
from version 0.4 onwards can be upgraded.

## Backwards compatibility

The internal numerics in mizer have changed over time, so there may be
small discrepancies between the results obtained with the upgraded
object in the new version and the original object in the old version. If
it is important for you to reproduce the exact results then you should
install the version of mizer with which you obtained the results. You
can do this with

    remotes::install_github("sizespectrum/mizer", ref = "v0.2")

where you should replace "v0.2" with the version number you require. You
can see the list of available releases at
<https://github.com/sizespectrum/mizer/tags>.

If you only have a serialised version of the old object, for example
created via [`saveRDS()`](https://rdrr.io/r/base/readRDS.html), and you
get an error when trying to read it in with
[`readRDS()`](https://rdrr.io/r/base/readRDS.html) then unfortunately
you will need to install the old version of mizer first to read the
params object into your workspace, then switch to the current version
and then call
[`validParams()`](https://sizespectrum.org/mizer/reference/validParams.md).
You can then save the new version again with
[`saveParams()`](https://sizespectrum.org/mizer/reference/saveParams.md).
