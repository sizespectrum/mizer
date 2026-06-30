# Validate MizerParams object and upgrade if necessary

Checks that the given MizerParams object is valid and upgrades it if
necessary.

## Usage

``` r
validParams(params, info_level = 3)
```

## Arguments

- params:

  The MizerParams object to validate

- info_level:

  Controls the amount of information messages and warnings that are
  shown. Higher levels lead to more messages.

## Value

A valid MizerParams object

## Details

It is possible to render a MizerParams object invalid by manually
changing its slots. This function checks that the object is valid and if
not it attempts to upgrade it to a valid object or gives an error
message. If the object is valid then it is returned unchanged. The
function reports an error if any of the rate arrays contain any
non-finite numbers (except for the maximum intake rate that is allowed
to be infinite).

Occasionally, during the development of new features for mizer, the
[MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
object gains extra slots. MizerParams objects created in older versions
of mizer are then no longer valid in the new version because of the
missing slots. You need to upgrade them with this function. It adds the
missing slots and fills them with default values. Any object from
version 0.4 onwards can be upgraded. Any old
[MizerSim](https://sizespectrum.org/mizer/reference/MizerSim-class.md)
objects should be similarly updated with
[`validSim()`](https://sizespectrum.org/mizer/reference/validSim.md).

This function uses
[`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.md)
to create a new MizerParams object using the parameters extracted from
the old MizerParams object.

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
and then call `validParams()`. You can then save the new version again
with
[`saveParams()`](https://sizespectrum.org/mizer/reference/saveParams.md).
