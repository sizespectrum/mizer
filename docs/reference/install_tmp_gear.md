# Add a gear that exerts the fishing mortality being scanned

Copies the selectivity of the first gear catching the species onto a new
gear called `gear_name` with catchability 1, so that the effort of that
gear is the fishing mortality it exerts, and switches off the
catchability of the gears it replaces. The fishing on every other
species is untouched. The effort of the new gear is left at whatever
[`gear_params<-()`](https://sizespectrum.org/mizer/reference/gear_params.md)
gives it; the caller sets it to the value being scanned.

## Usage

``` r
install_tmp_gear(params, species, gear = NULL, gear_name)
```

## Arguments

- params:

  A MizerParams object.

- species:

  The target species.

- gear:

  The gear whose mortality is replaced, or NULL for all of the gears
  catching the species.

- gear_name:

  The name to give the new gear.

## Value

The MizerParams object with the extra gear.
