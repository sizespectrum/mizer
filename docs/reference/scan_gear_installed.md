# Has a fishing-mortality scan already been installed in this model?

A gear of the right name proves nothing: it might be one the model
already had, and setting its effort would then leave the fishing the
scan is supposed to replace still switched on, so the scanned mortality
would be added to the existing mortality rather than replacing it. Nor
is it enough for the gear to *look* like the one
[`scanFishingMortality()`](https://sizespectrum.org/mizer/reference/scanEffort.md)
adds, for the same reason.

## Usage

``` r
scan_gear_installed(params, gear_name, species, gear = NULL)
```

## Arguments

- params:

  A MizerParams object.

- gear_name:

  The name of the gear to check.

- species:

  The target species.

- gear:

  The gear whose mortality the scan replaces, or NULL for all of the
  gears catching the species.

## Value

TRUE if this model already carries the installation.

## Details

What is checked is therefore the whole installation, exactly as
[`install_tmp_gear()`](https://sizespectrum.org/mizer/reference/install_tmp_gear.md)
leaves it: the name carried by exactly one row, which catches the target
species with catchability 1, **and** at least one original gear still
present with every gear it was supposed to replace switched off.
