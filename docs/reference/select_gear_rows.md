# The gear params rows whose fishing mortality is to be varied

The gear params rows whose fishing mortality is to be varied

## Usage

``` r
select_gear_rows(gp, species, gear = NULL)
```

## Arguments

- gp:

  The gear params data frame, with a character `gear` column.

- species:

  The target species.

- gear:

  The selected gear, or NULL for all gears catching the species.

## Value

An integer vector of row indices.
