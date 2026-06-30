# Check validity of gear parameters and set defaults

The function returns a valid gear parameter data frame that can be used
by
[`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md)
or it gives an error message.

## Usage

``` r
validGearParams(gear_params, species_params)
```

## Arguments

- gear_params:

  Gear parameter data frame

- species_params:

  Species parameter data frame

## Value

A valid gear parameter data frame

## Details

The gear_params data frame is allowed to have zero rows, but if it has
rows, then the following requirements apply:

- There must be columns `species` and `gear` and any species - gear pair
  is allowed to appear at most once. Any species that appears must also
  appear in the `species_params` data frame.

- There must be a `sel_func` column. If a selectivity function is not
  supplied, it will be set to "knife_edge".

- There must be a `catchability` column. If a catchability is not
  supplied, it will be set to 1.

- All the parameters required by the selectivity functions must be
  provided.

If gear_params is empty, then this function tries to find the necessary
information in the species_params data frame. This restricts each
species to be fished by only one gear. Defaults are used for information
that can not be found in the species_params dataframe, as follows:

- If there is no `gear` column or it is NA then a new gear named after
  the species is introduced.

- If there is no `sel_func` column or it is NA then `knife_edge` is
  used.

- If there is no `catchability` column or it is NA then this is set to
  1.

- If the selectivity function is `knife_edge` and no `knife_edge_size`
  is provided, it is set to `w_mat`.

The row names of the returned data frame are of the form "species,
gear".

When `gear_params` is `NULL` and there is no gear information in
`species_params`, then a gear called `knife_edge_gear` is set up with a
`knife_edge` selectivity for each species and a `knive_edge_size` equal
to `w_mat`. Catchability is set to 0.3 for all species.

## See also

[`gear_params()`](https://sizespectrum.org/mizer/reference/gear_params.md)
