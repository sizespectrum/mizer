# Get density-dependent reproduction rate during projection

S3 generic used by extension-aware projections to calculate the
density-dependent reproduction rate. The base method calls the selected
density-dependence function in `params@rates_funcs$RDD`.

## Usage

``` r
projectRDD(params, rdi, species_params = params@species_params, t = 0, ...)

# S3 method for class 'MizerParams'
projectRDD(params, rdi, species_params = params@species_params, t = 0, ...)
```

## Arguments

- params:

  A MizerParams object.

- rdi:

  Vector of density-independent reproduction rates \\R\_{di}\\ for all
  species.

- species_params:

  A species parameter dataframe. Must contain a column `R_max` holding
  the maximum reproduction rate \\R\_{max}\\ for each species.

- t:

  The time for which to do the calculation.

- ...:

  Unused

## Value

Vector of density-dependent reproduction rates.
