# Length-weight conversion

For each species, convert between length and weight using the
relationship \$\$w_i = a_i l_i^{b_i}\$\$ or \$\$l_i = (w_i /
a_i)^{1/b_i}\$\$ where `a` and `b` are taken from the species parameter
data frame and \\i\\ is the species index.

## Usage

``` r
l2w(l, species_params)

w2l(w, species_params)
```

## Arguments

- l:

  Lengths in cm. Either a single number used for all species or a vector
  with one number for each species.

- species_params:

  A species parameter data frame or a MizerParams object.

- w:

  Weights in grams. Either a single number used for all species or a
  vector with one number for each species.

## Value

A vector with one entry for each species. `l2w()` returns a vector of
weights in grams and `w2l()` returns a vector of lengths in cm.

## Details

This is useful for converting a length-based species parameter to a
weight-based species parameter.

If any `a` or `b` parameters are missing the default values `a = 0.01`
and `b = 3` are used for the missing values.
