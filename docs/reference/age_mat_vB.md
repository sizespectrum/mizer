# Calculate age at maturity from von Bertalanffy growth parameters

This is not a good way to determine the age at maturity because the von
Bertalanffy growth curve is not reliable for larvae and juveniles.
However this was used in previous versions of mizer and is supplied for
backwards compatibility.

## Usage

``` r
age_mat_vB(object, ...)
```

## Arguments

- object:

  A MizerParams object or a species_params data frame

- ...:

  Currently unused.

## Value

A named vector. The names are the species names and the values are the
ages at maturity.

## Details

Uses the age at maturity that is implied by the von Bertalanffy growth
curve specified by the `w_inf`, `k_vb`, `t0`, `a` and `b` parameters in
the species_params data frame.

If any of `k_vb` is missing for a species, the function returns NA for
that species. Default values of `b = 3` and `t0 = 0` are used if these
are missing. If `w_inf` is missing, `w_max` is used instead.
