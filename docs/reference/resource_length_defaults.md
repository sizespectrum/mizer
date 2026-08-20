# Default weight-length parameters for the resource

The resource is a composite of everything from bacteria to
macrozooplankton, so it has no taxonomic length-weight relationship. The
default is the geometric one that plankton ecology uses instead: the
**equivalent spherical diameter** of an organism with the density of
water, \$\$w = \frac{\pi}{6} l^3,\$\$ with \\w\\ in grams and \\l\\ in
centimetres. On a mizer size grid this puts the smallest resource sizes
at a fraction of a micrometre and a milligram organism at about a
millimetre, which is the right order for bacteria and copepods
respectively.

## Usage

``` r
resource_length_defaults
```

## Format

A list with entries `a` and `b`.

## Details

Note that this is a different convention from the one the species use: a
fish of a given weight is longer than a sphere of the same weight, by a
factor \\(a\_{fish}/a\_{resource})^{-1/3}\\, about 3.7 for the mizer
default `a = 0.01`. That difference is real rather than an artefact — a
1 mg copepod really is shorter than a 1 mg fish larva — but it does mean
the resource and the species sit on the plot at their own conventions.

## See also

[`resource_params()`](https://sizespectrum.org/mizer/reference/resource_params.md)
