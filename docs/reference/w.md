# Size bins

Functions to fetch information about the size bins used in the model
described by `params`.

## Usage

``` r
w(params)

w_full(params)

dw(params)

dw_full(params)
```

## Arguments

- params:

  A MizerParams object

## Value

`w()` returns a vector with the sizes at the start of each size bin of
the consumer spectrum.

`w_full()` returns a vector with the sizes at the start of each size bin
of the resource spectrum, which typically starts at smaller sizes than
the consumer spectrum.

`dw()` returns a vector with the widths of the size bins of the consumer
spectrum.

`dw_full()` returns a vector with the widths of the size bins of the
resource spectrum.

## Details

To represent the continuous size spectrum in the computer, the size
variable is discretized into a vector `w` of discrete weights, providing
a grid of sizes spanning the range from the smallest egg size to the
largest maximum size. These grid values divide the full size range into
a finite number of size bins. The size bins should be chosen small
enough to avoid the discretisation errors from becoming too big. You can
fetch this vector with `w()` and the vector of bin widths with `dw()`.

The weight grid is set up to be logarithmically spaced, so that
`w[j]=w[1]*10^(j*dx)` for some fixed `dx`. This means that the bin
widths increase with size: `dw[j] = w[j] * (10^dx - 1)`. This grid is
set up automatically when creating a MizerParams object.

Because the resource spectrum spans a larger range of sizes, these sizes
are discretized into a different vector of weights `w_full`. This
usually starts at a much smaller size than `w`, but also runs up to the
same largest size, so that the last entries of `w_full` have to coincide
with the entries of `w`. The logarithmic spacing for `w_full` is the
same as that for `w`, so that again `w_full[j]=w_full[1]*10^(j*dx)`. The
function `w_full()` gives the vector of sizes and `dw_full()` gives the
vector of bin widths.

You will need these vectors when converting number densities to numbers.
For example the size spectrum of a species is stored as a vector of
values that represent the *density* of fish in each size bin rather than
the *number* of fish. The number of fish in the size bin between `w[j]`
and `w[j+1]=w[j]+dw[j]` is obtained as `N[j]*dw[j]`.

The vector `w` can be used for example to convert the number of
individuals in a size bin into the biomass in the size bin. The biomass
in the `j`th bin is `biomass[j] = N[j] * dw[j] * w[j]`.

Of course all these calculations with discrete sizes and size bins are
only giving approximations to the continuous values, and these
approximations get better the smaller the size bins are, i.e., the more
size bins are used. However using more size bins also slows down the
calculations, so there is a trade-off. This is why the functions setting
up MizerParams objects allow you to choose the number of size bins
`no_w`.

## Examples

``` r
str(w(NS_params))
#>  num [1:100] 0.001 0.00119 0.00142 0.0017 0.00203 ...
str(dw(NS_params))
#>  num [1:100] 0.000193 0.000231 0.000275 0.000329 0.000392 ...
str(w_full(NS_params))
#>  num [1:226] 2.12e-13 2.53e-13 3.02e-13 3.61e-13 4.30e-13 ...
str(dw_full(NS_params))
#>  num [1:226] 4.10e-14 4.90e-14 5.84e-14 6.97e-14 8.32e-14 ...

# Calculating the biomass of Cod in each bin in the North Sea model
biomass <- initialN(NS_params)["Cod", ] * dw(NS_params) * w(NS_params)
# Summing to get total biomass
sum(biomass)
#> [1] 601591842546
```
