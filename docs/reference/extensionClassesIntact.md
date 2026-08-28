# Is the S4 marker class chain for a set of extensions intact?

Checks that every dispatch extension has both its params and its sim
marker class defined and that each one still extends the next class down
the chain.

## Usage

``` r
extensionClassesIntact(extensions)
```

## Arguments

- extensions:

  Named character vector of extensions (full chain or dispatch subset).

## Value

`TRUE` if nothing needs repairing.
