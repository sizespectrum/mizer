# Constructor for the `MizerSim` class

A constructor for the `MizerSim` class. This is used by
[`project()`](https://sizespectrum.org/mizer/reference/project.md) to
create `MizerSim` objects of the right dimensions. It is not necessary
for users to use this constructor.

## Usage

``` r
MizerSim(params, t_dimnames = NA, t_max = 100, t_save = 1)
```

## Arguments

- params:

  a
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
  object

- t_dimnames:

  Numeric vector that is used for the time dimensions of the slots.
  Default = NA.

- t_max:

  The maximum time step of the simulation. Only used if t_dimnames = NA.
  Default value = 100.

- t_save:

  How often should the results of the simulation be stored. Only used if
  t_dimnames = NA. Default value = 1.

## Value

An object of type
[MizerSim](https://sizespectrum.org/mizer/reference/MizerSim-class.md)
