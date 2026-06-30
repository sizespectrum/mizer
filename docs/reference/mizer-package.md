# mizer: Multi-species size-based modelling in R

The mizer package implements multi-species size-based modelling in R. It
has been designed for modelling marine ecosystems.

## Details

Using mizer is relatively simple. There are three main stages:

1.  *Setting the model parameters*. This is done by creating an object
    of class
    [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md).
    This includes model parameters such as the life history parameters
    of each species, and the range of the size spectrum. There are
    several setup functions that help to create a MizerParams objects
    for particular types of models:

    - [`newSingleSpeciesParams()`](https://sizespectrum.org/mizer/reference/newSingleSpeciesParams.md)

    - [`newCommunityParams()`](https://sizespectrum.org/mizer/reference/newCommunityParams.md)

    - [`newTraitParams()`](https://sizespectrum.org/mizer/reference/newTraitParams.md)

    - [`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.md)

2.  *Running a simulation*. This is done by calling the
    [`project()`](https://sizespectrum.org/mizer/reference/project.md)
    function with the model parameters. This produces an object of
    [MizerSim](https://sizespectrum.org/mizer/reference/MizerSim-class.md)
    that contains the results of the simulation.

3.  *Exploring results*. After a simulation has been run, the results
    can be explored using a range of
    [plotting_functions](https://sizespectrum.org/mizer/reference/plotting_functions.md),
    [summary_functions](https://sizespectrum.org/mizer/reference/summary_functions.md)
    and
    [indicator_functions](https://sizespectrum.org/mizer/reference/indicator_functions.md).

See the [mizer website](https://sizespectrum.org/mizer/) for full
details of the principles behind mizer and how the package can be used
to perform size-based modelling.

## See also

Useful links:

- <https://sizespectrum.org/mizer/>

- <https://github.com/sizespectrum/mizer>

- Report bugs at <https://github.com/sizespectrum/mizer/issues>

## Author

**Maintainer**: Gustav Delius <gustav.delius@york.ac.uk>
([ORCID](https://orcid.org/0000-0003-4092-8228)) \[copyright holder\]

Authors:

- Gustav Delius <gustav.delius@york.ac.uk>
  ([ORCID](https://orcid.org/0000-0003-4092-8228)) \[copyright holder\]

- Finlay Scott <drfinlayscott@gmail.com> \[copyright holder\]

- Julia Blanchard <julia.blanchard@utas.edu.au>
  ([ORCID](https://orcid.org/0000-0003-0532-4824)) \[copyright holder\]

- Ken Andersen <kha@aqua.dtu.dk>
  ([ORCID](https://orcid.org/0000-0002-8478-3430)) \[copyright holder\]

Other contributors:

- Richard Southwell <richard.southwell@york.ac.uk> \[contributor,
  copyright holder\]
