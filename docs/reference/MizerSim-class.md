# A class to hold the results of a simulation

A class that holds the results of projecting a
[MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
object through time using
[`project()`](https://sizespectrum.org/mizer/reference/project.md).

## Details

A new `MizerSim` object can be created with the
[`MizerSim()`](https://sizespectrum.org/mizer/reference/MizerSim.md)
constructor, but you will never have to do that because the object is
created automatically by
[`project()`](https://sizespectrum.org/mizer/reference/project.md) when
needed.

As a user you should never have to access the slots of a MizerSim object
directly. Instead there are a range of functions to extract the
information. [`N()`](https://sizespectrum.org/mizer/reference/N.md) and
[`NResource()`](https://sizespectrum.org/mizer/reference/N.md) return
arrays with the saved abundances of the species and the resource
population at size respectively.
[`getEffort()`](https://sizespectrum.org/mizer/reference/getEffort.md)
returns the fishing effort of each gear through time.
[`getTimes()`](https://sizespectrum.org/mizer/reference/getTimes.md)
returns the vector of times at which simulation results were stored and
[`idxFinalT()`](https://sizespectrum.org/mizer/reference/finalN.md)
returns the index with which to access specifically the value at the
final time in the arrays returned by the other functions.
[`getParams()`](https://sizespectrum.org/mizer/reference/getParams.md)
extracts the ecosystem state as a `MizerParams` object with initial
abundances set to values from the simulation;
[`finalParams()`](https://sizespectrum.org/mizer/reference/finalParams.md)
and
[`initialParams()`](https://sizespectrum.org/mizer/reference/initialParams.md)
are convenient shorthands for the final and initial time steps. There
are also several
[summary_functions](https://sizespectrum.org/mizer/reference/summary_functions.md)
and
[plotting_functions](https://sizespectrum.org/mizer/reference/plotting_functions.md)
available to explore the contents of a `MizerSim` object.

The arrays all have named dimensions. The names of the `time` dimension
denote the time in years. The names of the `w` dimension are weights in
grams rounded to three significant figures. The names of the `sp`
dimension are the same as the species name in the order specified in the
species_params data frame. The names of the `gear` dimension are the
names of the gears, in the same order as specified when setting up the
`MizerParams` object.

Extensions of mizer can use the `n_other` slot to store the abundances
of other ecosystem components and these extensions should provide their
own functions for accessing that information.

The `MizerSim` class has changed since previous versions of mizer. To
use a `MizerSim` object created by a previous version, you need to
upgrade it with
[`validSim()`](https://sizespectrum.org/mizer/reference/validSim.md).

## Slots

- `params`:

  An object of type
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md).
  If this params object uses extensions, the `MizerSim` object uses the
  same extension chain via `params@extensions`; `MizerSim` has no
  separate extension slot.

- `n`:

  Three-dimensional array (time x species x size) that stores the
  projected community number densities.

- `n_pp`:

  An array (time x size) that stores the projected resource number
  densities.

- `n_other`:

  A list array (time x component) that stores the projected values for
  other ecosystem components.

- `effort`:

  An array (time x gear) that stores the fishing effort by time and
  gear.

- `sim_params`:

  A named list of the parameters passed to
  [`project()`](https://sizespectrum.org/mizer/reference/project.md) or
  [`projectToSteady()`](https://sizespectrum.org/mizer/reference/projectToSteady.md)
  to produce this simulation, such as `method` and `dt`.
