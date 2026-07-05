# Description of indicator functions

Mizer provides a range of functions to calculate indicators from a
MizerSim or MizerParams object.

## Details

When called with a `MizerSim` object, these functions return a time
series of values. When called with a `MizerParams` object, they return a
single value calculated from the initial abundances stored in the params
object.

A list of available indicator functions is given in the table below

|  |  |  |
|----|----|----|
| Function | Returns | Description |
| [`getProportionOfLargeFish()`](https://sizespectrum.org/mizer/reference/getProportionOfLargeFish.md) | A vector with values at each time step (or a single value for MizerParams). | Calculates the proportion of large fish through time. The threshold value can be specified. It is possible to calculation the proportion of large fish based on either length or weight. |
| [`getMeanWeight()`](https://sizespectrum.org/mizer/reference/getMeanWeight.md) | A vector with values at each saved time step (or a single value for MizerParams). | The mean weight of the community through time. This is calculated as the total biomass of the community divided by the total abundance. |
| [`getMeanMaxWeight()`](https://sizespectrum.org/mizer/reference/getMeanMaxWeight.md) | Depends on the measure argument. If measure = “both” then you get a matrix with two columns, one with values by numbers, the other with values by biomass at each saved time step (or a named vector for MizerParams). If measure = “numbers” or “biomass” you get a vector of the respective values at each saved time step (or a single value for MizerParams). | The mean maximum weight of the community through time. This can be calculated by numbers or by biomass. See the help file for more details. |
| [`getCommunitySlope()`](https://sizespectrum.org/mizer/reference/getCommunitySlope.md) | A data.frame with four columns: time step, slope, intercept and the coefficient of determination (or a single-row data.frame for MizerParams). | Calculates the slope of the community abundance spectrum through time by performing a linear regression on the logged total numerical abundance and logged body size. |

## See also

[summary_functions](https://sizespectrum.org/mizer/reference/summary_functions.md),
[plotting_functions](https://sizespectrum.org/mizer/reference/plotting_functions.md)
