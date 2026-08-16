# Working with array wrapper classes (ArraySpeciesBySize, etc.)

Use this skill when writing or modifying functions that calculate, return, or
manipulate rate arrays and size/time-resolved model outputs in mizer.

## The wrapper classes

mizer provides lightweight S3 wrapper classes around multi-dimensional numeric arrays:

| Class | Underlying Shape | Key Functions Returning It |
|---|---|---|
| `ArraySpeciesBySize` | species × size | `getEncounter()`, `getFeedingLevel()`, `getPredRate()`, `getPredMort()`, `getMort()`, `getEReproAndGrowth()`, `getERepro()`, `getEGrowth()`, `getFMort()`, `getExtMort()`, `getExtEncounter()` |
| `ArrayTimeBySpecies` | time × species | `getBiomass()`, `getN()`, `getSSB()`, `getYield()`, `getYieldGear()` |
| `ArrayTimeBySpeciesBySize` | time × species × size | `N(sim)`, `getFeedingLevel(sim)`, `getFMort(sim)` |
| `ArrayResourceBySize` | 1 × size (or vector) | `getResourceMort()` |
| `ArrayTimeByResourceBySize` | time × 1 × size | `NResource(sim)` |

These classes attach metadata attributes (`value_name`, `units`, `type`, `params`, `representation`) so that outputs provide informative `print()`, `summary()`, `plot()`, and `as.data.frame()` methods.

## Core rules

### 1. Rate functions vs getter wrappers

- **Rate functions** (`mizerEncounter()`, `mizerMort()`, custom functions registered with `setRateFunction()`):
  Must return a **plain numeric matrix/array** with appropriate dimnames. They do **not** wrap their return value in `ArraySpeciesBySize`.
- **Public getter functions** (`getEncounter()`, `getMort()`, etc.):
  Are responsible for calling the underlying rate function and wrapping the numeric result with `ArraySpeciesBySize(...)`.
- **Standalone user-facing rate functions**:
  If writing a standalone function meant to be called directly by users (not dispatched through `getRateFunction()`), wrap the result using `ArraySpeciesBySize()`.

```r
# Inside a public getter wrapper:
getEncounter <- function(object, ...) {
    # ... compute plain numeric array `res` ...
    ArraySpeciesBySize(
        res,
        value_name = "Encounter rate",
        units = "g/year",
        type = "value",
        params = params
    )
}
```

### 2. Assignment to S4 slots: always use `slot[] <- value`

When storing an array (whether wrapped or plain) into a `MizerParams` or `MizerSim` slot, **always assign using `slot[] <- value`** (with empty brackets), never `slot <- value`:

```r
# CORRECT: Strips the S3 wrapper class and preserves existing slot dimensions/dimnames
params@metab[] <- value

# WRONG: Attaches S3 classes and attributes to the S4 slot, which violates class validity
params@metab <- value
```

### 3. Arithmetic strips the class; subsetting preserves it

- **Arithmetic (`+`, `-`, `*`, `/`, etc.)**:
  The group generic `Ops.ArraySpeciesBySize` (and sibling classes) **deliberately strips all wrapper attributes and class** using `unclass_rate()`, returning a plain matrix. This prevents attribute contamination (e.g. invalid `units` or `value_name`) during mathematical calculations.
- **Subsetting (`[`)**:
  Subsetting preserves the wrapper class and attributes as long as the result retains the required dimensionality (e.g. remaining a 2D matrix for `ArraySpeciesBySize`). If subsetting collapses a dimension, standard R matrix subsetting applies.

```r
enc <- getEncounter(NS_params)
is.ArraySpeciesBySize(enc)             # TRUE
is.ArraySpeciesBySize(enc["Cod", ])    # FALSE (collapsed to 1D vector)
is.ArraySpeciesBySize(enc["Cod", , drop = FALSE]) # TRUE (still 2D matrix)
is.ArraySpeciesBySize(enc * 2)         # FALSE (arithmetic returns plain matrix)
```

### 4. Attributes

- **`value_name`**: Human-readable label (e.g. `"Encounter rate"`, `"Biomass"`).
- **`units`**: Physical units string (e.g. `"g/year"`, `"1/year"`, `"g"`).
- **`type`**: `"value"` (default amount/rate), `"density"` (per gram body weight), or `"proportion"` (fraction 0–1). Tells `plot()` whether to apply the Jacobian when plotting against length (`size_axis = "l"`) or to use a fixed [0, 1] y-axis.
- **`representation`**: `"point"` (default, sampled at left bin edge) or `"average"` (finite-volume bin average). Bin averages are plotted at geometric bin midpoints when second-order quadrature (`second_order_w[["bin_average"]]`) is active.
- **`params`**: A reference to `MizerParams`. Used by `plot()` for species colours, linetypes, and natural size ranges (`w_min` to `w_max`).

### 5. `str()` safety

The `str()` method for these wrapper classes strips the `params` attribute before displaying the summary. This prevents `str(object)` from dumping all slots of the entire `MizerParams` model object to the console.

### 6. User-facing plotting and summaries

For how these wrapped arrays are plotted, summarized, and extracted in user workflows (including Jacobian adjustments, size axis conversions, and ggplot methods), consult `inst/skills/analyse-and-plot/SKILL.md`.
