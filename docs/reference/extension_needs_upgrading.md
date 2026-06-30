# Whether any extension recorded on an object needs upgrading

Returns TRUE if, for any extension recorded in the object's
`@extensions` slot whose package is installed, the recorded version
stamp is missing (`NA`) or older than the installed version of that
package. A missing stamp counts as needing an upgrade so that objects
created before extension-version tracking are brought up to date (and
stamped) on first use; this is safe because
[`runExtensionUpgrades()`](https://sizespectrum.org/mizer/reference/runExtensionUpgrades.md)
only calls an upgrade method if one is registered and such methods are
written to be idempotent.

## Usage

``` r
extension_needs_upgrading(params)
```

## Arguments

- params:

  A MizerParams object.

## Value

TRUE or FALSE.
