# Detect decoy/reverse identifiers (within-duplicate resolution)

Thin wrapper delegating to
[`prolfqua::is_decoy`](https://wolski.github.io/prolfqua/reference/is_decoy.html)
so annotation de-duplication and the quant layer share ONE detector:
built-in anchored default prefixes unioned with an optional configured
`pattern`; an empty / `NULL` / no-op (`"a^"`) pattern falls back to the
defaults only.

## Usage

``` r
.detect_decoy_ids(ids, pattern = NULL)
```

## Arguments

- ids:

  character vector of (prefixed) identifiers

- pattern:

  optional configured decoy regex

## Value

logical vector
