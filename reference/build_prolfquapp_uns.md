# Build the prolfquapp uns namespace for round-trip AnnData conversion

Serializes `AnalysisConfiguration` and `ProteinAnnotation` fields into a
nested list suitable for storage in the AnnData `uns` slot.

## Usage

``` r
build_prolfquapp_uns(
  config,
  protAnnot,
  layer_names,
  source_software = "unknown"
)
```

## Arguments

- config:

  AnalysisConfiguration object

- protAnnot:

  ProteinAnnotation object

- layer_names:

  character vector of layer names

- source_software:

  character, name of the source software

## Value

list to be assigned to `adata$uns`
